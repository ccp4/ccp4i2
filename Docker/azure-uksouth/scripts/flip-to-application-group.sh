#!/bin/bash
# Narrow an app registration's groups claim from "All" (or "SecurityGroup") to
# "ApplicationGroup", so JWT tokens only carry groups explicitly assigned to
# this enterprise app.
#
# Why: with groupMembershipClaims="All", users in many AAD groups get huge
# JWTs (10+ KB). When those JWTs are used as query-string tokens for file
# downloads, the total URL overflows the ingress header buffer and the
# browser sees HTTP 431 "Request Header Fields Too Large".
#
# Safety: before flipping the manifest, this script first assigns every
# group in ALLOWED_AZURE_AD_GROUPS to the enterprise app. If we flipped
# without doing that, the backend's group-membership check would reject
# every user (since the groups would no longer appear in the token) — a
# repeat of the 2026-04-12 outage.
#
# Usage:
#   ./flip-to-application-group.sh [--env <env-file>] [--dry-run] [--app-role-id <guid>]
#
#   --env <env-file>     env file under Docker/azure-uksouth/ (default: .env.deployment)
#   --dry-run            show what would change without making changes
#   --app-role-id <guid> app role to grant to each allowed group. If omitted:
#                        uses the zero GUID if the app has no custom roles,
#                        the sole custom role if it has one, otherwise errors
#                        and lists the available roles so you can pick.
#
# Prerequisites:
#   - `az login` as a user who owns the app registration (and can assign
#     principals to the enterprise app)
#   - env file sets AZURE_AD_CLIENT_ID and ALLOWED_AZURE_AD_GROUPS

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BICEP_DIR="$(dirname "$SCRIPT_DIR")"

ENV_FILE="$BICEP_DIR/.env.deployment"
DRY_RUN=false
APP_ROLE_ID=""  # empty = auto-select (zero GUID if app has no roles, sole role if it has one, else prompt)

while [ $# -gt 0 ]; do
    case "$1" in
        --env)          ENV_FILE="$BICEP_DIR/$2"; shift 2 ;;
        --dry-run)      DRY_RUN=true; shift ;;
        --app-role-id)  APP_ROLE_ID="$2"; shift 2 ;;
        -h|--help)
            sed -n '1,30p' "$0" | sed 's/^# \{0,1\}//'
            exit 0 ;;
        *)
            echo -e "${RED}Unknown argument: $1${NC}"
            exit 1 ;;
    esac
done

if [ ! -f "$ENV_FILE" ]; then
    echo -e "${RED}Env file not found: $ENV_FILE${NC}"
    exit 1
fi

echo -e "${YELLOW}Using env file: $(basename "$ENV_FILE")${NC}"
# shellcheck disable=SC1090
source "$ENV_FILE"

if [ -z "${AZURE_AD_CLIENT_ID:-}" ]; then
    echo -e "${RED}AZURE_AD_CLIENT_ID not set in $ENV_FILE${NC}"
    exit 1
fi
if [ -z "${ALLOWED_AZURE_AD_GROUPS:-}" ]; then
    echo -e "${RED}ALLOWED_AZURE_AD_GROUPS not set in $ENV_FILE${NC}"
    exit 1
fi

CLIENT_ID="$AZURE_AD_CLIENT_ID"

if [ "$DRY_RUN" = true ]; then
    echo -e "${BLUE}DRY RUN — no changes will be made${NC}"
fi
echo ""

# ── Resolve app registration (object id) and service principal (object id) ──

echo -e "${YELLOW}Resolving app registration $CLIENT_ID...${NC}"
APP_OBJECT_ID=$(az ad app show --id "$CLIENT_ID" --query id -o tsv 2>/dev/null || true)
if [ -z "$APP_OBJECT_ID" ]; then
    echo -e "${RED}  App registration not found. Are you logged in to the right tenant?${NC}"
    exit 1
fi
APP_DISPLAY_NAME=$(az ad app show --id "$CLIENT_ID" --query displayName -o tsv)
echo -e "  ${GREEN}✓ $APP_DISPLAY_NAME  (objectId=$APP_OBJECT_ID)${NC}"

echo -e "${YELLOW}Resolving service principal for $CLIENT_ID...${NC}"
SP_ID=$(az ad sp show --id "$CLIENT_ID" --query id -o tsv 2>/dev/null || true)
if [ -z "$SP_ID" ]; then
    echo -e "${RED}  No service principal found. Create one with:${NC}"
    echo -e "    az ad sp create --id $CLIENT_ID"
    exit 1
fi
echo -e "  ${GREEN}✓ servicePrincipal objectId=$SP_ID${NC}"

# ── Resolve which app role to grant ──

if [ -z "$APP_ROLE_ID" ]; then
    APP_ROLES_JSON=$(az ad app show --id "$CLIENT_ID" --query "appRoles[?isEnabled].{id:id, value:value, displayName:displayName}" -o json)
    ROLE_COUNT=$(echo "$APP_ROLES_JSON" | python3 -c "import sys,json; print(len(json.load(sys.stdin)))")

    if [ "$ROLE_COUNT" = "0" ]; then
        APP_ROLE_ID="00000000-0000-0000-0000-000000000000"
        echo -e "${YELLOW}App has no custom roles — using zero-GUID default.${NC}"
    elif [ "$ROLE_COUNT" = "1" ]; then
        APP_ROLE_ID=$(echo "$APP_ROLES_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin)[0]['id'])")
        ROLE_VALUE=$(echo "$APP_ROLES_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin)[0]['value'])")
        echo -e "${YELLOW}App has one custom role ($ROLE_VALUE) — using $APP_ROLE_ID.${NC}"
    else
        echo -e "${RED}App has multiple enabled roles. Pick one with --app-role-id:${NC}"
        echo "$APP_ROLES_JSON" | python3 -m json.tool | sed 's/^/  /'
        exit 1
    fi
fi

# ── Read current groupMembershipClaims ──

CURRENT_CLAIMS=$(az ad app show --id "$CLIENT_ID" --query groupMembershipClaims -o tsv)
echo ""
echo -e "${YELLOW}Current groupMembershipClaims:${NC} ${CURRENT_CLAIMS:-<null>}"
echo -e "${YELLOW}Target groupMembershipClaims:${NC}  ApplicationGroup"
echo ""

# ── Process each allowed group: verify + assign if needed ──

IFS=',' read -ra GROUP_IDS <<< "$ALLOWED_AZURE_AD_GROUPS"

for GROUP_ID in "${GROUP_IDS[@]}"; do
    GROUP_ID=$(echo "$GROUP_ID" | xargs)  # trim whitespace
    [ -z "$GROUP_ID" ] && continue

    echo -e "${YELLOW}Group $GROUP_ID${NC}"

    GROUP_JSON=$(az ad group show --group "$GROUP_ID" 2>/dev/null || true)
    if [ -z "$GROUP_JSON" ]; then
        echo -e "  ${RED}✗ Group not found in this tenant${NC}"
        exit 1
    fi

    GROUP_NAME=$(echo "$GROUP_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('displayName',''))")
    SEC_ENABLED=$(echo "$GROUP_JSON" | python3 -c "import sys,json; print(json.load(sys.stdin).get('securityEnabled', False))")
    GROUP_TYPES=$(echo "$GROUP_JSON" | python3 -c "import sys,json; print(','.join(json.load(sys.stdin).get('groupTypes') or []) or 'Security')")
    echo -e "  displayName:     $GROUP_NAME"
    echo -e "  securityEnabled: $SEC_ENABLED"
    echo -e "  groupTypes:      $GROUP_TYPES"

    # Check if already assigned to the enterprise app (client-side filter —
    # Graph's OData $filter on principalId is finicky with GUIDs).
    ALL_ASSIGNMENTS_JSON=$(az rest --method GET \
        --url "https://graph.microsoft.com/v1.0/servicePrincipals/$SP_ID/appRoleAssignedTo" \
        --query "value[?principalId=='$GROUP_ID'].{id:id, appRoleId:appRoleId}" -o json)
    EXISTING_ASSIGNMENT_ID=$(echo "$ALL_ASSIGNMENTS_JSON" | python3 -c "import sys,json; v=json.load(sys.stdin); print(v[0]['id'] if v else '')")
    EXISTING_ROLE_ID=$(echo "$ALL_ASSIGNMENTS_JSON" | python3 -c "import sys,json; v=json.load(sys.stdin); print(v[0]['appRoleId'] if v else '')")

    if [ -n "$EXISTING_ASSIGNMENT_ID" ]; then
        echo -e "  ${GREEN}✓ Already assigned (assignmentId=$EXISTING_ASSIGNMENT_ID, appRoleId=$EXISTING_ROLE_ID)${NC}"
    else
        if [ "$DRY_RUN" = true ]; then
            echo -e "  ${BLUE}→ Would POST appRoleAssignment (principalId=$GROUP_ID, appRoleId=$APP_ROLE_ID)${NC}"
        else
            echo -e "  ${YELLOW}→ Assigning group to enterprise app...${NC}"
            RESULT=$(az rest --method POST \
                --url "https://graph.microsoft.com/v1.0/servicePrincipals/$SP_ID/appRoleAssignments" \
                --headers "Content-Type=application/json" \
                --body "{\"principalId\": \"$GROUP_ID\", \"resourceId\": \"$SP_ID\", \"appRoleId\": \"$APP_ROLE_ID\"}" 2>&1) || {
                echo -e "  ${RED}✗ Assignment failed:${NC}"
                echo "$RESULT" | sed 's/^/    /'
                exit 1
            }
            NEW_ASSIGNMENT_ID=$(echo "$RESULT" | python3 -c "import sys,json; print(json.load(sys.stdin).get('id',''))" 2>/dev/null || echo "?")
            echo -e "  ${GREEN}✓ Assigned (assignmentId=$NEW_ASSIGNMENT_ID)${NC}"
        fi
    fi
    echo ""
done

# ── Flip the manifest ──

if [ "$CURRENT_CLAIMS" = "ApplicationGroup" ]; then
    echo -e "${GREEN}groupMembershipClaims is already ApplicationGroup — no manifest change needed.${NC}"
else
    if [ "$DRY_RUN" = true ]; then
        echo -e "${BLUE}→ Would PATCH application/$APP_OBJECT_ID: groupMembershipClaims=\"ApplicationGroup\"${NC}"
    else
        echo -e "${YELLOW}Patching app manifest: groupMembershipClaims → ApplicationGroup...${NC}"
        az rest --method PATCH \
            --url "https://graph.microsoft.com/v1.0/applications/$APP_OBJECT_ID" \
            --headers "Content-Type=application/json" \
            --body '{"groupMembershipClaims": "ApplicationGroup"}' > /dev/null
        echo -e "${GREEN}✓ Manifest patched${NC}"
    fi
fi
echo ""

# ── Verify ──

if [ "$DRY_RUN" = true ]; then
    echo -e "${BLUE}Dry run complete. Re-run without --dry-run to apply.${NC}"
    exit 0
fi

echo -e "${YELLOW}Reading back final state...${NC}"
FINAL_CLAIMS=$(az ad app show --id "$CLIENT_ID" --query groupMembershipClaims -o tsv)
echo -e "  groupMembershipClaims: ${GREEN}$FINAL_CLAIMS${NC}"

ASSIGNMENTS_JSON=$(az rest --method GET \
    --url "https://graph.microsoft.com/v1.0/servicePrincipals/$SP_ID/appRoleAssignedTo" \
    --query "value[].{principalId:principalId, principalType:principalType, principalDisplayName:principalDisplayName}")
echo -e "  appRoleAssignedTo:"
echo "$ASSIGNMENTS_JSON" | python3 -m json.tool | sed 's/^/    /'

echo ""
echo -e "${GREEN}Done.${NC}"
echo -e "${YELLOW}Note:${NC} AAD can take a few minutes to propagate. Users will need to"
echo -e "      re-login (or wait for token refresh) before their JWTs shrink."
echo -e "      After that, decode a token and confirm \`groups\` contains only:"
for GROUP_ID in "${GROUP_IDS[@]}"; do
    GROUP_ID=$(echo "$GROUP_ID" | xargs)
    [ -n "$GROUP_ID" ] && echo -e "        $GROUP_ID"
done
