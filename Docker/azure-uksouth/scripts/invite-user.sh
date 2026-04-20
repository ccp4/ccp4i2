#!/bin/bash
# Invite an external user to a CCP4i2 instance.
#
# Usage:
#   ./invite-user.sh <email-address>
#       Invite a single user (ad-hoc).
#
#   ./invite-user.sh --drain --server-url <https://...> [--group-id <GUID>]
#       Pull pending invites from the server's /pending-invites queue and
#       send each one, marking each row sent/failed. Intended for a
#       privileged operator to run periodically.
#
#   ./invite-user.sh --list --server-url <https://...>
#       Preview the pending queue without sending anything.
#
# Why this exists: Newcastle IT will not consent to Graph application
# permissions on the ccp4i2 app registrations, so the container apps cannot
# send B2B invitations directly. Instead, admins queue invites via the
# admin UI, and this script — run by a Newcastle staff user — drains the
# queue via their own delegated Graph permissions.
#
# Prerequisites:
# - `az login` as a user with permission to invite guests in the target
#   tenant and to add members to the access group.
# - The script assumes your Azure CLI session is already pointed at the
#   correct tenant (use `az login --tenant <id>` if you have multiple).

set -e

# Configuration (defaults target the ccp4i2-demo instance)
DEMO_GROUP_ID="4a3dbe2a-c700-4d29-a6a8-e34daf1aeec7"  # ccp4i2-demo-users
DEMO_URL="https://ccp4i2-demo-web.agreeablebeach-f023b4ee.uksouth.azurecontainerapps.io"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

usage() {
    echo "Usage:"
    echo "  $0 <email-address>                              Invite one user ad-hoc"
    echo "  $0 --drain --server-url <URL> [--group-id <GUID>]   Drain pending queue"
    echo "  $0 --list  --server-url <URL>                   Show pending queue"
    echo ""
    echo "  --group-id defaults to $DEMO_GROUP_ID"
    echo "  --server-url must be the instance's web URL (used as invite redirect too)"
    exit 1
}

# Send one invite + add to group. Usage: send_invite <email> <redirect_url> <group_id>
# Sets globals INVITE_STATUS (sent|failed) and INVITE_DETAIL.
send_invite() {
    local email="$1"
    local redirect="$2"
    local group="$3"

    INVITE_STATUS="failed"
    INVITE_DETAIL=""
    GUEST_OID=""

    local existing
    existing=$(az ad user list --filter "mail eq '${email}' or otherMails/any(m:m eq '${email}')" --query "[0].id" -o tsv 2>/dev/null || true)

    if [ -n "$existing" ]; then
        GUEST_OID="$existing"
    else
        local invite_result
        invite_result=$(az rest --method POST \
            --url "https://graph.microsoft.com/v1.0/invitations" \
            --headers "Content-Type=application/json" \
            --body "{
                \"invitedUserEmailAddress\": \"${email}\",
                \"inviteRedirectUrl\": \"${redirect}\",
                \"sendInvitationMessage\": true,
                \"invitedUserMessageInfo\": {
                    \"customizedMessageBody\": \"You have been invited to use CCP4i2, the Crystallographic Computing Environment. Click the link below to accept the invitation and access the application.\"
                }
            }" 2>&1) || {
            INVITE_DETAIL="$invite_result"
            return 1
        }

        GUEST_OID=$(echo "$invite_result" | python3 -c "import sys,json; print(json.load(sys.stdin)['invitedUser']['id'])" 2>/dev/null || true)
        if [ -z "$GUEST_OID" ]; then
            INVITE_DETAIL="Invitation response missing invitedUser.id: $invite_result"
            return 1
        fi
    fi

    local is_member
    is_member=$(az ad group member check --group "$group" --member-id "$GUEST_OID" --query "value" -o tsv 2>/dev/null || echo "false")

    if [ "$is_member" != "true" ]; then
        local tries=3
        local i=1
        while [ $i -le $tries ]; do
            if az ad group member add --group "$group" --member-id "$GUEST_OID" 2>&1 >/dev/null; then
                break
            fi
            if [ $i -eq $tries ]; then
                INVITE_DETAIL="Failed to add $GUEST_OID to group $group after $tries attempts"
                return 1
            fi
            sleep 5
            i=$((i+1))
        done
    fi

    INVITE_STATUS="sent"
    return 0
}

# ── Drain-queue and list-queue modes ────────────────────────────────────

DRAIN=false
LIST=false
SERVER_URL=""
GROUP_ID="$DEMO_GROUP_ID"

while [ $# -gt 0 ]; do
    case "$1" in
        --drain)       DRAIN=true; shift ;;
        --list)        LIST=true; shift ;;
        --server-url)  SERVER_URL="$2"; shift 2 ;;
        --group-id)    GROUP_ID="$2"; shift 2 ;;
        -h|--help)     usage ;;
        *)             POSITIONAL="$1"; shift ;;
    esac
done

if [ "$DRAIN" = true ] || [ "$LIST" = true ]; then
    if [ -z "$SERVER_URL" ]; then
        echo -e "${RED}--server-url is required for --drain/--list${NC}"
        usage
    fi

    # Acquire a user access token for the instance's app registration.
    # The token is sent to /api/proxy/users/pending-invites to authenticate.
    # Assumes the CLI session is in the same tenant as the app.
    echo -e "${YELLOW}Fetching pending invites from ${SERVER_URL}...${NC}"
    TOKEN=$(az account get-access-token --resource "https://graph.microsoft.com" --query accessToken -o tsv)
    QUEUE_JSON=$(curl -sS -H "Authorization: Bearer $TOKEN" "${SERVER_URL}/api/proxy/users/pending-invites/?status=pending")

    COUNT=$(echo "$QUEUE_JSON" | python3 -c "import sys,json; print(len(json.load(sys.stdin)))" 2>/dev/null || echo "0")
    echo -e "${YELLOW}${COUNT} pending invite(s)${NC}"

    if [ "$LIST" = true ]; then
        echo "$QUEUE_JSON" | python3 -m json.tool
        exit 0
    fi

    if [ "$COUNT" = "0" ]; then
        echo -e "${GREEN}Nothing to do.${NC}"
        exit 0
    fi

    # Iterate each invite
    echo "$QUEUE_JSON" | python3 -c "import sys,json; [print(f\"{r['id']}\t{r['email']}\") for r in json.load(sys.stdin)]" | while IFS=$'\t' read -r ID EMAIL; do
        echo ""
        echo -e "${YELLOW}→ Invite #${ID} ${EMAIL}${NC}"

        if send_invite "$EMAIL" "$SERVER_URL" "$GROUP_ID"; then
            echo -e "${GREEN}  ✓ Sent (guest ${GUEST_OID})${NC}"
            curl -sS -X POST -H "Authorization: Bearer $TOKEN" -H "Content-Type: application/json" \
                -d "{\"guest_object_id\": \"$GUEST_OID\"}" \
                "${SERVER_URL}/api/proxy/users/pending-invites/${ID}/mark_sent/" > /dev/null
        else
            echo -e "${RED}  ✗ Failed: ${INVITE_DETAIL}${NC}"
            REASON_JSON=$(python3 -c "import json,sys; print(json.dumps({'reason': sys.argv[1]}))" "$INVITE_DETAIL")
            curl -sS -X POST -H "Authorization: Bearer $TOKEN" -H "Content-Type: application/json" \
                -d "$REASON_JSON" \
                "${SERVER_URL}/api/proxy/users/pending-invites/${ID}/mark_failed/" > /dev/null
        fi
    done

    echo ""
    echo -e "${GREEN}Done.${NC}"
    exit 0
fi

# ── Ad-hoc mode (single email on argv) ──────────────────────────────────

EMAIL="${POSITIONAL:-}"
if [ -z "$EMAIL" ]; then
    usage
fi

echo -e "${YELLOW}Inviting ${EMAIL} to CCP4i2 demo...${NC}"
if send_invite "$EMAIL" "$DEMO_URL" "$DEMO_GROUP_ID"; then
    echo -e "${GREEN}Done!${NC}"
    echo -e "  Email:    ${EMAIL}"
    echo -e "  Guest ID: ${GUEST_OID}"
    echo -e "  Group:    ${GROUP_ID}"
    echo -e "  URL:      ${DEMO_URL}"
else
    echo -e "${RED}Invite failed: ${INVITE_DETAIL}${NC}"
    exit 1
fi
