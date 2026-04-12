#!/bin/bash
# Invite an external user to the CCP4i2 demo instance
#
# Usage: ./invite-user.sh <email-address>
#
# This script:
# 1. Sends a B2B guest invitation (user gets an email with a link)
# 2. Adds the guest to the ccp4i2-demo-users security group
#
# The user will be able to sign in from any Azure AD / Entra ID tenant,
# or via email one-time passcode if their org doesn't use Entra ID.
#
# Prerequisites:
# - az login (signed in as owner of the ccp4i2-demo app registration)
# - Permission to invite guests (Newcastle tenant allows this for all users)

set -e

# Configuration
DEMO_GROUP_ID="4a3dbe2a-c700-4d29-a6a8-e34daf1aeec7"  # ccp4i2-demo-users
DEMO_URL="https://ccp4i2-demo-web.agreeablebeach-f023b4ee.uksouth.azurecontainerapps.io"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

if [ -z "$1" ]; then
    echo -e "${RED}Usage: $0 <email-address>${NC}"
    echo ""
    echo "Examples:"
    echo "  $0 colleague@york.ac.uk"
    echo "  $0 researcher@diamond.ac.uk"
    echo "  $0 someone@gmail.com"
    exit 1
fi

EMAIL="$1"

echo -e "${YELLOW}Inviting ${EMAIL} to CCP4i2 demo...${NC}"

# Check if user already exists as a guest
EXISTING_USER=$(az ad user list --filter "mail eq '${EMAIL}' or otherMails/any(m:m eq '${EMAIL}')" --query "[0].id" -o tsv 2>/dev/null)

if [ -n "$EXISTING_USER" ]; then
    echo -e "${YELLOW}User already exists in directory (id: ${EXISTING_USER})${NC}"
    GUEST_OID="$EXISTING_USER"
else
    # Send invitation
    echo -e "${YELLOW}Sending invitation email...${NC}"
    INVITE_RESULT=$(az rest --method POST \
        --url "https://graph.microsoft.com/v1.0/invitations" \
        --headers "Content-Type=application/json" \
        --body "{
            \"invitedUserEmailAddress\": \"${EMAIL}\",
            \"inviteRedirectUrl\": \"${DEMO_URL}\",
            \"sendInvitationMessage\": true,
            \"invitedUserMessageInfo\": {
                \"customizedMessageBody\": \"You have been invited to use CCP4i2, the Crystallographic Computing Environment. Click the link below to accept the invitation and access the application.\"
            }
        }" 2>&1)

    if echo "$INVITE_RESULT" | grep -q "error"; then
        echo -e "${RED}Failed to send invitation:${NC}"
        echo "$INVITE_RESULT"
        exit 1
    fi

    GUEST_OID=$(echo "$INVITE_RESULT" | python3 -c "import sys,json; print(json.load(sys.stdin)['invitedUser']['id'])")
    echo -e "${GREEN}Invitation sent. Guest ID: ${GUEST_OID}${NC}"
fi

# Add to demo users group
echo -e "${YELLOW}Adding to ccp4i2-demo-users group...${NC}"

# Check if already a member
IS_MEMBER=$(az ad group member check --group "$DEMO_GROUP_ID" --member-id "$GUEST_OID" --query "value" -o tsv 2>/dev/null)

if [ "$IS_MEMBER" = "true" ]; then
    echo -e "${GREEN}Already a member of ccp4i2-demo-users${NC}"
else
    # Guest accounts can take a few seconds to propagate
    local RETRIES=3
    local DELAY=5
    for i in $(seq 1 $RETRIES); do
        az ad group member add --group "$DEMO_GROUP_ID" --member-id "$GUEST_OID" 2>&1 && break
        if [ $i -lt $RETRIES ]; then
            echo -e "${YELLOW}Waiting ${DELAY}s for guest account to propagate (attempt $i/$RETRIES)...${NC}"
            sleep $DELAY
        else
            echo -e "${RED}Failed to add to group after $RETRIES attempts. Run manually:${NC}"
            echo -e "  az ad group member add --group $DEMO_GROUP_ID --member-id $GUEST_OID"
            exit 1
        fi
    done
    echo -e "${GREEN}Added to ccp4i2-demo-users group${NC}"
fi

echo ""
echo -e "${GREEN}Done!${NC}"
echo -e "  Email:    ${EMAIL}"
echo -e "  Guest ID: ${GUEST_OID}"
echo -e "  Group:    ccp4i2-demo-users"
echo -e "  URL:      ${DEMO_URL}"
echo ""
echo -e "${YELLOW}The user will receive an invitation email. Once they accept,"
echo -e "they can sign in at the URL above.${NC}"
