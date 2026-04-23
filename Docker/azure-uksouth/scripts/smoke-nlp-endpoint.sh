#!/bin/bash
#
# Flow-B smoke test: hit POST /api/compounds/nlp/query/ on a live instance
# with an AAD-authenticated bearer token and print the response.
#
# This exercises the full production stack (i2remote MSAL token → web-app
# ingress → Next.js proxy → Django server view → parse_prompt via managed
# identity → Azure OpenAI → executor → scoped AnalysisResult walk → JSON
# response back). Complements the direct-to-Azure golden eval in
# apps/compounds/nlp/tests/test_golden.py, which bypasses the ccp4i2 server.
#
# Usage:
#   smoke-nlp-endpoint.sh                                         # default prompt
#   smoke-nlp-endpoint.sh "Show HTRF IC50 for ARd compounds"      # custom prompt
#   smoke-nlp-endpoint.sh -u https://other.host/api/proxy/compounds "..."
#
# Env / config:
#   CCP4I2_API_TOKEN  Overrides the token read from ~/.ccp4i2remote.json.
#                     Convenient for scripted runs — grab one with:
#                         export CCP4I2_API_TOKEN=$(az account get-access-token \
#                           --resource <aad-client-id> --query accessToken -o tsv)
#   CCP4I2_API_URL    Overrides the api_url from ~/.ccp4i2remote.json. The
#                     script swaps trailing `/ccp4i2` for `/compounds` so the
#                     i2remote-style URL works as-is.
#
# Expectations:
#   • The target instance has COMPOUNDS_NLP_ENABLED=true set on the server
#     container (see NLP_QUERY_PROPOSAL.md §16.2 / §16.3). Otherwise you get
#     a 404 with status=disabled — which still proves the endpoint routes
#     correctly, so it's a useful signal.
#   • The `ddu-openai` resource has customSubDomainName + the role
#     assignment applied (scripts/configure-openai-for-nlp.sh). Otherwise
#     the server returns 502 with status=parse (token acquisition failed).

set -e

if [ -t 1 ]; then
    RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'
    BLUE='\033[0;34m'; NC='\033[0m'
else
    RED=''; GREEN=''; YELLOW=''; BLUE=''; NC=''
fi

# -- Parse args --------------------------------------------------------------
PROMPT_DEFAULT="Show all HTRF IC50 values for ARd compounds"
PROMPT=""
API_URL_OVERRIDE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -u|--url) API_URL_OVERRIDE="$2"; shift 2 ;;
        -h|--help)
            sed -n '2,/^set -e/p' "$0" | sed 's/^# \{0,1\}//' | sed '$d'
            exit 0 ;;
        *) PROMPT="$1"; shift ;;
    esac
done
PROMPT="${PROMPT:-$PROMPT_DEFAULT}"

# -- Resolve token -----------------------------------------------------------
CONFIG_FILE="$HOME/.ccp4i2remote.json"
TOKEN="${CCP4I2_API_TOKEN:-}"

if [ -z "$TOKEN" ]; then
    if [ ! -f "$CONFIG_FILE" ]; then
        echo -e "${RED}No token: \$CCP4I2_API_TOKEN is unset and $CONFIG_FILE doesn't exist.${NC}"
        echo "Run \`i2remote login\` first (or export CCP4I2_API_TOKEN)."
        exit 1
    fi
    TOKEN=$(python3 -c "import json, sys; c=json.load(open('$CONFIG_FILE')); print(c.get('api_token') or '', end='')")
    if [ -z "$TOKEN" ]; then
        echo -e "${RED}No api_token in $CONFIG_FILE. Run \`i2remote login\` to refresh.${NC}"
        exit 1
    fi
fi

# -- Resolve API URL and derive compounds base -------------------------------
API_URL="${API_URL_OVERRIDE:-${CCP4I2_API_URL:-}}"
if [ -z "$API_URL" ] && [ -f "$CONFIG_FILE" ]; then
    API_URL=$(python3 -c "import json, sys; c=json.load(open('$CONFIG_FILE')); print(c.get('api_url') or '', end='')")
fi
if [ -z "$API_URL" ]; then
    echo -e "${RED}No API URL: pass -u <url>, set CCP4I2_API_URL, or run \`i2remote config set api_url ...\`.${NC}"
    exit 1
fi

# Accept either /api/proxy/ccp4i2 (i2remote convention) or a bare /api/proxy/compounds.
API_URL="${API_URL%/}"
if [[ "$API_URL" == */ccp4i2 ]]; then
    COMPOUNDS_BASE="${API_URL%/ccp4i2}/compounds"
elif [[ "$API_URL" == */compounds ]]; then
    COMPOUNDS_BASE="$API_URL"
else
    echo -e "${YELLOW}Warning: api_url doesn't end with /ccp4i2 or /compounds — using as-is.${NC}"
    COMPOUNDS_BASE="$API_URL"
fi

# The Next.js proxy appends the trailing slash when forwarding to Django, so
# clients MUST post without one (otherwise Next.js 308-redirects — see
# i2remote.py:470 for the same convention).
ENDPOINT="${COMPOUNDS_BASE}/nlp/query"

# -- POST --------------------------------------------------------------------
echo -e "${GREEN}→ POST${NC} $ENDPOINT"
echo -e "  prompt: ${BLUE}$PROMPT${NC}"
echo ""

HTTP_OUT=$(mktemp)
HTTP_CODE=$(curl -sS -w "%{http_code}" -o "$HTTP_OUT" \
    -X POST \
    -H "Authorization: Bearer $TOKEN" \
    -H "Content-Type: application/json" \
    --data "$(python3 -c 'import json,sys; print(json.dumps({"prompt": sys.argv[1]}))' "$PROMPT")" \
    "$ENDPOINT") || {
    echo -e "${RED}curl failed${NC}"
    rm -f "$HTTP_OUT"
    exit 1
}

# -- Report ------------------------------------------------------------------
case "$HTTP_CODE" in
    2??) echo -e "${GREEN}← HTTP $HTTP_CODE${NC}" ;;
    4??) echo -e "${YELLOW}← HTTP $HTTP_CODE${NC}" ;;
    *)   echo -e "${RED}← HTTP $HTTP_CODE${NC}" ;;
esac
echo ""

# Pretty-print JSON when possible; fall back to raw.
if python3 -c "import json,sys; json.load(open('$HTTP_OUT'))" 2>/dev/null; then
    python3 -m json.tool < "$HTTP_OUT"
else
    cat "$HTTP_OUT"
    echo ""
fi

rm -f "$HTTP_OUT"
case "$HTTP_CODE" in
    2??) exit 0 ;;
    *)   exit 1 ;;
esac
