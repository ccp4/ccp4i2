#!/bin/bash

# CCP4i2 API Client
# A CLI tool to interact with the CCP4i2 API for testing and debugging
#
# Usage:
#   ./api-client.sh [command] [options]
#
# Commands:
#   login          - Get Azure AD token (interactive)
#   projects       - List all projects
#   jobs           - List jobs for a project
#   job-status     - Get status of a specific job
#   clone-job      - Clone and optionally run a job
#   run-job        - Run an existing job
#   logs           - Get worker logs
#
# Examples:
#   ./api-client.sh login
#   ./api-client.sh projects
#   ./api-client.sh jobs --project <project-uuid>
#   ./api-client.sh clone-job --job <job-uuid> --run
#   ./api-client.sh job-status --job <job-uuid>

set -e

# Load environment
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_FILE="${SCRIPT_DIR}/../.env.deployment"

if [[ -f "$ENV_FILE" ]]; then
    source "$ENV_FILE"
fi

# Configuration - API_BASE_URL will be set after deployment, no default hardcoded
API_BASE_URL="${API_BASE_URL:-}"
TOKEN_CACHE_FILE="${SCRIPT_DIR}/.api-token-cache"
AAD_CLIENT_ID="${NEXT_PUBLIC_AAD_CLIENT_ID:-}"
AAD_TENANT_ID="${NEXT_PUBLIC_AAD_TENANT_ID:-}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Helper functions
log_info() { echo -e "${BLUE}ℹ️  $1${NC}"; }
log_success() { echo -e "${GREEN}✅ $1${NC}"; }
log_warn() { echo -e "${YELLOW}⚠️  $1${NC}"; }
log_error() { echo -e "${RED}❌ $1${NC}"; }

show_help() {
    cat << EOF
CCP4i2 API Client - Test and debug the CCP4i2 API

Usage: $(basename "$0") <command> [options]

Commands:
  login              Get Azure AD token (interactive browser login)
  token              Display current token (for debugging)
  projects           List all projects
  project            Get details of a specific project
  jobs               List jobs for a project
  job                Get details of a specific job
  job-status         Get status of a specific job
  clone-job          Clone a job (optionally run it)
  run-job            Run/submit an existing job
  logs               Show worker container logs

Options:
  --project, -p      Project UUID
  --job, -j          Job UUID
  --run, -r          Run job after cloning (for clone-job)
  --api-url          Override API base URL
  --help, -h         Show this help

Environment Variables:
  API_BASE_URL                  Base URL for the API
  NEXT_PUBLIC_AAD_CLIENT_ID     Azure AD Client ID
  NEXT_PUBLIC_AAD_TENANT_ID     Azure AD Tenant ID

Examples:
  # Login and cache token
  $(basename "$0") login

  # List projects
  $(basename "$0") projects

  # List jobs in a project
  $(basename "$0") jobs -p 5f9d8d47-dddc-4898-9048-95720ff10ae7

  # Clone and run a job
  $(basename "$0") clone-job -j 48141dfe-2a52-4d27-872a-dff86aeb9993 --run

  # Check job status
  $(basename "$0") job-status -j 48141dfe-2a52-4d27-872a-dff86aeb9993
EOF
}

# Get or refresh Azure AD token
get_token() {
    # Check if we have a cached token that's still valid
    if [[ -f "$TOKEN_CACHE_FILE" ]]; then
        local cached_token=$(cat "$TOKEN_CACHE_FILE" 2>/dev/null)
        local cached_expiry=$(echo "$cached_token" | jq -r '.expires_on // 0' 2>/dev/null)
        local current_time=$(date +%s)

        if [[ "$cached_expiry" -gt "$current_time" ]]; then
            echo "$cached_token" | jq -r '.access_token'
            return 0
        fi
    fi

    log_warn "No valid cached token. Run '$(basename "$0") login' first."
    return 1
}

do_login() {
    if [[ -z "$AAD_CLIENT_ID" ]] || [[ -z "$AAD_TENANT_ID" ]]; then
        log_error "AAD_CLIENT_ID and AAD_TENANT_ID must be set"
        log_info "Set NEXT_PUBLIC_AAD_CLIENT_ID and NEXT_PUBLIC_AAD_TENANT_ID in .env.deployment"
        exit 1
    fi

    log_info "Logging in via Azure AD..."
    log_info "Client ID: $AAD_CLIENT_ID"
    log_info "Tenant ID: $AAD_TENANT_ID"

    # Use Azure CLI to get a token for the application
    # The scope should be the API's app ID URI or the client ID with /.default
    local scope="api://${AAD_CLIENT_ID}/.default"

    log_info "Requesting token for scope: $scope"

    # Try using az cli to get token
    local token_response
    token_response=$(az account get-access-token \
        --resource "api://${AAD_CLIENT_ID}" \
        --tenant "$AAD_TENANT_ID" \
        2>/dev/null) || {
        log_warn "Could not get token via az cli resource method, trying alternative..."

        # Alternative: Use device code flow or browser login
        # For now, let's try a direct MSAL approach using Python
        token_response=$(python3 << EOF
import json
import sys
try:
    from msal import PublicClientApplication

    client_id = "$AAD_CLIENT_ID"
    tenant_id = "$AAD_TENANT_ID"
    authority = f"https://login.microsoftonline.com/{tenant_id}"
    scopes = [f"api://{client_id}/access_as_user"]

    app = PublicClientApplication(client_id, authority=authority)

    # Try interactive login
    result = app.acquire_token_interactive(scopes=scopes)

    if "access_token" in result:
        import time
        output = {
            "access_token": result["access_token"],
            "expires_on": int(time.time()) + result.get("expires_in", 3600)
        }
        print(json.dumps(output))
    else:
        print(json.dumps({"error": result.get("error_description", "Unknown error")}), file=sys.stderr)
        sys.exit(1)
except ImportError:
    print(json.dumps({"error": "msal not installed. Run: pip install msal"}), file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(json.dumps({"error": str(e)}), file=sys.stderr)
    sys.exit(1)
EOF
        ) || {
            log_error "Failed to get token"
            log_info "Make sure you have 'msal' installed: pip install msal"
            exit 1
        }
    }

    # Cache the token
    echo "$token_response" > "$TOKEN_CACHE_FILE"
    chmod 600 "$TOKEN_CACHE_FILE"

    log_success "Login successful! Token cached."
    log_info "Token expires at: $(echo "$token_response" | jq -r '.expires_on // .expiresOn' | xargs -I{} date -r {} 2>/dev/null || echo "unknown")"
}

show_token() {
    if [[ -f "$TOKEN_CACHE_FILE" ]]; then
        local token=$(cat "$TOKEN_CACHE_FILE")
        log_info "Cached token info:"
        echo "$token" | jq '{expires_on, token_type}'
        echo ""
        log_info "Access token (first 50 chars):"
        echo "$token" | jq -r '.access_token' | cut -c1-50
        echo "..."
    else
        log_warn "No cached token found. Run 'login' first."
    fi
}

api_call() {
    local method="${1:-GET}"
    local endpoint="$2"
    local data="$3"

    local token
    token=$(get_token) || {
        log_error "No valid token. Run 'login' first."
        exit 1
    }

    local url="${API_BASE_URL}${endpoint}"
    local curl_opts=(-s -w "\n%{http_code}" -H "Authorization: Bearer $token" -H "Content-Type: application/json")

    local response
    if [[ "$method" == "POST" ]] && [[ -n "$data" ]]; then
        response=$(curl "${curl_opts[@]}" -X POST -d "$data" "$url")
    elif [[ "$method" == "POST" ]]; then
        response=$(curl "${curl_opts[@]}" -X POST "$url")
    else
        response=$(curl "${curl_opts[@]}" "$url")
    fi

    local http_code=$(echo "$response" | tail -n1)
    local body=$(echo "$response" | sed '$d')

    if [[ "$http_code" -ge 200 ]] && [[ "$http_code" -lt 300 ]]; then
        echo "$body"
    else
        log_error "API request failed (HTTP $http_code)"
        echo "$body" | jq . 2>/dev/null || echo "$body"
        return 1
    fi
}

list_projects() {
    log_info "Fetching projects..."
    api_call GET "/api/projects/" | jq .
}

get_project() {
    local project_uuid="$1"
    if [[ -z "$project_uuid" ]]; then
        log_error "Project UUID required. Use --project <uuid>"
        exit 1
    fi
    log_info "Fetching project $project_uuid..."
    api_call GET "/api/projects/${project_uuid}/" | jq .
}

list_jobs() {
    local project_uuid="$1"
    if [[ -z "$project_uuid" ]]; then
        log_error "Project UUID required. Use --project <uuid>"
        exit 1
    fi
    log_info "Fetching jobs for project $project_uuid..."
    api_call GET "/api/projects/${project_uuid}/jobs/" | jq .
}

get_job() {
    local job_uuid="$1"
    if [[ -z "$job_uuid" ]]; then
        log_error "Job UUID required. Use --job <uuid>"
        exit 1
    fi
    log_info "Fetching job $job_uuid..."
    api_call GET "/api/jobs/${job_uuid}/" | jq .
}

get_job_status() {
    local job_uuid="$1"
    if [[ -z "$job_uuid" ]]; then
        log_error "Job UUID required. Use --job <uuid>"
        exit 1
    fi
    log_info "Fetching status for job $job_uuid..."
    local job=$(api_call GET "/api/jobs/${job_uuid}/")
    echo "$job" | jq '{uuid: .uuid, status: .status, task_name: .task_name, created: .created, modified: .modified}'
}

clone_job() {
    local job_uuid="$1"
    local run_after="$2"

    if [[ -z "$job_uuid" ]]; then
        log_error "Job UUID required. Use --job <uuid>"
        exit 1
    fi

    log_info "Cloning job $job_uuid..."
    local cloned=$(api_call POST "/api/jobs/${job_uuid}/clone/")

    if [[ $? -eq 0 ]]; then
        local new_uuid=$(echo "$cloned" | jq -r '.uuid')
        log_success "Job cloned successfully!"
        echo "$cloned" | jq '{uuid: .uuid, status: .status, task_name: .task_name}'

        if [[ "$run_after" == "true" ]]; then
            log_info "Running cloned job..."
            run_job "$new_uuid"
        fi
    else
        log_error "Failed to clone job"
    fi
}

run_job() {
    local job_uuid="$1"

    if [[ -z "$job_uuid" ]]; then
        log_error "Job UUID required. Use --job <uuid>"
        exit 1
    fi

    log_info "Submitting job $job_uuid for execution..."
    local result=$(api_call POST "/api/jobs/${job_uuid}/run/")

    if [[ $? -eq 0 ]]; then
        log_success "Job submitted!"
        echo "$result" | jq .
    else
        log_error "Failed to run job"
    fi
}

show_logs() {
    log_info "Fetching worker logs..."
    az containerapp logs show \
        --name ccp4i2-bicep-worker \
        --resource-group "${RESOURCE_GROUP:-ccp4i2-bicep-rg-uksouth}" \
        --type console \
        --tail 100 2>/dev/null | jq -r '.Log' 2>/dev/null || cat
}

# Parse arguments
COMMAND=""
PROJECT_UUID=""
JOB_UUID=""
RUN_AFTER="false"

while [[ $# -gt 0 ]]; do
    case "$1" in
        login|token|projects|project|jobs|job|job-status|clone-job|run-job|logs)
            COMMAND="$1"
            shift
            ;;
        --project|-p)
            PROJECT_UUID="$2"
            shift 2
            ;;
        --job|-j)
            JOB_UUID="$2"
            shift 2
            ;;
        --run|-r)
            RUN_AFTER="true"
            shift
            ;;
        --api-url)
            API_BASE_URL="$2"
            shift 2
            ;;
        --help|-h)
            show_help
            exit 0
            ;;
        *)
            log_error "Unknown option: $1"
            show_help
            exit 1
            ;;
    esac
done

# Execute command
case "$COMMAND" in
    login)
        do_login
        ;;
    token)
        show_token
        ;;
    projects)
        list_projects
        ;;
    project)
        get_project "$PROJECT_UUID"
        ;;
    jobs)
        list_jobs "$PROJECT_UUID"
        ;;
    job)
        get_job "$JOB_UUID"
        ;;
    job-status)
        get_job_status "$JOB_UUID"
        ;;
    clone-job)
        clone_job "$JOB_UUID" "$RUN_AFTER"
        ;;
    run-job)
        run_job "$JOB_UUID"
        ;;
    logs)
        show_logs
        ;;
    "")
        show_help
        ;;
    *)
        log_error "Unknown command: $COMMAND"
        show_help
        exit 1
        ;;
esac
