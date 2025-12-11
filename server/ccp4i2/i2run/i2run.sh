#!/usr/bin/env bash
# filepath: /Users/nmemn/Developer/ccp4i2-django/i2run.sh


args=()
dbfile=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --dbFile)
      dbfile="$2"
      shift 2
      ;;
    *)
      args+=("$1")
      shift
      ;;
  esac
done

if [[ -n "$dbfile" ]]; then
  export CCP4I2_DB_FILE="$dbfile"
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Running i2run with arguments: ${args[@]} in directory ${PWD}"

ccp4-python "$SCRIPT_DIR/../../manage.py" migrate

# Replace all --projectName with --project_name
# Ditto --projectPath => --project_path

for i in "${!args[@]}"; do
  if [[ "${args[$i]}" == "--projectName" ]]; then
    args[$i]="--project_name"
  fi
  if [[ "${args[$i]}" == "--projectPath" ]]; then
    args[$i]="--project_path"
  fi
done

exec ccp4-python "$SCRIPT_DIR/../../manage.py" i2run "${args[@]}"
