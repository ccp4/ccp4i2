#!/bin/bash
# The app the help pages are captured from: Django on 3421, Next on 3420,
# against a scratch home.
#
#   CCP4I2_HOME=<scratch home> docs/user/tools/devserver.sh start
#   ... stop | restart [django|next|all] | status
#   ... clear-reports <project> <job number>...   (then restart is not needed)
#
# Environment:
#   CCP4I2_HOME   the scratch home (required; a live home is refused)
#   CCP4_SETUP    ccp4.setup-sh for Django (default: ccp4-20260702 beside the repo)
#   SHELXDIR      passed to Django's jobs if set
#
# Django runs with --noreload, so a change to server code needs
# `restart django`; a report cached from before the change needs
# `clear-reports`. Next runs `next dev` and picks up client changes itself.
# Next must NOT see the CCP4 environment: its setup puts Node 20 first on
# the PATH, and the captures need Node 22+.
#
# Logs: $CCP4I2_HOME/django.log and $CCP4I2_HOME/next.log.

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
DJANGO_PORT=${DJANGO_PORT:-3421}
NEXT_PORT=${NEXT_PORT:-3420}
CCP4_SETUP=${CCP4_SETUP:-"$(dirname "$REPO")/ccp4-20260702/bin/ccp4.setup-sh"}

die() { echo "devserver: $*" >&2; exit 1; }

check_home() {
    [ -n "$CCP4I2_HOME" ] || die "set CCP4I2_HOME to a scratch home"
    local home live
    home=$(cd "$CCP4I2_HOME" 2>/dev/null && pwd -P) || die "no such directory: $CCP4I2_HOME"
    for live in .ccp4i2 .ccp4i2-django .ccp4i2x; do
        [ "$home" = "$(cd "$HOME/$live" 2>/dev/null && pwd -P)" ] && die "$CCP4I2_HOME is the live home $live"
    done
}

pids_on() { lsof -ti "tcp:$1" -sTCP:LISTEN 2>/dev/null; }
up() { [ "$(curl -s -o /dev/null -w '%{http_code}' "$1")" = "200" ]; }

wait_for() {  # url, name, seconds
    local i
    for i in $(seq 1 "$3"); do up "$1" && { echo "$2 up"; return 0; }; sleep 2; done
    die "$2 did not come up; see its log in $CCP4I2_HOME"
}

start_django() {
    [ -n "$(pids_on $DJANGO_PORT)" ] && { echo "django already on $DJANGO_PORT"; return; }
    [ -f "$CCP4_SETUP" ] || die "no CCP4 setup at $CCP4_SETUP (set CCP4_SETUP)"
    (
        cd "$REPO/server" || exit 1
        # shellcheck disable=SC1090
        source "$CCP4_SETUP" >/dev/null 2>&1
        env CCP4I2_HOME="$CCP4I2_HOME" DJANGO_SETTINGS_MODULE=ccp4i2.config.settings \
            ${SHELXDIR:+SHELXDIR="$SHELXDIR"} \
            nohup ccp4-python manage.py runserver "127.0.0.1:$DJANGO_PORT" --noreload \
            > "$CCP4I2_HOME/django.log" 2>&1 &
    )
    wait_for "http://127.0.0.1:$DJANGO_PORT/api/ccp4i2/projects/" django 30
}

start_next() {
    [ -n "$(pids_on $NEXT_PORT)" ] && { echo "next already on $NEXT_PORT"; return; }
    (
        cd "$REPO/client/renderer" || exit 1
        env BUILD_TARGET=web NEXT_PUBLIC_API_BASE_URL="http://127.0.0.1:$DJANGO_PORT" \
            nohup ../node_modules/.bin/next dev -p "$NEXT_PORT" \
            > "$CCP4I2_HOME/next.log" 2>&1 &
    )
    wait_for "http://127.0.0.1:$NEXT_PORT/" next 90
}

stop_port() {  # port, name
    local pids
    pids=$(pids_on "$1")
    [ -z "$pids" ] && { echo "$2 not running"; return; }
    # shellcheck disable=SC2086
    kill $pids
    local i
    for i in $(seq 1 15); do [ -z "$(pids_on "$1")" ] && { echo "$2 stopped"; return; }; sleep 1; done
    die "$2 on $1 would not stop"
}

clear_reports() {  # project, job numbers...
    local project=$1; shift
    [ $# -gt 0 ] || die "clear-reports needs a project and job numbers"
    local dir
    dir=$(cd "$REPO/server" && source "$CCP4_SETUP" >/dev/null 2>&1 && \
        env CCP4I2_HOME="$CCP4I2_HOME" DJANGO_SETTINGS_MODULE=ccp4i2.config.settings \
        ccp4-python manage.py shell -c \
        "from ccp4i2.db.models import Project; print('DIR=' + Project.objects.get(name='$project').directory)" \
        2>/dev/null | sed -n 's/^DIR=//p')
    [ -d "$dir" ] || die "no project $project"
    local n f cleared=0
    for n in "$@"; do
        f="$dir/CCP4_JOBS/job_${n//./\/job_}/report_xml.xml"
        [ -f "$f" ] && rm "$f" && cleared=$((cleared + 1))
    done
    echo "cleared $cleared cached report(s) in $project"
}

cmd=${1:-status}; shift
check_home
case "$cmd" in
    start) start_django; start_next ;;
    stop) stop_port $NEXT_PORT next; stop_port $DJANGO_PORT django ;;
    restart)
        case "${1:-django}" in
            django) stop_port $DJANGO_PORT django; start_django ;;
            next) stop_port $NEXT_PORT next; start_next ;;
            all) stop_port $NEXT_PORT next; stop_port $DJANGO_PORT django; start_django; start_next ;;
            *) die "restart django|next|all" ;;
        esac ;;
    status)
        up "http://127.0.0.1:$DJANGO_PORT/api/ccp4i2/projects/" && echo "django up on $DJANGO_PORT" || echo "django down"
        up "http://127.0.0.1:$NEXT_PORT/" && echo "next up on $NEXT_PORT" || echo "next down" ;;
    clear-reports) clear_reports "$@" ;;
    *) die "start | stop | restart [django|next|all] | status | clear-reports <project> <jobs>" ;;
esac
