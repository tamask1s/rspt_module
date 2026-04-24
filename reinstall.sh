#!/usr/bin/env bash
set -euo pipefail

# --- Colors ---
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# --- Icons ---
OK="✔"
FAIL="✖"
ARROW="➜"

LOG="/tmp/rspt_install.log"

# --- Spinner ---
spinner() {
    local pid=$!
    local delay=0.1
    local spinstr='|/-\'

    while kill -0 $pid 2>/dev/null; do
        local temp=${spinstr#?}
        printf " [%c]  " "$spinstr"
        spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "\b\b\b\b\b\b"
    done
    printf "      \b\b\b\b\b\b"
}

run_step() {
    local msg="$1"
    shift

    printf "${BLUE}${ARROW}${NC} %s..." "$msg"

    ("$@" >"$LOG" 2>&1) &
    spinner

    if wait $!; then
        printf "\r${GREEN}${OK}${NC} %s\n" "$msg"
    else
        printf "\r${RED}${FAIL}${NC} %s\n" "$msg"
        printf "${YELLOW}Log:${NC} $LOG\n"
        exit 1
    fi
}

echo ""
printf "${BLUE}=== RSPT MODULE INSTALL ===${NC}\n\n"

run_step "Uninstall old version" \
    python3 -m pip uninstall rspt-module -y

run_step "Clean build artifacts" \
    rm -rf build/ dist/ rspt_module.egg-info/ rspt_module/*.so

run_step "Install module" \
    python3 setup.py install --user

echo ""
printf "${GREEN}${OK} ALL DONE SUCCESSFULLY${NC}\n\n"