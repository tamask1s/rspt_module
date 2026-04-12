#!/usr/bin/env bash
set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color

LOG="/tmp/rspt_install.log"

python3 -m pip uninstall rspt-module -y >/dev/null 2>&1
rm -rf build/ dist/ rspt_module.egg-info/ rspt_module/*.so

if python3 setup.py install --user >"$LOG" 2>&1; then
    echo -e "${GREEN}*"
    echo -e "RSPT Module installation SUCCESSFUL"
    echo -e "* ${NC}"
else
    echo -e "${RED}*"
    echo -e "RSPT Module installation FAIL (log: $LOG)"
    echo -e "* ${NC}"
    exit 1
fi