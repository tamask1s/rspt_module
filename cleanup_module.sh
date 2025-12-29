#!/bin/bash
# Clean up corrupted pip package installations
rm -rf /home/tamas/.local/lib/python3.6/site-packages/*spt* /home/tamas/.local/lib/python3.6/site-packages/~* /home/tamas/.local/lib/python3.6/site-packages/-*
echo "Cleaned up corrupted package installations"
