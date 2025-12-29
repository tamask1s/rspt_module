#!/bin/bash
python3 -m pip uninstall rspt-module -y && rm -rf build/ dist/ rspt_module.egg-info/ rspt_module/*.so && python3 setup.py install --user
