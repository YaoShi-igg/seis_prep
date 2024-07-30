#!/usr/bin/sh
rm -rf build dist *.egg-info
python -m build
pip install dist/seis_prep-1.0-py3-none-any.whl
