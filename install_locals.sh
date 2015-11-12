#!/bin/bash

cd ~/work/code/hiveplotter
python setup.py install
cd ~/work/code/connectome/connectome_utils
python setup.py install
cd ~/work/code/bctpy
if ["$(git symbolic-ref --short -q HEAD)" != 'py3']
    then echo "WARNING: bctpy is not on branch py3"
    else python setup.py install
fi
