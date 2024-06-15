#!/bin/bash

# get the absolute path of this script
declare -r ChemScript_bashrc_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE}")"; pwd)"

export PATH="${PATH}:${ChemScript_bashrc_SCRIPT_DIR}/bash"
export PYTHONPATH="${PYTHONPATH}:${ChemScript_bashrc_SCRIPT_DIR}/pyg16"

# delete variables
unset ChemScript_bashrc_SCRIPT_DIR
