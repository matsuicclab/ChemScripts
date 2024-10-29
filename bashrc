#!/bin/bash

# get the absolute path of this script
ChemScript_bashrc_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE}")"; pwd)"

# set PATH
export PATH="${PATH}:${ChemScript_bashrc_SCRIPT_DIR}/bash"
export PYTHONPATH="${PYTHONPATH}:${ChemScript_bashrc_SCRIPT_DIR}/python"

# set complete
source "${ChemScript_bashrc_SCRIPT_DIR}/bash/utl/complete.sh"

# delete variables
unset ChemScript_bashrc_SCRIPT_DIR
