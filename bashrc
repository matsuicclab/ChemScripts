#!/bin/bash

# get the absolute path of this script
ChemScript_bashrc_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE}")"; pwd)"

# set PATH
export PATH="${PATH}:${ChemScript_bashrc_SCRIPT_DIR}/bash"
export PYTHONPATH="${PYTHONPATH}:${ChemScript_bashrc_SCRIPT_DIR}/python"

# set complete
for f in "${ChemScript_bashrc_SCRIPT_DIR}/bash/complete/"*.sh
do
	source "$f"
done


# delete variables
unset ChemScript_bashrc_SCRIPT_DIR
