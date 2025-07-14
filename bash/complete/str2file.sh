#!/bin/bash


###################################################################
###################################################################
# abstract
###################################################################
function __complete_str2file(){
	local cur prev cword
	_get_comp_words_by_ref -n : cur prev cword

	if [[ "$cur" =~ ^- ]]; then
		local totaloptionlist="${flagoptionlist[@]} ${valueoptionlist[@]} ${fileoptionlist[@]}"
		if [ "$cword" -eq 1 ]; then
			# add -O option
			totaloptionlist="-O ${totaloptionlist}"
		fi
		COMPREPLY=( $(compgen -W "$totaloptionlist" -- "$cur") )
	else
		if [ "$prev" = -O ] && [ "$cword" -eq 2 ]; then
			COMPREPLY=( $(compgen -f -- "$cur") )
		elif [[ " ${flagoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=()
		elif [[ " ${valueoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=()
		elif [[ " ${fileoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=( $(compgen -f -- "$cur") )
		else
			COMPREPLY=()
		fi
	fi
}
function __complete_str2file_no_option(){
	local flagoptionlist=(
		-h --help
	)
	local valueoptionlist=(
		-O
	)
	local fileoptionlist=()
	__complete_str2file
}



###################################################################
###################################################################
# implementation of str2file (alphabetical order)
###################################################################
complete -F __complete_str2file_no_option name2sdf


###################################################################
complete -F __complete_str2file_no_option pdbid2pdb



###################################################################
function __complete_smiles2png(){
	local flagoptionlist=(
		-h --help
		-p --add-proton
	)
	local valueoptionlist=(
		-O
	)
	local fileoptionlist=()

	__complete_str2file
}

complete -F __complete_smiles2png smiles2png


###################################################################
function __complete_smiles2xyz(){
	local flagoptionlist=(
		-h --help
		-n --num-threads
		-c --num-confs
	)
	local valueoptionlist=(
		-O
	)
	local fileoptionlist=()

	__complete_str2file
}

complete -F __complete_smiles2xyz smiles2xyz




