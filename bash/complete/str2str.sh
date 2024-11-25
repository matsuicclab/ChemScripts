#!/bin/bash


###################################################################
###################################################################
# abstract
###################################################################
function __complete_str2str(){
	local cur prev cword
	_get_comp_words_by_ref -n : cur prev cword

	if [[ "$cur" =~ ^- ]]; then
		local totaloptionlist="${flagoptionlist[@]} ${valueoptionlist[@]} ${fileoptionlist[@]}"
		COMPREPLY=( $(compgen -W "$totaloptionlist" -- "$cur") )
	else
		if [[ " ${flagoptionlist[@]} " =~ \ ${prev}\  ]]; then
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
function __complete_str2str_no_option(){
	local flagoptionlist=(
		-h --help
	)
	local valueoptionlist=()
	local fileoptionlist=()
	__complete_str2str
}


###################################################################
###################################################################
# implementation of str2str (alphabetical order)
###################################################################
function __complete_atomnum2symb(){
	local flagoptionlist=(
		-h --help
	)
	local valueoptionlist=(
		-d --delimiter
		-f --field
		-s --skip
	)
	local fileoptionlist=()

	__complete_str2str
}
complete -F __complete_atomnum2symb atomnum2symb


###################################################################
function __complete_chemblid2smiles(){
	local flagoptionlist=(
		-h --help
		-f --force
	)
	local valueoptionlist=(
		-d --default
	)
	local fileoptionlist=()

	__complete_str2str
}
complete -F __complete_chemblid2smiles chemblid2smiles


###################################################################
function __complete_energy2energy(){
	local flagoptionlist=(
		-h --help
		-w --whitespace
	)
	local valueoptionlist=(
		-d --delimiter
		-f --field
		-s --skip
		--format
		--from
		--to
	)
	local fileoptionlist=()

	__complete_str2str
}
complete -F __complete_energy2energy energy2energy


###################################################################
complete -F __complete_energy2energy len2len


###################################################################
complete -F __complete_chemblid2smiles name2smiles


###################################################################
complete -F __complete_chemblid2smiles pubchemid2smiles


###################################################################
complete -F __complete_atomnum2symb symb2atomnum


###################################################################
function __complete_symb2numbasis(){
	local flagoptionlist=(
		-h --help
		-t --total
	)
	local valueoptionlist=(
		-b --basis
		-g --gen
		--dtype
		--ftype
		--num-of
		-d --delimiter
		-f --field
		-s --skip
	)
	local fileoptionlist=()

	__complete_str2str
}
complete -F __complete_symb2numbasis symb2numbasis



