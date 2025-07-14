#!/bin/bash

###################################################################
###################################################################
# abstract
###################################################################
function __complete_file2seqfiles(){
	__complete_files2files
}
function __complete_file2seqfiles_no_option(){
	local flagoptionlist=(
		-h --help
	)
	local valueoptionlist=(
		-O
	)
	local fileoptionlist=()
	__complete_file2seqfiles
}


###################################################################
###################################################################
# implementation of abstract (alphabetical order)
###################################################################
function __complete_csvsmiles2png(){
	local flagoptionlist=(
		-h --help
		-p --add-proton
	)
	local valueoptionlist=(
		-O
		-d --delimiter
		-e --header
		-f --field
		-pc --per-column
		-pr --per-row
	)
	local fileoptionlist=()

	__complete_file2seqfiles
}
complete -F __complete_csvsmiles2png csvsmiles2png


###################################################################
###################################################################
# misc
###################################################################
function __complete_mkresp(){
	local cur prev cword
	_get_comp_words_by_ref -n : cur prev cword

	local flagoptionlist=(
		-h --help
	)
	local valueoptionlist=()
	local fileoptionlist=()

	if [[ "$cur" =~ ^- ]]; then
		local totaloptionlist="${flagoptionlist[@]} ${valueoptionlist[@]} ${fileoptionlist[@]}"
		COMPREPLY=( $(compgen -W "$totaloptionlist" -- "$cur") )
	else
		if [[ " ${flagoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=()
		elif [[ " ${valueoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=()
		elif [[ " ${fileoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=()
		else
			COMPREPLY=( $(compgen -f -- "$cur") )
		fi
	fi
}
complete -F __complete_mkresp mkresp


###################################################################
function __complete_plotmden(){
	local cur prev cword
	_get_comp_words_by_ref -n : cur prev cword

	local flagoptionlist=(
		-h
	)
	local valueoptionlist=(
		-d
	)
	local fileoptionlist=()

	if [[ "$cur" =~ ^- ]]; then
		local totaloptionlist="${flagoptionlist[@]} ${valueoptionlist[@]} ${fileoptionlist[@]}"
		COMPREPLY=( $(compgen -W "$totaloptionlist" -- "$cur") )
	else
		if [[ " ${flagoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=()
		elif [[ " ${valueoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=()
		elif [[ " ${fileoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=()
		else
			COMPREPLY=( $(compgen -f -- "$cur") )
		fi
	fi
}
complete -F __complete_plotmden plotmden


###################################################################
function __complete_rmexcept(){
	local cur prev cword
	_get_comp_words_by_ref -n : cur prev cword

	local flagoptionlist=(
		-h --help
		-r --recursive
	)
	local valueoptionlist=()
	local fileoptionlist=()

	if [[ "$cur" =~ ^- ]]; then
		local totaloptionlist="${flagoptionlist[@]} ${valueoptionlist[@]} ${fileoptionlist[@]}"
		COMPREPLY=( $(compgen -W "$totaloptionlist" -- "$cur") )
	else
		if [[ " ${flagoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=( $(compgen -f -- "$cur") )
		elif [[ " ${valueoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=()
		elif [[ " ${fileoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=()
		else
			COMPREPLY=( $(compgen -f -- "$cur") )
		fi
	fi
}
complete -F __complete_rmexcept rmexcept






