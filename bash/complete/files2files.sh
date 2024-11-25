#!/bin/bash


###################################################################
###################################################################
# abstract
###################################################################
function __complete_files2files(){
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
			COMPREPLY=( $(compgen -f -- "$cur") )
		elif [[ " ${valueoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=()
		elif [[ " ${fileoptionlist[@]} " =~ \ ${prev}\  ]]; then
			COMPREPLY=( $(compgen -f -- "$cur") )
		else
			COMPREPLY=( $(compgen -f -- "$cur") )
		fi
	fi
}
function __complete_files2files_no_option(){
	local flagoptionlist=(
		-h --help
	)
	local valueoptionlist=(
		-O
	)
	local fileoptionlist=()
	__complete_files2files
}


###################################################################
###################################################################
# implementation of files2files (alphabetical order)
###################################################################
complete -F __complete_files2files_no_option fchk2json


###################################################################
complete -F __complete_files2files_no_option gjf2xyz


###################################################################
complete -F __complete_files2files_no_option inpcrd2crd


###################################################################
complete -F __complete_files2files_no_option mden2csv


###################################################################
complete -F __complete_files2files_no_option mdout2csv


###################################################################
complete -F __complete_files2files_no_option mol22xyz


###################################################################
complete -F __complete_files2files_no_option mol2xyz


###################################################################
complete -F __complete_files2files_no_option pdb2csv


###################################################################
complete -F __complete_files2files_no_option pdb2xyz


###################################################################
function __complete_xyz2gjf(){
	local flagoptionlist=(
		-h --help
		--resp
	)
	local valueoptionlist=(
		-O
		--link0
		--route
		--title
		--charge
		--multi
		--gen
		--ecp
	)
	local fileoptionlist=()

	__complete_files2files
}
complete -F __complete_xyz2gjf xyz2gjf


###################################################################
function __complete_xyz2pdb(){
	local flagoptionlist=(
		-h --help
		-n --no-header
	)
	local valueoptionlist=(
		-O
		-r --residue
	)
	local fileoptionlist=()
	__complete_files2files
}
complete -F __complete_xyz2pdb xyz2pdb


###################################################################
function __complete_xyz2xyz(){
	local flagoptionlist=(
		-h --help
		--eliminate-header
		--translate-geometric-center
		--translate-mass-center
	)
	local valueoptionlist=(
		-O
	)
	local fileoptionlist=()
	__complete_files2files
}
complete -F __complete_xyz2xyz xyz2xyz

