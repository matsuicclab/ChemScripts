#!/bin/bash

###################################################################
function __complete_files2files(){
	local cur prev cword
	_get_comp_words_by_ref -n : cur prev cword

	local totaloptionlist="${flagoptionlist[@]} ${valueoptionlist[@]} ${fileoptionlist[@]}"
	if [[ "$cur" =~ ^- ]]; then
		COMPREPLY=( $(compgen -W "$totaloptionlist" -- "$cur") )
	else
		if [[ " ${flagoptionlist[@]} " =~ \ ${prev}\  ]]; then
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

function __complete_csvsmiles2png(){
	local flagoptionlist=(
		-h --help
		-p --add-proton
	)
	local valueoptionlist=(
		-d --delimiter
		-e --header
		-f --field
		-pc --per-column
		-pr --per-row
	)
	local fileoptionlist=()

	__complete_files2files
}
complete -F __complete_csvsmiles2png csvsmiles2png

#TODO How to handle the -o option
#TODO csvsmiles2png is not files2files
#TODO How to handle options with a single hyphen such as -help


###################################################################
function __complete_g16log2value(){
	local cur prev cword
	_get_comp_words_by_ref -n : cur prev cword
	local optionlist=(
		--help
		--num-atom
		--num-elec
		--num-a-elec
		--num-b-elec
		--charge
		--multiplicity
		--stoichiometry
		--basis
		--num-basis
		--num-primitive-basis
		--num-cartesian-basis
	    --theory
	    --solvation-model
	    --solvent
	    --SCF-energy
	    --total-energy
	    --ZPVE-energy
	    --total+ZPVE-energy
	    --internal-energy
	    --gibbs-energy
	    --temperature
	    --pressure
		--total-opt-step
		--geometry
		--opt-geometry
		--num-im-freq
		--scr-files-path
		--version
		--exit-status
		--cpu-time
		--total-cpu-time
		--elapsed-time
		--total-elapsed-time
	)
	optionlist="${optionlist[@]}"
	if [[ "$cur" =~ ^- ]]; then
		COMPREPLY=( $(compgen -W "$optionlist" -- "$cur") )
	elif [ "$prev" = --geometry ]; then
		COMPREPLY=()
	else
		COMPREPLY=( $(compgen -f -- "$cur") )
	fi
}

complete -F __complete_g16log2value g16log2value


###################################################################
function __complete_smiles2png(){
	local cur prev cword
	_get_comp_words_by_ref -n : cur prev cword
	local optionlist=(
		-h --help
		-p --add-proton
	)
	optionlist="${optionlist[@]}"
	if [[ "$cur" =~ ^- ]]; then
		COMPREPLY=( $(compgen -W "$optionlist" -- "$cur") )
	else
		case "$prev" in
		# previous option is the one to specify a flag
		-h | -help | --help)
		-p | -add-proton | --add-proton)
		*)
			COMPREPLY=( $(compgen -f -- "$cur") )
			;;
	fi
}

complete -F __complete_smiles2png smiles2png

###################################################################
function __complete_chemscript_file(){
	local cur prev cword
	_get_comp_words_by_ref -n : cur prev cword

	COMPREPLY=( $(compgen -f -- "$cur") )
}
complete -F __complete_chemscript_file fchk2json
complete -F __complete_chemscript_file inpcrd2crd
complete -F __complete_chemscript_file mden2csv
complete -F __complete_chemscript_file mdout2csv
complete -F __complete_chemscript_file mkresp
complete -F __complete_chemscript_file pdb2csv
complete -F __complete_chemscript_file xyz2pdb
complete -F __complete_chemscript_file rmexcept


###################################################################
function __complete_chemscript_void(){
	COMPREPLY=()
}
complete -F __complete_chemscript_void atomnum2symb
complete -F __complete_chemscript_void symb2atomnum



