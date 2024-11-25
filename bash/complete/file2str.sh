#!/bin/bash


###################################################################
###################################################################
# abstract
###################################################################
function __complete_file2str(){
	local cur prev cword
	_get_comp_words_by_ref -n : cur prev cword

	if [[ "$cur" =~ ^- ]]; then
		local totaloptionlist="${flagoptionlist[@]} ${valueoptionlist[@]} ${fileoptionlist[@]}"
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
function __complete_file2str_no_option(){
	local flagoptionlist=(
		-h --help
	)
	local valueoptionlist=()
	local fileoptionlist=()
	__complete_file2str
}


###################################################################
###################################################################
# implementation of file2str (alphabetical order)
###################################################################
function __complete_g16log2value(){
	local flagoptionlist=(
		-h --help
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
	    --is-unrestricted
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
	local valueoptionlist=(
		--alpha-orbital-energy
		--beta-orbital-energy
		--orbital-energy
		--geometry
	)
	local fileoptionlist=()

	__complete_file2str
}

complete -F __complete_g16log2value g16log2value



