#!/bin/bash

function _g16log2value(){
	local cur prev cword
	_get_comp_words_by_ref -n : cur prev cword
	local optionlist=(
		-h
		-help
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
	    --total-ZPVE-energy
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

complete -F _g16log2value g16log2value



