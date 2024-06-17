#!/bin/bash

###################################################################
function __complete_file2seqfiles(){
	__complete_files2files
}
function __complete_file2seqfiles_no_option(){
	local flagoptionlist=(
		-h --help
	)
	local valueoptionlist=()
	local fileoptionlist=()
	__complete_file2seqfiles
}

#################################
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

#################################
function __complete_files2files(){
	local cur prev cword
	_get_comp_words_by_ref -n : cur prev cword

	if [[ "$cur" =~ ^- ]]; then
		local totaloptionlist="${flagoptionlist[@]} ${valueoptionlist[@]} ${fileoptionlist[@]}"
		if [ "$cword" -eq 1 ]; then
			# add -o option
			totaloptionlist="-o ${totaloptionlist}"
		fi
		COMPREPLY=( $(compgen -W "$totaloptionlist" -- "$cur") )
	else
		if [ "$prev" = -o ] && [ "$cword" -eq 2 ]; then
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
	local valueoptionlist=()
	local fileoptionlist=()
	__complete_files2files
}

#################################
function __complete_str2file(){
	local cur prev cword
	_get_comp_words_by_ref -n : cur prev cword

	if [[ "$cur" =~ ^- ]]; then
		local totaloptionlist="${flagoptionlist[@]} ${valueoptionlist[@]} ${fileoptionlist[@]}"
		if [ "$cword" -eq 1 ]; then
			# add -o option
			totaloptionlist="-o ${totaloptionlist}"
		fi
		COMPREPLY=( $(compgen -W "$totaloptionlist" -- "$cur") )
	else
		if [ "$prev" = -o ] && [ "$cword" -eq 2 ]; then
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
	local valueoptionlist=()
	local fileoptionlist=()
	__complete_str2file
}

#################################
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
complete -F __complete_str2str_no_option atomnum2symb


###################################################################
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

	__complete_file2seqfiles
}
complete -F __complete_csvsmiles2png csvsmiles2png


###################################################################
complete -F __complete_files2files_no_option fchk2json


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
		--geometry
	)
	local fileoptionlist=()

	__complete_file2str
}

complete -F __complete_g16log2value g16log2value


###################################################################
complete -F __complete_files2files_no_option inpcrd2crd


###################################################################
complete -F __complete_files2files_no_option mden2csv


###################################################################
complete -F __complete_files2files_no_option mdout2csv


###################################################################
# no difine of mkresp


###################################################################
complete -F __complete_files2files_no_option pdb2csv


###################################################################
# no difine of plotmden


###################################################################
# no difine of rmexcept


###################################################################
function __complete_smiles2png(){
	local flagoptionlist=(
		-h --help
		-p --add-proton
	)
	local valueoptionlist=()
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
	local valueoptionlist=()
	local fileoptionlist=()

	__complete_str2file
}

complete -F __complete_smiles2xyz smiles2xyz


###################################################################
complete -F __complete_str2str_no_option symb2atomnum


###################################################################
complete -F __complete_files2files_no_option xyz2pdb



