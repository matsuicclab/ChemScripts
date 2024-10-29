#!/bin/bash

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

#################################
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
###################################################################
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
function __complete_energy2energy(){
	local flagoptionlist=(
		-h --help
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


###################################################################
complete -F __complete_files2files_no_option inpcrd2crd


###################################################################
complete -F __complete_energy2energy len2len


###################################################################
complete -F __complete_files2files_no_option mden2csv


###################################################################
complete -F __complete_files2files_no_option mdout2csv


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
complete -F __complete_files2files_no_option pdb2csv


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


###################################################################
complete -F __complete_atomnum2symb symb2atomnum


###################################################################
function __complete_xyz2gjf(){
	local flagoptionlist=(
		-h --help
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
complete -F __complete_files2files_no_option xyz2pdb



