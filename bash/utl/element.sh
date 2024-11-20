#!/bin/bash


function print_atomicnumber_elementsymbol_table(){
cat <<EOF
1 H
2 He
3 Li
4 Be
5 B
6 C
7 N
8 O
9 F
10 Ne
11 Na
12 Mg
13 Al
14 Si
15 P
16 S
17 Cl
18 Ar
19 K
20 Ca
21 Sc
22 Ti
23 V
24 Cr
25 Mn
26 Fe
27 Co
28 Ni
29 Cu
30 Zn
31 Ga
32 Ge
33 As
34 Se
35 Br
36 Kr
37 Rb
38 Sr
39 Y
40 Zr
41 Nb
42 Mo
43 Tc
44 Ru
45 Rh
46 Pd
47 Ag
48 Cd
49 In
50 Sn
51 Sb
52 Te
53 I
54 Xe
55 Cs
56 Ba
57 La
58 Ce
59 Pr
60 Nd
61 Pm
62 Sm
63 Eu
64 Gd
65 Tb
66 Dy
67 Ho
68 Er
69 Tm
70 Yb
71 Lu
72 Hf
73 Ta
74 W
75 Re
76 Os
77 Ir
78 Pt
79 Au
80 Hg
81 Tl
82 Pb
83 Bi
84 Po
85 At
86 Rn
87 Fr
88 Ra
89 Ac
90 Th
91 Pa
92 U
93 Np
94 Pu
95 Am
96 Cm
97 Bk
98 Cf
99 Es
100 Fm
101 Md
102 No
103 Lr
104 Rf
105 Db
106 Sg
107 Bh
108 Hs
109 Mt
110 Ds
111 Rg
112 Cn
113 Nh
114 Fl
115 Mc
116 Lv
117 Ts
118 Og
EOF
}

function print_elementsymbol_atomicnumber_table(){
	print_atomicnumber_elementsymbol_table |
		awk '{print $2, $1}'
}

function expand_symbol_list(){
	# $1: str
	# H-Ar,Ga-Kr,Fe,I
	# -> H,He,Li,Be,B,C,N,O,F,Ne,Na,Mg,Al,Si,P,S,Cl,Ar,Ga,Ge,As,Se,Br,Kr,Fe,I
	local str patt str2 numlist symblist

	# -> 'H-Ar Ga-Kr Fe I'
	str=`echo "$1" | sed -r 's/,/ /g'`

	# -> '%s-%s %s-%s %s %s'
	patt=`echo "$str" | sed -r 's/[A-z]+/%s/g'`

	# -> '1-18 31-36 26 53'
	str2=`echo "$str" |
			sed -r 's/([- ])/\n\1\n/g' |
			sed -r -n '/[A-z]/p' |
			symb2atomnum |
			tr '\n' ' ' |
			xargs printf "${patt}\n"`

	numlist=`echo "$str2" |
				sed -r 's/ /\n/g' |
				sed -r '/-/!s/^(.+)$/\1-\1/; /^-/s/-/1-/; /-$/s/-/-118/; /-/s/-/ /1' |
				xargs -L 1 seq`
	symblist=`echo "$numlist" |
				atomnum2symb |
				tr '\n' ','`

	echo "$symblist"
}

function append_to_genecp_dict(){
	# $1: gen / ecp
	# $2: 'symblist: basis'
	# results are saved in gendict (or ecpdict)
	# result (e.g.): 'H STO-3G\nHe STO-3G\n...'
	local symblist basis expandedsymblist newdict

	symblist=`echo "$2" | sed -r 's/:.*//' | sed -r 's/^ +//' | sed -r 's/ +$//'`
	basis=`echo "$2" | sed -r 's/^[^:]*://' | sed -r 's/^ +//' | sed -r 's/ +$//'`
	expandedsymblist=`expand_symbol_list "$symblist"`
	newdict=`echo "${expandedsymblist}" | sed -r 's/,/\n/g' | sed -r '/^$/d' | sed -r 's/$/ '"${basis}"'/'`
	if [ "$1" = 'gen' ]; then
		gendict=$(
			echo -e "${gendict-}\n${newdict}" | sed -r '/^$/d' # append to gendict. If gendict is undefined, assign newdict
		)
	elif [ "$1" = 'ecp' ]; then
		ecpdict=$(
			echo -e "${ecpdict-}\n${newdict}" | sed -r '/^$/d'
		)
	else
		error "append_to_genecp_dict: unknown option: '$1'"
	fi
}


