#!/bin/bash


function replace_tablecolumn_by_mappingdatabase(){
	# $1: source table string
	#     if $1 == '-', read from pipe input
	# $2: delimiter for source
	# $3: column number
	# $4: map database string
	# $5: delimiter for database
	# $6: default value
	#     if $6 == '-', do not replace
	#     if you want to set the default to '-', you can do it by escaping it.
	local sourcedata sourcedelimiter columnnum
	local mapdata mapdelimiter defaultvalue

	sourcedata="$1"
	sourcedelimiter="$2"
	columnnum="$3"
	mapdata="$4"
	mapdelimiter="$5"
	defaultvalue="$6"
	
	if [ "$sourcedata" = - ]; then
		sourcedata=`cat -`
	fi

	awkprocess=$(
		echo "$mapdata" |
			sed -r "s/${mapdelimiter}/\n/1" |
			awk -v cn="$columnnum" \
				-v dv="$defaultvalue" \
				'BEGIN  {
					printf "switch($%s){\n", cn
				}
				{
					if(NR%2==1){
						# $0 = match pattern
						printf "case \"%s\":\n", $0
					} else {
						# $0 = string after replacement
						if($0 ~ /^".+"$/){
							# If the string is enclosed in double quotes,
							# avoid enclosing it in double quotes.
							printf "$%s=%s; break\n", cn, $0
						}else{
							# If the string contains double quotes, escape them.
							gsub("\"","\\\"")
							printf "$%s=\"%s\"; break\n", cn, $0
						}
					}
				}
				END {
					# set the default value
					print "default:"
					if(dv=="-"){
						print "break"
					}else{
						printf "$%s=\"%s\"\n", cn, dv
					}
					print "}"
				}' |
			sed -r '1i {' |
			sed -r '$a print $0\n}'
	)

	echo "$sourcedata" |
		awk -v delimiter="${sourcedelimiter}" \
			'BEGIN {
				FS="["delimiter"]"; OFS=delimiter; # set delimiter
			}'"$awkprocess" # perform replacement
}

function replace_char_in_tablecolumn_with_char(){
	# $1: source table string
	#     if $1 == '-', read from pipe input
	# $2: delimiter for source
	# $3: column number
	# $4: chars before conversion
	# $5: chars after conversion
	local sourcedata sourcedelimiter columnnum
	local beforechars afterchars

	sourcedata="$1"
	sourcedelimiter="$2"
	columnnum="$3"
	beforechars="$4"
	afterchars="$5"

	if [ "$sourcedata" = - ]; then
		sourcedata=`cat -`
	fi

	if [ ! "${#beforechars}" -eq "${#afterchars}" ]; then
		error "number of characters before and after replacement does not match"
	fi
	
	echo "$sourcedata" |
		sed -r "s/^(([^${sourcedelimiter}]*${sourcedelimiter}){${columnnum}})/\1\n/" |
		sed -r "1~2s/([^${sourcedelimiter}]*${sourcedelimiter})$/\n\1/" |
		sed -r "2~3y/${beforechars}/${afterchars}/" |
		sed -r '{N;N;s/\n//g}'
	
}

