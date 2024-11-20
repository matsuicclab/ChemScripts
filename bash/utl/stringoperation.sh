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
			awk -v delimiter="$mapdelimiter" \
				-v cn="$columnnum" \
				-v dv="$defaultvalue" \
				'BEGIN  {
					FS="["delimiter"]"; OFS=delimiter; # set delimiter
					print "{"
					printf "switch($%s){\n", cn
				}
				(NF==2){
					# $1 = match pattern
					# $2 = string after replacement
					printf "case \"%s\":\n", $1

					if($2 ~ /^".+"$/){
						# If the string is enclosed in double quotes,
						# avoid enclosing it in double quotes.
						printf "$%s=%s; break\n", cn, $2
					}else{
						# If the string contains double quotes, escape them.
						gsub("\"","\\\"")
						printf "$%s=\"%s\"; break\n", cn, $2
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
					print "}" # end of switch
					print "print $0"
					print "}"
				}'
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
		sed -r "s/$/${sourcedelimiter}/" |
		sed -r "s/^(([^${sourcedelimiter}]*${sourcedelimiter}){${columnnum}})/\1\n/" |
		sed -r "1~2s/([^${sourcedelimiter}]*${sourcedelimiter})$/\n\1/" |
		sed -r "2~3y/${beforechars}/${afterchars}/" |
		sed -r '{N;N;s/\n//g}' |
		sed -r "s/${sourcedelimiter}$//"

}

function toupper_first_char_in_tablecolumn(){
	# $1: source table string
	#     if $1 == '-', read from pipe input
	# $2: delimiter for source
	# $3: column number
	local sourcedata sourcedelimiter columnnum

	sourcedata="$1"
	sourcedelimiter="$2"
	columnnum="$3"

	if [ "$sourcedata" = - ]; then
		sourcedata=`cat -`
	fi

	echo "$sourcedata" |
		awk -v delimiter="${sourcedelimiter}" \
			-v cn="$columnnum" \
			'function toupperfirst(s){
				len = length(s)
				return toupper(substr(s,1,1)) substr(s,2,len-1)
			}
			BEGIN{
				FS="["delimiter"]"; OFS=delimiter; # set delimiter
			}
			{
				$cn = toupperfirst($cn)
				print $0
			}'
}



