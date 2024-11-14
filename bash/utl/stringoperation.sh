#!/bin/bash


function replace_tablecolumn_by_mappingdatabase(){
	# $1: source table string
	# $2: delimiter for source
	# $3: column number
	# $4: map database string
	# $5: delimiter for database
	# $6: default value (if $6 == '--', not replace)
	local sourcedata sourcedelimiter columnnum
	local mapdata mapdelimiter defaultvalue

	sourcedata="$1"
	sourcedelimiter="$2"
	columnnum="$3"
	mapdata="$4"
	mapdelimiter="$5"
	defaultvalue="$6"

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
					if(dv=="--"){
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


