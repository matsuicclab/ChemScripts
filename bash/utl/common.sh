#!/bin/bash

function error(){
	help >&2
	while [ "$#" -ge 1 ]; do echo -e "\e[31m$1\e[m" >&2 ; shift; done
	exit 1
}

function warning(){
	while [ "$#" -ge 1 ]; do echo -e "\e[31m$1\e[m" >&2 ; shift; done
}



