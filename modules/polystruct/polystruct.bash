#!/bin/bash

prog_dir=$(cd $(dirname "${BASH_SOURCE[0]}"); pwd -P)

read line
while [ "${line}" != "" ]
do
    output=$(echo ${line} | timeout 300 ${prog_dir}/polystruct | head -c -1)
    if [ "${output}" != "" ]
    then
        echo "${output}"
        echo ""
    else
    	echo "None"
    fi
    sleep 0.01
    read line
done