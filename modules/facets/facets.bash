#!/bin/bash

read line
while [ "${line}" != "" ]
do
    polyid=$(echo ${line} | jq -r ".POLYID")
    nverts=$(echo ${line} | jq -r ".NVERTS")
    timeout 3600 ./facets "${polyid}" "${nverts}" 2>/dev/null
    echo "@"
    sleep 0.01
    read line
done