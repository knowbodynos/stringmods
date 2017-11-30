#!/bin/bash

read line
while [ "${line}" != "" ]
do
    nf=$(echo ${line} | jq -r ".NFORM")
    nf2skel=$(echo ${line} | jq -r ".NFORM2SKEL")
    nfinetriangs=$(echo ${nf2skel} | sed 's/{/[/g' | sed 's/}/]/g' | timeout 3600 points2nfinetriangs 2>/dev/null | head -c -1)
    echo "+FACET.{\"NFORM\":\"${nf}\"}>{\"FACETNREGTRIANG\":${nfinetriangs}}"
    echo "@"
    sleep 0.01
    read line
done