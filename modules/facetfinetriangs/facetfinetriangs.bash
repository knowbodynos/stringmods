#!/bin/bash

read line
while [ "${line}" != "" ]
do
    nf=$(echo ${line} | jq -r ".NFORM")
    nf2skel=$(echo ${line} | jq -r ".NFORM2SKEL")
    nfinetriangs=$(echo ${nf2skel} | sed 's/{/[/g' | sed 's/}/]/g' | timeout 60 points2nfinetriangs 2>/dev/null | head -c -1)
    if [ "${nfinetriangs}" != "" ]
    then
        echo "set FACET {\"NFORM\":\"${nf}\"} {\"FACETNREGTRIANG\":${nfinetriangs}}"
        #finetriangs=$(echo "{$(echo ${nf2skel} | sed 's/{/[/g' | sed 's/}/]/g' | timeout 3600 points2finetriangs --regular 2>/dev/null | rev | cut -d':' -f1 | rev | sed 's/];//g' | tr '\n' ',' | head -c -1)}")
        #nfinetriangs=$(($(echo ${finetriangs} | grep -o "}},{{" | wc -l)+1))
        #echo "+FACET.{\"NFORM\":\"${nf}\"}>{\"FACETNREGTRIANG\":${nfinetriangs},\"FACETREGTRIANG\":\"${finetriangs}\"}"
    fi
    echo ""
    sleep 0.01
    read line
done