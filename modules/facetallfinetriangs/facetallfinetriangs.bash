#!/bin/bash

read line
while [ "${line}" != "" ]
do
    nf=$(echo ${line} | jq -r ".NFORM")
    nf2skel=$(echo ${line} | jq -r ".NFORM2SKEL")
    nallfinetriangs=$(echo ${nf2skel} | sed 's/{/[/g' | sed 's/}/]/g' | points2nallfinetriangs 2>/dev/null | head -c -1)
    echo "+FACET.{\"NFORM\":\"${nf}\"}>{\"FACETNTRIANG\":${nallfinetriangs}}"
    #allfinetriangs=$(echo "{$(echo ${nf2skel} | sed 's/{/[/g' | sed 's/}/]/g' | timeout 3600 points2allfinetriangs 2>/dev/null | rev | cut -d':' -f1 | rev | sed 's/];//g' | tr '\n' ',' | head -c -1)}")
    #nallfinetriangs=$(($(echo ${allfinetriangs} | grep -o "}},{{" | wc -l)+1))
    #echo "+FACET.{\"NFORM\":\"${nf}\"}>{\"FACETNTRIANG\":${nallfinetriangs},\"FACETTRIANG\":\"${allfinetriangs}\"}"
    echo "@"
    sleep 0.01
    read line
done