#!/bin/bash

read line
while [ "${line}" != "" ]
do
    nf=$(echo ${line} | jq -r ".NFORM")
    nf2skel=$(echo ${line} | jq -r ".NFORM2SKEL")
    nallfinetriangs=$(echo ${nf2skel} | sed 's/{/[/g' | sed 's/}/]/g' | points2nallfinetriangs 2>/dev/null | head -c -1)
    echo "+FACET.{\"NFORM\":\"${nf}\"}>{\"FACETNTRIANG\":${nallfinetriangs}}"
    echo "@"
    sleep 0.01
    read line
done