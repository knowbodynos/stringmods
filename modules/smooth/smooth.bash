#!/bin/bash

progdir=$(cd $(dirname "${BASH_SOURCE[0]}"); pwd -P)

read line
while [ "${line}" != "" ]
do
    polyid=$(echo ${line} | jq -r ".POLYID")
    geomn=$(echo ${line} | jq -r ".GEOMN")
    triangn=$(echo ${line} | jq -r ".TRIANGN")
    involn=$(echo ${line} | jq -r ".INVOLN")
    dresverts=$(echo ${line} | jq -r ".DRESVERTS" | sed 's/[}{]//g')
    npoints=$(($(echo ${dresverts} | tr ',' '\n' | wc -l)/4))
    triangs=$(echo ${line} | jq -r ".TRIANG" | sed 's/[}{]//g')
    ntriangs=$(($(echo ${triangs} | tr ',' '\n' | wc -l)/4))
    symcypoly=$(echo ${line} | jq -r ".SYMCYPOLY | tostring" | sed 's/"//g' | sed 's/x\([0-9]*\)/x(\1)/g')
    smooth=$(timeout 3600 Singular -q -c "int ndim = 4; int npoints = ${npoints}; int ntriangs = ${ntriangs}; intmat points[npoints][ndim] = ${dresverts}; intmat triangs[ntriangs][ndim] = ${triangs}; ring r=1500450271, (x(1..npoints)), dp; vector pterms = ${symcypoly};" ${progdir}/smooth.sing 2>/dev/null | head -c -1)
    echo "+INVOL.{\"POLYID\":${polyid},\"GEOMN\":${geomn},\"TRIANGN\":${triangn},\"INVOLN\":${involn}}>{\"SMOOTH\":${smooth}}"
    echo "@"
    sleep 0.01
    read line
done