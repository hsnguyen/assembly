#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input.fastg>"
else
	awk -F '[:;]' -v q=\' 'BEGIN{flag=0;}/^>/{if(index($1,q)!=0) flag=0; else flag=1;}{if(flag==1) print $1;}' assembly_graph.fastg
fi

