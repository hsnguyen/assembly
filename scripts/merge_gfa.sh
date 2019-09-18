#!/bin/bash
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <1.gfa> <2.gfa>"
else
	awk -f gfa_merge.awk $1 $2
fi
