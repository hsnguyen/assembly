#!/bin/bash
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input.gfa>"
else
	awk '/^S/{print ">"$2; print $3;}' $1 | fold
fi

