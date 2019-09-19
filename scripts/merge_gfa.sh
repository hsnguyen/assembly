#!/bin/bash
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 <1.gfa> <2.gfa> ..."
else
	DIR=`dirname $0`
	ncount=0
	fcount=0
	for gfa in "$@"
	do
		fcount=$((fcount + 1))
		awk -v id="$ncount" -f ${DIR}/gfa_reindex.awk ${gfa} > reindexed_${fcount}.gfa
		sline=`cat reindexed_${gfa}|awk '/^S/{count++}END{print count}'`
		ncount=$((ncount + sline))
	done	
	#concatenate all reindexed gfa files together
	awk '{if(data[$1]=="") data[$1]=$$0; else data[$1]=data[$1]"\n"$0;}END{if(data["S"]!="") print data["S"]; if(data["L"]!="") print data["L"]; if(data["C"]!="") print data["C"]; if(data["P"]!="") print data["P"];}' reindexed_[1-${fcount}].gfa
fi