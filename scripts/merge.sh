#!/bin/bash
n="$#"
if [ $((n-1)) -lt 1 ] || [ $((n%2)) -eq 1 ]
then
    echo "Usage: $0 <1.gfa> <1.bin> <2.gfa> <2.bin>..."
else
	SDIR=`dirname $0`
	ncount=0 ##contigs count
	bcount=0 ##bin count
	fcount=1 ##file count
	while [ $((2*fcount-n)) -ge 0 ] 
	do
		awk -v cid=${ncount},fid=${fcount},bid=${bcount} -f ${SDIR}/reindex.awk $((2*fcount-1)) $((2*fcount)) 
		sline=`cat reindexed_${fcount}.gfa|awk '/^S/{count++}END{print count}'`
		ncount=$((ncount + sline))
		
		bmax=`cat reindexed_${fcount}.bin|awk '{if(max<$2){max=$2}}END{print max}'`
		bcount=$((bmax+1))
		fcount=$((fcount + 1))
	done	
	#concatenate all reindexed files together
	awk '{if(data[$1]=="") data[$1]=$$0; else data[$1]=data[$1]"\n"$0;}END{if(data["S"]!="") print data["S"]; if(data["L"]!="") print data["L"]; if(data["C"]!="") print data["C"]; if(data["P"]!="") print data["P"];}' reindexed_[1-${fcount}].gfa > combine.gfa
	cat reindexed_[1-${fcount}].bin > combine.bin
fi
