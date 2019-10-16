function repeat( str, n, rep, i )
{
    for( ; i<n; i++ )
        rep = rep str   
    return rep
}
BEGIN{
	CID=0;
	if(cid)#contig id
		CID=cid;
	BID=0;
	if(bid)#bin id
		BID=bid;
	FID=0;
	if(fid)#file id
		FID=fid;
} 
FNR==NR{ #Processing GFA file
	if(/^S/){#Segments must come first
		if(cidmap[$2]=="")
			cidmap[$2]=++CID;
		$2=cidmap[$2];
			
	}else if(/^[LC]/){#Links and Containment
		$2=cidmap[$2];
		$4=cidmap[$4];
	}else if(/^P/){#Paths
		#gsub(/[0-9]+/,ids["\&"],$3); didn't work: parameter must be direct value
		l=$3#need to do this to prevent re-match of the replaced id
		while(match(l, "[0-9]+")) {
    			fnd=substr(l, RSTART, RLENGTH);
    			$3=substr($3, 1, RSTART-1)""cidmap[fnd]""substr($3,RSTART+RLENGTH)
			l=substr(l, 1, RSTART-1)""repeat("x",length(cidmap[fnd]))""substr(l,RSTART+RLENGTH)
		}
	}else
		next;
	

	print >> "reindexed_"FID".gfa";
	
}
FNR!=NR{ #Processing bin file
	$1=cidmap[$1];
	if($2>0){
		if(bidmap[$2]="")
			bidmap[$2]=++BID;
		$2=bidmap[$2];
	}	
	print >> "reindexed_"FID".bin"
}
