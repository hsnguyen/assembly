function repeat( str, n, rep, i )
{
    for( ; i<n; i++ )
        rep = rep str   
    return rep
}
BEGIN{
	maxId=0;
} 
NR==FNR{
	if(/^S/){#Segments must come first
		if(ids1[$2]=="")
			ids1[$2]=++maxId;
		$2=ids1[$2];
			
	}else if(/^[LC]/){#Links and Containment
		$2=ids1[$2];
		$4=ids1[$4];
	}else if(/^P/){#Paths
		#gsub(/[0-9]+/,ids["\&"],$3); didn't work: parameter must be direct value
		l=$3#need to do this to prevent re-match of the replaced id
		while(match(l, "[0-9]+")) {
    			fnd=substr(l, RSTART, RLENGTH);
    			$3=substr($3, 1, RSTART-1)""ids1[fnd]""substr($3,RSTART+RLENGTH)
			l=substr(l, 1, RSTART-1)""repeat("x",length(ids1[fnd]))""substr(l,RSTART+RLENGTH)
		}
	}else
		next;
	

	if(data[$1]!="")
		data[$1]=data[$1]"\n"$0;
	else
		data[$1]=$0;

	
}
NR!=FNR{
	if(/^S/){#Segments must come first
		if(ids2[$2]=="")
			ids2[$2]=++maxId;
		$2=ids2[$2];
			
	}else if(/^[LC]/){#Links and Containment
		$2=ids2[$2];
		$4=ids2[$4];
	}else if(/^P/){#Paths
		#gsub(/[0-9]+/,ids["\&"],$3);
		l=$3#need to do this to prevent re-match of the replaced id
		while(match(l, "[0-9]+")) {
    			fnd=substr(l, RSTART, RLENGTH);
    			$3=substr($3, 1, RSTART-1)""ids2[fnd]""substr($3,RSTART+RLENGTH)
			l=substr(l, 1, RSTART-1)""repeat("x",length(ids2[fnd]))""substr(l,RSTART+RLENGTH)
		}
	}else
		next;
	
	if(data[$1]!="")
		data[$1]=data[$1]"\n"$0;
	else
		data[$1]=$0;
}
END{
	if("S" in data) 
		print data["S"];
	if("L" in data) 
		print data["L"];
	if("C" in data) 
		print data["C"];
	if("P" in data) 
		print data["P"];
}
