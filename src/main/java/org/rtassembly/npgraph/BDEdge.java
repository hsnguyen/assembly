package org.rtassembly.npgraph;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.graphstream.graph.implementations.AbstractEdge;
import org.graphstream.graph.implementations.AbstractNode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class BDEdge extends AbstractEdge{
	protected boolean dir0, dir1;//true: outward, false: inward
	//note that traversing direction (true: template, false: reverse complement) of destination node is opposite its defined direction (true: outward, false:inward) 
	
    private static final Logger LOG = LoggerFactory.getLogger(BDEdge.class);

	protected BDEdge(String id, AbstractNode node0, AbstractNode node1, boolean dir0, boolean dir1) {
		// id fuck off!!! we'll make one for ourselves
		this(node0,node1,dir0,dir1);
	}
	protected BDEdge(AbstractNode node0, AbstractNode node1, boolean dir0, boolean dir1) {
		//this(createID(src,dst,dir0,dir1), src, dst);

		super(createID(node0,node1,dir0,dir1),node0,node1,false);
		this.dir0=dir0;
		this.dir1=dir1;

	}
	
	
	/*
	 * To adapt Graph class from GraphStream library (called by newInstance())
	 */
	protected BDEdge(String id, AbstractNode source, AbstractNode dest){
		super(id, source, dest, false);
		
		String pattern = "^([0-9]*)([+-]),([0-9]*)([+-])$";
	    // Create a Pattern object
	    Pattern r = Pattern.compile(pattern);
	    // Now create matcher object.
	    Matcher m = r.matcher(id);
	    String 	leftID, rightID, 
	    		leftDir, rightDir;
	    if(m.find()){
	    	leftID=m.group(1);
	    	leftDir=m.group(2);
	    	rightID=m.group(3);
	    	rightDir=m.group(4);
	    	if(source.getId().equals(leftID)){
	    		dir0=(leftDir.equals("+")?true:false);
	    		dir1=(rightDir.equals("-")?true:false);
	    	}else if(source.getId().equals(rightID)){
	    		dir0=(rightDir.equals("-")?true:false);
	    		dir1=(leftDir.equals("+")?true:false);
	    	}else{
	        	System.err.println("ID does not match");
	        	System.exit(1);
	        }
	    } else{
	    	System.err.println("Illegal ID for a bidirected edge (id must have the form node_id[+/-],node_id[+/-])");
	    	System.exit(1);
	    }

	}
	
	//should have smt like this
	public static BDEdge makeEdge (AbstractNode src, AbstractNode dst, BDPath path) {
		return new BDEdge(path.getId(), src, dst);
	}


	
	@Override
	public String toString() {
		return String.format("%s:%s-%s-%s-%s", getId(), getNode0(), (dir0?">":"<"), (dir1?"<":">"), getNode1());
	}
	@Override
	public boolean equals(final Object obj){
        if (this == obj)
            return true; 

        if (obj == null || obj.getClass() != this.getClass()) { 
        	return false; 
    	}

        final BDEdge other = (BDEdge) obj;
        return other.getId().equals(getId());
	}
	
	public void setDir0(boolean dir){
		this.dir0=dir;
	}
	public void setDir1(boolean dir){
		this.dir1=dir;
	}
	public boolean getDir0(){
		return dir0;
	}
	public boolean getDir1(){
		return dir1;
	}

	public int getLength(){
		if(this.hasAttribute("path")) {
			BDPath tmp=(BDPath) this.getAttribute("path");
			return (int)(tmp.getLength()-tmp.getFirstNode().getNumber("len")-tmp.getLastNode().getNumber("len"));
		}else
			return -BDGraph.getKmerSize();
	}

	public Boolean getNodeDirection(AbstractNode node){
		Boolean retval=null;
		if(node==getNode0())
			retval=getDir0();
		if(node==getNode1()){
			if(retval!=null && retval!=getDir1())//infinity self-loop: cannot determine
				retval=null;
			else
				retval= getDir1();
		}
		return retval;
	}
	

	public static String createID(AbstractNode node0, AbstractNode node1, boolean dir0, boolean dir1){
		String 	srcDes = node0.getId(),
				dstDes = node1.getId();
		if(srcDes.compareTo(dstDes) == 0) { // loop
			if(dir0 == dir1)
				return String.format("%s%s,%s%s", srcDes, dir0?"+":"-", dstDes, dir1?"-":"+");
			else
				return String.format("%s%s,%s%s", srcDes, "+", dstDes, "+");
		}
		else if(srcDes.compareTo(dstDes) < 0)
			return String.format("%s%s,%s%s", srcDes, (dir0?"+":"-"), dstDes, (dir1?"-":"+"));
		else
			return String.format("%s%s,%s%s", dstDes, (dir1?"+":"-"), srcDes, (dir0?"-":"+"));
	}
}
