package org.rtassembly.npgraph;

public class BDEdgePrototype{
	BDNodeState n0,n1;

	public BDEdgePrototype(BDNode node, boolean dir){
		n0=new BDNodeState(node, dir);
	}
	public BDEdgePrototype(BDNode node1, BDNode node2, boolean dir1, boolean dir2){
		n0=new BDNodeState(node1, dir1);
		n1=new BDNodeState(node2,dir2);
	}
	public BDEdgePrototype(BDPath path) throws Exception {
		if(path==null || path.size()<=1)
			throw new Exception("Invalid path to make edge prototype!");
		else{
			BDNode 	start=(BDNode) path.getRoot(),
							end=(BDNode) path.peekNode();
			boolean startDir=((BDEdge) path.getEdgePath().get(0)).getDir(start),
					endDir=((BDEdge) path.peekEdge()).getDir(end);
			n0=new BDNodeState(start,startDir);
			n1=new BDNodeState(end,endDir);
			}
	}
	public BDEdgePrototype(Alignment from, Alignment to){
		if(from!=null)
			n0=new BDNodeState(from.node, from.strand);
		if(to!=null)
			n1=new BDNodeState(to.node, !to.strand);
		
	}
	
	public BDEdgePrototype reverse() {
		return new BDEdgePrototype(n1.getNode(), n0.getNode(), n1.getDir(), n0.getDir());
	}
	
	public BDNode getNode0(){
		return n0==null?null:n0.getNode();
	}
	public BDNode getNode1(){
		return n1==null?null:n1.getNode();
	}

	public boolean getDir0(){
		return n0.getDir();
	}
	public boolean getDir1(){
		return n1.getDir();
	}
	
	public int getEndingsNum(){
		int retval=0;
		if(n0!=null) retval++;
		if(n1!=null) retval++;
		
		return retval;
	}
	public String toString(){
		String retval=n0.toString();
		if(n1!=null)
			retval+=","+n1.toString();
		else
			retval+=",-";
		return retval;
	}
}

