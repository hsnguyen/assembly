package org.rtassembly.npgraph;
/*
 * A presentation of BDEdge with directed node (BDNodeState)
 * TODO: should be merged into one class?
 */
public class BDEdgePrototype{
	BDNodeState n0,n1;

	public BDEdgePrototype(BDNode node, boolean dir){
		n0=new BDNodeState(node, dir);
	}
	public BDEdgePrototype(BDNode node1, BDNode node2, boolean dir1, boolean dir2){
		n0=new BDNodeState(node1, dir1);
		n1=new BDNodeState(node2,dir2);
	}
	public BDEdgePrototype(BDPath path) {
		if(path==null || path.getEdgeCount()>0){
			n0=new BDNodeState(path.getFirstNode(),path.getFirstNodeDirection());
			n1=new BDNodeState(path.getLastNode(),path.getLastNodeDirection());
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
	public String getEdgeID(){
		if(n0==null||n1==null)
			return "";
		else
			return BDEdge.createID(n0.getNode(), n1.getNode(), n0.getDir(), n1.getDir());
	}
}

