package org.rtassembly.npgraph;

import org.graphstream.graph.implementations.AbstractNode;

public class BidirectedEdgePrototype{
	NodeDirection n0,n1;

	public BidirectedEdgePrototype(AbstractNode node, boolean dir){
		n0=new NodeDirection(node, dir);
	}
	public BidirectedEdgePrototype(AbstractNode node1, AbstractNode node2, boolean dir1, boolean dir2){
		n0=new NodeDirection(node1, dir1);
		n1=new NodeDirection(node2,dir2);
	}
	public BidirectedEdgePrototype(BidirectedPath path) throws Exception {
		if(path==null || path.size()<=1)
			throw new Exception("Invalid path to make edge prototype!");
		else{
			BidirectedNode 	start=(BidirectedNode) path.getRoot(),
							end=(BidirectedNode) path.peekNode();
			boolean startDir=((BidirectedEdge) path.getEdgePath().get(0)).getDir(start),
					endDir=((BidirectedEdge) path.peekEdge()).getDir(end);
			n0=new NodeDirection(start,startDir);
			n1=new NodeDirection(end,endDir);
			}
	}
	public BidirectedEdgePrototype(Alignment from, Alignment to){
		if(from!=null)
			n0=new NodeDirection(from.node, from.strand);
		if(to!=null)
			n1=new NodeDirection(to.node, !to.strand);
		
	}
	
	public BidirectedEdgePrototype reverse() {
		return new BidirectedEdgePrototype(n1.getNode(), n0.getNode(), n1.getDir(), n0.getDir());
	}
	
	public AbstractNode getNode0(){
		return n0==null?null:n0.getNode();
	}
	public AbstractNode getNode1(){
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

class NodeDirection implements Comparable<NodeDirection>{
	AbstractNode node;
	boolean dir;
	int weight;//for shortest path finding
	
	NodeDirection(AbstractNode node, boolean dir){
		this.node=node;
		this.dir=dir;
	}
	NodeDirection(AbstractNode node, boolean dir, int weight){
		this.node=node;
		this.dir=dir;
		this.weight=weight;
	}
	public AbstractNode getNode(){return node;}
	public boolean getDir(){
		assert node!=null: "Null node doesn't have direction!";
		return dir;
	}
	public int getWeight() {return  weight;}
	public void setWeight(int w) {weight=w;}
	
	public String toString(){
		return node.getId()+ (dir?"o":"i");
	}

	@Override
	public int compareTo(NodeDirection o) {
		return Integer.compare(weight, o.weight);
	}
    @Override
    public int hashCode() {
        return toString().hashCode();
    }

    @Override
    public boolean equals(final Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        final NodeDirection other = (NodeDirection) obj;
        if (this.toString()==other.toString())   
        	return true;
        else 
        	return false;
    }
}