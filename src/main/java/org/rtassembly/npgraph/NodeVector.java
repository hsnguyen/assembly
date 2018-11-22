package org.rtassembly.npgraph;

import java.util.Objects;

/*
 * Class represent a node together with its vector in relative to another root node.
 * Note: the direction here refers to the one of a ScaffoldVector, not bidirected-graph direction (in/out), neither sequence sense/antisense (+/-)
 * although we can translate into appropriate info
 */
public class NodeVector implements Comparable<NodeVector>{
	
	BidirectedNode node;
	ScaffoldVector vector;
	int nodeCover=1; //alignment score + number of occurences?
	
	NodeVector(){}
	NodeVector(BidirectedNode node, ScaffoldVector vector){
		this.node=node;
		this.vector=vector;
	}
	
	public NodeVector(BidirectedNode node, ScaffoldVector vector, int cov) {
		this(node, vector);
		this.nodeCover=cov;
	}
	public BidirectedNode getNode(){return node;}
	public ScaffoldVector getVector(){return vector;}
	
	public void setNode(BidirectedNode node){this.node=node;}
	public void setVector(ScaffoldVector vector){this.vector=vector;}
	
//	public PopBin getUniqueBin(){
//		return SimpleBinner.getUniqueBin(node);
//	}
	
	//get the direction of a node based on the root direction (not apply for the root itself!!!)
	public boolean getDirection(boolean rootDir){
		return rootDir!=(vector.getDirection()>0); //XOR: in-out not +/-
	}
	@Override
	public String toString(){
		return node.getId() + ":" + vector.toString() + ":" + nodeCover;
	}
	
	//compare in term of distance to the root node: for sortedset
	@Override
	public int compareTo(NodeVector o) {
		if(equals(o))
			return 0;
		else if(vector.isIdentity())
			return -1;
		else if(o.vector.isIdentity())
			return 1;
		else
			return Integer.compare(vector.relDistance(node), o.vector.relDistance(o.node));
	}
    @Override
    public int hashCode() {
    	String tmp=node.getId() + Math.signum(vector.direction);
        return tmp.hashCode();
    }

    @Override
    public boolean equals(final Object obj) {
        if (this == obj)
            return true; 

        if (obj == null || obj.getClass() != this.getClass()) { 
        	return false; 
    	}

        final NodeVector other = (NodeVector) obj;
        
        BidirectedNode 	thisNode=this.getNode(),
        				thatNode=other.getNode();
        
        if(!Objects.equals(thisNode, thatNode)){
        	return false;
        }else{ 
        	boolean retval;
        	if(SimpleBinner.getBinIfUnique(thatNode)!=null) {
//        		int dist=ScaffoldVector.composition(ScaffoldVector.reverse(other.getVector()), this.getVector()).distance(thatNode, thisNode);
//        		retval=(dist < -BidirectedGraph.getKmerSize());
        		retval= (Math.abs(this.getVector().relDistance(thisNode) - other.getVector().relDistance(thatNode)) < thisNode.getNumber("len") - BidirectedGraph.getKmerSize());
        	}
        	else {
        		retval=this.getVector().consistentWith(other.getVector());
        	}
    		return retval;

        }     

    }
	public boolean qc() {
		return nodeCover >= BidirectedGraph.MIN_COVER;
	}

	//Merging 2 equal NodeVectors
	public boolean merge(NodeVector nv){
		if(!this.equals(nv))
			return false;

		ScaffoldVector 	thisVector=this.getVector(),
						thatVector=nv.getVector();
		//update the vector
		if(!thisVector.isIdentity()){
			thisVector.setMagnitute( (thisVector.getMagnitute()*this.nodeCover+thatVector.getMagnitute()*nv.nodeCover)/(this.nodeCover+nv.nodeCover));
		}
		
		//update coverage
		this.nodeCover+=nv.nodeCover;
		return true;
	}
}
