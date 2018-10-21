package org.rtassembly.npgraph;

/*
 * Class represent a node together with its vector in relative to another root node.
 * Note: the direction here refers to the one of a ScaffoldVector, not bidirected-graph direction (in/out), neither sequence sense/antisense (+/-)
 * although we can translate into appropriate info
 */
public class NodeVector implements Comparable<NodeVector>{
	
	BidirectedNode node;
	ScaffoldVector vector;
	int score=1; //alignment score + number of occurences?
	
	NodeVector(){}
	NodeVector(BidirectedNode node, ScaffoldVector vector){
		this.node=node;
		this.vector=vector;
	}
	
	public NodeVector(BidirectedNode node, ScaffoldVector vector, int score) {
		this(node, vector);
		this.score=score;
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
		return node.getId() + ":" + vector.toString() + ":" + score;
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

        if (obj instanceof NodeVector){	        
	        final NodeVector other = (NodeVector) obj;
	        BidirectedNode n=other.getNode();
	        if(this.getNode()!=n){
	        	return false;
	        }else{ 
	        	if(SimpleBinner.getUniqueBin(n)!=null) {
	        		int dist=ScaffoldVector.composition(ScaffoldVector.reverse(other.getVector()), this.getVector()).distance(n, n);
	        		boolean retval=(dist < -BidirectedGraph.getKmerSize());
	        		return retval;
	        	}
	        	else {
	        		return this.getVector().consistentWith(other.getVector());
	        	}
	        }     
        }else 
        	return false;
    }
	public boolean qc() {
		return score >= BidirectedGraph.MIN_COVER;
	}

}
