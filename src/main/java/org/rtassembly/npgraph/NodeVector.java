package org.rtassembly.npgraph;

public class NodeVector implements Comparable<NodeVector>{
	
	BidirectedNode node;
	ScaffoldVector vector;
	boolean rootDir; //direction of the marker (out/in==+/-)
	int score; //alignment score + number of occurences
	
	NodeVector(){}
	NodeVector(BidirectedNode node, ScaffoldVector vector, boolean rootDir){
		this.node=node;
		this.vector=vector;
		this.rootDir=rootDir;
	}
	
	public BidirectedNode getNode(){return node;}
	public ScaffoldVector getVector(){return vector;}
	
	public void setNode(BidirectedNode node){this.node=node;}
	public void setVector(ScaffoldVector vector){this.vector=vector;}
	
//	public PopBin getUniqueBin(){
//		return SimpleBinner.getUniqueBin(node);
//	}
	
	//get the direction of a node based on the root direction (not apply for the root itself!!!)
	public boolean getDirection(){
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
		else if(vector.magnitude==0&&vector.direction==1)
			return -1;
		else if(o.vector.magnitude==0&&o.vector.direction==1)
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
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        final NodeVector other = (NodeVector) obj;
        if (this.getNode()!=other.getNode())   
        	return false;
        else{ 
        	return this.getVector().consistentWith(other.getVector());
        }
    }
	public boolean qc() {
		return score >= 2;
	}
}
