package org.rtassembly.npgraph;

import java.util.Objects;

/*
 * Class represent a node together with its vector in relative to another root node.
 * Note: the direction here refers to the one of a ScaffoldVector, not bidirected-graph direction (in/out), neither sequence sense/antisense (+/-)
 * although we can translate into appropriate info
 */
public class BDNodeVecState implements Comparable<BDNodeVecState>{
	BDNode node;
	ScaffoldVector vector;
	double nvsScore=0.0; //alignment score + number of occurences?
	
	BDNodeVecState(BDNode node, ScaffoldVector vector){
		this.node=node;
		this.vector=vector;
	}
	
	public BDNodeVecState(BDNode node, ScaffoldVector vector, double score) {
		this(node, vector);
		this.nvsScore=score;
	}
	public BDNode getNode(){return node;}
	public ScaffoldVector getVector(){return vector;}
	
	public void setNode(BDNode node){this.node=node;}
	public void setVector(ScaffoldVector vector){this.vector=vector;}
	
	
	//get the direction of a node based on the root direction (not apply for the root itself!!!)
	public boolean getDirection(boolean rootDir){
		return rootDir!=(vector.getDirection()>0); //XOR: in-out not +/-
	}
	@Override
	public String toString(){
		return node.getId() + ":" + vector.toString() + ":" + nvsScore;
	}
	
	//compare in term of distance to the root node: for sortedset
	@Override
	public int compareTo(BDNodeVecState o) {
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

        final BDNodeVecState other = (BDNodeVecState) obj;
        
        BDNode 	thisNode=this.getNode(),
				thatNode=other.getNode();
        
        if(!Objects.equals(thisNode, thatNode)){
        	return false;
        }else{ 
        	boolean retval;
        	if(SimpleBinner.getBinIfUnique(thatNode)!=null) {
//        		int dist=ScaffoldVector.composition(ScaffoldVector.reverse(other.getVector()), this.getVector()).distance(thatNode, thisNode);
//        		retval=(dist < -BDGraph.getKmerSize());
        		retval= (Math.abs(this.getVector().relDistance(thisNode) - other.getVector().relDistance(thatNode)) < thisNode.getNumber("len") - BDGraph.getKmerSize());
        	}
        	else {
        		retval=this.getVector().consistentWith(other.getVector());
        	}
    		return retval;

        }     

    }
	public boolean qc() {
		return nvsScore >= BDGraph.MIN_SCORE;
	}

	//Merging 2 equal NodeVectors
	public boolean merge(BDNodeVecState nv){
		if(!this.equals(nv))
			return false;

		ScaffoldVector 	thisVector=this.getVector(),
						thatVector=nv.getVector();
		//update the vector
		if(!thisVector.isIdentity()){
			thisVector.setMagnitute( (int)((thisVector.getMagnitute()*this.nvsScore+thatVector.getMagnitute()*nv.nvsScore)/(this.nvsScore+nv.nvsScore)));
		}
		
		//update coverage
		this.nvsScore+=nv.nvsScore;
		return true;
	}
}
