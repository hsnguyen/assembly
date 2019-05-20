package org.rtassembly.npgraph;

import java.util.List;
import java.util.Objects;
import java.util.TreeSet;

/*
 * Class represent a node together with its vector in relative to another root node.
 * Note: the direction here refers to the one of a ScaffoldVector, not bidirected-graph direction (in/out), neither sequence sense/antisense (+/-)
 * although we can translate into appropriate info
 */
public class BDNodeVecState implements Comparable<BDNodeVecState>{
	BDNode node;
	ScaffoldVector vector;
	int nvsScore=Alignment.MIN_QUAL; //alignment score + number of occurences?
	
	public BDNodeVecState(BDNode node, ScaffoldVector vector){
		this.node=node;
		this.vector=vector;
	}
	
	public BDNodeVecState(BDNode node, int qual, ScaffoldVector vector){
		this(node, vector);
		nvsScore=qual;
	}
	public BDNodeVecState(Alignment alg, ScaffoldVector vector) {
		this(alg.node, alg.quality, vector);
	}
	
	public BDNode getNode(){return node;}
	public ScaffoldVector getVector(){return vector;}
	
	public void setNode(BDNode node){this.node=node;}
	public void setVector(ScaffoldVector vector){this.vector=vector;}
	
	
	//get the direction of a node based on the root direction (not apply for the root itself!!!)
	public boolean getDirection(boolean rootDir){
		return rootDir!=(vector.getDirection()>0); //XOR: in-out not +/-
	}
	public int getScore(){
		return nvsScore;
	}
	public void setScore(int score){
		nvsScore=score;
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
//        	if(SimpleBinner.getBinIfUnique(thatNode)!=null) {
        		retval= (Math.abs(this.getVector().relDistance(thisNode) - other.getVector().relDistance(thatNode)) < thisNode.getNumber("len") - BDGraph.getKmerSize());
//        	}
//        	else {
//        		retval=this.getVector().consistentWith(other.getVector());
//        	}
    		return retval;

        }     

    }
    
	public boolean qc() {
		return nvsScore >= BDGraph.SAFE_COUNTS*Alignment.GOOD_QUAL; //60*2
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
	
	/*
	 * Needleman-Wunsch algorithm to return consensus steps
	 * by aligning two step-lists with identical first element
	 */
	public static TreeSet<BDNodeVecState> NWAlignment(TreeSet<BDNodeVecState> s0, TreeSet<BDNodeVecState> s1){
		TreeSet<BDNodeVecState> retval=new TreeSet<>();
		if(!s0.first().equals(s1.first()))
			return null;
		
		BDNodeVecState[] 	as0=s0.toArray( new BDNodeVecState[s0.size()]),
							as1=s1.toArray( new BDNodeVecState[s1.size()]);
		int[][] scoreTab = new int[as0.length][as1.length];
		for(int i=0; i<as0.length; i++)
			scoreTab[i][0] = -as0[i].getDistance();
		for(int i=0; i<as1.length; i++)
			scoreTab[0][i] = -as1[i].getDistance();
		int match, delete, insert;
		for(int i=1; i<as0.length; i++){
			for(int j=1; j<as1.length; j++){
				match=(as0[i].equals(as1[j])
						?scoreTab[i-1][j-1] + (int)as0[i].getNode().getNumber("len")-Math.abs(as0[i].getDistance()-as1[j].getDistance())
						:Integer.MIN_VALUE);
				delete=scoreTab[i-1][j] - as0[i].getDistance() + as0[i-1].getDistance();
				insert=scoreTab[i][j-1] -as1[j].getDistance() +as1[j-1].getDistance();
				scoreTab[i][j]= Math.max(delete>insert?delete:insert, match);
			}
		}
		
		System.out.println("Needleman-Wunch table: ");
		for(int i=0; i<as1.length; i++)
			System.out.printf("\t%s", as1[i].getNode().getId());
		System.out.println();
		for(int i=0; i<as0.length; i++){
			System.out.printf("%s\t", as0[i].getNode().getId());
			for(int j=0; j<as1.length; j++){
				System.out.printf("%d\t", scoreTab[i][j]);
			}
			System.out.println();
		}
		
		return retval;
	}


	public int getDistance(){
		return vector.relDistance(node);
	}
}
