package org.rtassembly.npgraph;

import java.util.ArrayList;
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
        	if(SimpleBinner.getBinIfUniqueNow(thatNode)!=null) {
        		retval= (Math.abs(this.getVector().relDistance(thisNode) - other.getVector().relDistance(thatNode)) < thisNode.getNumber("len") - BDGraph.getKmerSize());
        	}
        	else {
        		retval=this.getVector().consistentWith(other.getVector());
        	}
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
	//Another merge() that create the result
	public static BDNodeVecState merge(BDNodeVecState n0, BDNodeVecState n1){
		if(!n0.equals(n1))
			return null;
		
		ScaffoldVector 	v0=n0.getVector(),
						v1=n1.getVector();
		BDNodeVecState retval=new BDNodeVecState(n0.getNode(), n0.getScore()+n1.getScore(), new ScaffoldVector());

		//update the vector
		if(!v0.isIdentity()){
			retval.vector.setMagnitute( (int)((v0.getMagnitute()*n0.getScore()+v1.getMagnitute()*n1.getScore())/(n0.getScore()+n1.getScore())));
		}
		return retval;
	}
	
	/*
	 * Needleman-Wunsch algorithm to return consensus steps
	 * by aligning two step-lists with identical first element
	 */
	public static TreeSet<BDNodeVecState> NWAlignment(TreeSet<BDNodeVecState> s0, TreeSet<BDNodeVecState> s1){
		assert s0.first().equals(s1.first()):"First alignment must agree!";

		TreeSet<BDNodeVecState> retval=new TreeSet<>();
		
		BDNodeVecState[] 	as0=s0.toArray( new BDNodeVecState[s0.size()]),
							as1=s1.toArray( new BDNodeVecState[s1.size()]);
		int[][] scoreTab = new int[as0.length][as1.length];
		char[][] moveTab = new char[as0.length][as1.length];
		for(int i=0; i<as0.length; i++){
			scoreTab[i][0] = -as0[i].getDistance();
			moveTab[i][0] = '|';
		}
		for(int i=0; i<as1.length; i++){
			scoreTab[0][i] = -as1[i].getDistance();
			moveTab[0][i] = '_';
		}
		int match, delete, insert;
		for(int i=1; i<as0.length; i++){
			for(int j=1; j<as1.length; j++){
				match=(as0[i].equals(as1[j])
						?scoreTab[i-1][j-1] + (int)as0[i].getNode().getNumber("len")-Math.abs(as0[i].getDistance()-as1[j].getDistance())
						:Integer.MIN_VALUE);
				delete=scoreTab[i-1][j] - as0[i].getDistance() + as0[i-1].getDistance();
				insert=scoreTab[i][j-1] -as1[j].getDistance() +as1[j-1].getDistance();
				scoreTab[i][j]= Math.max(delete>insert?delete:insert, match);
				if(scoreTab[i][j]==match)
					moveTab[i][j]='\\';
				else if(scoreTab[i][j]==insert)
					moveTab[i][j]='_';
				else if(scoreTab[i][j]==delete)
					moveTab[i][j]='|';

			}
		}
		
		System.out.println("Needleman-Wunch table: ");
		for(int i=0; i<as1.length; i++)
			System.out.printf("\t%s", as1[i].getNode().getId());
		System.out.println();
		for(int i=0; i<as0.length; i++){
			System.out.printf("%s\t", as0[i].getNode().getId());
			for(int j=0; j<as1.length; j++){
				System.out.printf("%s%d\t", moveTab[i][j], scoreTab[i][j]);
			}
			System.out.println();
		}
		
		int i=as0.length-1, j=as1.length-1;
		ArrayList<Integer> 	m0 = new ArrayList<>(),
							m1 = new ArrayList<>();

		while(i>=0 && j>=0){
			switch(moveTab[i][j]){
				case '\\':
					m0.add(0, i--);
					m1.add(0, j--);
					break;
				case '_':
					j--;
					break;
				case '|':
					i--;
					break;
			}
		}
//		//and we always has a match at (0,0)
//		i=j=0;
//		
//		BDNodeVecState 	prevMatch=as0[0],
//						nextMatch=null, tmp=null; 
//		retval.add(prevMatch);
//		for(int idx=0;idx<m0.size();idx++){
//			int ii=m0.get(idx), jj=m1.get(idx);
//			nextMatch=merge(as0[ii], as1[jj]);
//			//add the vectors from first in-between
//			for(int i0=i+1;i0<ii;i0++){
//				tmp = new BDNodeVecState(as0[i0].getNode(), as0[i0].getScore(), as0[i0].getVector());
//				tmp.vector.setMagnitute(prevMatch.getVector().getMagnitute() + 
//										(as0[i0].getVector().getMagnitute() - as0[i].getVector().getMagnitute())
//										*(nextMatch.getVector().getMagnitute()-prevMatch.getVector().getMagnitute())
//										/(as0[ii].getVector().getMagnitute()-as0[i].getVector().getMagnitute()));
//				retval.add(tmp);
//			}
//			//add vectors from second in-between
//			for(int j0=j+1;j0<jj;j0++){
//				tmp = new BDNodeVecState(as1[j0].getNode(), as1[j0].getScore(), as1[j0].getVector());
//				tmp.vector.setMagnitute(prevMatch.getVector().getMagnitute() + 
//										(as1[j0].getVector().getMagnitute() - as1[j].getVector().getMagnitute())
//										*(nextMatch.getVector().getMagnitute()-prevMatch.getVector().getMagnitute())
//										/(as1[jj].getVector().getMagnitute()-as1[j].getVector().getMagnitute()));
//				retval.add(tmp);
//			}
//			
//			//add the vector of next match
//			prevMatch=nextMatch;
//			retval.add(prevMatch);
//			//move to next match coordinate
//			i=ii; j=jj;
//		}
//		
//		System.out.println("After NW merging: ");
//		for(BDNodeVecState nv:retval)
//			System.out.println(nv);
		
		return retval;
	}


	public int getDistance(){
		return vector.relDistance(node);
	}
}
