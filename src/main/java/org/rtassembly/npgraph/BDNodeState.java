package org.rtassembly.npgraph;

class BDNodeState implements Comparable<BDNodeState>{
	BDNode node;
	boolean dir;
	int weight;//for shortest path finding
	
	BDNodeState(BDNode node, boolean dir){
		this.node=node;
		this.dir=dir;
	}
	BDNodeState(BDNode node, boolean dir, int weight){
		this.node=node;
		this.dir=dir;
		this.weight=weight;
	}
	public BDNode getNode(){return node;}
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
	public int compareTo(BDNodeState o) {
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
        final BDNodeState other = (BDNodeState) obj;
        if (this.toString().equals(other.toString()))   
        	return true;
        else 
        	return false;
    }
}
