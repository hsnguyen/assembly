package org.rtassembly.npgraph;

public class POGBridge {
	public static volatile int SAFE_VOTE_DISTANCE=3;
	BidirectedGraph bridgeGraph; //partial order graph saving possible paths
	BidirectedNode source, sink;
	
	
	POGBridge(){
		bridgeGraph=new BidirectedGraph();
	}

	//This is for SPAdes path reader only
	POGBridge (BidirectedPath path) throws Exception{
		if(	path.size()<=1
		 || (SimpleBinner.getUniqueBin(path.getRoot())==null && SimpleBinner.getUniqueBin(path.peekNode())==null))
			throw new Exception("Invalid path to build bridge: " + path.getId());
		
		if(SimpleBinner.getUniqueBin(path.getRoot())!=null)
			source=(BidirectedNode) path.getRoot();
		if(SimpleBinner.getUniqueBin(path.peekNode())!=null)
			sink=(BidirectedNode) path.peekNode();
		
		//...
	}
	

	public void addHalfBridge(POGBridge brg){

	}
	
	public boolean append(Alignment alg) {
		return false;
		
	}
	//Return bridge status:
	//-1: half bridge, 0: complete bridge, unsolved; 1: solved
	public int getBridgeStatus(){
		int retval=-1;

		return retval;
	}
	
	
	public BidirectedPath getBestPath() {
		//return best path (unique path or best voted path); or null if undetermined
		return null;
	}
	
	public String getEndingsID() {
		return null;
		
		
	}
	public String getBridgeString() {
		String retval="<unknown>";

		return retval;
	}
	
	//Using another list of steps to rectify the current bridge's steps
	//This bridge HAS TO BE complete(2 unique ends), 
	//the reference need not be complete, instead one unique end is enough (half bridge)
	public void referencingTo(POGBridge brg){

	}
	//Merging with another half bridge, inheriting all info from it
	public void merging(POGBridge brg){

	}

	//Start to build bridge (i.e. identify possible paths) when 2 unique ends are determined
	public void bridging(BidirectedGraph graph, PopBin bin){

		
	}
		
	public String getAllPossiblePathsString(){
		String retval="";

		
		return retval;
	}
	

	//check if a bridge is worth to merge to this bridge
	public boolean checkIfMerge(POGBridge brg) {
		boolean retval=false;

		
		return retval;
	}
}
