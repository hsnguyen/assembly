package org.rtassembly.npgraph;

import java.util.ArrayList;
import java.util.HashMap;

import org.graphstream.graph.Node;
import org.rtassembly.scaffold.ScaffoldVector;

import japsa.seq.Sequence;

public class NewBridge {
	public static volatile int SAFE_VOTE_DISTANCE=3;
	BidirectedGraph graph; //partial order graph saving possible paths
	BidirectedNode source, sink; //must be unique nodes
	ArrayList<BridgeSegment> segments;
	HashMap<Node, ScaffoldVector> segmentSteps;
	
	NewBridge(BidirectedGraph graph){
		this.graph=graph;
	}

	//This is for SPAdes path reader only
	NewBridge (BidirectedPath path) throws Exception{
		if(	path.size()<=1
		 || (SimpleBinner.getUniqueBin(path.getRoot())==null && SimpleBinner.getUniqueBin(path.peekNode())==null))
			throw new Exception("Invalid path to build bridge: " + path.getId());
		
		if(SimpleBinner.getUniqueBin(path.getRoot())!=null)
			source=(BidirectedNode) path.getRoot();
		if(SimpleBinner.getUniqueBin(path.peekNode())!=null)
			sink=(BidirectedNode) path.peekNode();
		segments=new ArrayList<>();
		segments.add(new BridgeSegment(path));
	}
	

	public void addHalfBridge(NewBridge brg){

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
	public void referencingTo(NewBridge brg){

	}
	//Merging with another half bridge, inheriting all info from it
	public void merging(NewBridge brg){

	}

	//Start to build bridge (i.e. identify possible paths) when 2 unique ends are determined
	public void bridging(BidirectedGraph graph, PopBin bin){

		
	}
		
	public String getAllPossiblePathsString(){
		String retval="";

		
		return retval;
	}
	

	//check if a bridge is worth to merge to this bridge
	public boolean checkIfMerge(NewBridge brg) {
		boolean retval=false;

		
		return retval;
	}
	
	/*
	 * Class to represent a single segment of the whole bridge.
	 * A bridge consists of >=1 segments.
	 */
	public class BridgeSegment{
		ArrayList<BidirectedPath> connectedPaths;
		Node start, end;
		//TODO: scaffold vector??
		BridgeSegment(){}
		BridgeSegment(Alignment start, Alignment end, Sequence read){
			this.start=start.node;
			this.end=end.node;
			//invoke findPath()?
			
			int distance = start.readAlignmentStart()-end.readAlignmentEnd();
			
			if(distance<BidirectedGraph.D_LIMIT)
				connectedPaths = graph.getClosestPaths(start, end, distance);
			
		}
		BridgeSegment(BidirectedPath path) throws Exception{
			if(path==null || path.size()<=1)
				throw new Exception("Invalid path for a bridge segment: " + path.getId());
			else{
				start=path.getRoot();
				end=path.peekNode();
				connectedPaths=new ArrayList<>();
				connectedPaths.add(path);
			}
		}
		
		public int getNumberOfPaths(){
			if(connectedPaths==null)
				return 0;
			else 
				return connectedPaths.size();
		}
		public boolean isConnected(){
			return getNumberOfPaths()>=1;
		}
		public boolean isUnique(){
			return getNumberOfPaths()==1;
		}
	}
}
