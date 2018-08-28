package org.rtassembly.npgraph;

import java.util.ArrayList;
import java.util.HashMap;

import org.graphstream.graph.Node;
import org.rtassembly.scaffold.ScaffoldVector;

import japsa.seq.Sequence;

public class NewBridge {
	public static volatile int MAX_DIFF=3;
	private int status=-1;
	BidirectedGraph graph; //partial order graph saving possible paths
	BidirectedEdgePrototype pBridge; //must be unique nodes
	ArrayList<BridgeSegment> segments;
	HashMap<Node, ScaffoldVector> segmentSteps;
	
	NewBridge(BidirectedGraph graph){
		this.graph=graph;
	}

	//This is for SPAdes path reader only. The input path must be elementary unique path
	//(2 ending nodes are unique and not containing other unique path)
	NewBridge (BidirectedPath path) throws Exception{
		if(	path.size()<=1
		 || (SimpleBinner.getUniqueBin(path.getRoot())==null && SimpleBinner.getUniqueBin(path.peekNode())==null))
			throw new Exception("Invalid path to build bridge: " + path.getId());
		
		segments=new ArrayList<>();
		segments.add(new BridgeSegment(path));
		pBridge=segments.get(0).pSegment;
		
		status=1;
	}
	
	public NewBridge(AlignedRead bb) {
		// TODO Auto-generated constructor stub
	}

	/*
	 * Most important function: to read from a AlignedRead and
	 * try to build or complete the bridge
	 */
	public void buildFrom(AlignedRead alignedRead, int from, int to) {
		//1.scan the aligned read for the marker and direction to build
		ArrayList<Alignment> alignments = alignedRead.getAlignmentRecords();
		
		for(Alignment alg:alignments){
			
		}
		
		//2.compare the alignments to this bridge's steps and update & save nanopore read also if necessary
		
	
		
	}
	//Return bridge status:
	//-1: none or half bridge, 0: unsolved full bridge; 1: solved bridge
	public int getBridgeStatus(){
		return status;
	}
	
	public String getEndingsID() {
		if(pBridge==null)
			return "-,-";
		else
			return pBridge.toString();
	}
	
	@Override
	public String toString() {
		String retval="";
		//...print all steps
		
		return retval;
	}

	//Return the path on top of the possible list
	public BidirectedPath getBestPathPossible(){
		//TODO: implement it!
		return null;
	}
	

	public ArrayList<BidirectedPath> getAllPossiblePaths() {
		// TODO Auto-generated method stub
		return null;
	}

	public BidirectedPath getBestPath() {
		// TODO Auto-generated method stub
		return null;
	}

	public void compareTo(AlignedRead bb) {
		// TODO Auto-generated method stub
		
	}

	
	/************************************************************************************************
	 * Class to represent a single segment of the whole bridge.
	 * A bridge consists of >=1 segments.
	 ************************************************************************************************/
	public class BridgeSegment{
		ArrayList<Sequence> nnpReads; // to store nanopore data if needed
		ArrayList<BidirectedPath> connectedPaths;
		BidirectedEdgePrototype pSegment;
		//TODO: scaffold vector??
		BridgeSegment(){}
		BridgeSegment(Alignment start, Alignment end, Sequence read){
			pSegment=new BidirectedEdgePrototype(start.node, end.node, start.strand, !end.strand);
			//invoke findPath()?
			
			int distance = start.readAlignmentStart()-end.readAlignmentEnd();
			if(distance<BidirectedGraph.D_LIMIT)
				connectedPaths = graph.getClosestPaths(start, end, distance);
			
		}
		BridgeSegment(BidirectedPath path) throws Exception{
			pSegment=new BidirectedEdgePrototype(path);
			connectedPaths=new ArrayList<>();
			connectedPaths.add(path);
			
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
