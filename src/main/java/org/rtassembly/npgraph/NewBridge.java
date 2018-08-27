package org.rtassembly.npgraph;

import java.util.ArrayList;
import java.util.HashMap;

import org.graphstream.graph.Node;
import org.rtassembly.scaffold.ScaffoldVector;

import japsa.seq.Sequence;

public class NewBridge {
	public static volatile int MAX_DIFF=3;
	
	BidirectedGraph graph; //partial order graph saving possible paths
	BidirectedEdgePrototype pBridge; //must be unique nodes
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
		
		segments=new ArrayList<>();
		segments.add(new BridgeSegment(path));
		pBridge=segments.get(0).pSegment;
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
	//-1: half bridge, 0: complete bridge, unsolved; 1: solved
	public int getBridgeStatus(){
		int retval=-1;

		return retval;
	}
	
	public String getEndingsID() {
		if(pBridge==null)
			return "-,-";
		else
			return pBridge.toString();
	}

	//Return the path on top of the possible list
	public BidirectedPath getBestPathPossible(){
		//TODO: implement it!
		return null;
	}
	/*
	 * Class to represent a single segment of the whole bridge.
	 * A bridge consists of >=1 segments.
	 */
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
