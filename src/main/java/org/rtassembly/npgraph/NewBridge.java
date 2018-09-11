package org.rtassembly.npgraph;

import java.util.ArrayList;
import org.graphstream.graph.Node;
import japsa.seq.Sequence;

public class NewBridge {
	public static volatile int MAX_DIFF=3;
	private int status=-1;
	BidirectedGraph graph; //partial order graph saving possible paths
	BidirectedEdgePrototype pBridge; //note: fist node of the bridge is set unique
	ArrayList<BridgeSegment> segments;
	//HashMap<Node, ScaffoldVector> segmentSteps;
	
	NewBridge(BidirectedGraph graph){
		this.graph=graph;
	}

	//This is for SPAdes path reader only. The input path must be elementary unique path
	//(2 ending nodes are unique and not containing other unique path)
	NewBridge (BidirectedGraph graph, BidirectedPath path) throws Exception{
		this(graph);
		if(	path.size()<=1
		 || (SimpleBinner.getUniqueBin(path.getRoot())==null && SimpleBinner.getUniqueBin(path.peekNode())==null))
			throw new Exception("Invalid path to build bridge: " + path.getId());
		
		segments=new ArrayList<>();
		segments.add(new BridgeSegment(path));
		pBridge=segments.get(0).pSegment;
		
		status=1;
	}
	
	public NewBridge(BidirectedGraph graph, AlignedRead bb) {
		this(graph);
		buildFrom(bb);
	}

	/*
	 * Most important function: to read from a AlignedRead and
	 * try to build or complete the bridge
	 */
	public void buildFrom(AlignedRead alignedRead) {
		//1.scan the aligned read for the marker and direction to build
		//2.compare the alignments to this bridge's steps and update & save nanopore read also if necessary
		if(alignedRead.getAlignmentRecords().size() < 2)
			return;
		
		// Starting node of the aligned read must be unique
		if(	SimpleBinner.getUniqueBin(alignedRead.getFirstAlignment().node)==null 
				&& SimpleBinner.getUniqueBin(alignedRead.getFirstAlignment().node)!=null)
				alignedRead=alignedRead.reverse();
		
		if(status==-1){ // empty bridge: build from beginning	
			if(SimpleBinner.getUniqueBin(alignedRead.getFirstAlignment().node)==null)
				return;
			for(int i=0;i<alignedRead.getAlignmentRecords().size()-1;i++)
				segments.add(new BridgeSegment(alignedRead.getAlignmentRecords().get(i), alignedRead.getAlignmentRecords().get(i+1), alignedRead));
			status=0;	
			
		}else if(status==0){ // building in progress...
			 if(alignedRead.getLastAlignment().node == pBridge.getNode0())
				 alignedRead=alignedRead.reverse();
			 else if(alignedRead.getFirstAlignment().node != pBridge.getNode0()) {
				 System.err.println("Disagree first node of the bridge! Ignore.");
				 return;
			 }
			 boolean firstDirOnRead = alignedRead.getFirstAlignment().strand;
			 if(firstDirOnRead != pBridge.getDir0()) {
				 System.err.println("Disagree first node of the bridge! Ignore.");
				 return;
			 }
			 
			 //now we have an agreement between first node of alignedRead and this bridge (id and direction!)
			 int begCoord = alignedRead.getFirstAlignment().readAlignmentStart();
			 int foundCtgIdx=1; //save the index of found record that closest to the current alignment 
			 
			 for(int idx=1; idx<alignedRead.getAlignmentRecords().size(); idx++){
				 Alignment alg = alignedRead.getAlignmentRecords().get(idx);
				 //calculate distance to the first alignment
				 Range algRange = new Range(alg.readAlignmentStart()-begCoord, alg.readAlignmentEnd()-begCoord);
				 //1. locate the segment containing the corresponding alignment
				 while(segments.get(foundCtgIdx).getRangeOfLastNode().compareTo(algRange) < 0){
					 foundCtgIdx++;
					 if(foundCtgIdx>=segments.size())
						 break;
						 
				 }
				 if(foundCtgIdx < segments.size()){
					 BridgeSegment searchSegment = segments.get(foundCtgIdx);
					 int 	distance = searchSegment.getRangeOfLastNode().getLeft() - algRange.getRight(),
					 		numOfPathsToSearch = searchSegment.getNumberOfPaths();
					 
					 //check if exist path connecting this contig to the ends of this bridge segment
					 BidirectedNode fromNode = (BidirectedNode) searchSegment.pSegment.getNode1(),
							 		toNode = alg.node;
					 
					 boolean fromDir= searchSegment.pSegment.getDir1(),
							 toDir = !alg.strand;
					 ArrayList<BidirectedPath> paths = graph.getClosestPaths(fromNode, fromDir, toNode, toDir, distance);
					 //	
					 
					 //if not check the next segment (cause errors from nanopore data)
					 
					 

							
				 }
				 
				 
				 
				 
			 }
			 
		}else{ //completed bridge: do nothing?
			
		}
				
	
		
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


	
	/************************************************************************************************
	 * Class to represent a single segment of the whole bridge.
	 * A bridge consists of >=1 segments.
	 ************************************************************************************************/
	public class BridgeSegment{
		ArrayList<Sequence> nnpReads; // to store nanopore data if needed
		ArrayList<BidirectedPath> connectedPaths;
		BidirectedEdgePrototype pSegment;
		Range coverRange;
		//TODO: scaffold vector??
		BridgeSegment(){}
		BridgeSegment(Alignment start, Alignment end, AlignedRead read){
			pSegment=new BidirectedEdgePrototype(start.node, end.node, start.strand, !end.strand);
			int startCoordinate = read.getFirstAlignment().readAlignmentStart();
			coverRange=new Range(start.readAlignmentStart()-startCoordinate, end.readAlignmentEnd()-startCoordinate);
			//invoke findPath()?
			
			connectedPaths = graph.getClosestPaths(start, end);
			
		}
		
		BridgeSegment(BidirectedPath path) throws Exception{
			pSegment=new BidirectedEdgePrototype(path);
			connectedPaths=new ArrayList<>();
			connectedPaths.add(path);
			
			//only if path cover the whole bridge (1-segment bridge)
			coverRange=new Range(0,(int) path.getLength());
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
		
		public Range getRangeOfFirstNode(){
			return new Range(coverRange.getLeft(), (int) (coverRange.getLeft() + pSegment.getNode0().getNumber("len")-1));

		}
		
		public Range getRangeOfLastNode(){
			return new Range((int)(coverRange.getRight() - pSegment.getNode1().getNumber("len") +1), coverRange.getRight());

		}
		
		//take information about a node (position+direction) and reduce the number of possible paths 
		public boolean pathsReducing(Node node, int distance){
			//true if the hint node is found
			boolean retval=false;
			
			return retval;
		}
	}



}
