package org.rtassembly.npgraph;

import java.util.ArrayList;
import java.util.HashMap;

import org.graphstream.graph.implementations.AbstractNode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class NewBridge {
	private static final Logger LOG = LoggerFactory.getLogger(NewBridge.class);
	public static volatile int MAX_DIFF=3;
	BidirectedGraph graph; //partial order graph saving possible paths
	PopBin bin;
	BidirectedEdgePrototype pBridge; //note: fist node of the bridge is set unique
	ArrayList<BridgeSegment> segments;
	//HashMap<Node, ScaffoldVector> segmentSteps;
	
	NewBridge(BidirectedGraph graph, PopBin bin){
		this.graph=graph;
		segments=new ArrayList<>();
		this.bin=bin;
	}
	//This is for SPAdes path reader only. The input path must be elementary unique path
	//(2 ending nodes are unique and not containing other unique path)
	NewBridge (BidirectedGraph graph, BidirectedPath path){
		this(graph,path.getConsensusUniqueBinOfPath());
		if(	path.size()>1
		 && ((SimpleBinner.getUniqueBin(path.getRoot())!=null || SimpleBinner.getUniqueBin(path.peekNode())!=null))
		 ){
			segments.add(new BridgeSegment(path));
			pBridge=segments.get(0).pSegment;
		}else
			LOG.warn("Invalid path to build bridge: " + path.getId());
		

	}
	
	public NewBridge(BidirectedGraph graph, AlignedRead bb, PopBin b) {
		this(graph,b);
		if(SimpleBinner.getUniqueBin(bb.getFirstAlignment().node)==null)
			bb=bb.reverse();
		buildFrom(bb);
		if(!segments.isEmpty()){
			BridgeSegment 	firstSeg = segments.get(0),
							lastSeg = segments.get(segments.size()-1);
			pBridge=new BidirectedEdgePrototype(firstSeg.pSegment.getNode0(), lastSeg.pSegment.getNode1(), firstSeg.pSegment.getDir0(), lastSeg.pSegment.getDir1());
		}
	}

	public byte getNumberOfAnchors() {
		byte retval=0;
		if(pBridge!=null) {
			if(SimpleBinner.getUniqueBin(pBridge.getNode0())!=null) retval++;
			if(SimpleBinner.getUniqueBin(pBridge.getNode1())!=null) retval++;
		}
		return retval;
	}
	
	public boolean isComplete() {
		if(getNumberOfAnchors()!=2)
			return false;
		for(BridgeSegment seg:segments)
			if(seg.getNumberOfPaths()!=1)
				return false;
		return true;
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
		alignedRead.sortAlignment();
		
		// Starting node of the aligned read must be unique
		if(	SimpleBinner.getUniqueBin(alignedRead.getFirstAlignment().node)==null 
				&& SimpleBinner.getUniqueBin(alignedRead.getFirstAlignment().node)!=null)
				alignedRead=alignedRead.reverse();
		
		if(segments.isEmpty()){ // empty bridge: build from beginning	
			if(SimpleBinner.getUniqueBin(alignedRead.getFirstAlignment().node)==null)
				return;
			for(int i=0;i<alignedRead.getAlignmentRecords().size()-1;i++){
				BridgeSegment tmp=new BridgeSegment(alignedRead.getAlignmentRecords().get(i), alignedRead.getAlignmentRecords().get(i+1), alignedRead);
				if(tmp.isConnected())
					segments.add(tmp);
				else
					return;
			}
			
		}else{ // building on the existed one
			 if(alignedRead.getLastAlignment().node == pBridge.getNode0()) {
				 alignedRead=alignedRead.reverse();
			 }
			 else if(alignedRead.getFirstAlignment().node != pBridge.getNode0()) {
				 System.err.println("Disagree starting points of the alignment to the bridge! Ignored.");
				 return;
			 }
			 boolean firstDirOnRead = alignedRead.getFirstAlignment().strand;
			 if(firstDirOnRead != pBridge.getDir0()) {
				 System.err.println("Conflict between alignnment [" + alignedRead.getFirstAlignment().toString() + "] and bridge [" + pBridge.toString() + "]");
				 return;
			 }
			 
			 System.out.println("Applying aligned read on the bridge: " + alignedRead.getCompsString());

			 //now we have an agreement between first node of alignedRead and this bridge (id and direction!)
			 int startSearchingIdx=0; //save the last browsed segment index that closest to the current alignment 
			 int idx;
			 for(idx=1; idx<alignedRead.getAlignmentRecords().size(); idx++){
				 Alignment alg = alignedRead.getAlignmentRecords().get(idx);
				 System.out.println("Checking alignment of " + alg.node.getId());
				 ScaffoldVector algVec = alignedRead.getVector(alignedRead.getFirstAlignment(), alg),
						 		segEndVec = null, segStartVec = null;
				 //1. locate the segments approximately close by the corresponding aligned contig
				 int 	curMinDis=Integer.MAX_VALUE,
						tmpIdx=startSearchingIdx;
				 ArrayList<BidirectedPath> agreePaths = new ArrayList<>();
				 boolean foundSegmentCandidates = false;
				 HashMap<Integer, ArrayList<Alignment>> breakpoints = new HashMap<>(); 
				 while(tmpIdx<segments.size()){
					 BridgeSegment searchSegment = segments.get(tmpIdx);
					 segEndVec = searchSegment.getEndVector();
					 segStartVec = searchSegment.getStartVector();
					 System.out.print("-searching alignment on segment " + searchSegment.pSegment.toString());

					 if(	Math.signum(algVec.getMagnitute()-segStartVec.getMagnitute())!=Math.signum(algVec.getMagnitute()-segEndVec.getMagnitute()) 
						 || GraphUtil.approxCompare(segStartVec.getMagnitute(), algVec.getMagnitute()) == 0
						 || GraphUtil.approxCompare(segEndVec.getMagnitute(), algVec.getMagnitute()) == 0) {
						 foundSegmentCandidates=true;
						 System.out.println(" : close!");
						 ScaffoldVector diffVec=ScaffoldVector.composition(algVec, ScaffoldVector.reverse(segStartVec));//seg-0->alg
						 
						 if(searchSegment.getNumberOfPaths()==0) {
							 //identify alignments falling into this segment
						 }else {
							 for(BidirectedPath path:searchSegment.connectedPaths) {
								 if(!path.contains(alg.node))
									 continue;
								 int dist = diffVec.distance((BidirectedNode) searchSegment.pSegment.getNode0(), alg.node);
								 System.out.printf("...estimate distance between %s and %s is %d\n", searchSegment.pSegment.getNode0().getId(), alg.node.getId(), dist);
								 //FIXME: optimize this, use ScaffoldVector instead of distance for checkDistanceConsistency()
								 if(path.checkDistanceConsistency(searchSegment.pSegment.getNode0(), alg.node, searchSegment.pSegment.getDir0()==alg.strand, dist) ) {
									 //add this path to a candidate list for later consider...
									 if(!agreePaths.isEmpty() && GraphUtil.approxCompare(curMinDis, dist)!=0)
										 continue; 
									 else if(Math.abs(curMinDis)>Math.abs(dist)) {
										 agreePaths = new ArrayList<>();
										 curMinDis=dist;
										 startSearchingIdx=tmpIdx;//save its segment also?
									 } 
									 agreePaths.add(path);
									 
								 }
							 }
							 if(agreePaths.isEmpty())
								 startSearchingIdx++;
						 }
						 
					 }else {
						 System.out.println(": not close!");
						 if(foundSegmentCandidates){
							 System.out.println("Stop going further!");
							 break;
						 }
					 }
					 
					 
					 tmpIdx++;
					 
				 }
				 //reducing possible paths
				 if(!agreePaths.isEmpty()) {					 
					 //TODO: implement vote system here: 195i,76o
					 segments.get(startSearchingIdx).connectedPaths=agreePaths;
				 }else {
					 //TODO: how about pseudo edge (no path found)
					 
					 if(tmpIdx==segments.size()){ 			 
						 System.out.println("try extend instead!");
						 break;
					 }
				 }
			
			 }
			 if(idx < alignedRead.getAlignmentRecords().size()-1) {
				 //this aligned read is longer and have more information than this bridge
				 //first connect last node of bridge to first out-of-range node from alignedRead
				 Alignment lastBrowseAlg = alignedRead.getAlignmentRecords().get(idx);
				 BridgeSegment cnt = new BridgeSegment(
						 							pBridge.getNode1(), 
						 							lastBrowseAlg.node, 
						 							!pBridge.getDir1(), 
						 							!lastBrowseAlg.strand,
						 							segments.get(segments.size()-1).getEndVector(),
						 							alignedRead.getVector(alignedRead.getFirstAlignment(), lastBrowseAlg)
						 							);
				 if(cnt.getNumberOfPaths() > 0) {//to avoid error reads that give crap alignments???
					 segments.add(cnt);
					 //add others
					 for(int i=idx;i<alignedRead.getAlignmentRecords().size()-1;i++)
						 segments.add(new BridgeSegment(alignedRead.getAlignmentRecords().get(i), alignedRead.getAlignmentRecords().get(i+1), alignedRead));
				 }else
					 System.out.println("Cannot extend!");
		
			 }
		}
							
	
		
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
	

	public String getAllPossiblePaths() {
		String retval = "{\n";
		for(BridgeSegment seg:segments) {
			retval+=seg.pSegment.toString() + ": ";
			if(seg.getNumberOfPaths()==0)
				retval += "()";
			else
				for(BidirectedPath path:seg.connectedPaths)
					retval+="( "+path.getId()+" )";
			retval+="\n";
		}
		retval+="}";
		return retval;
	}

	public BidirectedPath getBestPath() {
		BidirectedPath retval=null;
		for(BridgeSegment seg:segments) {
			if(seg.getNumberOfPaths()>0){
				if(retval==null)
					retval=seg.connectedPaths.get(0);
				else
					retval=retval.join(seg.connectedPaths.get(0));//assuming first path has highest score!!!
			}else
				return null;
		}
		if(retval!=null)
			retval.setConsensusUniqueBinOfPath(bin);
		return retval;
	}

	
	/************************************************************************************************
	 * Class to represent a single segment of the whole bridge.
	 * A bridge consists of >=1 segments.
	 ************************************************************************************************/
	public class BridgeSegment{
//		ArrayList<Sequence> nnpReads; // to store nanopore data if needed
		ArrayList<BidirectedPath> connectedPaths;
		BidirectedEdgePrototype pSegment;
//		Range coverRange;
		ScaffoldVector startV, endV; // from bridge anchor (always +) to this segment's end
		//TODO: scaffold vector??
		BridgeSegment(){}
		BridgeSegment(Alignment start, Alignment end, AlignedRead read){
			pSegment=new BidirectedEdgePrototype(start.node, end.node, start.strand, !end.strand);
			//invoke findPath()?
			startV=read.getVector(read.getFirstAlignment(), start);
			endV=read.getVector(read.getFirstAlignment(), end);
			connectedPaths = graph.getClosestPaths(start, end);
			
		}
		BridgeSegment(AbstractNode node1, AbstractNode node2, boolean dir1, boolean dir2, ScaffoldVector v1, ScaffoldVector v2){
			assert Math.abs(v1.getMagnitute()) < Math.abs(v2.magnitude): "Illegal order of node position!";
			pSegment = new BidirectedEdgePrototype(node1,node2,dir1,dir2);
			startV = v1; endV = v2;
			int d = ScaffoldVector.composition(endV, ScaffoldVector.reverse(startV)).distance((BidirectedNode)node1, (BidirectedNode)node2);
			connectedPaths = graph.getClosestPaths((BidirectedNode)node1, dir1, (BidirectedNode)node2, dir2, d, false);
			
		}
		BridgeSegment(BidirectedPath path){
			try {
				pSegment=new BidirectedEdgePrototype(path);
				connectedPaths=new ArrayList<>();
				connectedPaths.add(path);
				
				//only if path cover the whole bridge (1-segment bridge)
				int dist=(int) (path.getLength()
						-(pSegment.getDir0()?0:pSegment.getNode0().getNumber("len"))
						-(pSegment.getDir1()?0:pSegment.getNode1().getNumber("len")));
				
				startV=new ScaffoldVector(0,1);
				endV=new ScaffoldVector(pSegment.getDir0()?dist:-dist, pSegment.getDir0()!=pSegment.getDir1()?1:-1);
			} catch (Exception e) {
				LOG.error("Cannot make bridge from this path!");
				e.printStackTrace();
			} 

		}
		
		public ArrayList<BridgeSegment> breakSegments(ArrayList<Alignment> algs){
			ArrayList<BridgeSegment> retval = new ArrayList<BridgeSegment>();
			//break and return list of new segments
			
			return retval;
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
		public ScaffoldVector getEndVector() {
			return endV;
		}
		public ScaffoldVector getStartVector() {
			return startV;
		}

	}





}
