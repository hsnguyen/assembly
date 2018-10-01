package org.rtassembly.npgraph;

import java.util.ArrayList;
import java.util.HashMap;

import org.graphstream.graph.implementations.AbstractNode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

//A bridge structure that the first node must be unique
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
	//This is for SPAdes path reader only. The input path must be have at least one end unique.
	//This unique ending node will be use as the first anchor of the bridge
	NewBridge (BidirectedGraph graph, BidirectedPath path){
		this(graph,path.getConsensusUniqueBinOfPath());
		if(	path.size()>1) {
			if(SimpleBinner.getUniqueBin(path.getRoot())!=null)
				addSegment(new BridgeSegment(path));
			else if(SimpleBinner.getUniqueBin(path.peekNode())!=null)
				addSegment(new BridgeSegment(path.reverse()));
			
		}
		
		

	}
	
	public NewBridge(BidirectedGraph graph, AlignedRead bb, PopBin b) {
		this(graph,b);

		buildFrom(bb);
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
	 * try to build or complete the bridge. Return true if the bridge
	 * is extended successfully to a new anchor
	 */
	public boolean buildFrom(AlignedRead alignedRead) {
		//1.scan the aligned read for the marker and direction to build
		//2.compare the alignments to this bridge's steps and update & save nanopore read also if necessary
		
		if(alignedRead.getAlignmentRecords().size() < 2)
			return false;
		alignedRead.sortAlignment();
		
		// Starting node of the aligned read must be unique
		if(	SimpleBinner.getUniqueBin(alignedRead.getFirstAlignment().node)==null 
				&& SimpleBinner.getUniqueBin(alignedRead.getLastAlignment().node)!=null)
				alignedRead=alignedRead.reverse();
		
		if(SimpleBinner.getUniqueBin(alignedRead.getFirstAlignment().node)==null)
			return false;
		
		BidirectedNode end=null;
		if(segments.isEmpty()){ // empty bridge: build from beginning	
			addAligments(alignedRead, 0, alignedRead.getAlignmentRecords().size()-1);			
		}else{ // building on the existed one
//			assert SimpleBinner.getUniqueBin(pBridge.getNode0())!=null:"First node must unique: illegal bridge "+ this.getEndingsID();
//			assert SimpleBinner.getUniqueBin(alignedRead.getFirstAlignment().node)!=null:"First node must unique: wrong read "+ alignedRead.getEndingsID();
			end=(BidirectedNode) pBridge.getNode1();
			
			if(pBridge.getNode0()!=alignedRead.getFirstAlignment().node) {
				if(pBridge.getNode0()==alignedRead.getLastAlignment().node) {
					alignedRead=alignedRead.reverse();
				}else{
					this.reverse();
					if(pBridge.getNode1()==alignedRead.getLastAlignment().node)
						alignedRead=alignedRead.reverse();
				}
			}
			 
			 System.out.println("Applying aligned read on the bridge: " + alignedRead.getCompsString());

			 //now we have an agreement between first node of alignedRead and this bridge (id and direction!)
			 int startSearchingIdx=0; //save the last browsed segment index that closest to the current alignment 
			 int idx;
			 for(idx=1; idx<alignedRead.getAlignmentRecords().size(); idx++){
				 Alignment alg = alignedRead.getAlignmentRecords().get(idx);
				 System.out.println("Checking alignment of " + alg.node.getId());
				 //0. Ignore if the node is unique: either this bridge's anchor or hidden anchor
//				 if(getNumberOfAnchors()==2 && SimpleBinner.getUniqueBin(alg.node)!=null) {
//					 System.out.println("Stop at unique node!");
//					 break;
//				 }
				 
				 ScaffoldVector algVec = alignedRead.getVector(alignedRead.getFirstAlignment(), alg),
						 		segEndVec = null, segStartVec = null;
				 //1. locate the segments approximately close by the corresponding aligned contig
				 int 	minDeviation=Integer.MAX_VALUE,
						tmpIdx=startSearchingIdx;
				 ArrayList<BidirectedPath> agreePaths = new ArrayList<>();
				 boolean foundSegmentCandidates = false;
				 while(tmpIdx<segments.size() && agreePaths.isEmpty()){
					 BridgeSegment searchSegment = segments.get(tmpIdx);
					 segEndVec = searchSegment.getEndVector();
					 segStartVec = searchSegment.getStartVector();
					 System.out.print("-searching alignment on segment " + searchSegment.pSegment.toString());

					 //TODO: review this!!
					 if(	Math.signum(algVec.getMagnitute()-segStartVec.getMagnitute())!=Math.signum(algVec.getMagnitute()-segEndVec.getMagnitute()) 
						 || GraphUtil.approxCompare(segStartVec.getMagnitute(), algVec.getMagnitute()) == 0
						 || GraphUtil.approxCompare(segEndVec.getMagnitute(), algVec.getMagnitute()) == 0) 
					 {
						 foundSegmentCandidates=true;
						 System.out.println(" : close!");
						 ScaffoldVector diffVec=ScaffoldVector.composition(algVec, ScaffoldVector.reverse(segStartVec));//seg-0->alg
						 
						 if(searchSegment.isConnected()) {
							 int score=0;
							 for(BidirectedPath path:searchSegment.connectedPaths) {
								 if(!path.contains(alg.node))
									 continue;
								 int dist = diffVec.distance((BidirectedNode) searchSegment.pSegment.getNode0(), alg.node);
								 System.out.printf("...estimate distance between %s and %s is %d\n", searchSegment.pSegment.getNode0().getId(), alg.node.getId(), dist);
								 int delta=path.checkDistanceConsistency(searchSegment.pSegment.getNode0(), alg.node, searchSegment.pSegment.getDir0()==alg.strand, dist);
								 if( delta != -1) {
									 //add this path to a candidate list for later consider...
									 if(!agreePaths.isEmpty() && GraphUtil.approxCompare(minDeviation, delta)!=0)
										 continue; 
									 else {
										 //FIXME: voting based on edit distance with the long read instead: deviation voting quite useless here!!!
										 if(Math.abs(minDeviation)>Math.abs(delta)) {
											 minDeviation=delta;
											 path.upVote(++score);
											 agreePaths.add(0,path);
										 }else if(Math.abs(minDeviation)==Math.abs(delta)) {
											 path.upVote(score);
											 agreePaths.add(0,path);
										 }else
											 agreePaths.add(path);
										 
										 startSearchingIdx=Math.max(startSearchingIdx, tmpIdx);
										 
									 }
								 }
							 }

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
					 //TODO: look at vote score and keep the best(s)
					 segments.get(startSearchingIdx).connectedPaths=agreePaths;
				 }else {
					 if(tmpIdx==segments.size()){ 			 
						 System.out.println("try extend instead!");
						
						 //Extending mode:
						 //this aligned read is longer and have more information than this bridge
						 //first connect last node of bridge to first out-of-range node from alignedRead
						 BridgeSegment cnt = new BridgeSegment(
								 							pBridge.getNode1(), 
								 							alg.node, 
								 							!pBridge.getDir1(), 
								 							!alg.strand,
								 							segments.get(segments.size()-1).getEndVector(),
								 							alignedRead.getVector(alignedRead.getFirstAlignment(), alg)
								 							);
						 if(cnt.isConnected()) {//to avoid error reads that give crap alignments???
							 addSegment(cnt);
							 //add others
							 addAligments(alignedRead, idx, alignedRead.getAlignmentRecords().size()-1);
							 System.out.println("Extend to have " + pBridge.toString());
							 return true;

						 }else{
							 System.out.println("Path not found: cannot extend " + cnt.pSegment.toString());
							 continue;
						 }
					
						 
						 
					 }
				 }
			
			 }

		}
		
		if(getNumberOfAnchors()==2 && pBridge.getNode1()!=end)
			return true;
		else 
			return false;
	}
							

	private void addAligments(AlignedRead read, int start, int end) {
		 for(int i=start;i<end;i++) {
			 BridgeSegment tmp = new BridgeSegment(read.getAlignmentRecords().get(i), read.getAlignmentRecords().get(i+1), read);
			 if(tmp.isConnected())
				 addSegment(tmp);
			 else return;
		 }
	}
	private void addSegment(BridgeSegment seg) {
		
		if(pBridge==null) {
			segments.add(seg);
			pBridge=new BidirectedEdgePrototype(seg.pSegment.getNode0(), seg.pSegment.getNode1(), seg.pSegment.getDir0(), seg.pSegment.getDir1());
		}
		else if(pBridge.getNode1()==seg.pSegment.getNode0() && pBridge.getDir1()!=seg.pSegment.getDir0()){
			segments.add(seg);
			pBridge.n1=seg.pSegment.n1;
		}
	}
	
	private void reverse() {
		ArrayList<BridgeSegment> tmp = segments;
		segments = new ArrayList<>();
		pBridge=null;
		if(!tmp.isEmpty()) {
			ScaffoldVector brgVector=tmp.get(tmp.size()-1).endV;
			for(int i=tmp.size()-1; i>=0 ;i--)
				addSegment(tmp.get(i).reverse(brgVector));
		}
		
	}
	
	public String getEndingsID() {
		if(pBridge==null)
			return "-,-";
		else
			return pBridge.toString();
	}
	

	public String getAllPossiblePaths() {
		String retval = "{\n";
		for(BridgeSegment seg:segments) {
			retval+=seg.pSegment.toString() + ": ";
			if(seg.getNumberOfPaths()==0)
				retval += "()";
			else
				for(BidirectedPath path:seg.connectedPaths)
					retval+="( "+path.getId()+ " : " + path.getVote() + " )";
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
			connectedPaths = graph.DFSAllPaths(start, end);
			
//			if(connectedPaths==null || connectedPaths.isEmpty())
//				connectedPaths = graph.getClosestPaths(start, end);


			
		}
		
		BridgeSegment(AbstractNode node1, AbstractNode node2, boolean dir1, boolean dir2, ScaffoldVector v1, ScaffoldVector v2){
			assert Math.abs(v1.getMagnitute()) < Math.abs(v2.magnitude): "Illegal order of node position!";
			pSegment = new BidirectedEdgePrototype(node1,node2,dir1,dir2);
			startV = v1; endV = v2;
			int d = ScaffoldVector.composition(endV, ScaffoldVector.reverse(startV)).distance((BidirectedNode)node1, (BidirectedNode)node2);
			connectedPaths = graph.DFSAllPaths((BidirectedNode)node1, (BidirectedNode)node2, dir1, dir2, d, false);
			
//			if(connectedPaths==null || connectedPaths.isEmpty())
//				connectedPaths = graph.getClosestPaths((BidirectedNode)node1, dir1, (BidirectedNode)node2, dir2, d, false);
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
		
		public BridgeSegment reverse(ScaffoldVector brgVector) {
			BridgeSegment retval = new BridgeSegment();
			retval.startV=ScaffoldVector.composition(endV, ScaffoldVector.reverse(brgVector));
			retval.endV=ScaffoldVector.composition(startV, ScaffoldVector.reverse(brgVector));
			retval.connectedPaths=new ArrayList<>();
			for(BidirectedPath p:connectedPaths)
				retval.connectedPaths.add(p.reverse());
			retval.pSegment=pSegment.reverse();
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
