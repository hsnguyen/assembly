package org.rtassembly.npgraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.SortedSet;
import java.util.TreeSet;

import org.graphstream.graph.implementations.AbstractNode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

//A bridge structure that the first node must be unique
public class GoInBetweenBridge {
	private static final Logger LOG = LoggerFactory.getLogger(GoInBetweenBridge.class);
	public static volatile int MAX_DIFF=3;
	BidirectedGraph graph; //partial order graph saving possible paths
	PopBin bin;
	BidirectedEdgePrototype pBridge; //note: fist node of the bridge is set unique
	ArrayList<BridgeSegment> segments;
	SortedSet<NodeVector> nodes;
	
	GoInBetweenBridge(BidirectedGraph graph, PopBin bin){
		this.graph=graph;
		segments=new ArrayList<>();
		nodes=new TreeSet<NodeVector>();
		this.bin=bin;
	}
	//This is for SPAdes path reader only. The input path must be have at least one end unique.
	//This unique ending node will be use as the first anchor of the bridge
	GoInBetweenBridge (BidirectedGraph graph, BidirectedPath path){
		this(graph,path.getConsensusUniqueBinOfPath());
		if(	path.size()>1) {
			if(SimpleBinner.getUniqueBin(path.getRoot())!=null)
				addSegment(new BridgeSegment(path),true);
			else if(SimpleBinner.getUniqueBin(path.peekNode())!=null)
				addSegment(new BridgeSegment(path.reverse()),true);
			
		}

	}
	
	public GoInBetweenBridge(BidirectedGraph graph, AlignedRead bb, PopBin b) {
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
		if(segments==null||segments.isEmpty())
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
		//building on existed one
		if(pBridge!=null){
			if(pBridge.getNode0()!=alignedRead.getFirstAlignment().node) {		//if the starting points don't agree
				if(pBridge.getNode0()==alignedRead.getLastAlignment().node) {	//the first unique contig actually located at the end of alignments list
					alignedRead=alignedRead.reverse();
				}else if(getNumberOfAnchors()==2){								//if not: this pBridge must have 2 anchors, the first anchor is not in the alignments list
					this.reverse();												//flip the first and second anchors: now the second anchor is used as shared maker between pBridge and alignedRead
					if(pBridge.getNode1()==alignedRead.getLastAlignment().node)	
						alignedRead=alignedRead.reverse();
				}else { 
					System.err.println("Error bridge prototype: " + pBridge.toString());
					return false;
				}
			}
		}
		 System.out.println("Applying aligned read on the bridge: " + alignedRead.getCompsString());
		 
		if(addAligments(alignedRead)){
			return connectBridge();
		}else
			return false;

	}
							

	// Return true if it reachs an unique nodes (good bridge)
	private boolean addAligments(AlignedRead read) {
		Alignment 	firstAlg = read.getFirstAlignment(),
					curAlg = null;
		boolean retval=false;

		nodes.add(new NodeVector(firstAlg.node, new ScaffoldVector(), firstAlg.strand));
		
		for(int i=1;i<read.getAlignmentRecords().size();i++) {
			curAlg = read.getAlignmentRecords().get(i);
			if(SimpleBinner.getUniqueBin(curAlg.node)!=null){
				retval=true;
			}
			NodeVector currentNV=new NodeVector(curAlg.node, read.getVector(firstAlg, curAlg),firstAlg.strand);
			if(nodes.contains(currentNV)){
//				System.out.println("Found node vector in this list: " + currentNV.toString());
				Iterator<NodeVector> ite = nodes.iterator();
				while(ite.hasNext()){
					NodeVector tmp=ite.next();
					if(tmp.equals(currentNV)){
						tmp.score++;
						break;
					}
				}
			}else
				nodes.add(currentNV);
		}
		if(getNumberOfAnchors()<2) {
			pBridge=new BidirectedEdgePrototype(firstAlg.node, nodes.last().getNode(), firstAlg.strand, nodes.last().getDirection());//getDirection doesn't work with self-vector
		}
		return retval;
	}
	//Try to make the continuous segments (those that have paths connected 2 ends). Return true if it possible
	private boolean connectBridge(){
		if(pBridge==null){
			System.err.println("pBridge must be set before using this function");
			return false;
		}
		if(nodes.size()<2 || SimpleBinner.getUniqueBin(nodes.last().getNode())==null || !nodes.last().qc())
			return false;
		System.out.println("Trying to connect bridge " + pBridge.toString() + ":\n" + getAllNodeVector());
		//First build shortest tree from the end
		HashMap<String,Integer> shortestMap = 
				graph.getShortestTreeFromNode(	(BidirectedNode)pBridge.getNode1(), 
												pBridge.getDir1(), 
												nodes.last().getVector().distance((BidirectedNode)pBridge.getNode0(), (BidirectedNode)pBridge.getNode1()));
		if(!shortestMap.containsKey(pBridge.getNode0().getId()+(pBridge.getDir0()?"i":"o"))){ //note the trick: direction BEFORE the path started
			System.err.println("Shortest tree couldn't reach to the other end!");
			return false;
		}
		Iterator<NodeVector> iterator = nodes.iterator();
		
		NodeVector 	prev=null, 
					current=null;
		boolean retval = false;
		while(iterator.hasNext()){
			if(current==null){
				prev=current=iterator.next();//first unique node, don't need to check quality
				continue;
			}
			
			current=iterator.next();
			//need a quality-checking here before including into a segment step
			if(!current.qc() || shortestMap.containsKey(current.getNode().getId())){
				continue;
			}else{
				BridgeSegment seg = new BridgeSegment(	prev.getNode(), current.getNode(), 
														prev.getDirection(), current.getDirection(), 
														prev.getVector(), current.getVector());
				System.out.print("...connecting " + prev.toString() + " to " + current.toString());
				if(seg.isConnected()){
					System.out.println(" :success!");
					addSegment(seg,false);
					//FIXME: direction change here!!!
					prev=current;
					if(current.equals(nodes.last()))
						retval=true;						
				}else
					System.out.println(" :fail!");

			}
				
			
		}
		
		return retval;
	}
	//adding new segment with/without update pBridge
	private void addSegment(BridgeSegment seg, boolean update) {
		
		if(pBridge==null) {
			segments.add(seg);
			if(update)
				pBridge=new BidirectedEdgePrototype(seg.pSegment.getNode0(), seg.pSegment.getNode1(), seg.pSegment.getDir0(), seg.pSegment.getDir1());
		}
		else if(pBridge.getNode1()==seg.pSegment.getNode0() && pBridge.getDir1()!=seg.pSegment.getDir0()){
			segments.add(seg);
			if(update)
				pBridge.n1=seg.pSegment.n1;
		}
	}

	private void reverse() {
		assert getNumberOfAnchors()==2:"Could only reverse a determined bridge (2 anchors)";
		ArrayList<BridgeSegment> tmp = segments;
		segments = new ArrayList<>();
		pBridge=new BidirectedEdgePrototype(pBridge.getNode1(), pBridge.getNode0(), pBridge.getDir1(), pBridge.getDir0());
		//reverse the segments
		if(!tmp.isEmpty()) {
			ScaffoldVector brgVector=tmp.get(tmp.size()-1).endV;
			for(int i=tmp.size()-1; i>=0 ;i--)
				addSegment(tmp.get(i).reverse(brgVector),false);
		}
		//reverse the nodes list
		SortedSet<NodeVector> reversedSet = new TreeSet<NodeVector>();
		ScaffoldVector rev=ScaffoldVector.reverse(nodes.last().getVector());//anchors number = 2 so there exist last()
		boolean rootDir=nodes.last().getDirection();
		for(NodeVector nv:nodes) {
			reversedSet.add(new NodeVector(nv.getNode(),ScaffoldVector.composition(nv.getVector(), rev),rootDir));;
		}
		nodes=reversedSet;
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
	
	public String getAllNodeVector(){
		String retval="";
		for(NodeVector nv:nodes)
			retval+=nv.toString() + "\n";
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
