package org.rtassembly.npgraph;

import java.io.File;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.TreeSet;

import org.graphstream.graph.Node;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import japsa.seq.Sequence;

//A bridge structure that the first node must be unique
public class GoInBetweenBridge {
    private static final Logger LOG = LoggerFactory.getLogger(GoInBetweenBridge.class);

	BDGraph graph;
	PopBin bin;
	BDEdgePrototype pBridge; //note: the ending nodes must be unique, or else omitted
	ArrayList<BridgeSegment> segments;
	BridgeSteps steps;
	
	public GoInBetweenBridge(BDGraph graph, AlignedRead bb, PopBin bin) {
		this.graph=graph;
		this.bin=bin;
		steps=new BridgeSteps(bb);
	}

	public int getNumberOfAnchors() {
		int retval=0;
		if(pBridge!=null) {
			if(pBridge.getNode0()!=null) retval++;
			if(pBridge.getNode1()!=null) retval++;
		}
		return retval;
	}
	
	public BDNodeVecState getLastExtendedTip() {
		BDNodeVecState retval= null;
		if(steps!=null) {
			retval=steps.start;
		}
		if(segments!=null && !segments.isEmpty())
			retval=segments.get(segments.size()-1).endNV;
		return retval;
	}
	/*
	 * How complete this bridge is built. Level of completion is as follow:
	 * 0: nothing, 1: one anchor, 2: two anchors determined, 3: bridge connected, 4: bridge completed
	 */
	
	public int getCompletionLevel() {
		int retval=getNumberOfAnchors();
		if(retval==2 && segments!=null && !segments.isEmpty()) {
			boolean isComplete=true, isConnected=true;
			for(BridgeSegment seg:segments) {
				if(!seg.isConnected())
					isConnected=false;;
				if(seg.getNumberOfPaths()!=1)
					isComplete=false;
			}
			if(isConnected)
				retval=3;
			if(isComplete)
				retval=4;
		}

		return retval;
	}

							
	//Merge 2 bridge (must share at least one same unique end-point) together
	//Return true if merging make the bridge reaching new anchor
	byte merge(GoInBetweenBridge qBridge, boolean toConnect) {
		byte retval=0b00; //unchanged state

		if(qBridge==null || qBridge.getCompletionLevel()==0 || getCompletionLevel()==4)
			return retval;
		
		BDNodeState 	sNode0=pBridge.n0,
						sNode1=pBridge.n1,
						qNode0=qBridge.pBridge.n0,
						qNode1=qBridge.pBridge.n1;

		if(HybridAssembler.VERBOSE)
			LOG.info("Merging Bridge {} lvl={} \n{}\nwith Bridge {} lvl={}\n{}\n", getEndingsID(), getCompletionLevel(), getAllNodeVector(), qBridge.getEndingsID(), qBridge.getCompletionLevel(), qBridge.getAllNodeVector());
//		System.out.println("sNode0=" + (sNode0==null?"null":sNode0.getId()));
//		System.out.println("sNode1=" + (sNode1==null?"null":sNode1.getId()));
//		System.out.println("qNode0=" + (qNode0==null?"null":qNode0.getId()));
//		System.out.println("qNode1=" + (qNode1==null?"null":qNode1.getId()));

		boolean found=true;
		if(!sNode0.equals(qNode0)) {			//if the starting points don't agree
			if(sNode0.equals(qNode1)) {			//the first unique contig actually located at the end of alignments list
				qBridge.reverse();
			}else if(sNode1!=null){				//if not: this pBridge must have 2 anchors, the first anchor is not in the alignments list
				reverse();						//flip the first and second anchors: now the second anchor is used as shared maker between subject and query bridge
				if(sNode1.equals(qNode1))	
					qBridge.reverse();
				else if(!sNode1.equals(qNode0))
					found=false;
			}else { 
				found=false;
			}
		}
		if(!found) {
//			System.err.println("Not found common end point between merging bridges");
			return retval;
		}
		
		BridgeSteps qSteps=qBridge.steps;
		
		Iterator<BDNodeVecState> iterator = qSteps.nodes.iterator();		
		BDNodeVecState 	current=null;
		int lastIdx=0, numberOfAnchorsBefore=getNumberOfAnchors();
		Set<BridgeSegment> changedSegments = new HashSet<>();
		
		if(getCompletionLevel() < 3){
			BDNodeVecState.NWAlignment(steps.nodes, qSteps.nodes);
			scanForAnEnd(false);
		}else{
			while(iterator.hasNext()){
				current=iterator.next();
				if(segments!=null){//&& !steps.nodes.contains(current)???
					int prev=-1, cur;
					for(int i=lastIdx; i<segments.size();i++) {
						BridgeSegment curSeg=segments.get(i);
						cur=curSeg.locateAndVote(current);
						//dont go pass the estimated coordinates
						if(cur<prev)
							break;
						else {
							prev=cur;
							lastIdx=i;
							if(cur>0)
								changedSegments.add(curSeg);
						}
					}
					
				}
				
			}	
		}
		for(BridgeSegment sg:changedSegments)
			if(sg.removeUnlikelyPaths())
				retval=0b01;//code representing number of paths of at least 1 segment reduced to 1
		
		if(toConnect && getCompletionLevel()<3)
			steps.connectBridgeSteps(false);	
		
		if(getNumberOfAnchors()>numberOfAnchorsBefore)
			retval=(byte)(retval+0b10);
		
		return retval;
	}

	//return true if there is change in the number of anchor from the updated bridge

	byte merge(AlignedRead read, boolean toConnect) {
		byte retval=0b00; //unchanged state
		if(read==null || read.getAlignmentRecords().size() < 2 || getCompletionLevel()==4)
			return retval;
		
		BDNodeState 	sNode0=pBridge.n0,
						sNode1=pBridge.n1,
						qNode0=new BDNodeState(read.getFirstAlignment().node, read.getFirstAlignment().strand),
						qNode1=new BDNodeState(read.getLastAlignment().node, !read.getLastAlignment().strand);

		if(HybridAssembler.VERBOSE)
			LOG.info("Merging Bridge {} lvl={} \n{}\nwith AlignedRead {}\n{}\n", getEndingsID(), getCompletionLevel(), getAllNodeVector(), read.getEndingsID(), read.getCompsString());
//		System.out.println("sNode0=" + (sNode0==null?"null":sNode0.getId()));
//		System.out.println("sNode1=" + (sNode1==null?"null":sNode1.getId()));
//		System.out.println("qNode0=" + (qNode0==null?"null":qNode0.getId()));
//		System.out.println("qNode1=" + (qNode1==null?"null":qNode1.getId()));

		boolean found=true;
		if(!sNode0.equals(qNode0)) {	//if the starting points don't agree
			if(sNode0.equals(qNode1)) {	//the first unique contig actually located at the end of alignments list
				read.reverse();
			}else if(sNode1!=null){		//if not: this pBridge must have 2 anchors, the first anchor is not in the alignments list
				reverse();				//flip the first and second anchors: now the second anchor is used as shared maker between subject and query bridge
				if(sNode1.equals(qNode1))	
					read.reverse();
				else if(!sNode1.equals(qNode0))
					found=false;
			}else { 
				found=false;
			}
		}
		if(!found) {
//			System.err.println("Not found common end point to merge");
			return retval;
		}
		
		Alignment start=read.getFirstAlignment();
		BDNodeVecState 	current=null;
		int lastIdx=0, numOfAnchorsBefore=getNumberOfAnchors();
		Set<BridgeSegment> changedSegments = new HashSet<>();
		
		if(getCompletionLevel()<3){
			BDNodeVecState.NWAlignment(steps.nodes, new BridgeSteps(read).nodes);	
			scanForAnEnd(false);
		}else{
			for(Alignment alg:read.getAlignmentRecords()) {
				current=new BDNodeVecState(start, alg, read.getVector(start,alg));
				if(segments!=null){ //&& !steps.nodes.contains(current)???
					int prev=-1, cur;
					for(int i=lastIdx; i<segments.size();i++) {
						BridgeSegment curSeg=segments.get(i);
						cur=curSeg.locateAndVote(current);
						//dont go pass the estimated coordinates
						if(cur<prev)
							break;
						else {
							prev=cur;
							lastIdx=i;
							if(cur>0)
								changedSegments.add(curSeg);
						}
					}
					
				}
			}
		}
		
		for(BridgeSegment sg:changedSegments)
			if(sg.removeUnlikelyPaths())
				retval=0b01;//code representing number of paths of at least 1 segment reduced to 1
		
		if(toConnect && getCompletionLevel()<3)
			steps.connectBridgeSteps(false);		
				
		
		if(HybridAssembler.VERBOSE){
			LOG.info("After merging: {}\nstart={}; end={}; completion level={}", pBridge.toString(), (steps.start==null?"null":steps.start.toString()), (steps.end==null?"null":steps.end.toString()), getCompletionLevel());
			LOG.info(getAllPossiblePaths());
		}

		if(getNumberOfAnchors()>numOfAnchorsBefore)
			retval=(byte)(retval+0b10);
				
		return retval;
	
	}
	
	private void addSegment(BridgeSegment seg) {
		if(segments==null)
			segments=new ArrayList<>();
		segments.add(seg);
	}
	
	private void reverse() {
		assert getNumberOfAnchors()==2:"Could only reverse a determined bridge (2 anchors)";
		if(HybridAssembler.VERBOSE)
			LOG.info("Reversing the bridge {} (start={} end={}):\n{}\n", getEndingsID(), steps.start.toString(), steps.end.toString(), getAllNodeVector());

		pBridge=pBridge.reverse();
		int direction=pBridge.getDir0()?1:-1;
		//reverse the segments
		if(segments!=null){
			ArrayList<BridgeSegment> tmp = segments;
			segments = new ArrayList<>();
			if(!tmp.isEmpty()) {
				ScaffoldVector brgVector=tmp.get(tmp.size()-1).getEndVector();
				for(int i=tmp.size()-1; i>=0 ;i--)
					addSegment(tmp.get(i).reverse(steps.end.getNode(), brgVector));
			}
		}
		//reverse the nodes list
		TreeSet<BDNodeVecState> reversedSet = new TreeSet<BDNodeVecState>();
		ScaffoldVector rev=ScaffoldVector.reverse(steps.end.getVector());//anchors number = 2 so there exist end node
		BDNodeVecState tmp = null;
		
		//need to do this to re-sort the changed elements
		for(BDNodeVecState nv:steps.subSet(steps.start, steps.end)) {
			nv.setVector(ScaffoldVector.composition(nv.getVector(), rev));
			nv.setRoot(steps.end.getNode());
			if(nv.getVector().isIdentity() || (nv.getVector().getMagnitute()-BDGraph.A_TOL)*direction>0)
				reversedSet.add(nv);

		}
		//re-assign start and end
		tmp=steps.start;
		steps.start=steps.end;
		steps.end=tmp;
		
		steps.nodes=reversedSet;
		
		if(HybridAssembler.VERBOSE)
			LOG.info("Reversed bridge = {} (start={} end={}):\n{}\n", getEndingsID(), steps.start.toString(), steps.end.toString(), getAllNodeVector());

	}
	
	public String getEndingsID() {
		if(pBridge==null)
			return "-,-";
		else
			return pBridge.toString();
	}
	

	public String getAllPossiblePaths() {
		if(segments==null)
			return "{none}";
		
		String retval = "{\n";
		for(BridgeSegment seg:segments) {
			retval+=seg.pSegment.toString() + ": ";
			if(seg.getNumberOfPaths()==0)
				retval += "()";
			else
				for(BDPath path:seg.connectedPaths)
					retval+="( "+path.getId()+  ": deviation=" + path.getDeviation() + " score=" + path.getPathEstats() + " )";
			retval+="\n";
		}
		retval+="}";
		return retval;
	}
	
	public String getAllNodeVector(){
		if(steps==null)
			return "empty list of steps";
		return steps.toString();
	}

	//Get the most probable path out of the candidates based on deviation and likelihood
	public BDPath getBestPath(Node startFrom, Node endAt) { //the markers must be (transformed) unique
		if(HybridAssembler.VERBOSE)
			LOG.info("Finding best path from " + startFrom.getId() + " to " + endAt.getId() + " among: \n" + getAllPossiblePaths());
		BDPath retval=null;

		if(segments==null || segments.isEmpty())
			return null;
		//1. Find all paths with the same deviation score
		ArrayList<BDPath> 	candidates = new ArrayList<>(),
							tmpList = new ArrayList<>();
		boolean found=false;
		int tolerance=0, keep_max=1; //increase these number for better accuracy (?) but much slower
		
		for(BridgeSegment seg:segments) {
			//finding the start node
			if(!found){
				if(seg.pSegment.getNode0()!=startFrom){
					continue;
				}
				else
					found=true;
			}

			
			if(seg.getNumberOfPaths()>0){			
				int bestDeviation=Math.abs(seg.connectedPaths.get(0).getDeviation());
				seg.connectedPaths.removeIf(p->(Math.abs(p.getDeviation())>bestDeviation+tolerance));
				tmpList = new ArrayList<>();
				if(candidates.isEmpty()){
					candidates.addAll(seg.connectedPaths);
					continue;
				}
				
				for(BDPath p1:seg.connectedPaths)
					for(BDPath p0:candidates){
							BDPath p01=p0.join(p1);
							if(p01==null)
								return null;
							else
								tmpList.add(p01);
						}
				
				tmpList.sort(Comparator.comparing(BDPath::getPathEstats));
				candidates=new ArrayList<>(tmpList.subList(0, keep_max>tmpList.size()?tmpList.size():keep_max));
				
			}else
				return null;
			
			//escape if end node is reached
			if(seg.pSegment.getNode1()==endAt)
				break;
		}
		
		//2. Pick the best one with highest likelihood score
		candidates.sort(Comparator.comparing(BDPath::getPathEstats));
		retval=candidates.get(0);
		
		if(retval!=null) {
			retval.setConsensusUniqueBinOfPath(bin);
			if(HybridAssembler.VERBOSE)
				LOG.info("...best path found: " + retval.getId());

		}
		return retval;
	}
	
	public int countPathsBetween(Node startFrom, Node endAt){ //the markers must be (transformed) unique
		int retval=0;
		boolean found=false;
		for(BridgeSegment seg:segments) {
			//finding the start node
			if(!found){
				if(seg.pSegment.getNode0()!=startFrom)
					continue;
				else
					found=true;
			}
			if(!seg.isConnected())
				return 0;
			if(retval!=0)
				retval*=seg.connectedPaths.size();
			else
				retval=seg.connectedPaths.size();
			//escape if end node is reached
			if(seg.pSegment.getNode1()==endAt)
				break;
		}
		return retval;
	}
	
	//combination of getBestPath() + countPathsBetween() + path.chopAtAnchors()
	//return new unique path to reduce in a building bridge
	public List<BDPath> scanForNewUniquePaths(){
		List<BDPath> retval = new ArrayList<>();
		if(segments==null || segments.isEmpty())
			return retval;
		if(HybridAssembler.VERBOSE)
			LOG.info("Scanning on bridge with segments:\n" + getAllPossiblePaths());		
		BDPath curPath=null;
		PopBin sbin=null;
		for(BridgeSegment seg:segments){
			if(seg.isConnected() && seg.connectedPaths.size()==1){
				sbin=SimpleBinner.getBinIfUniqueNow(seg.startNV.getNode());
				if(PopBin.isCloseBins(sbin,bin)){
					if(curPath!=null){ 
						graph.getNewSubPathsToReduce(curPath).stream().forEach(p->retval.add(p)); //in case there're anchors within segment
					}
					
					curPath=seg.connectedPaths.get(0);
					curPath.setConsensusUniqueBinOfPath(bin);

				}else if(curPath!=null)
					curPath=curPath.join(seg.connectedPaths.get(0)); //don't need to calculate likelihood of unique path 

			}else{
				curPath=null;
			}
		}
		
		if(curPath!=null && curPath.getEdgeCount()>0) 
			graph.getNewSubPathsToReduce(curPath).stream().forEach(p->retval.add(p)); //in case there're anchors within segment
		
		
		return retval;
	}
	
	//Try to look for an ending unique node of a unidentifiable bridge
	//by using isUniqueNow()
	public boolean scanForAnEnd(boolean greedy){
		if(steps==null)
			return false;
		Iterator<BDNodeVecState> ite = steps.nodes.descendingIterator();//reversed order
		while(ite.hasNext()){
			BDNodeVecState tmp=ite.next();
			if(tmp==getLastExtendedTip()) //there is no (new) end detected
				return false;
			
			PopBin b=SimpleBinner.getBinIfUniqueNow(tmp.dest);
			if(PopBin.isCloseBins(b,bin)){
				if(steps.end==null || greedy || steps.end.nvsScore < tmp.nvsScore){
					if(steps.end!=tmp){
						steps.end=tmp;
						if(HybridAssembler.VERBOSE)
							LOG.info("FOUND NEW END: " + steps.end.getNode().getId());
						return true;
					}
				}
			}
				
		}
		return false;
		
	}

	
	//TODO: check if it's enough to run MSA for long reads consensus
	//need estimation of gap (depend on level of completion) versus number of spanning reads
	public boolean checkMSACoverage(){
		return false;
	}

	/************************************************************************************************
	 ************************************************************************************************
	 * Class to represent a single segment of the whole bridge.
	 * A bridge consists of >=1 segments.
	 ************************************************************************************************/
	class BridgeSegment{
//		ArrayList<Sequence> nnpReads; // to store nanopore data if needed
		ArrayList<BDPath> connectedPaths=null;
		BDEdgePrototype pSegment;
		BDNodeVecState startNV, endNV;
		BridgeSegment(){}
		
		BridgeSegment(BDNodeVecState nv1, BDNodeVecState nv2, boolean dir1, boolean dir2, boolean greedy){
			assert Math.abs(nv1.getVector().getMagnitute()) < Math.abs(nv2.getVector().magnitude): "Illegal order of node position!";
			BDNode srcNode=nv1.getNode(), dstNode=nv2.getNode();
			pSegment = new BDEdgePrototype(srcNode,dstNode,dir1,dir2);
			startNV = nv1; endNV = nv2;
			int d = ScaffoldVector.composition(nv2.getVector(), ScaffoldVector.reverse(nv1.getVector())).distance(srcNode, dstNode);
			if(d<=(greedy?10:1)*BDGraph.D_LIMIT) //greedy search will tolerate 10 times longer gap
	    		connectedPaths = graph.DFSAllPaths(srcNode, dstNode, dir1, dir2, d);

			//call consensus when time come!
			if(	(connectedPaths==null || connectedPaths.isEmpty())){				
				//TODO: more anchors connecting!!! Here just connect unique dead-end unique nodes
				if(	PopBin.isCloseBins(SimpleBinner.getBinIfUnique(srcNode), SimpleBinner.getBinIfUnique(dstNode)) 
					&& !graph.isConflictBridge(pSegment))
				{					
					connectedPaths=new ArrayList<>();
					BDPath path = new BDPath(srcNode);
					String id = BDEdge.createID(srcNode, dstNode, dir1, dir2);				
					try{
						//Option 1: hide long-read consensus from the graph
						BDNode n=new BDNode(graph, "000"+AlignedRead.PSEUDO_ID++);
						Sequence seq=graph.consensus.getConsensus(id, greedy);
						if(seq==null)
							return;
						
						n.setAttribute("seq", seq);
						n.setAttribute("len", seq.length());
						n.setAttribute("cov",SimpleBinner.getBinIfUnique(srcNode).estCov);
						boolean disagreement=srcNode.getId().compareTo(dstNode.getId()) > 0 
								|| (srcNode==dstNode && dir1==false && dir2==true);
						BDEdge 	e0=new BDEdge(srcNode, n, dir1, disagreement),
								e1=new BDEdge(n,dstNode,!disagreement,dir2);	
						BDPath p = new BDPath(srcNode);					
						p.add(e0);
						p.add(e1);
						BDEdge pseudoEdge=graph.addEdge(srcNode, dstNode, dir1, dir2);
						pseudoEdge.setAttribute("path", p);
						pseudoEdge.setAttribute("consensus");
						path.add(pseudoEdge);

//							//Option 2: or using this to create&display a "pseudo node"
//							BDNode n=(BDNode) graph.addNode("000"+AlignedRead.PSEUDO_ID++);
//							Sequence seq=GraphUtil.consensusSequence(AlignedRead.tmpFolder+File.separator+id+".fasta", distance, id, "poa");
//							if(seq==null)
//								return null;
//							n.setAttribute("seq", seq);
//							n.setAttribute("len", seq.length());
//							n.setAttribute("cov",SimpleBinner.getBinIfUnique(srcNode).estCov);
//							n.setGUI("red", "diamond");
//							n.setAttribute("unique", SimpleBinner.getBinIfUnique(srcNode));
//							boolean disagreement=srcNode.getId().compareTo(dstNode.getId()) > 0 
//									|| (srcNode==dstNode && srcDir==false && dstDir==true);
//							Edge 	e0=addEdge(srcNode, n, srcDir, disagreement),
//									e1=addEdge(n,dstNode,!disagreement,dstDir);	
	//
//							path.add(e0);
//							path.add(e1);
						
						connectedPaths.add(path);
					}catch(Exception e){
						System.err.println("Failed to make consensus sequence for " + id +"!");
						System.err.println("Reason: "+ e.getMessage());
					}	
	    		}

				
			}
		}
		
//		BridgeSegment(BDPath path){
//			try {
//				pSegment=new BDEdgePrototype(path);
//				connectedPaths=new ArrayList<>();
//				connectedPaths.add(path);
//				
//				//only if path cover the whole bridge (1-segment bridge)
//				int dist=(int) (path.getLength()
//						-(pSegment.getDir0()?0:pSegment.getNode0().getNumber("len"))
//						-(pSegment.getDir1()?0:pSegment.getNode1().getNumber("len")));
//				
//				startNV=new BDNodeVecState(pSegment.getNode0(),pSegment.getNode0(),new ScaffoldVector());
//				endNV=new BDNodeVecState(pSegment.getNode1(), new ScaffoldVector(pSegment.getDir0()?dist:-dist, pSegment.getDir0()!=pSegment.getDir1()?1:-1));
//			} catch (Exception e) {
//				System.err.println("Cannot make bridge from this path!");
//				e.printStackTrace();
//			} 
//
//		}
		
		BridgeSegment reverse(BDNode newRoot, ScaffoldVector brgVector) {
			BridgeSegment retval = new BridgeSegment();
			retval.pSegment=pSegment.reverse();
			retval.startNV=new BDNodeVecState(newRoot, endNV.getNode(), endNV.getScore(), ScaffoldVector.composition(getEndVector(), ScaffoldVector.reverse(brgVector)));
			retval.endNV=new BDNodeVecState(newRoot, startNV.getNode(), startNV.getScore(),  ScaffoldVector.composition(getStartVector(), ScaffoldVector.reverse(brgVector)));

			retval.connectedPaths=new ArrayList<>();
			for(BDPath p:connectedPaths)
				retval.connectedPaths.add(p.reverse());
			
			return retval;
		}
		
		/*
		 * Find a node vector if it appear in one of the paths, only then update the scores
		 * Return retval=-1 if not located within the range of this segment
		 * Return retval>0 if found in x paths
		 * Return retval=0 if located in the interval of this segment but not found in any path 
		 */
		int locateAndVote(BDNodeVecState nv) {
			int retval=-1;
			if(!isConnected())
				return retval;
			HashMap<BDPath,Integer> scores = new HashMap<>();
			if(nv.compareTo(startNV)*nv.compareTo(endNV)<=0) {
				retval=0;
				ScaffoldVector start2nv=ScaffoldVector.composition(nv.getVector(), ScaffoldVector.reverse(getStartVector()));
				int d=start2nv.distance((BDNode) pSegment.getNode0(), nv.getNode());
				for(BDPath p:connectedPaths) {
					//update the divergence no matter what it's found or not. 
					int bd = p.getClosestDistance(pSegment.getNode0(), nv.getNode(), start2nv.direction>0, d);
					if(bd<Integer.MAX_VALUE){
						retval++;
					}else
						bd=Math.max(BDGraph.A_TOL, (int)(BDGraph.R_TOL*d));
					scores.put(p, bd);

				}
			}

			if(retval>0){
				//apply score changes and eliminate paths with too much deviations
				for(BDPath p:scores.keySet())
					p.updatePathDeviation(scores.get(p));
			}
		
			return retval;
		}
		//return true iff there is only one candidate left <=> result found!
		public boolean removeUnlikelyPaths(){
			if(connectedPaths==null || connectedPaths.isEmpty())
				return false;
			connectedPaths.sort((a,b)->Integer.compare(Math.abs(a.getDeviation()), Math.abs(b.getDeviation())));
			int bestDiff = 	connectedPaths.get(0).getDeviation();
			connectedPaths.removeIf(p->(Math.abs(p.getDeviation()) > Math.abs(bestDiff)+BDGraph.A_TOL));
			if(connectedPaths.size()==1)
				return true;
			else 
				return false;
		}
		
		int getNumberOfPaths(){
			if(connectedPaths==null)
				return 0;
			else 
				return connectedPaths.size();
		}
		boolean isConnected(){
			return getNumberOfPaths()>=1;
		}
		boolean isUnique(){
			return getNumberOfPaths()==1;
		}
		ScaffoldVector getEndVector() {
			return endNV.getVector();
		}
		ScaffoldVector getStartVector() {
			return startNV.getVector();
		}
		String getId() {
			return pSegment.toString();
		}
	}
	
	/*
	 * Class representing list of found nodes and their position (depicted as vectors)
	 * in-between two ends of a bridge. 
	 * Used for determine segments 
	 */
	class BridgeSteps{
		TreeSet<BDNodeVecState> nodes;
		//save the two unique nodes at 2 endings of the bridge. When end is covered more than threshold, pBridge is set
		BDNodeVecState start, end; 
		
		BridgeSteps(){
			nodes=new TreeSet<BDNodeVecState>();
		}
		
		BridgeSteps(AlignedRead read){
			this();
			if(read.getAlignmentRecords().size()<2)
				return;
			if(SimpleBinner.getBinIfUnique(read.getFirstAlignment().node)==null)
				if(SimpleBinner.getBinIfUnique(read.getLastAlignment().node)!=null)
					read.reverse();
				else
					return;
					
					
			Alignment 	firstAlg = read.getFirstAlignment(),
						curAlg = null;				
			//First alignment must be from an unique node to use as start point of the bridge
			start=new BDNodeVecState(firstAlg.node, firstAlg.node, firstAlg.quality, new ScaffoldVector());
			nodes.add(start);
				
			for(int i=1;i<read.getAlignmentRecords().size();i++) {
				curAlg = read.getAlignmentRecords().get(i);
				BDNodeVecState tmp=new BDNodeVecState(firstAlg, curAlg, read.getVector(firstAlg, curAlg));
				nodes.add(tmp);
				
				if(SimpleBinner.getBinIfUnique(curAlg.node)!=null)
					end=tmp;

			}
			//update the pBridge accordingly
			if(end!=null && end.qc()) {
				pBridge=new BDEdgePrototype(firstAlg.node, end.getNode(), firstAlg.strand, end.getDirection(firstAlg.strand));//getDirection doesn't work with self-vector
			}else if(pBridge==null)
				pBridge=new BDEdgePrototype(firstAlg.node, firstAlg.strand);
		}
		
		//////////////////////////////////////////////////////////////////////////////////////
		/*
		 * To connect the bridge based on the steps identified inbetween
		 */
		//////////////////////////////////////////////////////////////////////////////////////
		
		//need this because BDNodeVecState.compareTo() not fit for this purpose
		private TreeSet<BDNodeVecState> subSet(BDNodeVecState left, BDNodeVecState right){
			TreeSet<BDNodeVecState> retval = new TreeSet<>();
			boolean found=false;
			for(BDNodeVecState nvs:nodes){
				if(found){
					if(nvs==right)
						break;
					retval.add(nvs);
				}else if(nvs==left)
					found=true;
				
			}	
			return retval;
		}
		private PriorityQueue<BDNodeVecState> getInBetweenSteps(BDNodeVecState left, BDNodeVecState right){
			PriorityQueue<BDNodeVecState> retval = new PriorityQueue<>(nodes.size(), Comparator.comparing(BDNodeVecState::getScore).reversed());
			retval.addAll(subSet(left, right));
			return retval;
		}
		

		private ArrayList<BridgeSegment> greedyConnect(HashMap<String, ArrayList<BridgeSegment>> memory, BDNodeVecState left, BDNodeVecState right, boolean greedy){
			if(HybridAssembler.VERBOSE)
				LOG.info("\nConnecting " + left + " to " + right + "...");
			if(left==right)
				return null;
			String pairKey=left.toString()+right.toString();
			ArrayList<BridgeSegment> retval;

			if(memory.containsKey(pairKey)){
				if(memory.get(pairKey)==null)
					retval=null;
				else
					retval=new ArrayList<>(memory.get(pairKey));
				
				if(HybridAssembler.VERBOSE)
					LOG.info("\nConnecting " + left + " to " + right + ": " + (retval==null?"unreachable!":retval.size()+"-reachable!"));
				return retval;
			}
			
			PriorityQueue<BDNodeVecState> inBetween=getInBetweenSteps(left, right);
			
			ArrayList<BridgeSegment> lBridge, rBridge;
			while(!inBetween.isEmpty()){
				BDNodeVecState mid=inBetween.poll();
				if(HybridAssembler.VERBOSE)
					LOG.info("..cut at mid="+mid.toString());
				
				lBridge = greedyConnect(memory, left, mid, greedy);
				if(lBridge==null||lBridge.isEmpty()){
					continue;
				}
				else{
					retval=new ArrayList<>(lBridge);
				}
				
				rBridge = greedyConnect(memory, mid, right, greedy);
				if(rBridge==null||rBridge.isEmpty()){
					continue;
				}
				else
					retval.addAll(rBridge);
				
				memory.put(pairKey, retval);
				return retval;
			}
			
			//If there are no intermediate step left in between 
			BridgeSegment seg=new BridgeSegment( 
									left, 
									right, 
									left.getVector().isIdentity()?pBridge.getDir0():!left.getDirection(pBridge.getDir0()), 
									right.getDirection(pBridge.getDir0()), 
									greedy);
			if(seg.isConnected()){
				retval=new ArrayList<>();
				retval.add(seg);
			}
			else{
				retval=null;	
			}
			
			memory.put(pairKey, retval);	
			return retval; 
						
		}
		//Try to make the continuous segments (those that have paths connected 2 ends). Return true if it possible
		//starting from a marker (unique or transformed unique node)
		boolean connectBridgeSteps(boolean force){
			if(isIdentifiable()){
				if(!connectable() && !force){
					if(HybridAssembler.VERBOSE)
						LOG.info("Bridge {} is not qualified to connect yet: start={} end={}", pBridge.toString(), (start==null?"null":start.toString()), (end==null?"null":end.toString()));
					return false;
				}
			}else{
				if(HybridAssembler.VERBOSE)
					LOG.info("Bridge {} is not complete to connect!\n", pBridge.toString());
				return false;
			}	
			
			BDNodeVecState startFrom=getLastExtendedTip(),
							endAt=end;

			
			if(HybridAssembler.VERBOSE)
				LOG.info("Trying to connect bridge {} from {} to {}:\n{}\n", pBridge.toString(), startFrom.getNode().getId(), endAt.getNode().getId(), getAllNodeVector());
				
			HashMap<String, ArrayList<BridgeSegment>> memory = new HashMap<>();
			ArrayList<BridgeSegment> segs=greedyConnect(memory, startFrom, endAt, force);
			
			if(segs==null||segs.isEmpty()){
				segments=null;
				if(HybridAssembler.VERBOSE)
					LOG.info("Failed to connect " + pBridge.toString());
				return false;
			}else{
				segs.stream().forEach(segment->addSegment(segment));
				//set pBridge end iff endAt node is original unique node
				if(SimpleBinner.getBinIfUnique(endAt.dest)!=null){
					pBridge.n1=new BDNodeState(endAt.dest, endAt.getDirection(pBridge.getDir0()));
					if(HybridAssembler.VERBOSE)
						LOG.info("Success to finish " + pBridge.toString());
				} else
					if(HybridAssembler.VERBOSE)
						LOG.info("Success to extend " + pBridge.toString() + " to " + endAt.dest.getId());

				return true;
			}
		}
		boolean isIdentifiable(){
			return start!=null && end!=null;
		}
		
		boolean connectable(){
			if(!isIdentifiable())
				return false;
			
			return start.qc() && end.qc();
		}

		public void setNodes(TreeSet<BDNodeVecState> nodes){
			this.nodes=nodes;
		}
		public void setStart(BDNodeVecState start){
			this.start=start;
		}
		public void setEnd(BDNodeVecState end){
			this.end=end;
		}
		
		public String toString(){
			String retval="";
			for(BDNodeVecState nv:nodes)
				retval+=nv.toString() + "\n";
			return retval;
		}
	}


}
