package org.rtassembly.npgraph;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import org.graphstream.graph.implementations.AbstractNode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

//A bridge structure that the first node must be unique
public class GoInBetweenBridge {
	private static final Logger LOG = LoggerFactory.getLogger(GoInBetweenBridge.class);
	
	BidirectedGraph graph; //partial order graph saving possible paths
	PopBin bin;
	BidirectedEdgePrototype pBridge; //note: the ending nodes must be unique, or else omitted
	ArrayList<BridgeSegment> segments;
	BridgeSteps steps;
	
	GoInBetweenBridge(BidirectedGraph graph, PopBin bin){
		this.graph=graph;
//		segments=new ArrayList<>();
//		steps=new BridgeSteps();
		this.bin=bin;
	}
	//This is for SPAdes path reader only. The input path must have at least one end unique.
	//This unique ending node will be use as the first anchor of the bridge
	GoInBetweenBridge (BidirectedGraph graph, BidirectedPath path){
		this(graph,path.getConsensusUniqueBinOfPath());
		if(	path.size()>1) {
			if(SimpleBinner.getUniqueBin(path.getRoot())!=null) {
				addSegment(new BridgeSegment(path));
			}
			else if(SimpleBinner.getUniqueBin(path.peekNode())!=null) {
				addSegment(new BridgeSegment(path.reverse()));
			}
			
			try {
				pBridge=new BidirectedEdgePrototype(path.getFirstNode(), path.getFirstNodeDirection());
			} catch (Exception e) {
				System.err.println("Illegal path to construct pBridge: " + path.getId());
				e.printStackTrace();
			}

		}

	}
	
	public GoInBetweenBridge(BidirectedGraph graph, AlignedRead bb, PopBin b) {
		this(graph,b);
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
	boolean merge(GoInBetweenBridge qBridge) {
	
		if(qBridge==null || qBridge.getCompletionLevel()==0 || getCompletionLevel()==4)
			return false;
		
		NodeDirection 	sNode0=pBridge.n0,
						sNode1=pBridge.n1,
						qNode0=qBridge.pBridge.n0,
						qNode1=qBridge.pBridge.n1;

		System.out.printf("Merging Bridge %s lvl=%d \n%s\nwith Bridge %s lvl=%d\n%s\n", getEndingsID(), getCompletionLevel(), getAllNodeVector(), qBridge.getEndingsID(), qBridge.getCompletionLevel(), qBridge.getAllNodeVector());
//		System.out.println("sNode0=" + (sNode0==null?"null":sNode0.getId()));
//		System.out.println("sNode1=" + (sNode1==null?"null":sNode1.getId()));
//		System.out.println("qNode0=" + (qNode0==null?"null":qNode0.getId()));
//		System.out.println("qNode1=" + (qNode1==null?"null":qNode1.getId()));

		boolean found=true;
		if(!sNode0.equals(qNode0)) {							//if the starting points don't agree
			if(sNode0.equals(qNode1)) {						//the first unique contig actually located at the end of alignments list
				qBridge.reverse();
			}else if(sNode1!=null){						//if not: this pBridge must have 2 anchors, the first anchor is not in the alignments list
				reverse();							//flip the first and second anchors: now the second anchor is used as shared maker between subject and query bridge
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
			return false;
		}
		
		BridgeSteps sSteps=steps, qSteps=qBridge.steps;
		
		Iterator<NodeVector> iterator = qSteps.nodes.iterator();		
		NodeVector 	current=null;
		int lastIdx=0, numberOfAnchorsBefore=getNumberOfAnchors();
		while(iterator.hasNext()){
			current=iterator.next();
			if(getCompletionLevel()==3){
				if(!sSteps.nodes.contains(current)){
					int prev=-1, cur;
					for(int i=lastIdx; i<segments.size();i++) {
						BridgeSegment curSeg=segments.get(i);
						cur=curSeg.locateAndVote(current);
						if(cur<prev)
							break;
						else {
							prev=cur;
							lastIdx=i;
						}
					}
				}
			}else{
				//just adding
				steps.addNode(current);
			}	
			
		}		
		if(getCompletionLevel()<3)
			steps.connectBridgeSteps(false);	
		return getNumberOfAnchors()>numberOfAnchorsBefore;
	}

	//return true if there is change in the number of anchor from the updated bridge

	boolean merge(AlignedRead read) {
		
		if(read==null || read.getAlignmentRecords().size() < 2 || getCompletionLevel()==4)
			return false;
		
		NodeDirection 	sNode0=pBridge.n0,
						sNode1=pBridge.n1,
						qNode0=new NodeDirection(read.getFirstAlignment().node, read.getFirstAlignment().strand),
						qNode1=new NodeDirection(read.getLastAlignment().node, !read.getLastAlignment().strand);

		System.out.printf("Merging Bridge %s lvl=%d \n%s\nwith AlignedRead %s\n%s\n", getEndingsID(), getCompletionLevel(), getAllNodeVector(), read.getEndingsID(), read.getCompsString());
//		System.out.println("sNode0=" + (sNode0==null?"null":sNode0.getId()));
//		System.out.println("sNode1=" + (sNode1==null?"null":sNode1.getId()));
//		System.out.println("qNode0=" + (qNode0==null?"null":qNode0.getId()));
//		System.out.println("qNode1=" + (qNode1==null?"null":qNode1.getId()));

		boolean found=true;
		if(!sNode0.equals(qNode0)) {							//if the starting points don't agree
			if(sNode0.equals(qNode1)) {						//the first unique contig actually located at the end of alignments list
				read.reverse();
			}else if(sNode1!=null){						//if not: this pBridge must have 2 anchors, the first anchor is not in the alignments list
				reverse();							//flip the first and second anchors: now the second anchor is used as shared maker between subject and query bridge
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
			return false;
		}
		
		Alignment start=read.getFirstAlignment();
		NodeVector 	current=null;
		int lastIdx=0, numOfAnchorsBefore=getNumberOfAnchors();
		for(Alignment alg:read.getAlignmentRecords()) {
			current=new NodeVector(alg.node, read.getVector(start,alg));
			if(getCompletionLevel()==3){
				if(!steps.nodes.contains(current)){
					int prev=-1, cur;
					for(int i=lastIdx; i<segments.size();i++) {
						BridgeSegment curSeg=segments.get(i);
						cur=curSeg.locateAndVote(current);
						if(cur<prev)
							break;
						else {
							prev=cur;
							lastIdx=i;
						}
					}
				}
			}else{
				//just adding
				steps.addNode(current);
			}
		}
		if(getCompletionLevel()<3)
			steps.connectBridgeSteps(false);		
				
		
		System.out.printf("After merging: %s\nstart=%s\nend=%s\n", pBridge.toString(), (steps.start==null?"null":steps.start.toString()), (steps.end==null?"null":steps.end.toString()));
		return getNumberOfAnchors()>numOfAnchorsBefore;
	
	}
	
	//adding new segment with/without update pBridge
	private void addSegment(BridgeSegment seg) {
		if(segments==null)
			segments=new ArrayList<>();
//		if(pBridge==null) {
//			segments.add(seg);
//			if(update)
//				pBridge=new BidirectedEdgePrototype(seg.pSegment.getNode0(), seg.pSegment.getNode1(), seg.pSegment.getDir0(), seg.pSegment.getDir1());
//		}
//		else if(pBridge.getNode1()==seg.pSegment.getNode0() && pBridge.getDir1()!=seg.pSegment.getDir0()){
//			segments.add(seg);
//			if(update)
//				pBridge.n1=seg.pSegment.n1;
//		}
		segments.add(seg);
	}


	
	
	private void reverse() {
		assert getNumberOfAnchors()==2:"Could only reverse a determined bridge (2 anchors)";
//		System.out.printf("Reversing the bridge %s:\n%s\n", getEndingsID(), getAllNodeVector());

		pBridge=new BidirectedEdgePrototype(pBridge.getNode1(), pBridge.getNode0(), pBridge.getDir1(), pBridge.getDir0());
		//reverse the segments
		if(segments!=null){
			ArrayList<BridgeSegment> tmp = segments;
			segments = new ArrayList<>();
			if(!tmp.isEmpty()) {
				ScaffoldVector brgVector=tmp.get(tmp.size()-1).endV;
				for(int i=tmp.size()-1; i>=0 ;i--)
					addSegment(tmp.get(i).reverse(brgVector));
			}
		}
		//reverse the nodes list
		steps.reverse();
		
//		System.out.printf("Reversed bridge %s:\n%s\n", getEndingsID(), getAllNodeVector());

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
				for(BidirectedPath path:seg.connectedPaths)
					retval+="( "+path.getId()+ " : vote=" + path.getVote() + " deviation=" + path.getDeviation() + " )";
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

	public BidirectedPath getBestPath() {
		BidirectedPath best=null,retval=null;
		for(BridgeSegment seg:segments) {
			if(seg.getNumberOfPaths()>0){
				if(seg.getNumberOfPaths()>1)
					//FIXME: double-check this
					seg.connectedPaths.sort(Comparator.comparing(BidirectedPath::getVote, Comparator.reverseOrder()).thenComparing(BidirectedPath::getDeviation));	
				
				best=seg.connectedPaths.get(0);
				if(retval==null)
					retval=best;
				else
					retval=retval.join(best);//assuming first path has highest score!!!
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
	class BridgeSegment{
//		ArrayList<Sequence> nnpReads; // to store nanopore data if needed
		ArrayList<BidirectedPath> connectedPaths;
		BidirectedEdgePrototype pSegment;
		ScaffoldVector startV, endV; // from bridge anchor (always +) to this segment's end
		int bestElections=0;
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
		
		BridgeSegment reverse(ScaffoldVector brgVector) {
			BridgeSegment retval = new BridgeSegment();
			retval.startV=ScaffoldVector.composition(endV, ScaffoldVector.reverse(brgVector));
			retval.endV=ScaffoldVector.composition(startV, ScaffoldVector.reverse(brgVector));
			retval.connectedPaths=new ArrayList<>();
			for(BidirectedPath p:connectedPaths)
				retval.connectedPaths.add(p.reverse());
			retval.pSegment=pSegment.reverse();
			return retval;
		}
		
		/*
		 * Find a node vector if it appear in one of the paths
		 * Return retval=-1 if not located within the range of this segment
		 * Return retval>0 if found in x paths
		 * Return retval=0 if located in the interval of this segment but not found in any path 
		 */
		int locateAndVote(NodeVector nv) {
			int retval=-1;
			if(!isConnected())
				return retval;
			
			NodeVector 	start=new NodeVector((BidirectedNode) pSegment.getNode0(), startV),
						end=new NodeVector((BidirectedNode) pSegment.getNode1(),endV);
			List<BidirectedPath> tobeRemoved=new ArrayList<>();
			if(nv.compareTo(start)*nv.compareTo(end)<=0) {
				retval=0;
				ScaffoldVector start2nv=ScaffoldVector.composition(nv.getVector(), ScaffoldVector.reverse(startV));
				int d=start2nv.distance((BidirectedNode) pSegment.getNode0(), nv.getNode());
				for(BidirectedPath p:connectedPaths) {
					if(p.checkDistanceConsistency(pSegment.getNode0(), nv.getNode(), start2nv.direction>0, d) >= 0) {
						p.upVote(1);
						if(p.getVote() > bestElections)
							bestElections=p.getVote();
						retval++;
					}
					else{
						p.downVote(1);
						if(p.getVote() < bestElections-BidirectedGraph.MAX_DIFF)
							tobeRemoved.add(p);
					}
				}
			}
			if(retval>0 && retval<connectedPaths.size())
				connectedPaths.removeAll(tobeRemoved);

			
			return retval;
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
			return endV;
		}
		ScaffoldVector getStartVector() {
			return startV;
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
		SortedSet<NodeVector> nodes;
		//save the two unique nodes at 2 endings of the bridge. When end is covered more than threshold, pBridge is set
		NodeVector start, end; 
		
		BridgeSteps(){
			nodes=new TreeSet<NodeVector>();
		}
		
		BridgeSteps(AlignedRead read){
			this();
			if(read.getAlignmentRecords().size()<2)
				return;
			if(SimpleBinner.getUniqueBin(read.getFirstAlignment().node)==null)
				if(SimpleBinner.getUniqueBin(read.getLastAlignment().node)!=null)
					read.reverse();
				else
					return;
					
					
			Alignment 	firstAlg = read.getFirstAlignment(),
						curAlg = null;				
					
			start=new NodeVector(firstAlg.node, new ScaffoldVector());
			nodes.add(start);
				
			for(int i=1;i<read.getAlignmentRecords().size();i++) {
				curAlg = read.getAlignmentRecords().get(i);
				NodeVector tmp=new NodeVector(curAlg.node, read.getVector(firstAlg, curAlg));
				nodes.add(tmp);
				
				if(SimpleBinner.getUniqueBin(curAlg.node)!=null)
					end=tmp;

			}
			//update the pBridge accordingly
			if(end!=null && end.qc()) {
				pBridge=new BidirectedEdgePrototype(firstAlg.node, end.getNode(), firstAlg.strand, end.getDirection(firstAlg.strand));//getDirection doesn't work with self-vector
			}else if(pBridge==null)
				pBridge=new BidirectedEdgePrototype(firstAlg.node, firstAlg.strand);
		}
		
	
		
		void addNode(NodeVector nv) {
			//TODO: optimize it! hashmap??? not linear searching like this!
			//FIXME: acting weird, not call equals() properly for 141o,20o: Because TreeSet used compareTo() instead of equals()!!! 
			//(https://dzone.com/articles/the-hidden-contract-between-equals-and-comparable)
			
			Iterator<NodeVector> ite = nodes.iterator();
			boolean found=false;
			while(ite.hasNext()){
				NodeVector tmp=ite.next();
				if(tmp.merge(nv)){
					//re-assign end node if there is another unique node with higher score
					if(SimpleBinner.getUniqueBin(tmp.getNode())!=null && !tmp.getVector().isIdentity()) {
						if(!end.equals(tmp)){
							if(end.nodeCover < tmp.nodeCover)
								end=tmp;
							else
								System.out.println("Conflict detected on end node: " + tmp.toString() + " to " + end.toString());
						}
						
					}
					//nv is already in the set
					found=true;
					break;

				}
			}
			
			if(!found) {
				nodes.add(nv);
				if(SimpleBinner.getUniqueBin(nv.getNode())!=null) {
					if(nv.getVector().isIdentity()){
						start=nv;	
					}else if(end==null){
						end=nv;
					}
				}
			}

			
//			if(nodes.contains(nv)){
//				System.out.println("Found " + nv);
//				Iterator<NodeVector> ite = nodes.iterator();
//				while(ite.hasNext()){
//					NodeVector tmp=ite.next();
//					if(tmp.equals(nv)){
//						tmp.score+=nv.score;
//						//assign end node
//						if(SimpleBinner.getUniqueBin(tmp.getNode())!=null && !tmp.getVector().isIdentity()) {
//							if(end==null)
//								end=tmp;
//							else if(!end.equals(tmp)){
//								if(end.score < tmp.score)
//									end=tmp;
//								else
//									System.out.println("Conflict detected on end node: " + tmp.toString() + " to " + end.toString());
//							}
//							
//						}
//						break;
//
//					}
//				}
//				
//			}else {
//				System.out.println("Not Found " + nv);
//				//assign start node
//				if(SimpleBinner.getUniqueBin(nv.getNode())!=null && nv.getVector().isIdentity())
//					start=nv;
//				nodes.add(nv);
//			}
		}
			
		
		//Try to make the continuous segments (those that have paths connected 2 ends). Return true if it possible
		boolean connectBridgeSteps(boolean force){
			if(isIdentifiable()){
				if(connectable() || force){
					pBridge.n1=new NodeDirection(end.node, end.getDirection(pBridge.getDir0()));
				}else{
					System.err.printf("Bridge %s is not qualified to connect yet: start=%s end=%s\n", pBridge.toString(), (start==null?"null":start.toString()), (end==null?"null":end.toString()));
					return false;
				}

			}else{
				System.err.printf("Bridge %s is not complete to connect!\n", pBridge.toString());
				return false;
			}	
			
			System.out.println("Trying to connect bridge " + pBridge.toString() + ":\n" + getAllNodeVector());
			//First build shortest tree from the end
			//TODO: optimize finding path by using this list
			HashMap<String,Integer> shortestMap = 
					graph.getShortestTreeFromNode(	(BidirectedNode)pBridge.getNode1(), 
													pBridge.getDir1(), 
													end.getVector().distance((BidirectedNode)pBridge.getNode0(), (BidirectedNode)pBridge.getNode1()));
			String key=pBridge.getNode0().getId()+(pBridge.getDir0()?"i":"o");
			if(!shortestMap.containsKey(key)){ //note the trick: direction BEFORE the path started
				System.out.printf("Shortest tree couldn't reach to the other end: ");

				if(!force){
					pBridge=new BidirectedEdgePrototype(pBridge.getNode0(), pBridge.getDir0());	
					System.out.println(" ignored!");
					return false;
				}else{
					System.out.println(" proceed for the last attempt!");

				}
			}else{
				System.out.printf("Shortest tree contain the other end: %s=%d\n", key, shortestMap.get(key));
			}
			
			Iterator<NodeVector> iterator = nodes.iterator();
			
			NodeVector 	prev=null, 
						current=null;
			boolean retval = false;
			ArrayList<NodeVector> inbetween = new ArrayList<>();
			while(iterator.hasNext()){
				if(current==null){
					prev=current=iterator.next();//first unique node, don't need to check quality
					continue;
				}
				
				current=iterator.next();
				key=current.getNode().getId() + (current.getDirection(pBridge.getDir0())?"o":"i");
				//need a quality-checking here before including into a segment step
				if(shortestMap.containsKey(key)){			 
					if((force&&current==end) || current.qc()){
						BridgeSegment seg = null;
						if(prev.getVector().isIdentity())
							seg=new BridgeSegment(	prev.getNode(), current.getNode(), 
													pBridge.getDir0(), 
													current.getDirection(pBridge.getDir0()), 
													prev.getVector(), current.getVector());
						else
							seg=new BridgeSegment(	prev.getNode(), current.getNode(), 
									!prev.getDirection(pBridge.getDir0()), 
									current.getDirection(pBridge.getDir0()), 
									prev.getVector(), current.getVector());
						
						System.out.print("...connecting " + seg.getId());
						if(seg.isConnected()){
							System.out.println(" :success!");		
							prev=current;
							//use qc-failed nodes to vote
							if(seg.getNumberOfPaths()>1)
								for(NodeVector nv:inbetween)
									seg.locateAndVote(nv);
							addSegment(seg);
							if(current.equals(end)) {
								retval=true;	
								break;
							}
							inbetween=new ArrayList<>();
							
						}else {
							System.out.println(" :skip!");
						}
					}else
						inbetween.add(current);
				}
					
			}
			
			if(!retval) {		
				segments=null;
				System.out.println("Failed to connect + " + pBridge.toString());
				
			}
			
			return retval;
		}
		
		boolean isIdentifiable(){
			return start!=null && end!=null;
		}
		
		boolean connectable(){
			if(!isIdentifiable())
				return false;
			
			return start.qc() && end.qc();
		}
		void reverse() {
			//reverse the nodes list
			SortedSet<NodeVector> reversedSet = new TreeSet<NodeVector>();
			ScaffoldVector rev=ScaffoldVector.reverse(end.getVector());//anchors number = 2 so there exist end node
			NodeVector tmp = null;
			for(NodeVector nv:nodes) {
				tmp=new NodeVector(nv.getNode(), ScaffoldVector.composition(nv.getVector(), rev), nv.nodeCover);
				reversedSet.add(tmp);

			}
			nodes=reversedSet;
			//re-assign start and end
			tmp=start;
			start=end;
			end=tmp;
			//change the vector also
			start.setVector(new ScaffoldVector());
			end.setVector(rev);
		}
		
		public String toString(){
			String retval="";
			for(NodeVector nv:nodes)
				retval+=nv.toString() + "\n";
			return retval;
		}
	}





}
