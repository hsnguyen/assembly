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
	BidirectedEdgePrototype pBridge; //note: the ending nodes must be unique, or else omitted
	ArrayList<BridgeSegment> segments;
	BridgeSteps steps;
	
	GoInBetweenBridge(BidirectedGraph graph, PopBin bin){
		this.graph=graph;
//		segments=new ArrayList<>();
//		steps=new BridgeSteps();
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
		boolean isComplete=false;
		if(retval==2 && segments!=null) {
			retval++;
			isComplete=true;
			for(BridgeSegment seg:segments)
				if(seg.getNumberOfPaths()!=1)
					isComplete=false;
		}
		if(isComplete)
			retval++;
		return retval;
	}

							
	//Merge 2 bridge (must share at least one same unique end-point) together
	GoInBetweenBridge merge(GoInBetweenBridge qBridge) {
		if(qBridge==null)
			return null;
		
		GoInBetweenBridge 	subject = this, 
							query = qBridge;
		if(subject.getCompletionLevel() < query.getCompletionLevel()) {
			subject=qBridge;
			query=this;
		}
		if(query.getCompletionLevel()==0 || subject.getCompletionLevel()==4)
			return subject;
		
		BidirectedNode 	sNode0=(BidirectedNode) subject.pBridge.getNode0(),
						sNode1=(BidirectedNode) subject.pBridge.getNode1(),
						qNode0=(BidirectedNode) query.pBridge.getNode0(),
						qNode1=(BidirectedNode) query.pBridge.getNode1();

		boolean found=true;
		if(sNode0!=qNode0) {							//if the starting points don't agree
			if(sNode0==qNode1) {						//the first unique contig actually located at the end of alignments list
				query.reverse();
			}else if(sNode1!=null){						//if not: this pBridge must have 2 anchors, the first anchor is not in the alignments list
				subject.reverse();							//flip the first and second anchors: now the second anchor is used as shared maker between subject and query bridge
				if(sNode1==qNode1)	
					query.reverse();
				else if(sNode1!=qNode0)
					found=false;
			}else { 
				found=false;
			}
		}
		if(!found) {
			System.err.println("Not found common end point between merging bridges");
			return subject;
		}
		
		SortedSet<NodeVector> sSteps=getStepsList(), qSteps=qBridge.getStepsList();
		
		Iterator<NodeVector> iterator = qSteps.iterator();		
		NodeVector 	current=null;
		int lastIdx=0;
		while(iterator.hasNext()){
			current=iterator.next();
			if(subject.getCompletionLevel()==3){
				if(!sSteps.contains(current)){
					int prev=-1, cur;
					for(int i=lastIdx; i<subject.segments.size();i++) {
						BridgeSegment curSeg=subject.segments.get(i);
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
				sSteps.add(current);
			}	
			
		}		
			
		return subject;
	}

	SortedSet<NodeVector> getStepsList(){
		return steps.nodes;
	}
	
	//adding new segment with/without update pBridge
	private void addSegment(BridgeSegment seg, boolean update) {
		if(segments==null)
			segments=new ArrayList<>();
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
//		System.out.printf("Reversing the bridge %s:\n%s\n", getEndingsID(), getAllNodeVector());

		pBridge=new BidirectedEdgePrototype(pBridge.getNode1(), pBridge.getNode0(), pBridge.getDir1(), pBridge.getDir0());
		//reverse the segments
		if(segments!=null){
			ArrayList<BridgeSegment> tmp = segments;
			segments = new ArrayList<>();
			if(!tmp.isEmpty()) {
				ScaffoldVector brgVector=tmp.get(tmp.size()-1).endV;
				for(int i=tmp.size()-1; i>=0 ;i--)
					addSegment(tmp.get(i).reverse(brgVector),false);
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
					retval+="( "+path.getId()+ " : " + path.getVote() + " )";
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
	class BridgeSegment{
//		ArrayList<Sequence> nnpReads; // to store nanopore data if needed
		ArrayList<BidirectedPath> connectedPaths;
		BidirectedEdgePrototype pSegment;
//		Range coverRange;
		ScaffoldVector startV, endV; // from bridge anchor (always +) to this segment's end
		
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
		 * Return -1 if not located within the range of this segment
		 * Return 1 if found in one of the paths
		 * Return 0 otherwise
		 */
		int locateAndVote(NodeVector nv) {
			int retval=-1;
			NodeVector 	start=new NodeVector((BidirectedNode) pSegment.getNode0(), startV),
						end=new NodeVector((BidirectedNode) pSegment.getNode1(),endV);
			if(nv.compareTo(start)*nv.compareTo(end)<=0) {
				retval++;
				ScaffoldVector start2nv=ScaffoldVector.composition(nv.getVector(), ScaffoldVector.reverse(startV));
				int d=start2nv.distance((BidirectedNode) pSegment.getNode0(), nv.getNode());
				for(BidirectedPath p:connectedPaths) {
					if(p.checkDistanceConsistency(pSegment.getNode0(), nv.getNode(), start2nv.direction>0, d) >= 0) {
						p.upVote(1);
						retval++;
					}
					else
						p.downVote(1);
				}
			}
			//TODO: filter out low score paths + check completeness here???
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
					read=read.reverse();
				else
					return;
					
					
			Alignment 	firstAlg = read.getFirstAlignment(),
						curAlg = null;				
					
			start=new NodeVector(firstAlg.node, new ScaffoldVector());
			addNode(start);;
			
			for(int i=1;i<read.getAlignmentRecords().size();i++) {
				curAlg = read.getAlignmentRecords().get(i);
				addNode(new NodeVector(curAlg.node, read.getVector(firstAlg, curAlg)));
				
			}
			//update the pBridge accordingly
			if(end!=null) {
				pBridge=new BidirectedEdgePrototype(firstAlg.node, end.getNode(), firstAlg.strand, end.getDirection(firstAlg.strand));//getDirection doesn't work with self-vector
			}else if(pBridge==null)
				pBridge=new BidirectedEdgePrototype(firstAlg.node, firstAlg.strand);
		}
		
		void addNode(NodeVector nv) {
			if(nodes.contains(nv)){
				Iterator<NodeVector> ite = nodes.iterator();
				while(ite.hasNext()){
					NodeVector tmp=ite.next();
					if(tmp.equals(nv)){
						tmp.score+=nv.score;
						break;
					}
				}
			}else {
				//assign end nodes by checking its uniqueness 
				if(SimpleBinner.getUniqueBin(nv.getNode())!=null && end==null){
					if(nv.getVector().isIdentity())
						start=nv;
					else
						end=nv;
				}
				
				if(!isComplete() || nv.compareTo(start)*nv.compareTo(end)<=0) //opt:loss information here
					nodes.add(nv);
			}
		}
			
		
		//Try to make the continuous segments (those that have paths connected 2 ends). Return true if it possible
		boolean connectBridgeSteps(){
			if(!connectable()){
				System.err.println("Bridge is not qualified to connect yet!");
				return false;
			}
			System.out.println("Trying to connect bridge " + pBridge.toString() + ":\n" + getAllNodeVector());
			//First build shortest tree from the end
			//TODO: optimize finding path by using this list
			HashMap<String,Integer> shortestMap = 
					graph.getShortestTreeFromNode(	(BidirectedNode)pBridge.getNode1(), 
													pBridge.getDir1(), 
													nodes.last().getVector().distance((BidirectedNode)pBridge.getNode0(), (BidirectedNode)pBridge.getNode1()));
			if(!shortestMap.containsKey(pBridge.getNode0().getId()+(pBridge.getDir0()?"i":"o"))){ //note the trick: direction BEFORE the path started
				System.err.println("Shortest tree couldn't reach to the other end!");
				//remove the the end node since it's unreachable
				nodes.remove(end);
				end=null;
				pBridge=new BidirectedEdgePrototype(pBridge.getNode0(), pBridge.getDir0());
				
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
				if((current!=end && !current.qc()) || shortestMap.containsKey(current.getNode().getId())){
					continue;
				}else{
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
						addSegment(seg,true);
						prev=current;
						if(current.equals(end))
							retval=true;						
					}else
						System.out.println(" :fail!");

				}
					
				
			}
			
			return retval;
		}
		
		boolean isComplete(){
			return start!=null && end!=null;
		}
		
		boolean connectable(){
			if(!isComplete())
				return false;
			
			boolean retval=true;
			//TODO: check NodeVector.qc() & density???
			
			return retval;
		}
		void reverse() {
			//reverse the nodes list
			SortedSet<NodeVector> reversedSet = new TreeSet<NodeVector>();
			ScaffoldVector rev=ScaffoldVector.reverse(end.getVector());//anchors number = 2 so there exist last()
			NodeVector tmp = null;
			for(NodeVector nv:nodes) {
				tmp=new NodeVector(nv.getNode(), ScaffoldVector.composition(nv.getVector(), rev), nv.score);
				reversedSet.add(tmp);

			}
			nodes=reversedSet;
			//re-assign start and end
			start=nodes.first();
			end=nodes.last();

		}
		
		public String toString(){
			String retval="";
			for(NodeVector nv:nodes)
				retval+=nv.toString() + "\n";
			return retval;
		}
	}




}
