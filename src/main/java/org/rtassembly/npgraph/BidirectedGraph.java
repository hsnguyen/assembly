package org.rtassembly.npgraph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Optional;
import java.util.PriorityQueue;
import java.util.Stack;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.graphstream.graph.*;
import org.graphstream.graph.implementations.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;


public class BidirectedGraph extends MultiGraph{
    static int KMER=127;
    static double RCOV=0.0;
	SimpleBinner binner;

    volatile static double ILLUMINA_READ_LENGTH=300; //Illumina MiSeq
    
    static final double R_TOL=.3;// relative tolerate: can be interpreted as long read error rate (10-25%)
    static final int A_TOL=200;// absolute tolerate: can be interpreted as long read absolute error bases (200bp)

    static final int D_LIMIT=2105; //distance bigger than this will be ignored
    public static int S_LIMIT=125;// maximum number of DFS steps
    
    //provide mapping from unique directed node to its corresponding bridge
    //E.g: 103-: <103-82-> also 82+:<82+103+>
    private HashMap<String, NewBridge> bridgesMap; 
    
    private static final Logger LOG = LoggerFactory.getLogger(BidirectedGraph.class);

    // *** Constructors ***
	/**
	 * Creates an empty graph.
	 * 
	 * @param id
	 *            Unique identifier of the graph.
	 * @param strictChecking
	 *            If true any non-fatal error throws an exception.
	 * @param autoCreate
	 *            If true (and strict checking is false), nodes are
	 *            automatically created when referenced when creating a edge,
	 *            even if not yet inserted in the graph.
	 * @param initialNodeCapacity
	 *            Initial capacity of the node storage data structures. Use this
	 *            if you know the approximate maximum number of nodes of the
	 *            graph. The graph can grow beyond this limit, but storage
	 *            reallocation is expensive operation.
	 * @param initialEdgeCapacity
	 *            Initial capacity of the edge storage data structures. Use this
	 *            if you know the approximate maximum number of edges of the
	 *            graph. The graph can grow beyond this limit, but storage
	 *            reallocation is expensive operation.
	 */
	public BidirectedGraph(String id, boolean strictChecking, boolean autoCreate,
			int initialNodeCapacity, int initialEdgeCapacity) {
		super(id, strictChecking, autoCreate);
		bridgesMap=new HashMap<String, NewBridge>(initialNodeCapacity*2);
		
		// All we need to do is to change the node & edge factory
		setNodeFactory(new NodeFactory<BidirectedNode>() {
			public BidirectedNode newInstance(String id, Graph graph) {
				return new BidirectedNode((AbstractGraph) graph, id);
			}
		});

		setEdgeFactory(new EdgeFactory<BidirectedEdge>() {
			public BidirectedEdge newInstance(String id, Node src, Node dst, boolean directed) { //stupid??
				return new BidirectedEdge(id, (AbstractNode)src, (AbstractNode)dst);
			}
		});
		
	}
	public HashSet<NewBridge> getUnsolvedBridges(){
		HashSet<NewBridge> retval = new HashSet<NewBridge>();
		for(NewBridge brg:bridgesMap.values()){
			if(!brg.isComplete())
				retval.add(brg);
		}
		
		return retval;
	}

	/**
	 * Creates an empty graph with default edge and node capacity.
	 * 
	 * @param id
	 *            Unique identifier of the graph.
	 * @param strictChecking
	 *            If true any non-fatal error throws an exception.
	 * @param autoCreate
	 *            If true (and strict checking is false), nodes are
	 *            automatically created when referenced when creating a edge,
	 *            even if not yet inserted in the graph.
	 */
	public BidirectedGraph(String id, boolean strictChecking, boolean autoCreate) {
		this(id, strictChecking, autoCreate, DEFAULT_NODE_CAPACITY,
				DEFAULT_EDGE_CAPACITY);
	}

	/**
	 * Creates an empty graph with strict checking and without auto-creation.
	 * 
	 * @param id
	 *            Unique identifier of the graph.
	 */
	public BidirectedGraph(String id) {
		this(id, true, false, 1000, 10000);
	}
	
    public BidirectedGraph(){
    	this("Assembly graph",true,false, 1000, 10000);
        setKmerSize(127);//default kmer size used by SPAdes to assembly MiSeq data
    }
	
	protected BidirectedEdge addEdge(AbstractNode src, AbstractNode dst, boolean dir0, boolean dir1){
		BidirectedEdge tmp = (BidirectedEdge) addEdge(BidirectedEdge.createID(src, dst, dir0, dir1), src, dst);
		return tmp;
	}
	
	
	public String printEdgesOfNode(BidirectedNode node){
		Stream<Edge> 	ins = getNode(node.getId()).enteringEdges(),
						outs = getNode(node.getId()).leavingEdges();
		String retval=node.getId() + ": IN={ ";
//		ins.map(e->(retval += e.getId() + " "));
		retval+=ins.map(e -> (e.getId()+" ")).reduce((a,b)->(a+b));
		retval+="}; OUT={ ";
		retval+=outs.map(e -> (e.getId()+" ")).reduce((a,b)->(a+b));
		retval+="}";
		return retval;		
	}
	
	public void outputFASTA(String fileName) throws IOException {
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(fileName);
		
		for(Node node:this) {
			Sequence seq=(Sequence) node.getAttribute("seq");
			if((node.getDegree()==0 && (seq.length() < SimpleBinner.SIG_CTG_LEN) || node.getNumber("cov") < RCOV*0.1))//noises
				continue;
			seq.writeFasta(out);
		}
		
		out.close();
	}
	public void outputGFA(String fileName) throws IOException {
		//TODO: implement it
	}
    
	
    public static int getKmerSize(){
    	return BidirectedGraph.KMER;
    }
    public static void setKmerSize(int kmer){
    	BidirectedGraph.KMER=kmer;
    }
    
    /**************************************************************************************************
     ********************** utility functions to serve the assembly algo ****************************** 
     * ***********************************************************************************************/
    
    synchronized public NewBridge getBridgeFromMap(Node unqNode, boolean direction){
    	String key=unqNode.getId()+(direction?"o":"i"); //true:going outward, false:going inward
    	return bridgesMap.get(key);
    }
    
    // when this unique node actually contained by a bridge
    synchronized public void updateBridgesMap(Node unqNode, NewBridge bidirectedBridge){
    	bridgesMap.put(unqNode.getId()+"o", bidirectedBridge);
    	bridgesMap.put(unqNode.getId()+"i", bidirectedBridge);
    }
    // when there is new unique bridge 
    synchronized protected void updateBridgesMap(NewBridge bidirectedBridge) {
    	if(bidirectedBridge==null || bidirectedBridge.getNumberOfAnchors()==0)
    		return;

		try {
	    	Node 	startNode=bidirectedBridge.pBridge.getNode0(),
	    			endNode=bidirectedBridge.pBridge.getNode1();
	    	boolean startNodeDir, endNodeDir;
			startNodeDir = bidirectedBridge.pBridge.getDir0();
			endNodeDir = bidirectedBridge.pBridge.getDir1();
			

	    	if(SimpleBinner.getUniqueBin(startNode)!=null)
	    		bridgesMap.put(startNode.getId()+(startNodeDir?"o":"i"), bidirectedBridge);
	    	
	    	if(SimpleBinner.getUniqueBin(endNode)!=null)
	    		bridgesMap.put(endNode.getId()+(endNodeDir?"o":"i"), bidirectedBridge);
		} catch (Exception e) {
			LOG.error("Invalid bridge to add to bridge map: " + bidirectedBridge.toString());
			e.printStackTrace();
		}
    	
    }
    
    //when there is a path that could represent a bridge (half or full)
    synchronized protected void updateBridgesMap(BidirectedPath path){
    	if(path==null || path.size() < 2)
    		return;
    	try{
	    	BidirectedNode 	startNode=path.getFirstNode(),
	    					endNode=path.getLastNode();
	    	boolean startNodeDir=path.getFirstNodeDirection(),
	    			endNodeDir=path.getLastNodeDirection();
	    	if(SimpleBinner.getUniqueBin(startNode)!=null){
	    		bridgesMap.put(startNode.getId()+(startNodeDir?"o":"i"), new NewBridge(this,path));
	    	}
	    	if(SimpleBinner.getUniqueBin(endNode)!=null){
	    		bridgesMap.put(endNode.getId()+(endNodeDir?"o":"i"), new NewBridge(this,path));
	    	}
		} catch (Exception e) {
			LOG.error("Invalid path to add to bridge map: " + path.getId());
			e.printStackTrace();
		}
    	
    }
    
    //Return bridge in the map (if any) that share the same bases (unique end) 
    synchronized public NewBridge getHomoBridgeFromMap(AlignedRead algRead){
    	NewBridge retval = null, tmp = null;
    	if(algRead!=null){
	    	Node 	startNode=algRead.getFirstAlignment().node,
	    			endNode=algRead.getLastAlignment().node;
	    	boolean startNodeDir=algRead.getFirstAlignment().strand,
	    			endNodeDir=!algRead.getLastAlignment().strand;
	    	
	    	if(SimpleBinner.getUniqueBin(startNode)!=null){
	    		tmp=bridgesMap.get(startNode.getId()+(startNodeDir?"o":"i"));
	    		if(tmp!=null)
	    			retval=tmp;
	    	}
	    	
	    	if(SimpleBinner.getUniqueBin(endNode)!=null){
	    		tmp=bridgesMap.get(endNode.getId()+(endNodeDir?"o":"i"));
	    		if(tmp!=null){
	    			if(retval==null || retval.getNumberOfAnchors()<tmp.getNumberOfAnchors())
	    				retval=tmp;
	    			
	    		}
	    	}
	    	
    	}
    	
    	return retval;
    }
    
    synchronized public void removeBridgesFromNode(Node node){
    	
    }
    
    //If the node was wrongly identified as unique before, do things...
    synchronized public void destroyFalseBridges(Node node){
//    	node.removeAttribute("unique");
//    	binner.node2BinMap.remove(node);
//    	
//    	NewBridge 	falseBridgeFrom = bridgesMap.get(node.getId()+ "o"),
//					falseBridgeTo = bridgesMap.get(node.getId() + "i");
//    	if(falseBridgeFrom!=null){
//    		if(falseBridgeFrom.getBridgeStatus()>=0){
//    			if(falseBridgeFrom.getBestPath()!=null)
//    				falseBridgeFrom.getBestPath().revert();
//    			  		
//    			falseBridgeFrom.halfPaths.addAll(falseBridgeFrom.fullPaths);
//    			falseBridgeFrom.fullPaths=new ArrayList<BidirectedPath>();
//    			//todo: remove 2 entries from bridgeMap here...
//    			
//    		}else{    		
//    			//destroy the half bridge
//        		bridgesMap.remove(node.getId()+ "o");
//    		}
//
//    	}
//    	
//    	if(falseBridgeTo!=null){
//    		if(falseBridgeTo.getBridgeStatus()>=0){
//    			if(falseBridgeTo.getBestPath()!=null)
//    				falseBridgeTo.getBestPath().revert();
//    			  		
//    			falseBridgeTo.halfPaths.addAll(falseBridgeTo.fullPaths);
//    			falseBridgeTo.fullPaths=new ArrayList<BidirectedPath>();
//    			//todo: remove 2 entries from bridgeMap here...
//
//    		}else{    		
//    			//destroy the half bridge
//        		bridgesMap.remove(node.getId()+ "i");
//    		}
//    	}
    		
    }
    
    synchronized public void binning() {
    	binner=new SimpleBinner(this);
    	binner.estimatePathsByCoverage();
    }

    
//    /*
//     * This function deduces a full path in this graph between 2 nodes aligned with a long read
//     * 
//     */
//    synchronized protected ArrayList<BidirectedPath> getClosestPaths(Alignment from, Alignment to){
//    	assert from.readID==to.readID && to.compareTo(from)>=0:"Illegal alignment pair to find path!";
//    	
//    	int distance=to.readAlignmentStart()-from.readAlignmentEnd();
//    	if(distance>BidirectedGraph.D_LIMIT)
//    		return null;
//    	
//    	BidirectedNode srcNode = from.node,
//    					dstNode = to.node;
//    	System.out.println("Looking for path between " + srcNode.getId() + " to " + dstNode.getId() + " with distance " + distance);
//    	BidirectedPath 	tmp = new BidirectedPath();
//    	ArrayList<BidirectedPath>	possiblePaths = new ArrayList<BidirectedPath>(),
//    								retval = new ArrayList<BidirectedPath>();
//    	tmp.setRoot(srcNode);  	
//   
//    	traverse(tmp, dstNode, possiblePaths, distance, distance>A_TOL?(int) (R_TOL*distance):A_TOL, from.strand, to.strand, 0);
////    	traverse(tmp, dstNode, possiblePaths, distance, 500, from.strand, to.strand, 0);
//
//    	/**************************************************************************************
//    	 * To cover the missing edges due to big k-mer of DBG
//    	 **************************************************************************************/
//
//    	if(possiblePaths.isEmpty()){
//			if(SimpleBinner.getUniqueBin(srcNode)!=null && SimpleBinner.getUniqueBin(dstNode)!=null && srcNode.getDegree() == 1 && dstNode.getDegree()==1 &&
//				Math.min(from.quality, to.quality) >= Alignment.GOOD_QUAL)
//    		{
////    			BidirectedEdge pseudoEdge = new BidirectedEdge(srcNode, dstNode, from.strand, to.strand);
//    			//save the corresponding content of long reads to this edge
//	    		//FIXME: save nanopore reads into this pseudo edge to run consensus later
//    			BidirectedEdge pseudoEdge = addEdge(srcNode, dstNode, from.strand, !to.strand);
//    			pseudoEdge.setAttribute("dist", distance);
//    			tmp.add(pseudoEdge);
//    			retval.add(tmp);
//    			System.out.println("pseudo path from " + srcNode.getId() + " to " + dstNode.getId() + " distance=" + distance);
//    			
////    			HybridAssembler.promptEnterKey();
//    			return retval;
//    		}else
//    			return null;
//    	}
//    	/*****************************************************************************************/
//    	
////    	double bestScore=possiblePaths.get(0).getDeviation();
////    	for(int i=0;i<possiblePaths.size();i++){
////    		BidirectedPath p = possiblePaths.get(i);
////    		if(p.getDeviation()>bestScore*(1+R_TOL))
////    			break;
////    		retval.add(p);
////    	}
////    	
////    	return retval;
//    	//FIXME: reduce the number of returned paths here (based on the score?)
//    	return possiblePaths;
//    	
//    }
//    
//    synchronized protected ArrayList<BidirectedPath> getClosestPaths(BidirectedNode srcNode, boolean srcDir, BidirectedNode dstNode, boolean dstDir, int distance, boolean force){
//    	if(distance>BidirectedGraph.D_LIMIT)
//    		return null;
//    	System.out.println("Looking for path between " + srcNode.getId() + " to " + dstNode.getId() + " with distance " + distance);
//		BidirectedPath 	tmp = new BidirectedPath();
//		ArrayList<BidirectedPath>	possiblePaths = new ArrayList<BidirectedPath>(),
//									retval = new ArrayList<BidirectedPath>();
//		tmp.setRoot(srcNode);  		
////		DFSAllPaths(srcNode, dstNode, srcDir, dstDir, distance, true);
//		traverse(tmp, dstNode, possiblePaths, distance, distance>A_TOL?(int) (R_TOL*distance):A_TOL, srcDir, !dstDir, 0);
//		if(possiblePaths.isEmpty()){
//			if(SimpleBinner.getUniqueBin(srcNode)!=null && SimpleBinner.getUniqueBin(dstNode)!=null && srcNode.getDegree() == 1 && dstNode.getDegree()==1 && force){
//				//save the corresponding content of long reads to this edge
//				//FIXME: save nanopore reads into this pseudo edge to run consensus later
//				BidirectedEdge pseudoEdge = addEdge(srcNode, dstNode, srcDir, dstDir);
//				pseudoEdge.setAttribute("dist", distance);
//				tmp.add(pseudoEdge);
//				retval.add(tmp);
//				System.out.println("pseudo path from " + srcNode.getId() + " to " + dstNode.getId() + " distance=" + distance);
//				
//				return retval;
//    		}else
//    			return null;
//
//		}
//		//double bestScore=possiblePaths.get(0).getDeviation();
//		//for(int i=0;i<possiblePaths.size();i++){
//		//	BidirectedPath p = possiblePaths.get(i);
//		//	if(p.getDeviation()>bestScore*(1+R_TOL))
//		//		break;
//		//	retval.add(p);
//		//}
//		//
//		//return retval;
//		//FIXME: reduce the number of returned paths here (calculate edit distance with nanopore read: dynamic programming?)
//		return possiblePaths;
//    }
//    
//    synchronized private void traverse(	BidirectedPath path, BidirectedNode dst, ArrayList<BidirectedPath> curResult, 
//    						int distance, int tolerance, boolean srcDir, boolean dstDir, int stepCount)
//    {
//    	//stop if it's going too deep!
//    	if(stepCount >= S_LIMIT) {
//			System.out.println("Stop going too deep with  "+stepCount+" levels already!");
//    		return;
//    	}
//    	
//    	BidirectedNode currentNode=(BidirectedNode) path.peekNode();
//    	BidirectedEdge currentEdge;
//    	boolean curDir;//direction to the next node, = ! previous'
//    	
//    	Iterator<Edge> ite;
//    	if(path.size() <= 1) //only root
//			curDir=srcDir;//re-check
//		else{
//			currentEdge = (BidirectedEdge) path.peekEdge();
//			curDir = !((BidirectedEdge) currentEdge).getDir(currentNode);
//		}
//		ite=curDir?currentNode.leavingEdges().iterator():currentNode.enteringEdges().iterator();
//    	while(ite.hasNext()){
//    		BidirectedEdge e = (BidirectedEdge) ite.next();
//			path.add(e);
//			
//			int delta=Math.abs(distance-e.getLength());
//			//note that traversing direction (true: template, false: reverse complement) of destination node is opposite its defined direction (true: outward, false:inward) 
//			if(e.getOpposite(currentNode)==dst && e.getDir(dst)!=dstDir && delta < tolerance){ 
//		    	BidirectedPath 	curPath=curResult.isEmpty()?new BidirectedPath():curResult.get(0), //the best path saved among all possible paths from the list curResult
//		    					tmpPath=new BidirectedPath(path);
//		    	tmpPath.setDeviation(delta);
//		    	if(	delta < curPath.getDeviation() )
//		    		curResult.add(0, tmpPath);
//		    	else
//		    		curResult.add(tmpPath);
//				
//				System.out.println("Hit added: "+path.getId()+"(candidate deviation: "+delta+")");
//			}else{
//				int newDistance = distance - ((Sequence) e.getOpposite(currentNode).getAttribute("seq")).length() - e.getLength();
////				System.out.println("adding edge: " + e.getId() + " length=" + e.getLength() +" -> distance=" + newDistance);
//				if (newDistance - e.getLength() < -tolerance){
//					System.out.println("Stop go to edge " + e.getId() + " from path with distance "+newDistance+" already! : "+path.getId());
//				}else
//					traverse(path, dst, curResult, newDistance, tolerance, srcDir, dstDir, stepCount+1);
//
//			}
//			path.popNode();
//    	
//    	}
//    }
    
    synchronized ArrayList<BidirectedPath> DFSAllPaths(Alignment from, Alignment to){
    	assert from.readID==to.readID && to.compareTo(from)>=0:"Illegal alignment pair to find path!"; 	
    	int distance=to.readAlignmentStart()-from.readAlignmentEnd();
    	BidirectedNode srcNode = from.node,
						dstNode = to.node;
    	boolean srcDir = from.strand, dstDir = !to.strand;
    	return DFSAllPaths(srcNode, dstNode, srcDir, dstDir, distance, Math.min(from.quality, to.quality) >= Alignment.GOOD_QUAL);
    }
    
	synchronized ArrayList<BidirectedPath> DFSAllPaths(BidirectedNode srcNode, BidirectedNode dstNode, boolean srcDir, boolean dstDir, int distance, boolean assureFlag)
	{
    	if(distance>BidirectedGraph.D_LIMIT)
    		return null;
    	System.out.println("Looking for DFS path between " + srcNode.getId() + " to " + dstNode.getId() + " with distance " + distance);
		ArrayList<BidirectedPath> possiblePaths = new ArrayList<BidirectedPath>(), 
									retval=new ArrayList<BidirectedPath>();
		//1. First build shortest tree from dstNode 		
		HashMap<String,Integer> shortestMap = getShortestTreeFromNode(dstNode, dstDir, distance);
		BidirectedPath 	path = new BidirectedPath();
		path.setRoot(srcNode);  		

		//2. DFS from srcNode with the distance info above
		NodeState curNodeState = new NodeState(srcNode, !srcDir); // first node is special
		if(shortestMap.containsKey(curNodeState.toString())) {
			
			Stack<List<Edge>> stack = new Stack<>();
			List<Edge> curList = (curNodeState.getDir()?curNodeState.getNode().enteringEdges():curNodeState.getNode().leavingEdges()).collect(Collectors.toList());
			stack.push(curList);
			
			int shortestDist2Dest = shortestMap.get(curNodeState.toString());
			int tolerance = A_TOL, 
				delta;
			BidirectedEdge curEdge = null;
			System.out.println("Found " + curNodeState.toString() + " with shortest distance=" + shortestDist2Dest);
			
			while(true) {
//				System.out.println("\nCurrent stack: ");
//				for(List<Edge> l:stack) {
//					System.out.print("[");
//					for(Edge e:l)
//							System.out.printf("%s; ", e.getId());
//					System.out.println("]");
//				}
//				System.out.println("Current path: " + path.getId());
//				System.out.println("Current distance: " + distance);

				curList=stack.peek();
				
				if(curList.isEmpty()) {
					if(path.size() <= 1)
						break;
					stack.pop();
//					System.out.print("removing edge " + path.peekEdge().getId());
					distance += (int)path.peekNode().getNumber("len") + ((BidirectedEdge) path.popEdge()).getLength();
//					System.out.println(" -> distance = " + distance);
				}else {
					curEdge=(BidirectedEdge) curList.remove(0);
					BidirectedNode from = (BidirectedNode) path.peekNode(),
									to = (BidirectedNode) curEdge.getOpposite(from);
					boolean dir = curEdge.getDir(to);

					AtomicInteger limit = new AtomicInteger(distance + tolerance);
			    	stack.push(
			    			(curEdge.getDir(to)?to.enteringEdges():to.leavingEdges())
			    			.filter(e->{
			    				BidirectedNode n=(BidirectedNode) e.getOpposite(to);
			    				NodeState ns = new NodeState(n, ((BidirectedEdge) e).getDir(n));
			    				
			    				return shortestMap.containsKey(ns.toString()) && shortestMap.get(ns.toString()) < limit.get();
			    			})
			    			.collect(Collectors.toList())			    			
			    			);
			    	
					path.add(curEdge);

					delta=Math.abs(distance-curEdge.getLength());
					//note that traversing direction (true: template, false: reverse complement) of destination node is opposite its defined direction (true: outward, false:inward) 
					if(to==dstNode && dir==dstDir && delta < tolerance){ 
				    	BidirectedPath 	curPath=possiblePaths.isEmpty()?new BidirectedPath():possiblePaths.get(0), //the best path saved among all possible paths from the list curResult
				    					tmpPath=new BidirectedPath(path);
				    	tmpPath.setDeviation(delta);
				    	if(	delta < curPath.getDeviation() )
				    		possiblePaths.add(0, tmpPath);
				    	else
				    		possiblePaths.add(tmpPath);
						
						System.out.println("Hit added: "+path.getId()+"(candidate deviation: "+delta + "; depth: " + path.size()+")");
						
						if(possiblePaths.size() > S_LIMIT) //not go too far
							break;
					}
					
					distance -= (int)to.getNumber("len") + curEdge.getLength();
//					System.out.println("adding edge: " + curEdge.getId() + " length=" + (int)to.getNumber("len") +" -> distance=" + distance);

					
				}
				
			}
		} 
		
//		traverse(tmp, dstNode, possiblePaths, distance, distance>A_TOL?(int) (R_TOL*distance):A_TOL, srcDir, !dstDir, 0);
		
		if(possiblePaths.isEmpty()){
			if(SimpleBinner.getUniqueBin(srcNode)!=null && SimpleBinner.getUniqueBin(dstNode)!=null && srcNode.getDegree() == 1 && dstNode.getDegree()==1 && assureFlag){
				//save the corresponding content of long reads to this edge
				//FIXME: save nanopore reads into this pseudo edge to run consensus later
				BidirectedEdge pseudoEdge = addEdge(srcNode, dstNode, srcDir, dstDir);
				pseudoEdge.setAttribute("dist", distance);
				path.add(pseudoEdge);
				possiblePaths.add(path);
				System.out.println("pseudo path from " + srcNode.getId() + " to " + dstNode.getId() + " distance=" + distance);
				
				return possiblePaths;
    		}else
    			return null;

		}
		
		double bestScore=possiblePaths.get(0).getDeviation();
		int keepMax = 10;//only keep this many possible paths 
		for(int i=0;i<possiblePaths.size();i++){
			BidirectedPath p = possiblePaths.get(i);
			if(p.getDeviation()>bestScore+Math.abs(distance+getKmerSize())*R_TOL || i>=keepMax)
				break;
			retval.add(p);
		}
		
		//FIXME: reduce the number of returned paths here (calculate edit distance with nanopore read: dynamic programming?)
		return retval;
	}    
    
    /*
     * Get a map showing shortest distances from surrounding nodes to a *rootNode* expanding to a *direction*, within a *distance*
     * based on Dijkstra algorithm
     */
    private HashMap<String, Integer> getShortestTreeFromNode(BidirectedNode rootNode, boolean expDir, int distance){
		PriorityQueue<NodeState> pq = new PriorityQueue<>();
		HashMap<String,Integer> retval = new HashMap<>();
		
		int curDistance=(int) -rootNode.getNumber("len"), newDistance=0;
		NodeState curNS = new NodeState(rootNode, expDir, curDistance);
		pq.add(curNS);
		
		System.out.println("Building shortest tree for " + rootNode.getId());
		retval.put(curNS.toString(), curDistance); // direction from the point of srcNode
		while(!pq.isEmpty()) {
//			System.out.println("Current queue: ");
//			for(NodeState n:pq) {
//				System.out.printf("%s:%d\t", n.toString(), n.getWeight());
//			}
//			System.out.println();
			curNS=pq.poll();
			curDistance=curNS.getWeight();

			Iterator<Edge> ite=curNS.getDir()?curNS.getNode().leavingEdges().iterator():curNS.getNode().enteringEdges().iterator();
			while(ite.hasNext()) {
	    		BidirectedEdge edge = (BidirectedEdge) ite.next();
	    		BidirectedNode nextNode = (BidirectedNode) edge.getOpposite(curNS.getNode());
	    		boolean direction =  !edge.getDir(nextNode);
	    		newDistance=curDistance+edge.getLength()+(int)curNS.getNode().getNumber("len");
	    		if(newDistance > distance+A_TOL)
	    			continue;
	    		
	    		NodeState nextNS = new NodeState(nextNode, direction, newDistance);
	    		String key=nextNS.toString();
	    		if(retval.containsKey(key)) {
	    			if(retval.get(key) > newDistance) {
	    				//update the map
	    				retval.put(key, newDistance);
	    				//update the queue
	    				pq.remove(nextNS);
	    				pq.add(nextNS);
	    			}
	    		}else {
	    			retval.put(key, newDistance);
    				pq.add(nextNS);
	    		}

			}
			
		}
		
		System.out.println("Shortest tree from " + rootNode.getId() + (expDir?"o":"i"));
		for(String str:retval.keySet()) {
			System.out.printf("%s:%d; ", str, retval.get(str));
		}
		System.out.println();
		return retval;
    }
    
    
    /*
     * Find bridges based on list of Alignments.
     * Return list of bridges with endings as markers and alignments of non-markers in-between.
     */
    synchronized protected List<BidirectedPath> uniqueBridgesFinding(Sequence nnpRead, ArrayList<Alignment> alignments) {
		if(nnpRead==null || alignments.size()<=1)
			return null;
		
		System.out.println("=================================================");
		for(Alignment alg:alignments)
			System.out.println("\t"+alg.toString());
		System.out.println("=================================================");
		//First bin the alignments into different overlap regions			
		//only considering useful alignments
		HashMap<Range,Alignment> allAlignments = new HashMap<Range,Alignment>();
	
		for(Alignment alg:alignments){
			if(alg.useful){
				Range range = new Range(alg.readAlignmentStart(),alg.readAlignmentEnd());
				allAlignments.put(range, alg);
			}
		}
		//now get all the bin in order
		List<Range> baseRanges=new ArrayList<Range>(allAlignments.keySet());
		List<List<Range>> rangeGroups = MetaRange.getOverlappingGroups(baseRanges);
		
		if(rangeGroups.size() < 2)
			return null;
		
		List<Range> stepRanges=new ArrayList<>(rangeGroups.size());
		
		System.out.println("Step ranges: ");
	    for(List<Range> group : rangeGroups){
	    	int maxscore=-1;
	    	Range rangeOfBest=null;
	    	for(Range range:group) { 
	    		System.out.print(allAlignments.get(range).node.getId() + " "+ binner.getBinsOfNode(allAlignments.get(range).node) + ": " + range + "; ");	    
	    		if(allAlignments.get(range).quality > maxscore)
	    			rangeOfBest=range;
	    	}
	  	
	    	stepRanges.add(rangeOfBest);
	    	System.out.println();
	    }

	    PopBin leastAbundancePop;
		Optional<PopBin> aBin=
				stepRanges.stream()
	    			.map(r->SimpleBinner.getUniqueBin(allAlignments.get(r).node))
	    			.filter(b->b!=null)
	    			.reduce((a,b)->(a.estCov > b.estCov?b:a));
		if(aBin.isPresent())
			leastAbundancePop=aBin.get();
		else 
			return null;
	    
		Range curRange = stepRanges.get(0);

		Alignment 	curAlignment =allAlignments.get(curRange),
					nextAlignment;
		ArrayList<AlignedRead> allBuildingBlocks = new ArrayList<>();
		AlignedRead	curBuildingBlocks=new AlignedRead(nnpRead, curAlignment);
		PopBin tmpBin=null;
		HashMap<PopBin, Long> bins2Length = new HashMap<PopBin,Long>();

		
		tmpBin=SimpleBinner.getUniqueBin(curAlignment.node);
		if( tmpBin != null){
			bins2Length.put(tmpBin, (long)curAlignment.node.getNumber("len"));
		}
				
		for(int i=1; i<stepRanges.size();i++){
			Range nextRanges = stepRanges.get(i);
			nextAlignment = allAlignments.get(nextRanges);
			tmpBin=SimpleBinner.getUniqueBin(nextAlignment.node);
			if(tmpBin!=null){				
//				System.out.print("Node " + nextAlignment.node.getId() + " has unique bin " + tmpBin.binID + ": ");
				if(tmpBin.isCloseTo(leastAbundancePop)){
//					System.out.println(" close to " + leastAbundancePop.binID);
					curBuildingBlocks.append(nextAlignment);
					allBuildingBlocks.add(curBuildingBlocks);
					
					curBuildingBlocks=new AlignedRead(nnpRead, nextAlignment);
					
					if(bins2Length.containsKey(tmpBin)){
						long newval=bins2Length.get(tmpBin)+(long)nextAlignment.node.getNumber("len");
						bins2Length.replace(tmpBin, newval);
					}else
						bins2Length.put(tmpBin, (long)nextAlignment.node.getNumber("len"));
				}else{
					//revert its bridges! (correct later when traverse through induced path)
//					System.out.println(" not close to " + leastAbundancePop.binID);
					curBuildingBlocks.append(nextAlignment);
					continue;
				}	
			}else
				curBuildingBlocks.append(nextAlignment);
				
				
			
		}
		if(curBuildingBlocks.getAlignmentRecords().size() > 1)
			allBuildingBlocks.add(curBuildingBlocks);
		
		
		//determine the global unique bin of the whole path
		long ltmp=0;
		for(PopBin b:bins2Length.keySet())
			if(bins2Length.get(b)>ltmp){
				ltmp=bins2Length.get(b);
				tmpBin=b;
			}
		
		ArrayList<BidirectedPath> retrievedPaths = new ArrayList<>();

		
		// Now we got all possible bridges from chopping the alignments at unique nodes
		System.out.println("\n=> bridges list: ");
		for(AlignedRead bb:allBuildingBlocks) {
			NewBridge storedBridge=getHomoBridgeFromMap(bb);
			System.out.printf("+++%s <=> %s\n", bb.getEndingsID(), storedBridge==null?"null":storedBridge.getEndingsID());
			
			if(storedBridge!=null) {
				if(storedBridge.isComplete()){
					System.out.println(storedBridge.getEndingsID() + ": already solved: ignore!");
					continue;
				}else{
					System.out.println(storedBridge.getEndingsID() + ": already built: fortify!");
					System.out.println(storedBridge.getAllPossiblePaths());
					if(storedBridge.buildFrom(bb))
						updateBridgesMap(storedBridge);

				}			

			}else{
				storedBridge=new NewBridge(this,bb,tmpBin);
				updateBridgesMap(storedBridge);
				
			}
			
			//storedBridge.bridging();
			if(storedBridge.isComplete()){
				retrievedPaths.add(storedBridge.getBestPath());
			}
			
		
		}
		return retrievedPaths;
	}
 	

    /**
     * Another reduce that doesn't remove the unique nodes
     * Instead redundant edges are removed on a path way
     * @param path: unique path to simplify the graph (from origGraph)
     */

    
    //This assuming path is surely unique!!!
    public boolean reduceUniquePath(BidirectedPath path){
    	//do nothing if the path has only one node
    	if(path.getEdgeCount()<2)
    		return false;
    	else
    		System.out.println("Reducing path: " + path.getId());
    	//loop over the edges of path (like spelling())
    	BidirectedNode 	startNode = (BidirectedNode) path.getRoot(),
    					endNode = (BidirectedNode) path.peekNode();
	
    	boolean startDir=((BidirectedEdge) path.getEdgePath().get(0)).getDir(startNode),
    			endDir=((BidirectedEdge) path.peekEdge()).getDir(endNode);

    	ArrayList<BidirectedEdge> 	potentialRemovedEdges = binner.walkAlongUniquePath(path);
		HashMap<PopBin, Integer> oneBin = new HashMap<>();
		oneBin.put(path.getConsensusUniqueBinOfPath(), 1);
    	
    	if(potentialRemovedEdges!=null && potentialRemovedEdges.size()>1){
	    	//remove appropriate edges
	    	for(BidirectedEdge e:potentialRemovedEdges){
	    		LOG.info("REMOVING EDGE " + e.getId() + " from " + e.getNode0().getGraph().getId() + "-" + e.getNode1().getGraph().getId());
	    		LOG.info("before: \n\t" + printEdgesOfNode((BidirectedNode) e.getNode0()) + "\n\t" + printEdgesOfNode((BidirectedNode) e.getNode1()));
	    		removeEdge(e.getId());
	    		LOG.info("after: \n\t" + printEdgesOfNode((BidirectedNode) e.getNode0()) + "\n\t" + printEdgesOfNode((BidirectedNode) e.getNode1()));
	    	}
	    	
	    	//add appropriate edges

    		BidirectedEdge reducedEdge = addEdge(startNode,endNode,startDir,endDir);
    		LOG.info("ADDING EDGE " + reducedEdge.getId()+ " from " + reducedEdge.getNode0().getGraph().getId() + "-" + reducedEdge.getNode1().getGraph().getId());
    		
			if(reducedEdge!=null){
//				reducedEdge.setAttribute("ui.label", path.getId());
//				reducedEdge.setAttribute("ui.style", "text-offset: -10; text-alignment: along;"); 
//				reducedEdge.setAttribute("isReducedEdge", true);
//				reducedEdge.setAttribute("ui.class", "marked");
//				reducedEdge.addAttribute("layout.weight", 10);
				reducedEdge.setAttribute("path", path);
				binner.edge2BinMap.put(reducedEdge, oneBin);
//				updateGraphMap(reducedEdge, path);
			}
			
	    	return true;
    	}else {
    		LOG.info("Path {} has failed reduce operation!", path.getId());
    		return false;
    	}

    }
    
    // old reducing. Now use for SPAdes paths only
    public boolean reduceFromSPAdesPath(BidirectedPath path){
    	//do nothing if the path has only one node
    	if(path==null || path.getEdgeCount()<1)
    		return false;
    	else
    		System.out.println("Input SPAdes path: " + path.getId());
    	//loop over the edges of path (like spelling())
    	BidirectedNode 	markerNode = null,
    			curNodeFromSimGraph = (BidirectedNode) path.getRoot();
	
    	boolean markerDir=true, curDir;
    	PopBin 	curUniqueBin = SimpleBinner.getUniqueBin(curNodeFromSimGraph);
    	if(curUniqueBin!=null){
    		markerNode=curNodeFromSimGraph;
    		markerDir=((BidirectedEdge) path.getEdgePath().get(0)).getDir(markerNode);
    	}
    	BidirectedPath curPath = new BidirectedPath();
		curPath.setRoot(curNodeFromSimGraph);

    	//search for an unique node as the marker. 
    	ArrayList<BidirectedEdge> 	tobeRemoved = new ArrayList<BidirectedEdge>(),
    								tobeAdded = new ArrayList<BidirectedEdge>();
    	for(Edge edge:path.getEdgePath()){
    		curPath.add(edge);	
    		curNodeFromSimGraph=(BidirectedNode) edge.getOpposite(curNodeFromSimGraph);
    		   		
    		curDir=((BidirectedEdge) edge).getDir(curNodeFromSimGraph);
    		
    		curUniqueBin = SimpleBinner.getUniqueBin(curNodeFromSimGraph);
    		if(curUniqueBin!=null){//only when reach the end of path
				curPath.setConsensusUniqueBinOfPath(curUniqueBin);
				if(markerNode!=null){
					//create an edge connect markerNode to curNode with curPath
					BidirectedEdge reducedEdge = new BidirectedEdge(markerNode, curNodeFromSimGraph, markerDir, curDir);
					reducedEdge.setAttribute("path", curPath);

					tobeAdded.add(reducedEdge);
					updateBridgesMap(curPath);
					
					ArrayList<BidirectedEdge> potentialRemovedEdges = binner.walkAlongUniquePath(curPath);
					if(potentialRemovedEdges!=null)
						tobeRemoved.addAll(potentialRemovedEdges);
					
					HashMap<PopBin, Integer> oneBin = new HashMap<>();
					oneBin.put(curUniqueBin, 1);
					binner.edge2BinMap.put(reducedEdge, oneBin);

				}else{
					if(curPath.size()>2)
						updateBridgesMap(curPath);
				}
				
				
				markerNode=curNodeFromSimGraph;
        		markerDir=!curDir; //in-out, out-in
				curPath= new BidirectedPath();
				curPath.setRoot(curNodeFromSimGraph);
				curPath.setConsensusUniqueBinOfPath(curUniqueBin);

				
    		}
    		
    		
		}
    	if(curPath.size() > 2)
			updateBridgesMap(curPath);

    	
    	if(tobeRemoved.size()>0){
	    	//remove appropriate edges
	    	for(BidirectedEdge e:tobeRemoved){
	    		LOG.info("REMOVING EDGE " + e.getId() + " from " + e.getNode0().getGraph().getId() + "-" + e.getNode1().getGraph().getId());
	    		LOG.info("before: \n\t" + printEdgesOfNode((BidirectedNode) e.getNode0()) + "\n\t" + printEdgesOfNode((BidirectedNode) e.getNode1()));
	    		removeEdge(e.getId());
	    		LOG.info("after: \n\t" + printEdgesOfNode((BidirectedNode) e.getNode0()) + "\n\t" + printEdgesOfNode((BidirectedNode) e.getNode1()));
	    	}
	    	
	    	//add appropriate edges
	    	for(BidirectedEdge e:tobeAdded){
	    		LOG.info("ADDING EDGE " + e.getId()+ " from " + e.getNode0().getGraph().getId() + "-" + e.getNode1().getGraph().getId());
	    		LOG.info("before: \n\t" + printEdgesOfNode((BidirectedNode) e.getNode0()) + "\n\t" + printEdgesOfNode((BidirectedNode) e.getNode1()));
	    		
	    		BidirectedEdge reducedEdge = addEdge((BidirectedNode)e.getSourceNode(),(BidirectedNode)e.getTargetNode(),e.getDir0(),e.getDir1());
	    		
				if(reducedEdge!=null){
	//				reducedEdge.addAttribute("ui.label", reducedEdge.getId());
	//				reducedEdge.setAttribute("ui.style", "text-offset: -10; text-alignment: along;"); 
//					reducedEdge.setAttribute("isReducedEdge", true);
//					reducedEdge.setAttribute("ui.class", "marked");
	//				reducedEdge.addAttribute("layout.weight", 10);
//					reducedEdge.setAttribute("path", path);
				}
	    		LOG.info("after: \n\t" + printEdgesOfNode((BidirectedNode) e.getNode0()) + "\n\t" + printEdgesOfNode((BidirectedNode) e.getNode1()));
	
	    	}
	    	return true;
    	}else
    		return false;

    }
    //return path in the graph that contain only unique nodes
    synchronized protected BidirectedPath getLongestLinearPathFromNode(BidirectedNode startNode, boolean direction){
//    	assert (direction?startNode.getOutDegree()<=1:startNode.getInDegree()<=1):" Node " + startNode.getId() + "has more than one possible extending way!";
    	BidirectedPath retval = new BidirectedPath();
    	retval.setRoot(startNode);
    	BidirectedNode currentNode = startNode;
    	boolean curDirection=direction;
    	while(curDirection?currentNode.getOutDegree()==1:currentNode.getInDegree()==1){
    		BidirectedEdge curEdge=curDirection?currentNode.leavingEdges().toArray(BidirectedEdge[]::new)[0]
    										:currentNode.enteringEdges().toArray(BidirectedEdge[]::new)[0];
    		retval.add(curEdge);
    		currentNode=(BidirectedNode) curEdge.getOpposite(currentNode);
    		curDirection=!curEdge.getDir(currentNode);
    	}
    	
    	return retval;
    }
}
