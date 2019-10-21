package org.rtassembly.npgraph;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Stack;
import java.util.TreeSet;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.graphstream.graph.*;
import org.graphstream.graph.implementations.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.util.concurrent.AtomicDouble;

import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;


public class BDGraph extends MultiGraph{
    private static final Logger LOG = LoggerFactory.getLogger(BDGraph.class);

    static int KMER=127;
	static boolean 	circular=true; // circular or linear preference for the output genome
    static double RCOV=0.0;
    SimpleBinner binner;
	//not gonna change these parameters in other thread
    public static final double R_TOL=.25;// relative tolerate: can be interpreted as long read error rate (10-25%)
    public static final int A_TOL=300;// absolute tolerate: can be interpreted as long read absolute error bases (100bp)

    public static final int GOOD_SUPPORT=20; //number of minimum spanning reads for an affirmative long-reads-based bridge.
	public static final double ALPHA=.5; //coverage less than alpha*bin_cov will be considered noise
    public static final int D_LIMIT=5000; //distance bigger than this will be ignored
    public static int S_LIMIT=300;// maximum number of graph traversing steps
    public static int T_LIMIT=21588;// maximum time allowed to run DFS (milliseconds)
    public static int MAX_LISTING=100; //maximum number of whatever paths
    
    //these should be changed in another thread, e.g. settings from GUI
	public static volatile double ILLUMINA_READ_LENGTH=300; //Illumina MiSeq
	public static volatile int MIN_SUPPORT=3; //minimal support reads for bridging

	
    //provide mapping from unique directed node to its corresponding bridge
    //E.g: 103-: <103-82-> also 82+:<82+103+>
    private HashMap<String, GoInBetweenBridge> bridgesMap; 
    ConsensusCaller consensus;

    //mapping long but unknown contigs to unique successors and predecessor 
    //with number of supported reads
    private static HashMap<String, Set<BDNodeState>> unknownBinMap=new HashMap<>();
    
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
	public BDGraph(String id, boolean strictChecking, boolean autoCreate,
			int initialNodeCapacity, int initialEdgeCapacity) {
		super(id, strictChecking, autoCreate);
		
		bridgesMap=new HashMap<String, GoInBetweenBridge>(initialNodeCapacity*2);
		consensus = new ConsensusCaller();
//		adjacencyMap=new HashMap<Node, Set<Node>>();
		// All we need to do is to change the node & edge factory
		setNodeFactory(new NodeFactory<BDNode>() {
			public BDNode newInstance(String id, Graph graph) {
				return new BDNode((AbstractGraph) graph, id);
			}
		});

		setEdgeFactory(new EdgeFactory<BDEdge>() {
			public BDEdge newInstance(String id, Node src, Node dst, boolean directed) { //stupid??
				return new BDEdge(id, (AbstractNode)src, (AbstractNode)dst);
			}
		});
		
	}
	public HashSet<GoInBetweenBridge> getUnsolvedBridges(){
		HashSet<GoInBetweenBridge> retval = new HashSet<GoInBetweenBridge>();
		for(GoInBetweenBridge brg:bridgesMap.values()){
			if(brg.getCompletionLevel()<4)
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
	public BDGraph(String id, boolean strictChecking, boolean autoCreate) {
		this(id, strictChecking, autoCreate, DEFAULT_NODE_CAPACITY,
				DEFAULT_EDGE_CAPACITY);
	}

	/**
	 * Creates an empty graph with strict checking and without auto-creation.
	 * 
	 * @param id
	 *            Unique identifier of the graph.
	 */
	public BDGraph(String id) {
		this(id, true, false, 10000, 10000);
	}
	
    public BDGraph(){
    	this("Assembly graph",true,false, 10000, 10000);
    }
	
	protected BDEdge addEdge(AbstractNode src, AbstractNode dst, boolean dir0, boolean dir1){
		BDEdge tmp = (BDEdge) addEdge(BDEdge.createID(src, dst, dir0, dir1), src, dst);
		return tmp;
	}
	
	@Override
	public Node removeNode(Node node){
//		System.out.println("Remove node " + node.getAttribute("name"));
		Node retval=super.removeNode(node);
		//clear bridge map
//		removeNodeFromBridgesMap(node);
		return retval;
	}
	
	public String printEdgesOfNode(BDNode node){
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

    public static int getKmerSize(){
    	return BDGraph.KMER;
    }
    public static void setKmerSize(int kmer){
    	BDGraph.KMER=kmer;
    }
        	
    
    //Several functionalities based on APSP
    // ...
    
    /**************************************************************************************************
     ********************** utility functions to serve the assembly algo ****************************** 
     * ***********************************************************************************************/   
    
    // when this unique node actually contained by a bridge
    synchronized public void removeNodeFromBridgesMap(Node unqNode){
    	bridgesMap.remove(unqNode.getId()+"o");
    	bridgesMap.remove(unqNode.getId()+"i");
    }
    // when there is new unique bridge 
    synchronized protected void updateBridgesMap(GoInBetweenBridge bidirectedBridge) {
    	if(bidirectedBridge==null || bidirectedBridge.pBridge==null)
    		return;

		if(bidirectedBridge.getNumberOfAnchors()==1)
    		bridgesMap.put(bidirectedBridge.pBridge.n0.toString(), bidirectedBridge);
    	
		if(bidirectedBridge.getNumberOfAnchors()==2){
			String 	end0=bidirectedBridge.pBridge.n0.toString(),
					end1=bidirectedBridge.pBridge.n1.toString();
			GoInBetweenBridge 	brg0=bridgesMap.get(end0),		
								brg1=bridgesMap.get(end1);
			if(brg0!=null&&brg0!=bidirectedBridge) 
				bidirectedBridge.merge(brg0,true);			
			
			if(brg1!=null&&brg1!=bidirectedBridge) 
				bidirectedBridge.merge(brg1,true);

			
			bridgesMap.put(end0, bidirectedBridge);
    		bridgesMap.put(end1, bidirectedBridge);

		}

    	
    }
    
    
    //Return bridge in the map (if any) that share the same bases (unique end) 
    synchronized public GoInBetweenBridge getBridgeFromMap(AlignedRead algRead){
    	GoInBetweenBridge retval = null, tmp = null;
    	if(algRead!=null){
	    	Node 	startNode=algRead.getFirstAlignment().node,
	    			endNode=algRead.getLastAlignment().node;
	    	boolean startNodeDir=algRead.getFirstAlignment().strand,
	    			endNodeDir=!algRead.getLastAlignment().strand;
	    	
	    	if(SimpleBinner.getBinIfUnique(startNode)!=null){
	    		tmp=bridgesMap.get(startNode.getId()+(startNodeDir?"o":"i"));
	    		if(tmp!=null)
	    			retval=tmp;
	    	}
	    	
	    	if(SimpleBinner.getBinIfUnique(endNode)!=null){
	    		tmp=bridgesMap.get(endNode.getId()+(endNodeDir?"o":"i"));
	    		if(tmp!=null){
	    			if(retval==null||retval.getCompletionLevel()<tmp.getCompletionLevel())
	    				retval=tmp;
	    			
	    		}
	    	}
	    	
    	}
    	
    	return retval;
    }
    
    synchronized public boolean isConflictBridge(BDEdgePrototype brg){
    	GoInBetweenBridge 	brg0=bridgesMap.get(brg.n0.toString()),
    						brg1=bridgesMap.get(brg.n1.toString());
    	return (brg0!=null&&brg0.getCompletionLevel()>=3) || (brg1!=null&&brg1.getCompletionLevel()>=3);
    }

    /*****************************************************************
     * Utility functions for real-time binning based on long reads
     * for unknown long contigs from initial binning step (SimpleBinner or metabat)
     *****************************************************************/
    //check if node is significant but not enough evidence to assign uniqueness
    public static boolean isSuspectedNode(BDNode node){
    	return unknownBinMap.containsKey(node.getId()+"o")&&unknownBinMap.containsKey(node.getId()+"i");
    }
    public void addUnknownNodes(BDNode node){
    	unknownBinMap.put(node.getId()+"o", new TreeSet<BDNodeState>());
    	unknownBinMap.put(node.getId()+"i", new TreeSet<BDNodeState>());
    }
    
    //add successors&predecessors of unknown nodes based on AlignedRead
    //read must have either start or end as unique and no more
    public void addReadsToUnknowMap(AlignedRead read){
    	if(read.getEFlag() < 1)
    		return;
    	Alignment 	a0=read.getFirstAlignment(),
    				a1=read.getLastAlignment();
    	BDNodeState ns0=null,ns1=null;
    	if(SimpleBinner.getBinIfUnique(a0.node)!=null)
    		ns0=new BDNodeState(a0.node, a0.strand,1);
    	if(SimpleBinner.getBinIfUnique(a1.node)!=null)
    		ns1=new BDNodeState(a1.node, !a1.strand,1);
    	
    	for(Alignment alg:read.getAlignmentRecords()){
    		if(SimpleBinner.getBinIfUnique(alg.node)!=null)
    			continue;
    		String key;
    		if(ns0!=null){
    			key=alg.node.getId()+(alg.strand?"i":"o");
    			Set<BDNodeState> values=unknownBinMap.get(key);
    			if(values!=null){
    				Iterator<BDNodeState> iterator=values.iterator();
    				boolean found=false;
    				while(iterator.hasNext()){
    					BDNodeState ns=iterator.next();
    					if(ns.equals(ns0)){
    						ns.weight++;
    						found=true;
    						break;
    					}
    				}
    				if(!found)
    					values.add(ns0);
    			}
    		}
    		
    		if(ns1!=null){
      			key=alg.node.getId()+(alg.strand?"o":"i");
    			Set<BDNodeState> values=unknownBinMap.get(key);
    			if(values!=null){
    				Iterator<BDNodeState> iterator=values.iterator();
    				boolean found=false;
    				while(iterator.hasNext()){
    					BDNodeState ns=iterator.next();
    					if(ns.equals(ns1)){
    						ns.weight++;
    						found=true;
    						break;
    					}
    				}
    				if(!found)
    					values.add(ns1);
    			}
    		}
    		
    	}
    	
    }
    
    //get multiplicity after a while (20X???)...
    //NOPE: can only know if unique or repetitive based on the previous&next unique nodes
    synchronized public static PopBin getUniqueBinFromLongReads(BDNode node){
    	String 	ko=node.getId()+"o",
    			ki=node.getId()+"i";
    	Set<BDNodeState> 	successors=unknownBinMap.get(ko),
							predecessors=unknownBinMap.get(ki);
    	int c=0, co=0, ci=0;

    	if(successors==null && predecessors==null)
    		return null;
    	
    	//FIXME: checking stats from GraphWatcher instead???
    	if(successors!=null)
    		c+=successors.stream().mapToInt(ns->ns.getWeight()).sum();
    	if(predecessors!=null)
    		c+=predecessors.stream().mapToInt(ns->ns.getWeight()).sum();
    	
    	if(c<GOOD_SUPPORT)//just need a number
    		return null;
    	
    	//now validate successors and predecessors based on number of support reads
    	PopBin retval=null;

    	for(BDNodeState ns:successors){
    		if(ns.getWeight()<.1*BDGraph.GOOD_SUPPORT)//less than 10% will be ignored
    			continue;
    		co++;
    		if(co>1)
    			return null;
    		retval=SimpleBinner.getBinIfUnique(ns.getNode());//check??
    	} 
    	for(BDNodeState ns:predecessors){
    		if(ns.getWeight()<.1*BDGraph.GOOD_SUPPORT)
    			continue;
    		ci++;
    		if(ci>1)
    			return null;
    		if(retval!=null&&!PopBin.isCloseBins(retval, SimpleBinner.getBinIfUnique(ns.getNode())))
				return null;
    	} 
    	
    	//reove this from unknowmap
    	unknownBinMap.remove(ko);
    	unknownBinMap.remove(ki);
    	node.setAttribute("unique", retval);
    	if(HybridAssembler.VERBOSE)
    		LOG.info("FOUND NEW UNIQUE CONTIG BY LONG READS: ID=" + node.getId() + " degree= " + node.getDegree() + " out=" + co + " in=" + ci + " read count="+c);
    	return retval;
    }
    
    //If the node was wrongly identified as unique before, do things...
//    synchronized public void destroyFalseBridges(Node node){
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
//    			falseBridgeFrom.fullPaths=new ArrayList<BDPath>();
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
//    			falseBridgeTo.fullPaths=new ArrayList<BDPath>();
//    			//todo: remove 2 entries from bridgeMap here...
//
//    		}else{    		
//    			//destroy the half bridge
//        		bridgesMap.remove(node.getId()+ "i");
//    		}
//    	}
//    		
//    }
    
    synchronized public void binning(String binFileName) {   		
    	binner=new SimpleBinner(this, binFileName);
    	binner.estimatePathsByCoverage();
    	initGraphComponents();
    }

    //Scanning for shorter overlaps (<k) in a DBG graph. Not making difference??! 
    synchronized public void fixDeadEnds(){
    	List<BDNodeState> weirdNodes = new ArrayList<>();
    	for(Node node:this){

    		if(node.getNumber("len") < SimpleBinner.UNIQUE_CTG_LEN)
    			continue;
    		
    		double 	inCov=node.enteringEdges().map(e->e.getOpposite(node).getNumber("cov")).mapToDouble(Double::doubleValue).sum(),
    				outCov=node.leavingEdges().map(e->e.getOpposite(node).getNumber("cov")).mapToDouble(Double::doubleValue).sum();
    		double compare=GraphUtil.approxCompare(inCov, outCov);
    		if(inCov*outCov==0){
    			if(inCov==0)
        			weirdNodes.add(new BDNodeState((BDNode) node, false));
    			if(outCov==0)
        			weirdNodes.add(new BDNodeState((BDNode) node, true));
    		}else if(compare > 0)
    			weirdNodes.add(new BDNodeState((BDNode) node, true));
    		else if(compare < 0)
    			weirdNodes.add(new BDNodeState((BDNode) node, false));
			else
				continue;
    		
    		if(HybridAssembler.VERBOSE)
    			LOG.info("Found one weird node {} inCov={}/{}, outCov={}/{}",(String)node.getAttribute("name"), inCov, node.getInDegree(), outCov, node.getOutDegree());

    	}
    	BDNodeState n0,n1;
    	Sequence seq0,seq1;
    	for(int i=0;i<weirdNodes.size();i++){
    		n0 = weirdNodes.get(i);
    		seq0=(Sequence)(n0.getNode().getAttribute("seq"));
    		if(!n0.getDir())
    			seq0=Alphabet.DNA.complement(seq0);
    		
    		for(int j=i+1; j<weirdNodes.size();j++){
    			n1 = weirdNodes.get(j);
        		seq1=(Sequence)(n1.getNode().getAttribute("seq"));
        		if(n1.getDir())
        			seq1=Alphabet.DNA.complement(seq0);
        		
        		int overlap=GraphUtil.overlap(seq0,seq1);
        		if(overlap > 55){
        			BDEdge e=addEdge(n0.getNode(), n1.getNode(), n0.getDir(), n1.getDir());
        			if(HybridAssembler.VERBOSE)
            			LOG.info("...adding potential edge {} length={}", e.getId(), overlap);
        		}
    		}
    	}
    }
	
    
    //Depth First Search strategy
	synchronized ArrayList<BDPath> DFSAllPaths(BDNode srcNode, BDNode dstNode, boolean srcDir, boolean dstDir, int distance)
	{
    	if(HybridAssembler.VERBOSE)
    		LOG.info("Looking for DFS path between {}{} to {}{} with distance={}",srcNode.getId(), srcDir?"o":"i", dstNode.getId(), dstDir?"o":"i" ,distance);
		ArrayList<BDPath> possiblePaths = new ArrayList<BDPath>(), 
									retval=new ArrayList<BDPath>();

		//1. First build shortest tree from dstNode 		
		HashMap<String,Integer> shortestMap = getShortestTreeFromNode(dstNode, dstDir, distance);
		// The good thing is that we need only 1 temporary path variable
		BDPath path = new BDPath(srcNode);

		//2. DFS from srcNode with the distance info above
		BDNodeState curNodeState = new BDNodeState(srcNode, !srcDir); // first node is special
		if(shortestMap.containsKey(curNodeState.toString())) {
			
			Stack<List<Edge>> stack = new Stack<>();
			final List<Edge> tmpList = new ArrayList<>(); 
			List<Edge> curList = null;
			int tolerance = A_TOL, 
					delta;
			BDEdge curEdge = null;

			if(HybridAssembler.VERBOSE)
	    		LOG.info("Found " + curNodeState.toString() + " with shortest distance=" + shortestMap.get(curNodeState.toString()));

			AtomicDouble limit = new AtomicDouble(distance+tolerance);
			
			(curNodeState.getDir()?curNodeState.getNode().enteringEdges():curNodeState.getNode().leavingEdges())
				.forEach(e->{
					BDNode n=(BDNode) e.getOpposite(srcNode);

					BDNodeState ns = new BDNodeState(n, ((BDEdge) e).getNodeDirection(n)!=null?
														((BDEdge) e).getNodeDirection(n):!srcDir);
					
    				if(	shortestMap.containsKey(ns.toString()) 
						&& shortestMap.get(ns.toString()) < limit.get()
						)
	    					tmpList.add(e);
	    				
					
				});
		
			stack.push(new ArrayList<>(tmpList));
			Stack<Double> nodeScores = new Stack<>();
			double pathScore=0.0;
			long startTime = System.currentTimeMillis();
			while(true) {
				curList=stack.peek();
				if(curList.isEmpty()) {
					if(path.size() <= 1)
						break;
					stack.pop();
					distance += (int)path.peekNode().getNumber("len") + ((BDEdge) path.popEdge()).getLength();
					pathScore-=nodeScores.pop();
					
				}else {
					curEdge=(BDEdge) curList.remove(0);
					BDNode 	from = (BDNode) path.peekNode(),
							to = (BDNode) curEdge.getOpposite(from);
					
					//Important: an anchor is not allowed in the result path
					if(SimpleBinner.getBinIfUnique(to)!=null && to!=dstNode)
						continue;
					
					boolean dir = 	curEdge.getNodeDirection(to)!=null?
									curEdge.getNodeDirection(to):
									path.getLastNodeDirection();
					
					/*
					 * Update path and if terminated condition is met, add the candidate path
					 */
					nodeScores.push(path.getExtendLikelihood(to));
					pathScore+=nodeScores.peek();
					path.add(curEdge);
					
					delta=distance-curEdge.getLength();
					//note that traversing direction (true: template, false: reverse complement) of destination node 
					//is opposite its defined direction (true: outward, false:inward) 
					if(to==dstNode && dir==dstDir && Math.abs(delta) < tolerance){ 

				    	BDPath 	tmpPath=new BDPath(path);
				    	tmpPath.setDeviation(delta);
				    	tmpPath.setPathEstats(pathScore);

				    	//insert to the list with sorting
				    	if(possiblePaths.isEmpty())
				    		possiblePaths.add(tmpPath);
				    	else{
				    		int idx=0;
				    		for(BDPath p:possiblePaths){
				    			if(Math.abs(delta)>Math.abs(p.getDeviation()))
				    				idx++;
				    			else
				    				break;
				    		}
				    		possiblePaths.add(idx,tmpPath);
				    	}
						
						if(possiblePaths.size() > S_LIMIT || System.currentTimeMillis() > startTime+T_LIMIT) //not go too far
							break;
					}
					/*
					 * Looking for next candidate set of edges to traverse
					 */
					limit.set(distance + tolerance);
			    	//get possible next edges to traverse
					tmpList.clear(); 

	    			(dir?to.enteringEdges():to.leavingEdges())
	    			.forEach(e->{
	    				BDNode n=(BDNode) e.getOpposite(to);
	    				BDNodeState ns = new BDNodeState(n, ((BDEdge) e).getNodeDirection(n)!=null?
															((BDEdge) e).getNodeDirection(n):dir);
	    				
	    				if(shortestMap.containsKey(ns.toString()) 
    						&& shortestMap.get(ns.toString()) < limit.get()
							)
	    	    					tmpList.add(e);
	    				
	    			});
//	    			
			    	stack.push(new ArrayList<>(tmpList));
			    						
					distance -= (int)to.getNumber("len") + curEdge.getLength();
					
				}
				
			}
		} 
		if(HybridAssembler.VERBOSE)
			LOG.info("select from list of " + possiblePaths.size() + " DFS paths:");
		
		if(possiblePaths.isEmpty())
			return null;
		
		double closestDist=Math.abs(possiblePaths.get(0).getDeviation());
		int keepMax = MAX_LISTING;//only keep this many possible paths 
		for(int i=0;i<possiblePaths.size();i++){
			BDPath p = possiblePaths.get(i);
			if(Math.abs(p.getDeviation())>closestDist+Math.abs(distance+getKmerSize())*R_TOL || i>=keepMax)
				break;
			retval.add(p);
			if(HybridAssembler.VERBOSE)
	    		LOG.info("Hit added: {} deviation={}; depth={}; likelihood score={}", p.getId(), p.getDeviation(), p.size(), p.getPathEstats());

		}
		
		return retval;
	}    
	
	
    /*
     * Get a map showing shortest distances from surrounding nodes to a *rootNode* expanding to a *direction*, within a *distance*
     * based on Dijkstra algorithm
     */
    public static HashMap<String, Integer> getShortestTreeFromNode(BDNode rootNode, boolean expDir, int distance){
		PriorityQueue<BDNodeState> pq = new PriorityQueue<>();
		HashMap<String,Integer> retval = new HashMap<>();
		
		int curDistance=(int) -rootNode.getNumber("len"), newDistance=0;
		BDNodeState curND = new BDNodeState(rootNode, expDir, curDistance);
		pq.add(curND);
		
//		if(HybridAssembler.VERBOSE)
//    		LOG.info("Building shortest tree for " + rootNode.getId() + " with distance=" + distance);

		retval.put(curND.toString(), curDistance); // direction from the point of srcNode
		
		boolean direction=expDir;
		while(!pq.isEmpty()) {
			curND=pq.poll();
			curDistance=curND.getWeight();

			Iterator<Edge> ite=curND.getDir()?curND.getNode().leavingEdges().iterator():curND.getNode().enteringEdges().iterator();
			while(ite.hasNext()) {
	    		BDEdge edge = (BDEdge) ite.next();
	    		BDNode nextNode = (BDNode) edge.getOpposite(curND.getNode());
	    		
	    		if(edge.getNodeDirection(nextNode)!=null)	    			
	    			direction = !edge.getNodeDirection(nextNode);
	    		
	    		newDistance=curDistance+edge.getLength()+(int)curND.getNode().getNumber("len");
    			if(newDistance-distance > BDGraph.A_TOL && GraphUtil.approxCompare(newDistance, distance)>0)
	    			continue;
	    		
	    		BDNodeState nextND = new BDNodeState(nextNode, direction, newDistance);
	    		String key=nextND.toString();
	    		if(retval.containsKey(key)) {
	    			if(retval.get(key) > newDistance) {
	    				//update the map
	    				retval.put(key, newDistance);
	    				//update the queue
	    				pq.remove(nextND);
	    				pq.add(nextND);
	    			}
	    		}else {
	    			retval.put(key, newDistance);
    				pq.add(nextND);
	    		}

			}
			
		}
		return retval;
    }
    
    
    /*
     * Find bridges based on list of Alignments.
     * Return list of bridges with endings as markers and alignments of non-markers in-between.
     */ 
    synchronized protected List<BDPath> uniqueBridgesFinding(Sequence nnpRead, ArrayList<Alignment> alignments) {
 		if(nnpRead==null || alignments.size()<=1)
 			return null;
 		
 		if(HybridAssembler.VERBOSE) {
	 		LOG.info("=================================================");
	 		for(Alignment alg:alignments)
	 			LOG.info("\t"+alg.toString());
	 		LOG.info("=================================================");
 		}
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

 		ArrayList<BDPath> retrievedPaths = new ArrayList<>();

 		List<Range> curRanges=rangeGroups.get(0);
 		AlignedRead	curBuildingBlocks; //building blocks for a bridge, taken from alignments with unique end(s)
 		PopBin curBin=null, prevUnqBin=null;


 		int flag=0; //+1 for unique startAlignment.node, +2 for unique endAlignment.node flag={0,1,2,3}
 		if(curRanges.size()==1){
 			curBin=SimpleBinner.getBinIfUnique(allAlignments.get(curRanges.get(0)).node);
 			if( curBin!= null){
 				flag+=1;
 				prevUnqBin=curBin;
 			}
 			curBuildingBlocks=new AlignedRead(nnpRead, allAlignments.get(curRanges.get(0)));
 		}else{
 			curBuildingBlocks=new AlignedRead(nnpRead, (ArrayList<Alignment>) curRanges.stream().map(g->allAlignments.get(g)).collect(Collectors.toList()) );
 		
 		}
 		
 		if(HybridAssembler.VERBOSE) {
 			LOG.info("Step ranges: ");
 			String log="";
	    	for(Range range:curRanges) { 
	    		log+=(allAlignments.get(range).node.getId() + " "+ binner.getBinsOfNode(allAlignments.get(range).node) + ": " + range + "; ");	  
	    	}
	    		LOG.info(log);
 		}
 		
 		Alignment curAlg=null;
 	    for(int i=1; i<rangeGroups.size(); i++){
 	    	curRanges = rangeGroups.get(i);
 	    	
 	 		if(HybridAssembler.VERBOSE) {
 	 			String log="";
 		    	for(Range range:curRanges) { 
 		    		log+=(allAlignments.get(range).node.getId() + " "+ binner.getBinsOfNode(allAlignments.get(range).node) + ": " + range + "; ");	  
 		    	}
 		    		LOG.info(log);
 	 		}
 	      	
 	    	if(curRanges.size()==1){
 	    		curAlg=allAlignments.get(curRanges.get(0));
 	    		curBin=SimpleBinner.getBinIfUnique(curAlg.node);	    		
    			curBuildingBlocks.append(curAlg);	    		
 	    	}else{
 	    		curBin=null;
 	    		curBuildingBlocks.appendAll((ArrayList<Alignment>) curRanges.stream().map(g->allAlignments.get(g)).collect(Collectors.toList()));
 	    	}
 	    	
  			
    		if(curBin!=null){
    			if(PopBin.isCloseBins(curBin, prevUnqBin))
    				flag+=2;
    			else if(prevUnqBin!=null)
    				continue;
    			
    			curBuildingBlocks.setEFlag(flag);
    			////////////////////////////////////////////////////////////////////////////////////
    			
    			retrievedPaths.addAll(buildBridge(curBuildingBlocks, PopBin.getDominateBin(curBin,prevUnqBin)));
    			   					
    			////////////////////////////////////////////////////////////////////////////////////
    			//start new building block
				curBuildingBlocks=new AlignedRead(nnpRead, curAlg);
				flag=1;
				prevUnqBin=curBin;
    		}
	    		
 	    }
 	    
 	    //last time iteration: shouldn't be any new paths found, just merging information
 	    
 	    if(curBuildingBlocks.alignments.size() > 1){
 	    	curBuildingBlocks.setEFlag(1);
 	    	retrievedPaths.addAll(buildBridge(curBuildingBlocks, PopBin.getDominateBin(curBin,prevUnqBin)));
 	    }

 	    return retrievedPaths;
 	}
  	
    
    private List<BDPath> buildBridge(AlignedRead read, PopBin bin){
    	List<BDPath> retval=new ArrayList<BDPath>();

		GoInBetweenBridge 	storedBridge=getBridgeFromMap(read), reversedBridge;
		if(HybridAssembler.VERBOSE) 
			LOG.info("+++{} <=> {}\n", read.getEndingsID(), storedBridge==null?"null":storedBridge.getEndingsID());
		//long read give info about unique successor and predecessor
		addReadsToUnknowMap(read);
		if(storedBridge!=null) {
			if(storedBridge.getCompletionLevel()==4){//important since it will ignore the wrong transformed unique nodes here!!!
				if(HybridAssembler.VERBOSE) 
					LOG.info(storedBridge.getEndingsID() + ": already solved and reduced: ignore!");
				return retval;

			}else{
				if(HybridAssembler.VERBOSE) 
					LOG.info(storedBridge.getEndingsID() + ": already built: fortify!");
				//update available bridge using alignments
				byte state=storedBridge.merge(read,true);
				
				if((state&0b10)>0)//number of anchors has changed after merging
					updateBridgesMap(storedBridge);
				
				//scan for transformed unique nodes
				boolean extend=false;
				if(storedBridge.getCompletionLevel()==1){
					if(storedBridge.scanForAnEnd(false)){
						if(HybridAssembler.VERBOSE) 
							LOG.info("FOUND NEW TRANSFORMED END: " + storedBridge.steps.end.getNode().getId());
						extend=storedBridge.steps.connectBridgeSteps(false);
					}					
				}
				if(storedBridge.getCompletionLevel()==4 || (state&0b01)>0 || extend)
					retval.addAll(storedBridge.scanForNewUniquePaths());
					
//				//also update the reversed bridge: important e.g. Shigella_dysenteriae_Sd197. WHY??? (already updated and merged 2 homo bridges)
				if(read.getEFlag()==3) {
					read.reverse();
					reversedBridge = getBridgeFromMap(read);
					if(reversedBridge!=storedBridge) {
						byte anotherState=reversedBridge.merge(read, true);
						if((anotherState&0b10)>0)//number of anchors has changed after merging
							updateBridgesMap(reversedBridge);											
						
						if(reversedBridge.getCompletionLevel()==4 || (anotherState&0b01)>0)
							retval.addAll(reversedBridge.scanForNewUniquePaths());

					}

				}

			}			
			

		}else{
			storedBridge=new GoInBetweenBridge(this,read, bin);
			updateBridgesMap(storedBridge);		
			
			read.reverse();
			reversedBridge = new GoInBetweenBridge(this,read, bin);
			updateBridgesMap(reversedBridge);
			
		}
		consensus.saveBridgingReadsFromAlignments(read);
		
		return retval;
    }

    /**
     * Another reduce that doesn't remove the unique nodes
     * Instead redundant edges are removed on a path way
     * @param path: unique path to simplify the graph (from origGraph)
     */
    //This assuming path is surely unique!!!
    public boolean reduceUniquePath(BDPath path){
    	//do nothing if the path has only one node
    	if(path==null||path.getEdgeCount()<1)
    		return false;

    	else if(HybridAssembler.VERBOSE) 
    			LOG.info("Reducing path: " + path.getId());


    	Set<Edge> 	potentialRemovedEdges = binner.walkAlongUniquePath(path);
		Multiplicity oneBin = new Multiplicity(path.getConsensusUniqueBinOfPath(), 1);
    	
    	if(potentialRemovedEdges!=null && potentialRemovedEdges.size()>1){
	    	//remove appropriate edges
    		
	    	for(Edge e:potentialRemovedEdges){
	    		if(HybridAssembler.VERBOSE) {
		    		LOG.info("REMOVING EDGE " + e.getId() + " from " + e.getNode0().getGraph().getId() + "-" + e.getNode1().getGraph().getId());
		    		LOG.info("before: \n\t" + printEdgesOfNode((BDNode) e.getNode0()) + "\n\t" + printEdgesOfNode((BDNode) e.getNode1()));
	    		}
	    		removeEdge(e.getId());
	    		if(HybridAssembler.VERBOSE)
	    			LOG.info("after: \n\t" + printEdgesOfNode((BDNode) e.getNode0()) + "\n\t" + printEdgesOfNode((BDNode) e.getNode1()));

	    	}
	    	
	    	//add appropriate edges
			BDEdge reducedEdge = addEdge(path.getFirstNode(), path.getLastNode(), path.getFirstNodeDirection(), path.getLastNodeDirection());
    		if(HybridAssembler.VERBOSE)
    			LOG.info("ADDING EDGE " + reducedEdge.getId()+ " from " + reducedEdge.getNode0().getGraph().getId() + "-" + reducedEdge.getNode1().getGraph().getId());

			if(reducedEdge!=null){
				if(path.getPrimitivePath().getEdgeCount()>1)
					reducedEdge.setAttribute("path", path);
				binner.edge2BinMap.put(reducedEdge, oneBin);
			}
			if(HybridAssembler.VERBOSE)
				LOG.info("after adding: \n\t" + printEdgesOfNode((BDNode) reducedEdge.getNode0()) + "\n\t" + printEdgesOfNode((BDNode) reducedEdge.getNode1()));
	    	return true;
    	}else {
    		if(HybridAssembler.VERBOSE)
    			LOG.info("Path {} doesn't need to be reduced!", path.getId());
    		return false;
    	}

    }
    
    // old reducing. Now use for SPAdes paths only
    public boolean reduceFromSPAdesPath(BDPath path){
    	boolean retval=false;
    	//do nothing if the path has only one node
    	if(path==null || path.getEdgeCount()<1)
    		return retval;
    	else if(HybridAssembler.VERBOSE)
    		LOG.info("Input SPAdes path: " + path.getId());
    	//loop over the edges of path (like spelling())
    	for(BDPath p:getNewSubPathsToReduce(path)){
    		if(reduceUniquePath(p))
    			retval=true;
    	}
    	return retval;

    }

    
	//Only call for the final reduce path with 2 unique ends: if path containing other unique nodes than 2 ends then we have list of paths to reduce
    //Exclude already-reduced path (reduced edges would have composite "path" attribute,
    //but consensus pseudo edge is special case)
	synchronized public ArrayList<BDPath> getNewSubPathsToReduce(BDPath path){
		ArrayList<BDPath> retval=new ArrayList<>();
		if(path!=null){
			BDPath curPath = new BDPath(path.getRoot(), path.getConsensusUniqueBinOfPath());
			BDNode curNode = (BDNode) path.getRoot(), nextNode=null;
			String id=null;
			for(Edge e:path.getEdgePath()) {
				nextNode=(BDNode) e.getOpposite(curNode);
				curPath.add(e);
				if(SimpleBinner.getBinIfUniqueNow(nextNode)!=null) {
					id=curPath.getEndingID();
					if(id!=null){
						Edge rdEdge = getEdge(id);
						if(rdEdge==null || !rdEdge.hasAttribute("path") || rdEdge.hasAttribute("consensus"))
							retval.add(curPath);		
					}
					curPath=new BDPath(nextNode, path.getConsensusUniqueBinOfPath());
				}
				curNode=nextNode;
			}
			//return intact if no extra unique node has been found
			if(retval.isEmpty()) {
				id=curPath.getEndingID();
				if(id!=null&&getEdge(id)==null){
					retval.add(curPath);
				}
			}
		}
		return retval;
	}
	
	
    /***********************************************************************
     * For graph status reporting
     * Compliance with with GraphWatcher.update()
     **********************************************************************/
    int n50, n75, maxl; //in Kbp
    int numOfCtgs, numOfCircularCtgs;
    public synchronized void updateStats() {
    	numOfCtgs=getNodeCount();
		int count=0;
		numOfCircularCtgs=0;
		maxl=0;
		int [] lengths = new int[numOfCtgs];
		double sum = 0;		
    	for (Node node : this) {
    		/*
    		 * Re-calculate stats
    		 */
			int nlen=(int)node.getNumber("len"); 
			if(maxl < nlen)
				maxl=nlen;
			if(node.hasAttribute("circular"))
				numOfCircularCtgs++;
			
			lengths[count++]=nlen; 
			sum+=nlen;	
    	}
    	
    	/*
    	 * Calculate N50, N75
    	 */
		Arrays.sort(lengths);

		int i50 = lengths.length,
			i75 = lengths.length;
		double contains = 0;
		while (true){
			if(contains < .5*sum)
				i50 --;
			if(contains < .75*sum) {
				i75--;
			}else
				break;
			contains += lengths[i75];
		}
		n50=lengths[i50];
		n75=lengths[i75];
    }
    
    private void initGraphComponents() {	
    	List<String> pallette=getUniqueColors(binner.binList.size());

    	HashMap<Integer, String> bin2color = new HashMap<Integer, String>();
    	for(int i=0;i<pallette.size();i++){
    		bin2color.put(binner.binList.get(i).getId(), pallette.get(i));
    	}
    	
    	for (Node node : this) {		
			String color="rgb(255,255,255)";
			PopBin bin=SimpleBinner.getBinIfUnique(node);
			if(bin!=null){
				color=bin2color.get(bin.getId());
			}
			
			((BDNode)node).setGUI(color, "circle");
    	}

    }
	private List<String> getUniqueColors(int amount) {
        final int lowerLimit = 0x10;
        final int upperLimit = 0xE0;    
        final int colorStep = (int) ((upperLimit-lowerLimit)/Math.pow(amount,1f/3));

        final List<String> colors = new ArrayList<String>(amount);

        for (int R = lowerLimit;R <= upperLimit; R+=colorStep)
            for (int G = lowerLimit;G <= upperLimit; G+=colorStep)
                for (int B = lowerLimit;B <= upperLimit; B+=colorStep) {
                	if(R==G && G==B){
                    	continue;
                    }
                    if (colors.size() >= amount) { //The calculated step is not very precise, so this safeguard is appropriate
                        return colors;
                    } 
                    else {
                        colors.add("rgb("+R+","+G+","+B+")");
                    }               
                }
        return colors;
    }
	public void outputFASTA(String fileName) throws IOException {
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(fileName);
		for(Node node:this) {
			Sequence seq=(Sequence) node.getAttribute("seq");
//			if( (node.getDegree()==0 && (seq.length() < SimpleBinner.ANCHOR_CTG_LEN)) 
//				|| node.getNumber("cov") < 10.0 )	//not display <10% abundance pops
//				continue;
			seq.writeFasta(out);
		}
		
		out.close();
	}
	public void outputJAPSA(String fileName) throws IOException {
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(fileName);
		JapsaAnnotation annotation;
		for(Node node:this) {
			annotation=(JapsaAnnotation) node.getAttribute("annotation");
			annotation.write(out);
		}	
		out.close();
	}
	
	public void outputGFA(String fileName) throws IOException {
	    PrintWriter printWriter = new PrintWriter(new FileWriter(fileName));
	    
	    //Differentiate composite edges and normal edges
	    List<Edge> compositeEdges = new ArrayList<Edge>();
	    List<BDEdgePrototype> normalEdges = new ArrayList<BDEdgePrototype>();
	    edges().forEach(e->{
	    					if(e.hasAttribute("path")) 
	    						compositeEdges.add(e); 
	    					else 
	    						normalEdges.add(new BDEdgePrototype((BDNode)e.getNode0(), (BDNode)e.getNode1(), ((BDEdge)e).getDir0(), ((BDEdge)e).getDir1()));
	    					}
	    );
	    
	    Set<String> addedNodes = new HashSet<String>();
	    //Print S (Segment)
	    for(Node node:this){
	    	Sequence seq=(Sequence) node.getAttribute("seq");
	    	int kmer_count=(int)(GraphUtil.getRealCoverage(node.getNumber("cov"))*(BDGraph.ILLUMINA_READ_LENGTH-BDGraph.getKmerSize())/BDGraph.ILLUMINA_READ_LENGTH);
	    	printWriter.printf("S\t%s\t%s\tKC:i:%d\n", node.getId(), seq.toString(),kmer_count);
	    	addedNodes.add(node.getId());
	    }	    
	    for(Edge ce:compositeEdges){
	    	BDPath p = ((BDPath)ce.getAttribute("path")).getPrimitivePath();

	    	BDNode curNode=(BDNode) p.getRoot(), nextNode=null;
	    	String curID=curNode.getId(), nextID=null;
	    	boolean curDir=p.getFirstNodeDirection(), nextDir=!curDir;

	    	for(Edge e:p.getEdgePath()){
	    		nextNode=(BDNode) e.getOpposite(curNode);
	    		nextID=nextNode.getId();
	    		if(((BDEdge)e).getNodeDirection(nextNode)!=null)
	    			nextDir=((BDEdge)e).getNodeDirection(nextNode);
	    		
	    		if(nextNode!=p.peekNode()){
		    		//create ID
		    		int count=1;
		    		String tmpID=nextID;
		    		while(addedNodes.contains(tmpID)){
		    			tmpID=nextID+"."+(count++);
		    		}
		    		nextID=tmpID;
		    		addedNodes.add(nextID);
		    		
		    		Sequence seq=(Sequence) nextNode.getAttribute("seq");
			    	int kmer_count=(int)(GraphUtil.getRealCoverage(nextNode.getNumber("cov"))*(BDGraph.ILLUMINA_READ_LENGTH-BDGraph.getKmerSize())/BDGraph.ILLUMINA_READ_LENGTH);
			    	printWriter.printf("S\t%s\t%s\tKC:i:%d\n", nextID, seq.toString(),kmer_count);

	    		}
		    	
	    		normalEdges.add(new BDEdgePrototype(new BDNode(this, curID), new BDNode(this, nextID), curDir, nextDir));
	    		
	    		curNode=nextNode;
	    		curID=nextID;
	    		curDir=!nextDir;
	    	}
	    }
	    
	    //Print L (Links)
	    for(BDEdgePrototype e:normalEdges){
	    	printWriter.printf("L\t%s\t%s\t%s\t%s\t%dM\n", 
	    						e.getNode0().getId(), 
	    						e.getDir0()?"+":"-",
								e.getNode1().getId(), 
	    						e.getDir1()?"-":"+",
								BDGraph.getKmerSize());
	    }

	    
	    printWriter.close();
	}
    
}
