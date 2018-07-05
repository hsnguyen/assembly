package org.rtassembly.npgraph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.graphstream.graph.*;
import org.graphstream.graph.implementations.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;


public class BidirectedGraph extends MultiGraph{
    static int KMER=127;
    static double RCOV=0.0;
	SimpleBinner binner;

    volatile static double ILLUMINA_READ_LENGTH=300; //Illumina MiSeq
    
    static final double TOLERATE=.3;//can be interpreted as long read error rate (10-25%)
    static final int D_LIMIT=5000; //distance bigger than this will be ignored
    static final int S_LIMIT=15;// maximum number of DFS steps
    
    //provide mapping from unique directed node to its corresponding bridge
    //E.g: 103-: <103-82-> also 82+:<82+103+>
    private HashMap<String, BidirectedBridge> bridgesMap; 
    
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
		bridgesMap=new HashMap<String, BidirectedBridge>(initialNodeCapacity*2);
		
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
	public HashSet<BidirectedBridge> getUnsolvedBridges(){
		HashSet<BidirectedBridge> retval = new HashSet<BidirectedBridge>();
		for(BidirectedBridge brg:bridgesMap.values()){
			System.out.printf("Bridge %s : status=%d \n", brg.getBridgeString(), brg.getBridgeStatus());
			if(brg.getBridgeStatus()<1)
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
    
    synchronized public BidirectedBridge getBridgeFromMap(Node unqNode, boolean direction){
    	String key=unqNode.getId()+(direction?"o":"i"); //true:going outward, false:going inward
    	return bridgesMap.get(key);
    }
    
    synchronized protected void updateBridgesMap(BidirectedBridge bridge) {
    	if(bridge==null)
    		return;
    	Node 	startNode=bridge.getStartAlignment().node,
    			endNode=bridge.getEndAlignment().node;
    	boolean startNodeDir=bridge.getStartAlignment().strand,
    			endNodeDir=!bridge.getEndAlignment().strand;
    	BidirectedBridge tmp = null;
    	//FIXME: check if worth it!
//    	bridge.addHalfBridge(bridge);
    	if(SimpleBinner.getUniqueBin(startNode)!=null){
    		tmp=getBridgeFromMap(startNode,startNodeDir);
    		bridge.merging(tmp);
    		bridgesMap.put(startNode.getId()+(startNodeDir?"o":"i"), bridge);
    	}
    	if(SimpleBinner.getUniqueBin(endNode)!=null){
    		tmp=getBridgeFromMap(endNode,endNodeDir);
    		bridge.merging(tmp);
    		bridgesMap.put(endNode.getId()+(endNodeDir?"o":"i"), bridge);
    	}
    	
    }
    
    synchronized protected void updateBridgesMap(BidirectedPath path, boolean isFullPath) {
    	if(path==null || path.size() < 2)
    		return;
    	BidirectedNode 	startNode=(BidirectedNode) path.getRoot(),
    					endNode=(BidirectedNode) path.peekNode();
    	boolean startNodeDir=((BidirectedEdge)path.getEdgePath().get(0)).getDir(startNode),
    			endNodeDir=((BidirectedEdge)path.peekEdge()).getDir(endNode);
    	BidirectedBridge brg = new BidirectedBridge(path,isFullPath);
    	if(SimpleBinner.getUniqueBin(startNode)!=null){
    		bridgesMap.put(startNode.getId()+(startNodeDir?"o":"i"), brg);
    	}
    	if(SimpleBinner.getUniqueBin(endNode)!=null){
    		bridgesMap.put(endNode.getId()+(endNodeDir?"o":"i"), brg);
    	}
    	
    }
    
    //Return bridge in the map (if any) that share the same bases (unique end) 
    synchronized public BidirectedBridge getHomoBridgeFromMap(BidirectedBridge bridge){
    	BidirectedBridge retval = null;
    	if(bridge!=null){
	    	Node 	startNode=bridge.getStartAlignment().node,
	    			endNode=bridge.getEndAlignment().node;
	    	boolean startNodeDir=bridge.getStartAlignment().strand,
	    			endNodeDir=!bridge.getEndAlignment().strand;
	    	
	    	if(SimpleBinner.getUniqueBin(startNode)!=null){
	    		retval=bridgesMap.get(startNode.getId()+(startNodeDir?"o":"i"));
	    		if(retval!=null)
	    			return retval;
	    	}
	    	
	    	if(SimpleBinner.getUniqueBin(endNode)!=null){
	    		retval=bridgesMap.get(endNode.getId()+(endNodeDir?"o":"i"));
	    	}
	    	
    	}
    	
    	return retval;
    }
    
//    synchronized protected  BidirectedBridge updateBridgesMap(BidirectedEdge e, BidirectedBridge bridge) {
//    	//Get the ending 21-mer of each nodes to find hidden potential edges? 
//    	//NOPE, only do this if a suspicious alignment appeared!
//    	if(e==null || bridge==null)
//    		return null;
//    	else
//    		return bridgesMap.put(e.getId(), bridge);
//    }
    
    synchronized public void binning() {
    	binner=new SimpleBinner(this);
    	binner.estimatePathsByCoverage();
    }

    
    /*
     * This function deduces a full path in this graph between 2 nodes aligned with a long read
     * 
     */
    synchronized protected ArrayList<BidirectedPath> getClosestPaths(Alignment from, Alignment to, int distance){
    	BidirectedNode srcNode = from.node,
    					dstNode = to.node;
    	System.out.println("Looking for path between " + srcNode.getId() + " to " + dstNode.getId() + " with distance " + distance);
    	BidirectedPath 	tmp = new BidirectedPath();
    	ArrayList<BidirectedPath>	possiblePaths = new ArrayList<BidirectedPath>(),
    								retval = new ArrayList<BidirectedPath>();
    	tmp.setRoot(srcNode);  	
    	
    	//TODO: set tolerance dynamically based on distance
//    	traverse(tmp, dstNode, possiblePaths, distance, distance>200?(int) (TOLERATE*distance):200, from.strand, to.strand, 0);
    	traverse(tmp, dstNode, possiblePaths, distance, 500, from.strand, to.strand, 0);

    	/**************************************************************************************
    	 * To cover the missing edges due to big k-mer of DBG
    	 **************************************************************************************/
    	if(possiblePaths.isEmpty()){
    		//try to find an overlap < kmer size (dead-end)
    		if(distance < 0) {
    			Sequence seq1, seq2;
    			BidirectedNode node1, node2;
    			if(Math.max(from.readStart, from.readEnd) < Math.max(to.readStart, to.readEnd)) {
    				node1=srcNode;
    				node2=dstNode;
        			seq1=(Sequence)(node1.getAttribute("seq"));
        			seq2=(Sequence)(node2.getAttribute("seq"));
        			if(!from.strand)
        				seq1=Alphabet.DNA.complement(seq1);
        			if(!to.strand)
        				seq2=Alphabet.DNA.complement(seq2);
    			}else {
    				node1=dstNode;
    				node2=srcNode;
        			seq1=(Sequence)(node1.getAttribute("seq"));
        			seq2=(Sequence)(node2.getAttribute("seq"));
        			if(!from.strand)
        				seq2=Alphabet.DNA.complement(seq1);
        			if(!to.strand)
        				seq1=Alphabet.DNA.complement(seq2);
        			
    			}
    			// stupid scanning for overlap < k
    			String prev=seq1.toString(), next=seq2.toString();
    			int index=-1;
    			for(int i=BidirectedGraph.getKmerSize()-1; i >30; i--) {
    				if(prev.substring(prev.length()-i, prev.length()).compareTo(next.substring(0,i-1))==0) {
    					index=i;
    					break;
    				}
    			}
    			if(index>0) {
        			BidirectedEdge overlapEdge = new BidirectedEdge(srcNode, dstNode, from.strand, to.strand);
        			//TODO: save the corresponding content of long reads to this edge
        			overlapEdge.setAttribute("dist", index);
        			tmp.add(overlapEdge);
        			retval.add(tmp);
        			System.out.println("Overlap path from " + srcNode.getId() + " to " + dstNode.getId() + " d=-" + index);
        			HybridAssembler.promptEnterKey();
        			return retval;
    			}else
    				return null;

    		}
    		//if a path couldn't be found between 2 dead-ends but alignments quality are insanely high
    		//FIXME: return a pseudo path having an nanopore edge
    		else if(SimpleBinner.getUniqueBin(srcNode)!=null && SimpleBinner.getUniqueBin(dstNode)!=null && srcNode.getDegree() == 1 && dstNode.getDegree()==1 &&
				Math.min(from.quality, to.quality) >= Alignment.GOOD_QUAL)
    		{
    			BidirectedEdge pseudoEdge = new BidirectedEdge(srcNode, dstNode, from.strand, to.strand);
    			//save the corresponding content of long reads to this edge
    			pseudoEdge.setAttribute("dist", distance);
    			tmp.add(pseudoEdge);
    			retval.add(tmp);
    			System.out.println("pseudo path from " + srcNode.getId() + " to " + dstNode.getId());
//    			HybridAssembler.promptEnterKey();
    			return retval;
    		}else
    			return null;
    	}
    	/*****************************************************************************************/
    	
//    	double bestScore=possiblePaths.get(0).getDeviation();
//    	for(int i=0;i<possiblePaths.size();i++){
//    		BidirectedPath p = possiblePaths.get(i);
//    		if(p.getDeviation()>bestScore*(1+TOLERATE))
//    			break;
//    		retval.add(p);
//    	}
//    	
//    	return retval;
    	//FIXME: reduce the number of returned paths here (based on the score?)
    	return possiblePaths;
    	
    }
    
    synchronized private void traverse(	BidirectedPath path, BidirectedNode dst, ArrayList<BidirectedPath> curResult, 
    						int distance, int tolerance, boolean srcDir, boolean dstDir, int stepCount)
    {
    	//stop if it's going too deep!
    	if(stepCount >= S_LIMIT)
    		return;
    	
    	BidirectedNode currentNode=(BidirectedNode) path.peekNode();
    	BidirectedEdge currentEdge;
    	boolean curDir;//direction to the next node, = ! previous'
    	
    	Iterator<Edge> ite;
    	if(path.size() <= 1) //only root
			curDir=srcDir;//re-check
		else{
			currentEdge = (BidirectedEdge) path.peekEdge();
			curDir = !((BidirectedEdge) currentEdge).getDir(currentNode);
		}
		ite=curDir?currentNode.leavingEdges().iterator():currentNode.enteringEdges().iterator();
    	while(ite.hasNext()){
    		BidirectedEdge e = (BidirectedEdge) ite.next();
			path.add(e);
			
			int delta=Math.abs(distance-e.getLength());
			if(e.getOpposite(currentNode)==dst && e.getDir(dst)!=dstDir && delta < tolerance){
		    	BidirectedPath 	curPath=curResult.isEmpty()?new BidirectedPath():curResult.get(0), //the best path saved among all possible paths from the list curResult
		    					tmpPath=new BidirectedPath(path);
		    	tmpPath.setDeviation(delta);
		    	if(	delta < curPath.getDeviation() )
		    		curResult.add(0, tmpPath);
		    	else
		    		curResult.add(tmpPath);
				
				System.out.println("Hit added: "+path.getId()+"(candidate deviation: "+delta+")");
			}else{
				int newDistance = distance - ((Sequence) e.getOpposite(currentNode).getAttribute("seq")).length() - e.getLength();
//				System.out.println("adding edge: " + e.getId() + " length=" + e.getLength() +" -> distance=" + newDistance);
				if (newDistance - e.getLength() < -tolerance){
					System.out.println("Stop go to edge " + e.getId() + " from path with distance "+newDistance+" already! : "+path.getId());
				}else
					traverse(path, dst, curResult, newDistance, tolerance, srcDir, dstDir, stepCount++);
			}
			path.popNode();
    	
    	}
    }
    
    /*
     * Find bridges based on list of Alignments.
     * Return list of bridges with endings as markers and alignments of non-markers in-between.
     */
    synchronized protected List<BidirectedPath> uniqueBridgesFinding(ArrayList<Alignment> alignments) {
		if(alignments.size()<=1)
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

		Range curRange = stepRanges.get(0);

		Alignment 	curAlignment =allAlignments.get(curRange),
					nextAlignment;
		ArrayList<BidirectedBridge> bridges = new ArrayList<>();
		BidirectedBridge 	curBridge=new BidirectedBridge(curAlignment);
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
			if(tmpBin!=null) {				
				
				curBridge.append(nextAlignment);
				bridges.add(curBridge);
				
				curBridge=new BidirectedBridge(nextAlignment);
				
				if(bins2Length.containsKey(tmpBin)){
					long newval=bins2Length.get(tmpBin)+(long)nextAlignment.node.getNumber("len");
					bins2Length.replace(tmpBin, newval);
				}else
					bins2Length.put(tmpBin, (long)nextAlignment.node.getNumber("len"));

					
				
			}else{
				curBridge.append(nextAlignment);
				
			}	
			
		}
		if(curBridge.steps.size() > 1)
			bridges.add(curBridge);
		
		
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
		for(BidirectedBridge brg:bridges) {
			System.out.printf("+++%s <=> ", brg.getEndingsID());
			BidirectedBridge storedBridge=getHomoBridgeFromMap(brg);
			if(storedBridge!=null) {
				if(storedBridge.getBridgeStatus()==1){
					System.out.println(storedBridge.getEndingsID() + ": already solved: ignore!");
					break;
				}else if(storedBridge.getBridgeStatus()==0){
					System.out.println(storedBridge.getEndingsID() + ": already built: fortify!");
					System.out.println(storedBridge.getAllPossiblePathsString());
					
//					BidirectedNode startNode=brg.getStartAlignment().node;
//					if(binner.getUniqueBin(startNode)==null)
//						startNode=brg.getEndAlignment().node;
					
//					if(storedBridge.checkIfMerge(brg))
						storedBridge.referencingTo(brg);
					if(storedBridge.getBridgeStatus()==1){
//						storedBridge.getBestPath().setConsensusUniqueBinOfPath(tmp);
						retrievedPaths.add(storedBridge.getBestPath());
					}
					continue;
				}
				

			}
			System.out.println();
			//check if brg is complete or not (only bridging complete bridge)
			if(checkCompleteBridge(brg)){				
				brg.bridging(this, tmpBin);
				if(!brg.fullPaths.isEmpty()){
					updateBridgesMap(brg);//must be here

					BidirectedPath bestPath=brg.getBestPath();
					if(bestPath!=null){
	//					brg.getBestPath().setConsensusUniqueBinOfPath(tmp);
						retrievedPaths.add(brg.getBestPath());
					}
				}
			}else {//a half bridge is already in the map... merge it!
				System.out.println(brg.getBridgeString() + ": half built: storing!");
				updateBridgesMap(brg);
			}
			
			
		}
		return retrievedPaths;
	}
 	
    //simple check if a bridge connecting 2 unique nodes
    //todo: combine info from bridgesMap also?
    private boolean checkCompleteBridge(BidirectedBridge brg){
    	if(brg==null)
    		return false;
    	else return (SimpleBinner.getUniqueBin(brg.getStartAlignment().node) != null) 
    				&& (SimpleBinner.getUniqueBin(brg.getEndAlignment().node) != null);
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

    	ArrayList<BidirectedEdge> 	potentialRemovedEdges = binner.reducedUniquePath(path);
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
    		System.out.println("Reducing path: " + path.getId());
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
        		
				if(markerNode!=null){
					//create an edge connect markerNode to curNode with curPath
					BidirectedEdge reducedEdge = new BidirectedEdge(markerNode, curNodeFromSimGraph, markerDir, curDir);
				
					tobeAdded.add(reducedEdge);
					updateBridgesMap(curPath,true);
					
					curPath.setConsensusUniqueBinOfPath(curUniqueBin);
					ArrayList<BidirectedEdge> potentialRemovedEdges = binner.reducedUniquePath(curPath);
					if(potentialRemovedEdges!=null)
						tobeRemoved.addAll(potentialRemovedEdges);
					
					HashMap<PopBin, Integer> oneBin = new HashMap<>();
					oneBin.put(curUniqueBin, 1);
					binner.edge2BinMap.put(reducedEdge, oneBin);

				}else{
					if(curPath.size()>2)
						updateBridgesMap(curPath,false);
				}
				
				
				markerNode=curNodeFromSimGraph;
        		markerDir=!curDir; //in-out, out-in
				curPath= new BidirectedPath();
				curPath.setRoot(curNodeFromSimGraph);
				
    		}
    		
    		
		}
    	if(curPath.size() > 2)
			updateBridgesMap(curPath,false);

    	
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
