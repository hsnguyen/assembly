package org.rtassembly.npgraph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
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
	private SimpleBinner binner;

    volatile static double ILLUMINA_READ_LENGTH=300; //Illumina MiSeq
    
    static final double TOLERATE=.3;//can be interpreted as long read error rate (10-25%)
    static final int D_LIMIT=5000; //distance bigger than this will be ignored
    static final int S_LIMIT=100;// maximum number of DFS steps
    
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
	
	public void printNodeSequencesToFile(String fileName) throws IOException {
		SequenceOutputStream out = SequenceOutputStream.makeOutputStream(fileName);
		
		for(Node node:this) {
			Sequence seq=(Sequence) node.getAttribute("seq");
			seq.writeFasta(out);
		}
		
		out.close();
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
    /*
     * Adding bridge to the map, return false if entry already exists
     * TODO: conflict inductions of bridges should be solved here
     */
    synchronized protected boolean updateBridgesMap(String brgID, BidirectedBridge bridge) {
    	if(brgID==null || bridge==null)
    		return false;
    	else{
    		//FIXME: should check for uniqueness here???
    		String[] keys=GraphUtil.getHeadOfPathString(brgID);
    		return bridgesMap.put(keys[0], bridge)==null && bridgesMap.put(keys[1], bridge)==null;
    	}
    } 
    synchronized public BidirectedBridge getBridgeFromMap(String brgID){
    	BidirectedBridge retval=null;
		String[] keys=GraphUtil.getHeadOfPathString(brgID);
		if(bridgesMap.get(keys[0]) != null)
			retval=bridgesMap.get(keys[0]);
		
		if(bridgesMap.get(keys[1]) != null)
			retval=bridgesMap.get(keys[1]);	
		
    	return retval;
    }
    synchronized public BidirectedBridge getBridgeFromMap(Node unqNode, boolean direction){
    	String key=unqNode.getId()+(direction?"+":"-");
    	return bridgesMap.get(key);
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
    		else if(binner.getUniqueBin(srcNode)!=null && binner.getUniqueBin(dstNode)!=null && srcNode.getDegree() == 1 && dstNode.getDegree()==1 &&
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
    	//stop if it's going too far!
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
		//FIXME: use stream directly?
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
	    	int maxscore=0;
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
		BidirectedBridge curBridge=new BidirectedBridge(curAlignment);
		PopBin tmp=null;
		
		HashMap<PopBin, Long> bins2Length = new HashMap<PopBin,Long>();
		
		tmp=binner.getUniqueBin(curAlignment.node);
		if( tmp != null){
			bins2Length.put(tmp, (long)curAlignment.node.getNumber("len"));
		}
				
		for(int i=1; i<stepRanges.size();i++){
			Range nextRanges = stepRanges.get(i);
			nextAlignment = allAlignments.get(nextRanges);
			tmp=binner.getUniqueBin(nextAlignment.node);
			if(tmp!=null) {				
				curBridge.append(nextAlignment);
				bridges.add(curBridge);
				
				curBridge=new BidirectedBridge(nextAlignment);
				
				if(bins2Length.containsKey(tmp)){
					long newval=bins2Length.get(tmp)+(long)nextAlignment.node.getNumber("len");
					bins2Length.replace(tmp, newval);
				}else
					bins2Length.put(tmp, (long)nextAlignment.node.getNumber("len"));

					
				
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
				tmp=b;
			}
		
		ArrayList<BidirectedPath> retrievedPaths = new ArrayList<>();
		// Now we got all possible bridges from chopping the alignments at unique nodes
		System.out.println("\n=> bridges list: ");
		for(BidirectedBridge brg:bridges) {
			System.out.printf("+++%s <=> ", brg.getEndingsID());
			BidirectedBridge storedBridge=getBridgeFromMap(brg.getEndingsID());
			if(storedBridge!=null) {
				if(storedBridge.isSolved){
					System.out.println(storedBridge.getEndingsID() + ": already processed: ignore!");
					continue;
				}else{
					System.out.println(storedBridge.getEndingsID() + ": already processed: fortify!");
					System.out.println(storedBridge.getAllPossiblePathsString());
					
					BidirectedNode startNode=brg.getStartAlignment().node;
					if(binner.getUniqueBin(startNode)==null)
						startNode=brg.getEndAlignment().node;
					storedBridge.merging(brg, startNode);
					if(storedBridge.isSolved){
						storedBridge.getBestPath().setConsensusUniqueBinOfPath(tmp);
						retrievedPaths.add(storedBridge.getBestPath());
						continue;
					}
						
				}
			}
			System.out.println();
			
			//check if brg is unique or not (only bridging unique bridge)
			if(checkUniqueBridge(brg)){				
				brg.bridging(this);
				if(!brg.paths.isEmpty())
					updateBridgesMap(brg.getEndingsID(), brg);//must be here

				BidirectedPath bestPath=brg.getBestPath();
				if(bestPath!=null){
					brg.getBestPath().setConsensusUniqueBinOfPath(tmp);
					retrievedPaths.add(brg.getBestPath());
				}
			}
		}
		return retrievedPaths;
	}
 	
    //simple check if a bridge connecting 2 unique nodes
    //todo: combine info from bridgesMap also?
    private boolean checkUniqueBridge(BidirectedBridge brg){
    	if(brg==null)
    		return false;
    	else return (binner.getUniqueBin(brg.getStartAlignment().node) != null) 
    				&& (binner.getUniqueBin(brg.getEndAlignment().node) != null);
    }
    /**
     * Another reduce that doesn't remove the unique nodes
     * Instead redundant edges are removed on a path way
     * @param path: unique path to simplify the graph (from origGraph)
     */

    
    //This assuming path is surely unique!!!
    public boolean reduce(BidirectedPath path){
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

    	//search for an unique node as the marker. 
    	ArrayList<BidirectedEdge> 	tobeRemoved = binner.reducedUniquePath(path);
		HashMap<PopBin, Integer> oneBin = new HashMap<>();
		oneBin.put(path.getConsensusUniqueBinOfPath(), 1);
    	
    	if(tobeRemoved!=null && tobeRemoved.size()>1){
	    	//remove appropriate edges
	    	for(BidirectedEdge e:tobeRemoved){
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
    
    
}
