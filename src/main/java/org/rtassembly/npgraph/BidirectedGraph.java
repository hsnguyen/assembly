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
    
    static final int TOLERATE=500;
    static final int D_LIMIT=10000; //distance bigger than this will be ignored
    static final int S_LIMIT=50;
    
    //provide dynamic state of particular pair of nodes with directions (e.g. A+B- : state)
    //state of edge: -1=removed, 0=connect, 1=new connect(reduce edge)
    private HashMap<String, BidirectedPath> graphMap; 
    
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
		graphMap=new HashMap<String, BidirectedPath>(initialNodeCapacity*(initialNodeCapacity+1)/2);
		
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
    synchronized protected  BidirectedPath updateGraphMap(String id, BidirectedPath path) {
    	//Get the ending 21-mer of each nodes to find hidden potential edges? 
    	//NOPE, only do this if a suspicious alignment appeared!
    	if(id==null || path==null)
    		return null;
    	else
    		return graphMap.put(id, path);
    } 
    synchronized protected  BidirectedPath updateGraphMap(BidirectedEdge e, BidirectedPath path) {
    	//Get the ending 21-mer of each nodes to find hidden potential edges? 
    	//NOPE, only do this if a suspicious alignment appeared!
    	if(e==null || path==null)
    		return null;
    	else
    		return graphMap.put(e.getId(), path);
    }
    
    synchronized public void binning() {
    	binner=new SimpleBinner(this);
    	binner.estimatePathsByCoverage();
    }

    
    /*
     * This function deduces a full path in this graph between 2 nodes aligned with a long read
     * 
     */
    synchronized protected ArrayList<BidirectedPath> getClosestPath(Alignment from, Alignment to, int distance){
    	BidirectedNode srcNode = from.node,
    					dstNode = to.node;
    	System.out.println("Looking for path between " + srcNode.getId() + " to " + dstNode.getId() + " with distance " + distance);
    	BidirectedPath 	tmp = new BidirectedPath();
    	ArrayList<BidirectedPath>	possiblePaths = new ArrayList<BidirectedPath>(),
    								retval = new ArrayList<BidirectedPath>();
    	tmp.setRoot(srcNode);  	
    	
    	traverse(tmp, dstNode, possiblePaths, distance, from.strand, to.strand, 0);
    	//only get the best ones
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
    			//TODO: save the corresponding content of long reads to this edge
    			pseudoEdge.setAttribute("dist", distance);
    			tmp.add(pseudoEdge);
    			retval.add(tmp);
    			System.out.println("pseudo path from " + srcNode.getId() + " to " + dstNode.getId());
//    			HybridAssembler.promptEnterKey();
    			return retval;
    		}else
    			return null;
    	}
    	double bestScore=possiblePaths.get(0).getDeviation();
    	for(int i=0;i<possiblePaths.size();i++){
    		BidirectedPath p = possiblePaths.get(i);
    		if(p.getDeviation() != bestScore)
    			break;
    		retval.add(p);
    	}
    	
    	return retval;
    	
    }
    
    synchronized private void traverse(	BidirectedPath path, BidirectedNode dst, ArrayList<BidirectedPath> curResult, 
    						int distance, boolean srcDir, boolean dstDir, int stepCount)
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
			
			int toTarget=Math.abs(distance-e.getLength());
			if(e.getOpposite(currentNode)==dst && e.getDir(dst)!=dstDir && toTarget < TOLERATE){
		    	BidirectedPath 	curPath=curResult.isEmpty()?new BidirectedPath():curResult.get(0), //the best path saved among all possible paths from the list curResult
		    					tmpPath=new BidirectedPath(path);
		    	tmpPath.setDeviation(toTarget);
		    	if(	toTarget < curPath.getDeviation() )
		    		curResult.add(0, tmpPath);
		    	else
		    		curResult.add(tmpPath);
				
				System.out.println("Hit added: "+path.getId()+"(candidate deviation: "+toTarget+")");
			}else{
				int newDistance = distance - ((Sequence) e.getOpposite(currentNode).getAttribute("seq")).length() - e.getLength();
//				System.out.println("adding edge: " + e.getId() + " length=" + e.getLength() +" -> distance=" + newDistance);
				if (newDistance - e.getLength() < -TOLERATE){
					System.out.println("Stop go to edge " + e.getId() + " from path with distance "+newDistance+" already! : "+path.getId());
				}else
					traverse(path, dst, curResult, newDistance, srcDir, dstDir, stepCount++);
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
		BidirectedBridge curBridge=null;
		PopBin tmp=null;
		HashMap<PopBin, Long> bins2Length = new HashMap<PopBin,Long>();
		
		tmp=binner.getUniqueBin(curAlignment.node);
		if( tmp != null){
			curBridge=new BidirectedBridge(curAlignment);
			bins2Length.put(tmp, (long)curAlignment.node.getNumber("len"));
		}
				
		for(int i=1; i<stepRanges.size();i++){
			Range nextRanges = stepRanges.get(i);
			nextAlignment = allAlignments.get(nextRanges);
			tmp=binner.getUniqueBin(nextAlignment.node);
			if(tmp!=null) {
				if(curBridge!=null) {
					curBridge.append(nextAlignment);
					bridges.add(curBridge);
				}
				curBridge=new BidirectedBridge(nextAlignment);
				
				if(bins2Length.containsKey(tmp)){
					long newval=bins2Length.get(tmp)+(long)nextAlignment.node.getNumber("len");
					bins2Length.replace(tmp, newval);
				}else
					bins2Length.put(tmp, (long)nextAlignment.node.getNumber("len"));

					
				
			}else if(curBridge!=null){
				curBridge.append(nextAlignment);
				
			}	
			
		}
		//determine the global unique bin of the whole path
		long ltmp=0;
		for(PopBin b:bins2Length.keySet())
			if(bins2Length.get(b)>ltmp){
				ltmp=bins2Length.get(b);
				tmp=b;
			}
		
		ArrayList<BidirectedPath> retrievedPaths = new ArrayList<>();
		// Now we got all possible unique bridges from the alignments, do smt with them:
		System.out.println("\n=> bridges list: ");
		for(BidirectedBridge brg:bridges) {
			//TODO: check already-found path here: both ends must have reasonable bin!
			System.out.printf("...%s ", brg.getEndingsID());
			if(brg==null||graphMap.get(brg.getEndingsID())!=null) {
				System.out.println(": ignored, already processed!");
				continue;
			}
			System.out.println();
			brg.bridging(this);
			if(brg.getPath()==null)
				continue;
			brg.getPath().setConsensusUniqueBinOfPath(tmp);
			retrievedPaths.add(brg.getPath());
		}
		return retrievedPaths;
	}
 	
    /**
     * Another reduce that doesn't remove the unique nodes
     * Instead redundant edges are removed on a path way
     * @param path: unique path to simplify the graph (from origGraph)
     */
    public boolean reduce(BidirectedPath path){
    	//do nothing if the path has only one node
    	if(path==null || path.getEdgeCount()<1)
    		return false;
    	else
    		System.out.println("Reducing path: " + path.getId());
    	//loop over the edges of path (like spelling())
    	BidirectedNode 	markerNode = null,
    			curNodeFromSimGraph = (BidirectedNode) path.getRoot();
	
    	BidirectedPath curPath= null;
    	boolean markerDir=true, curDir;
    	PopBin 	curUniqueBin = binner.getUniqueBin(curNodeFromSimGraph);
    	if(curUniqueBin!=null){
    		markerNode=curNodeFromSimGraph;
    		markerDir=((BidirectedEdge) path.getEdgePath().get(0)).getDir(markerNode);
    		curPath = new BidirectedPath();
    		curPath.setRoot(curNodeFromSimGraph);
    	}
    	

    	//search for an unique node as the marker. 
    	ArrayList<BidirectedEdge> 	tobeRemoved = new ArrayList<BidirectedEdge>(),
    								tobeAdded = new ArrayList<BidirectedEdge>();
    	for(Edge edge:path.getEdgePath()){
    			
    		curNodeFromSimGraph=(BidirectedNode) edge.getOpposite(curNodeFromSimGraph);
    		   		
//    		curNodeFromSimGraph = simGraph.getNode(curNodeFromOrigGraph.getId()); //change back to Node belong to simGraph (instead of origGraph)
    		curDir=((BidirectedEdge) edge).getDir(curNodeFromSimGraph);
    		
    		curUniqueBin = binner.getUniqueBin(curNodeFromSimGraph);//TODO: check consistency
    		if(curUniqueBin!=null){//only when reach the end of path
        		
				if(markerNode!=null){
					//this is when we have 1 jumping path (both ends are markers)
					curPath.add(edge);	
//					LOG.info("Processing path {} with marker {}:{}:{} and curNode {}:{}:{}", curPath.getId(), markerNode.getId(), markerDir?"out":"in", markerNode.getGraph().getId(), curNodeFromSimGraph.getId(), curDir?"out":"in", curNodeFromSimGraph.getGraph().getId());
					//create an edge connect markerNode to curNode with curPath
					BidirectedEdge reducedEdge = new BidirectedEdge(markerNode, curNodeFromSimGraph, markerDir, curDir);

//					if(reducedEdge!=null)
//						reducedEdge.addAttribute("path", new BidirectedPath(curPath));
				
					tobeAdded.add(reducedEdge);
					updateGraphMap(reducedEdge, curPath);
					
					curPath.setConsensusUniqueBinOfPath(path.getConsensusUniqueBinOfPath());
					tobeRemoved=binner.reducedUniquePath(curPath);
					
					HashMap<PopBin, Integer> oneBin = new HashMap<>();
					oneBin.put(curUniqueBin, 1);
					binner.edge2BinMap.put(reducedEdge, oneBin);

				}
				
				
				markerNode=curNodeFromSimGraph;
        		markerDir=!curDir; //in-out, out-in
				curPath= new BidirectedPath();
				curPath.setRoot(curNodeFromSimGraph);
				
    		}
    		else{
    			if(markerNode!=null){
    				curPath.add(edge);
    			}
    		}
    		
		}
    	
    	if(tobeRemoved!=null && tobeRemoved.size()>0){
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
					reducedEdge.setAttribute("isReducedEdge", true);
					reducedEdge.setAttribute("ui.class", "marked");
	//				reducedEdge.addAttribute("layout.weight", 10);
				}
	    		LOG.info("after: \n\t" + printEdgesOfNode((BidirectedNode) e.getNode0()) + "\n\t" + printEdgesOfNode((BidirectedNode) e.getNode1()));
	
	    	}
	    	return true;
    	}else {
    		LOG.info("Nothing to remove!");
    		return false;
    	}

    }
    
    //This assuming path is surely unique!!!
    public boolean reduce2(BidirectedPath path){
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
				reducedEdge.setAttribute("ui.label", path.getId());
				reducedEdge.setAttribute("ui.style", "text-offset: -10; text-alignment: along;"); 
				reducedEdge.setAttribute("isReducedEdge", true);
				reducedEdge.setAttribute("ui.class", "marked");
//				reducedEdge.addAttribute("layout.weight", 10);
				binner.edge2BinMap.put(reducedEdge, oneBin);
				updateGraphMap(reducedEdge, path);
			}
			
	    	return true;
    	}else {
    		LOG.info("Path {} has failed reduce operation!", path.getId());
    		return false;
    	}

    }
    /*
     * Important function: determine if a node is able to be removed or not
     * TODO: re-implement it based on statistics of coverage also
     * 1. pick the least coverage ones among a path as the base
     * 2. global base
     */
//    synchronized public static boolean isMarker(Node node){
//    	boolean res = false;
//    	
//    	if(node.getDegree()<=2){ // not always true, e.g. unique node in a repetitive component
//    		Sequence seq = (Sequence) node.getAttribute("seq");
////    		if(seq.length() > 1000 && node.getNumber("astats") > 10 && node.getNumber("cov")/RCOV > .3)
//    		if(seq.length() > 300 && node.getNumber("astats") > 10 && node.getNumber("cov")/RCOV > .3)
//    			res=true;
//    	}
//    	
////    	if(res)
////    		LOG.info(node.getAttribute("name") + " with coverage " + node.getNumber("cov") + " is a marker!");
////    	else
////    		LOG.info(node.getAttribute("name") + " with coverage " + node.getNumber("cov") + " is NOT a marker!");
//
//    	return res;
//    }
    
    /*
     * Print adjacent edges' coverage of a node
     */
    public void printEdgesCov(Node node) {
    	Iterator<Edge> in = node.enteringEdges().iterator(), out = node.leavingEdges().iterator();
    	while(in.hasNext()) {
    		BidirectedEdge e = (BidirectedEdge) in.next();
    		System.out.printf("[%s]=%.2f; ", e.getId(), e.getNumber("cov"));
    	}
    	System.out.print("\n\tOUT: ");
    	while(out.hasNext()) {
    		BidirectedEdge e = (BidirectedEdge) out.next();
    		System.out.printf("[%s]=%.2f; ", e.getId(), e.getNumber("cov"));
    	}
    	System.out.println();
    }
    
}
