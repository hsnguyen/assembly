package org.rtassembly.npgraph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Comparator;
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
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;


public class BidirectedGraph extends MultiGraph{
    static int KMER=127;
    static double RCOV=0.0;
    
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
    
	//just to make AbstractGraph.removeEdge(AbstractEdge, boolean, boolean, boolean) visible
	protected void removeEdgeDup(AbstractEdge edge, boolean graphCallback,
			boolean sourceCallback, boolean targetCallback) {
		this.removeEdge(edge, graphCallback, sourceCallback, targetCallback);
	
	}
	
	protected BidirectedEdge addEdge(AbstractNode src, AbstractNode dst, boolean dir0, boolean dir1){
		BidirectedEdge tmp = (BidirectedEdge) addEdge(BidirectedEdge.createID(src, dst, dir0, dir1), src, dst);
//		String s1=tmp.toString();
//		//tmp.setDir0(dir0);
//		//tmp.setDir1(dir1);
//		if(!s1.equals(tmp.toString()))
//			System.out.println(s1 + " ---> " + tmp);
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
    

   
//    /**
//     * 
//     * @param p Path to be grouped as a virtually vertex
//     */
//    public AbstractNode reduce(BidirectedPath p){
//    	//do nothing if the path has only one node
//    	if(p==null || p.getEdgeCount()<1)
//    		return null;
//    	
//    	//now only work with path containing more than 2 unique nodes
//    	int uniqueCount=0;
//    	for(Node n:p.getEachNode()){
//    		if(isUnique(n))
//    			uniqueCount++;
//    	}
//    	if(uniqueCount < 2)
//    	{
//    		System.out.println("ignore path with less than 1 unique contig!");
//    		return null;
//    	}
//    	//add the new composite Node to the graph
//    	//compare id from sense & anti-sense to get the unique one
//    	AbstractNode comp = addNode(p.getId().compareTo(p.getReversedComplemented().getId())>0?
//    								p.getReversedComplemented().getId():p.getId());
//    	
//    	comp.addAttribute("path", p);
//    	comp.addAttribute("seq", p.spelling());
//        comp.addAttribute("ui.label", comp.getId());
//        comp.setAttribute("ui.style", "text-offset: -10;"); 
//        comp.setAttribute("ui.class", "marked");
//        try { Thread.sleep(100); } catch (Exception e) {}
//
//    	//store unique nodes on p for removing
//    	ArrayList<String> tobeRemoved=new ArrayList<String>();
//    	for(Node n:p.getEachNode()){
//    		if(isUnique(n))
//    			tobeRemoved.add(n.getId());
//    	}
//    	BidirectedNode 	start = (BidirectedNode) p.getRoot(),
//    					end = (BidirectedNode) p.peekNode();
//    	boolean startDir = ((BidirectedEdge) p.getEdgePath().get(0)).getDir(start), 
//    			endDir = ((BidirectedEdge) p.peekEdge()).getDir(end);
//    	//set neighbors of the composite Node
//    	Iterator<Edge> startEdges = startDir?start.getEnteringEdgeIterator():start.getLeavingEdgeIterator(),
//    					endEdges = endDir?end.getEnteringEdgeIterator():end.getLeavingEdgeIterator();
//    	while(startEdges.hasNext()){
//    		BidirectedEdge e = (BidirectedEdge) startEdges.next();
//    		BidirectedNode opNode = e.getOpposite(start);
//    		boolean opDir = e.getDir(opNode);
//    		//Edge tmp=
//    		addEdge(BidirectedEdge.createID(comp, opNode, false, opDir), comp, opNode);//always into start node
//    		//System.out.println("From " + start.getId() + ": " + tmp.getId() + " added!");
//    	}
//    	
//    	while(endEdges.hasNext()){
//    		BidirectedEdge e = (BidirectedEdge) endEdges.next();
//    		BidirectedNode opNode = e.getOpposite(end);
//    		boolean opDir = e.getDir(opNode);
//    		//Edge tmp=
//    		addEdge(BidirectedEdge.createID(comp, opNode, true, opDir), comp, opNode);//always out of end node
//    	
//    		//System.out.println("From " + end.getId() + ": " + tmp.getId() + " added!");
//
//    	}
//
//    	for(String nLabel:tobeRemoved){
//    		//System.out.println("About to remove " + nLabel);
//    		removeNode(nLabel);
//    	}
//    		
//    	//todo: remove bubbles...
//    	return comp;
//    }
//    
//    
//    /**
//     * 
//     * @param v Node to be reverted (1-level reverting)
//     */
//    public void revert(AbstractNode v){
//    	System.out.println("Reverting...");
//    	Path p=v.getAttribute("path");
//    	if(p==null) return;
//    	
//    	BidirectedNode 	start = (BidirectedNode) p.getRoot(),
//    					end = (BidirectedNode) p.peekNode();
//    	boolean startDir = ((BidirectedEdge) p.getEdgePath().get(0)).getDir(start), 
//    			endDir = ((BidirectedEdge) p.peekEdge()).getDir(end);
//    	
//    	//add back all neighbor edges of this composite vertex
//    	Iterator<Edge> 	startEdges = v.getEnteringEdgeIterator(),
//						endEdges = v.getLeavingEdgeIterator();
//    	//add back all nodes from the path
//		for(Node n:p.getNodeSet()){
//			if(getNode(n.getId())!=null)
//				continue;
//			Node tmp = addNode(n.getId());
//			tmp.addAttribute("seq", (japsa.seq.Sequence)n.getAttribute("seq"));
//			tmp.addAttribute("name", (String)n.getAttribute("name"));
//			tmp.addAttribute("path", (BidirectedPath)n.getAttribute("path"));
//
//			//System.out.println("Adding back edge "+tmp.getId());
//		}
//		while(startEdges.hasNext()){
//			BidirectedEdge e = (BidirectedEdge) startEdges.next();
//			BidirectedNode opNode = e.getOpposite(v);
//			boolean opDir = e.getDir(opNode);
//			//Edge tmp = 
//			addEdge(BidirectedEdge.createID(start, opNode, !startDir, opDir), start, opNode);
//			//System.out.println("Adding back edge "+tmp.getId());
//		}
//		
//		while(endEdges.hasNext()){
//			BidirectedEdge e = (BidirectedEdge) endEdges.next();
//			BidirectedNode opNode = e.getOpposite(v);
//			boolean opDir = e.getDir(opNode);
//			//Edge tmp = 
//			addEdge(BidirectedEdge.createID(end, opNode, !endDir, opDir), end, opNode);
//			//System.out.println("Adding back edge "+tmp.getId());
//		}
//
//    	//add back all edges from the path
//		for(Edge e:p.getEdgeSet()){
//			//Edge tmp = 
//			addEdge(e.getId(), e.getSourceNode().getId(), e.getTargetNode().getId());
//			//System.out.println("Adding back edge "+tmp.getId());
//		}
//    	//finally remove the composite node
//    	removeNode(v);
//	}
    
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
    	
    	//traverse(tmp, dest, retval, distance+source.getSeq().length()+dest.getSeq().length());
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
    		else if(isMarker(srcNode) && isMarker(dstNode) && srcNode.getDegree() == 1 && dstNode.getDegree()==1 &&
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
		
		System.out.println("Binning ranges: ");
	    for(List<Range> group : rangeGroups){
	    	int maxscore=0;
	    	Range rangeOfBest=null;
	    	for(Range range:group) { 
	    		System.out.print(allAlignments.get(range).node.getId() + ": " + range + "; ");	    
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
		
		if(isMarker(curAlignment.node))
			curBridge=new BidirectedBridge(curAlignment);
				
		for(int i=1; i<stepRanges.size();i++){
			Range nextRanges = stepRanges.get(i);
			nextAlignment = allAlignments.get(nextRanges);
			if(isMarker(nextAlignment.node)) {
				if(curBridge!=null) {
					curBridge.append(nextAlignment);
					bridges.add(curBridge);
				}
				curBridge=new BidirectedBridge(nextAlignment);
				
			}else if(curBridge!=null){
				curBridge.append(nextAlignment);
				
			}	
			
		}
		ArrayList<BidirectedPath> retrievedPaths = new ArrayList<>();
		// Now we got all possible unique bridges from the alignments, do smt with them:
		for(BidirectedBridge brg:bridges) {
			//TODO: check already-found path here! (increase score, conflict resolve???)
			if(graphMap.get(brg.getEndingsID())!=null)
				continue;
			brg.bridging(this);
			retrievedPaths.add(brg.getPath());
		}
		
		return retrievedPaths;
	}
    /*
     * Find a path based on list of Alignments
     * Only best alignment in each group are chosen to find the path to
     */
    synchronized public List<BidirectedPath> pathFinding(ArrayList<Alignment> alignments) {
		if(alignments.size()<=1)
			return null;
		
		System.out.println("=================================================");
		for(Alignment alg:alignments)
			System.out.println("\t"+alg.toString());
		System.out.println("=================================================");
		//First bin the alignments into different overlap regions			
		//only considering useful alignments
		HashMap<Range,Alignment> allAlignments = new HashMap<Range,Alignment>();
		ArrayList<BidirectedPath> joinPaths = new ArrayList<>();;
		
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
		
		System.out.println("Binning ranges: ");
	    for(List<Range> group : rangeGroups){
	    	int maxscore=0;
	    	Range rangeOfBest=null;
	    	for(Range range:group) { 
	    		System.out.print(allAlignments.get(range).node.getId() + ": " + range + "; ");	    
	    		if(allAlignments.get(range).quality > maxscore)
	    			rangeOfBest=range;
	    	}
	  	
	    	stepRanges.add(rangeOfBest);
	    	System.out.println();
	    }

		Range curRange = stepRanges.get(0);

		Alignment 	curAlignment =allAlignments.get(curRange),
					nextAlignment;

		BidirectedPath curPath=null;
		for(int i=1; i<stepRanges.size();i++){
			Range nextRanges = stepRanges.get(i);
			nextAlignment = allAlignments.get(nextRanges);
			
			BidirectedPath curBridge=null;
			
			

			int distance = nextAlignment.readAlignmentStart()-curAlignment.readAlignmentEnd();
			if(distance<D_LIMIT){
				ArrayList<BidirectedPath> bridges = getClosestPath(curAlignment, nextAlignment, distance);
				if(bridges!=null)
					curBridge = bridges.get(0);
			}else 
				continue;

			//join all paths from previous to the new ones
			//TODO:optimize it
			System.out.println("====Before curPath: " + curPath + "; curBridge:" + curBridge);

			if(curPath==null)
				curPath=curBridge;
			else{
				if(!curPath.join(curBridge)) {
					joinPaths.add(curPath);		
					curPath=curBridge;
				}
			}
			curAlignment=nextAlignment;
			
			System.out.println("====After curPath: " + curPath + "; curBridge:" + curBridge);
		}
		if(curPath!=null)
			joinPaths.add(curPath);
		

		if(joinPaths.isEmpty())
			return null;
		else {
			for(BidirectedPath path:joinPaths)
				System.out.println("A member Path: " + path.toString() + " deviation: " + path.getDeviation());
			
			return joinPaths;
		}
	}
	
    /**
     * Another reduce that doesn't remove the unique nodes
     * Instead redundant edges are removed on a path way
     * @param path Path to simplify the graph (from origGraph)
     * @param target Subjected graph for the simplification
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
    	
    	if(BidirectedGraph.isMarker(curNodeFromSimGraph)){
    		markerNode=curNodeFromSimGraph;
    		markerDir=((BidirectedEdge) path.getEdgePath().get(0)).getDir(markerNode);
    		curPath = new BidirectedPath();
    		curPath.setRoot(curNodeFromSimGraph);
    	}
    	

    	//search for an unique node as the marker. 
    	ArrayList<BidirectedEdge> 	tobeRemoved = new ArrayList<BidirectedEdge>(),
    								tobeAdded = new ArrayList<BidirectedEdge>();
    	double aveCov=0;
    	for(Edge edge:path.getEdgePath()){
    			
    		curNodeFromSimGraph=(BidirectedNode) edge.getOpposite(curNodeFromSimGraph);
    		   		
//    		curNodeFromSimGraph = simGraph.getNode(curNodeFromOrigGraph.getId()); //change back to Node belong to simGraph (instead of origGraph)
    		curDir=((BidirectedEdge) edge).getDir(curNodeFromSimGraph);
    		
    		if(BidirectedGraph.isMarker(curNodeFromSimGraph)){
        		
				if(markerNode!=null){
					//this is when we have 1 jumping path (both ends are markers)
					curPath.add(edge);	
//					LOG.info("Processing path {} with marker {}:{}:{} and curNode {}:{}:{}", curPath.getId(), markerNode.getId(), markerDir?"out":"in", markerNode.getGraph().getId(), curNodeFromSimGraph.getId(), curDir?"out":"in", curNodeFromSimGraph.getGraph().getId());
					//create an edge connect markerNode to curNode with curPath
					//Edge reducedEdge = simGraph.addEdge(markerNode, curNodeFromSimGraph, markerDir, curDir);
					BidirectedEdge reducedEdge = new BidirectedEdge(markerNode, curNodeFromSimGraph, markerDir, curDir);

//					if(reducedEdge!=null)
//						reducedEdge.addAttribute("path", new BidirectedPath(curPath));
				
					tobeAdded.add(reducedEdge);
					updateGraphMap(reducedEdge, curPath);
					//loop over curPath to find out edges needed to be removed
//					Node  	n0 = curPath.getRoot(),
//							n1 = null;
					
					Node  	curNode = curPath.getRoot();
					aveCov=curPath.averageCov();
					for(Edge ep:curPath.getEdgePath()){
						LOG.info("--edge {} coverage:{} to {}",ep.getId(),ep.getNumber("cov"),ep.getNumber("cov") - aveCov);
						ep.setAttribute("cov", ep.getNumber("cov") - aveCov);	

						if(ep.getNumber("cov")/aveCov < .5 && (BidirectedGraph.isMarker(ep.getSourceNode()) || BidirectedGraph.isMarker(ep.getTargetNode())) ) //plasmid coverage is different!!!
							tobeRemoved.add((BidirectedEdge) ep);
						
						if(!BidirectedGraph.isMarker(curNode))
							curNode.setAttribute("cov", curNode.getNumber("cov")>aveCov?curNode.getNumber("cov")-aveCov:0);
						curNode=ep.getOpposite(curNode);
//						n1 = ep.getOpposite(n0);
//						if(!BidirectedGraph.isUnique(n0) == BidirectedGraph.isUnique(n1)){
//			    			tobeRemoved.add((BidirectedEdge)ep);
////			    			if(BidirectedGraph.isUnique(n0))
////			    				n1.setAttribute("cov", Math.max(n1.getNumber("cov")-n0.getNumber("cov"),n0.getNumber("cov")));   
////			    			if(BidirectedGraph.isUnique(n1))
////			    				n0.setAttribute("cov", Math.max(n0.getNumber("cov")-n1.getNumber("cov"),n1.getNumber("cov")));   
//
//						}
//						
//						
//			    		n0=n1;
					}

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
    	
    	boolean retval = tobeRemoved.size()>0?true:false;
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
    	
    	return retval;

    }
    
    /*
     * Important function: determine if a node is able to be removed or not
     * TODO: re-implement it based on statistics of coverage also
     * 1. pick the least coverage ones among a path as the base
     * 2. global base
     */
    synchronized public static boolean isMarker(Node node){
    	boolean res = false;
    	
    	if(node.getDegree()<=2){ // not always true, e.g. unique node in a repetitive component
    		Sequence seq = (Sequence) node.getAttribute("seq");
//    		if(seq.length() > 1000 && node.getNumber("astats") > 10 && node.getNumber("cov")/RCOV > .3)
    		if(seq.length() > 300 && node.getNumber("astats") > 10 && node.getNumber("cov")/RCOV > .3)
    			res=true;
    	}
    	
//    	if(res)
//    		LOG.info(node.getAttribute("name") + " with coverage " + node.getNumber("cov") + " is a marker!");
//    	else
//    		LOG.info(node.getAttribute("name") + " with coverage " + node.getNumber("cov") + " is NOT a marker!");

    	return res;
    }

    	
    /*
     * Traverse the graph and assign weights to every edges based on coverage of ending nodes
     * and update a node's coverage if possible (later)
     */
    synchronized public void balancing() {
    	List<Edge> unknownEdges=this.edges().collect(Collectors.toList());
    	ArrayList<Edge> newUnknownEdges = new ArrayList<Edge>();
    	while(true) {
	    	for(Edge e:unknownEdges) {
	    		if(!Double.isNaN(e.getNumber("cov")))
	    			continue;
	    			
	    		BidirectedNode n0 = (BidirectedNode) e.getNode0(), n1=(BidirectedNode) e.getNode1();
	    		boolean dir0 = ((BidirectedEdge) e).getDir0(), dir1 = ((BidirectedEdge) e).getDir1();
	    		Iterator<Edge> 	in0 = n0.enteringEdges().iterator(), out0 = n0.leavingEdges().iterator(),
	    						in1 = n1.enteringEdges().iterator(), out1 = n1.leavingEdges().iterator();
	    		int unknwIn0 = 0, unknwOut0 = 0, unknwIn1  = 0, unknwOut1 = 0;
	    		double inWeight0 = 0, outWeight0 = 0, inWeight1 = 0, outWeight1 = 0;
	    		while(in0.hasNext()) {
	    			BidirectedEdge tmp = (BidirectedEdge) in0.next();
	    			double tmpW = tmp.getNumber("cov");
	    			if(!Double.isNaN(tmpW))
	    				inWeight0+=tmpW;
	    			else
	    				unknwIn0++;
	    		}
	    		while(out0.hasNext()) {
	    			BidirectedEdge tmp = (BidirectedEdge) out0.next();
	    			double tmpW = tmp.getNumber("cov");
	    			if(!Double.isNaN(tmpW))
	    				outWeight0+=tmpW;
	    			else
	    				unknwOut0++;
	    		}
	    		
	    		while(in1.hasNext()) {
	    			BidirectedEdge tmp = (BidirectedEdge) in1.next();
	    			double tmpW = tmp.getNumber("cov");
	    			if(!Double.isNaN(tmpW))
	    				inWeight1+=tmpW;
	    			else
	    				unknwIn1++;
	    		}
	    		while(out1.hasNext()) {
	    			BidirectedEdge tmp = (BidirectedEdge) out1.next();
	    			double tmpW = tmp.getNumber("cov");
	    			if(!Double.isNaN(tmpW))
	    				outWeight1+=tmpW;
	    			else
	    				unknwOut1++;
	    		}
	    		
	    		//calculate sumOfCovToN0, sumOfCovFromN0; sumOfCovToN1, sumOfCovFromN1;
	    		//should be: sumOfCovToN0==sumOfCovFromN0==covOfN0; sumOfCovToN1==sumOfCovFromN1==covOfN1 ??? 
	    		double covInferredFromN0 = Double.NaN, covInferredFromN1=Double.NaN, estimatedCov=Double.NaN;
	    		//int confidence = 0; //maximum=sum of length of 2 end nodes
	    		double aveCov0=n0.getNumber("cov"), aveCov1=n1.getNumber("cov");
	    		if(dir0) {
	    			if(unknwOut0 == 1) {
	    				if(unknwIn0==0) {
	    					if(approxCompare(inWeight0, aveCov0)!=0)
	    						System.out.println("Node " + n0.getAttribute("name") + ": sum of entering edges = " + inWeight0);
	    					aveCov0=calibrate(aveCov0,inWeight0);
	    				}
	    				covInferredFromN0 = aveCov0-outWeight0;
	    			}
	    		}else {
	    			if(unknwIn0 == 1) {
	    				if(unknwOut0==0) {
	    					if(approxCompare(outWeight0, aveCov0)!=0)
	    						System.out.println("Node " + n0.getAttribute("name") + ": sum of out-going edges = " + outWeight0);
	    					aveCov0=calibrate(aveCov0,outWeight0);
	    				}
	    				covInferredFromN0 = aveCov0-inWeight0;
	    			}
	    		}
	    		
	    		if(dir1) {
	    			if(unknwOut1 == 1) {
	    				if(unknwIn1==0) {
	    					if(approxCompare(inWeight1, aveCov1)!=0)
	    						System.out.println("Node " + n1.getAttribute("name") + ": sum of entering edges = " + inWeight1);
	    					aveCov1=calibrate(aveCov1,inWeight1);
	    				}
	    				covInferredFromN1 = aveCov1-outWeight1;
	    			}
	    		}else {
	    			if(unknwIn1 == 1) {
	    				if(unknwOut1==0) {
	    					if(approxCompare(outWeight1, aveCov1)!=0)
	    						System.out.println("Node " + n1.getAttribute("name") + ": sum of out-going edges = " + outWeight1);
	    					aveCov1=calibrate(aveCov1,outWeight1);
	    				}
	    				covInferredFromN1 = aveCov1-inWeight1;
	    				
	    			}
	    		}
	    		
	    		
	    		if(Double.isNaN(covInferredFromN0) && Double.isNaN(covInferredFromN1)) { //both are unknown
	    			System.out.println("==> Unable to resolve: " + e.getId() + " : " + (dir0?unknwOut0:unknwIn0) + "(" + n0.getId()+ ") --- " + (dir1?unknwOut1:unknwIn1) + "(" + n1.getId() + ")");
	    		}else if (!Double.isNaN(covInferredFromN0) && !Double.isNaN(covInferredFromN1)){
					if(approxCompare(covInferredFromN0, covInferredFromN1)!=0) {
						System.out.println("==> Infer from node " + n0.getAttribute("name") + " cov="+ n0.getNumber("cov") + ": " + covInferredFromN0 
								+ " and from node " + n1.getAttribute("name") + " cov="+ n1.getNumber("cov") + ": " + covInferredFromN1);
						printEdgesCov(n0);
						printEdgesCov(n1);
					}
	    			
					//estimatedCov=calibrate(covInferredFromN0,covInferredFromN1);
					estimatedCov=calibrateWithWeight(	covInferredFromN0, ((Sequence)(n0.getAttribute("seq"))).length(), 
														covInferredFromN1, ((Sequence)(n1.getAttribute("seq"))).length());
	    		}else {
	    			estimatedCov=Double.isNaN(covInferredFromN0)?covInferredFromN1:covInferredFromN0;
	    		}
	    		
	    		if(Double.isNaN(estimatedCov))
	    			newUnknownEdges.add(e);
	    		else {
	    			System.out.println("Final estimated coverage for " + e.getId() + " is " + estimatedCov);

	    			e.setAttribute("cov", estimatedCov);
	    		}
					
	    	}

	    	unknownEdges = newUnknownEdges;
	    	newUnknownEdges = new ArrayList<Edge>();
	    	if(unknownEdges.isEmpty())
	    		break;
	    	else {
	    		
	    		System.out.println("Next traversal round on "+unknownEdges.size() + " edges");
	    		continue;
	    	}
    	}
    	
    }
    /*
     * Compare 2 double values x,y
     * Return 0 if x~=y, 1 if x>>y, -1 if x<<y
     */
    public static int approxCompare(double x, double y) {
    	int retval=0;
    	double ratio=Math.abs(x-y)/(Math.max(Math.abs(x), Math.abs(y)));
    	if(ratio > .3)
    		retval=x>y?1:-1;
    	
    	return retval;
    }
    
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
    
    /*
     * Balance function
     */
    private double calibrate(double x, double y) {
    	return Math.max(x, y);
//    	return (a+b)/2;
    }
    private double calibrateWithWeight(double x, int xw, double y, int yw){
    	return (x*xw+y*yw)/(xw+yw);
    }
    
}
