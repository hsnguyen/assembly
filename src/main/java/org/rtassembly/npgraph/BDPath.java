package org.rtassembly.npgraph;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;

import java.util.HashMap;
import java.util.List;

import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.graphstream.graph.Path;
import org.graphstream.graph.implementations.AbstractNode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class BDPath extends Path{
	private int diff; //deviation from path to aligned read (now just length): the more the worse
	private double pathEstats=Double.MIN_NORMAL; //sum of estats of inbetween nodes
    private long len=0;
	private static final Logger LOG = LoggerFactory.getLogger(BDPath.class);
    private PopBin uniqueBin;//the unique population bin that this path belongs to (can only be set from outside)
    @Override
    public void add(Edge edge) {
    	super.add(edge);
    	Node lastNode = peekNode();    	
    	len+=((long)lastNode.getNumber("len"))+((BDEdge)edge).getLength();
    	  	
    }
    @Override
    public void setRoot(Node root) {
    	super.setRoot(root);
    	len=(long) root.getNumber("len");
    }
    
	public BDPath(Node root){
		this(root, SimpleBinner.getBinIfUnique(root));
	}
	
	public BDPath(Node root, PopBin bin) {
		super();
		setRoot(root);
		uniqueBin=bin;
	}
	
	public BDPath(BDPath p){
		super();
		if(p!=null && !p.empty()){
			setRoot(p.getRoot());
			for(Edge e:p.getEdgePath())
				add(e);
		}
		diff=p.diff;
		len=p.len;
		pathEstats=p.pathEstats;
		uniqueBin=p.uniqueBin;
	}
	
	//This constructor is used e.g. to load in contigs.path from SPAdes
	//So no recursive path here (path contains all primitive edges)
	public BDPath(BDGraph graph, String paths){
		super();
//		paths=paths.replace(";", ","); //not yet process path with gap!!!
		String[] comps = paths.split(",");
		if(comps.length<1)
			return;
		String curID = comps[0], nextID;
		boolean curDir = curID.contains("+")?true:false,
				nextDir;
		BDNode curNode = (BDNode) graph.getNode(curID.substring(0,curID.length()-1)),
						nextNode;
		if(curNode==null){
			LOG.info("Node {} already removed from graph! Stop reading path here!",curID.substring(0,curID.length()-1));
			return;
		}
		setRoot(curNode);
		HashMap<PopBin,Long> bin2length= new HashMap<PopBin,Long>();
		PopBin curBin=null, pathBin=null;
		for(int i=1; i<comps.length; i++){
			nextID = comps[i];
			nextDir = nextID.contains("+")?true:false;
			nextNode = (BDNode) graph.getNode(nextID.substring(0,nextID.length()-1));		
			if(nextNode==null){
				LOG.info("Node {} already removed from graph! Stop reading path here!",nextID.substring(0,nextID.length()-1));
				return;
			}
			curBin=SimpleBinner.getBinIfUnique(curNode);
			if(curBin!=null){
				if(!bin2length.containsKey(curBin))
					bin2length.put(curBin, (long) (curNode.getNumber("len")));
				else{
					bin2length.replace(curBin, bin2length.get(curBin) + (long) (curNode.getNumber("len")));
				}
			}
			String edgeID=BDEdge.createID(curNode, nextNode, curDir, !nextDir);
			BDEdge curEdge=(BDEdge) graph.getEdge(edgeID);
			if(curEdge!=null){
				add(curEdge);
				curDir=nextDir;
				curNode=nextNode;
			}else{
				if(HybridAssembler.VERBOSE)
					LOG.error("Graph " + graph.getId() + " doesn't contain " + edgeID);
				break;
			}
		}
		//last time
		curBin=SimpleBinner.getBinIfUnique(curNode);
		if(curBin!=null){
			if(!bin2length.containsKey(curBin))
				bin2length.put(curBin, (long) (curNode.getNumber("len")));
			else{
				bin2length.replace(curBin, bin2length.get(curBin) + (long) (curNode.getNumber("len")));
			}
		}
		long tmpLen=0;
		for(PopBin b:bin2length.keySet()){
			if(bin2length.get(b)>tmpLen)
				pathBin=b;
		}
		uniqueBin=pathBin;
	}
	public BDPath reverse(){
		BDPath rcPath = new BDPath(this.peekNode(), uniqueBin);
		List<Edge> edges = this.getEdgePath();
		for(int i = edges.size()-1; i>=0; i--)
			rcPath.add(edges.get(i));
		rcPath.pathEstats=pathEstats;
		return rcPath;
	}
	//It is not really ID because Path doesn't need an ID
	public String getId(){
		//need to make the Id unique for both sense and antisense spelling???
		BDNode curNode = (BDNode) getRoot();
		if(getEdgeCount()<1)
			return curNode.getId();

		String 	retval=curNode.getId(),
				curDir=((BDEdge) getEdgePath().get(0)).getDir(curNode)?"+":"-";
		retval+=curDir;
		for(Edge e:getEdgePath()){
			curNode=(BDNode) e.getOpposite(curNode);
			retval+=","+curNode.getId();
			curDir=((BDEdge) e).getDir(curNode)?"-":"+"; //note that curNode is target node
			retval+=curDir;
		}

		return retval.trim();
	}
	
	public String toString(){
		return "path:(" + getId() + ")";
	}
	//Get the norm path of edges and nodes from a recursive path 
	public BDPath getPrimitivePath() {	
		BDNode curNode = (BDNode) getRoot(), nextNode=null;
		BDPath retval = new BDPath(curNode);		
		for(Edge e:getEdgePath()) {
			nextNode=(BDNode) e.getOpposite(curNode);
			if(e.hasAttribute("path")) {
				BDPath subPath=(BDPath) e.getAttribute("path");
				if(subPath.getRoot()!=curNode)
					subPath=subPath.reverse();
//				retval=retval.join(subPath);
				retval=retval.join(subPath.getPrimitivePath());
			}
			else
				retval.add(e);
			
			curNode=nextNode;
		}
		
		return retval;
	}
	/*
	 * Return DNA sequence of a path (recursively if need)
	 */
	public Sequence spelling(){
		BDPath realPath=getPrimitivePath();
		
		BDNode curNode = (BDNode) realPath.getRoot();
		Sequence curSeq = (Sequence) curNode.getAttribute("seq");
		if(realPath.getEdgeCount()==0)
			return curSeq;
		
		SequenceBuilder seq = new SequenceBuilder(Alphabet.DNA5(), 1024*1024, toString());
		seq.setDesc(realPath.toString());
		boolean curDir=((BDEdge) realPath.getEdgePath().get(0)).getDir(curNode);
		curSeq = curDir?curSeq:Alphabet.DNA.complement(curSeq);
		//If path is circular: don't need to duplicate the closing node
		if(realPath.getRoot() != realPath.peekNode())
			seq.append(curSeq);
		BDNode nextNode=null;
		for(Edge e:realPath.getEdgePath()){
			nextNode=(BDNode) e.getOpposite(curNode);

			curSeq= (Sequence) nextNode.getAttribute("seq");
			curDir=!((BDEdge) e).getDir(nextNode);
			curSeq = curDir?curSeq:Alphabet.DNA.complement(curSeq);
			

			int overlap=((BDEdge) e).getLength();
			//if length of edge > 0: should add NNNN...NN to seq (in case there are gaps in NGS assembly graph)
			if(overlap < 0)
				seq.append(curSeq.subSequence(-overlap, curSeq.length())); 
			else {
				String filler=new String(new char[overlap]).replace("\0", "N");
				Sequence fillerSeq=new Sequence(Alphabet.DNA5(), filler, "gap");
				seq.append(fillerSeq.concatenate(curSeq));				
				if(HybridAssembler.VERBOSE)
					LOG.error("Edge {} has length={} > 0: filled with Ns", e.getId(), overlap);

			}
			
			curNode=nextNode;
			
		}
		return seq.toSequence();
	}
	
	/*
	 * Get a path with two ends removed
	 */
	public BDPath trimEndingNodes(){
		BDPath realPath=getPrimitivePath();
		if(realPath == null || realPath.getNodeCount() < 3)
			return null;
		List<Edge> edges = realPath.getEdgePath();
		BDPath retval=new BDPath(edges.get(0).getOpposite(getRoot()), uniqueBin);
		for(int i = 1; i<edges.size(); i++){
			retval.add(edges.get(i));
		}
		return retval;
	}
	 /*
	  * Add a path to the current path. The path to be added must start with the last node
	  * of the current path. Return TRUE if the joining valid and succeed.
	  */
	public BDPath join(BDPath newPath) {
		
		if(newPath==null || newPath.getNodeCount() <=1 ){
			return new BDPath(this); 
		}else if(this.size() <=1 && this.getRoot() == newPath.getRoot())
			return new BDPath(newPath);
		
		if(newPath.getRoot() != peekNode()){
			if(HybridAssembler.VERBOSE)
				LOG.error("Cannot join path {} to path {} with disagreed first node: {} != {}", newPath.getId(), this.getId(), newPath.getRoot().getId() ,peekNode().getId());
			return null;
		}
		if(((BDEdge) newPath.getEdgePath().get(0)).getDir((AbstractNode) newPath.getRoot())
			== ((BDEdge) peekEdge()).getDir((AbstractNode) peekNode())){
			if(HybridAssembler.VERBOSE)
				LOG.error("Conflict direction from the first node " + newPath.getRoot().getId());
			return null;
		}
		BDPath retval=new BDPath(this);
		//TODO: need a way to check coverage consistent (or make change less significant bin?)			
		for(Edge e:newPath.getEdgePath()){
			retval.add(e);
		}
		diff+=newPath.diff;
		return retval;
	}
	
	public int getDeviation(){
		return this.diff;
	}
	public void setDeviation(int deviation){
		this.diff=deviation;
	}
	public void updatePathDeviation(int diff){
		this.diff+=diff;
	}
	//set path score actively
	public void setPathEstats(double pathAstat) {
		this.pathEstats=pathAstat;
	}
	
	public void updatePathEstats(double nodeAstat){
		this.pathEstats+=nodeAstat;
	}
	
	public double getPathEstats() {	
		return pathEstats;
	}
	

	public long getLength() {
		return len;
	}
	
	public void setConsensusUniqueBinOfPath(PopBin bin){
		uniqueBin=bin;
	}
	public PopBin getConsensusUniqueBinOfPath(){
		return uniqueBin;
	}
	/*
	 * Scan for the closest distance from a node (to) to an end (from) 
	 * with regards to a predefined value (distance) 
	 */
	public int getClosestDistance(Node from, Node to, boolean direction, int distance){
		int retval=Integer.MAX_VALUE;
		boolean dirOfFrom, dirOfTo;
		BDPath ref=null;
		
		if(from==getRoot()){
			ref=this;
		}else if(from==peekNode()){
			ref=this.reverse();
		}else{
			if(HybridAssembler.VERBOSE)
				LOG.warn("Node {} couldn't be found as one of the end node in path {}!", from.getId(), getId());
			return retval;
		}
		int curDistance=0;
		dirOfFrom = ((BDEdge) ref.getEdgePath().get(0)).getDir((BDNode)from);

		BDNode curNode= (BDNode)from;
		for(Edge e:ref.getEdgePath()){
			curNode=(BDNode) e.getOpposite(curNode);
			curDistance+=((BDEdge) e).getLength();
			if(curNode==to){
				dirOfTo=!((BDEdge) e).getDir((BDNode) curNode);
				if((dirOfFrom == dirOfTo) == direction){
					if(Math.abs(curDistance-distance) < BDGraph.A_TOL || GraphUtil.approxCompare(curDistance, distance)==0){
						if(HybridAssembler.VERBOSE)
							LOG.info("|-> agree distance between node %s, node %s: %d and given distance %d\n",
								from.getId(), to.getId(), curDistance, distance);
						if(Math.abs(retval) > Math.abs(curDistance-distance))
							retval=curDistance-distance;
					}
					else if(HybridAssembler.VERBOSE)
						LOG.info("!-> inconsistence direction between node {}:{}, node {}:{} and given direction {}",
								from.getId(), dirOfFrom?"+":"-", to.getId(), dirOfTo?"+":"-", direction);
				}

			}
			curDistance+=curNode.getNumber("len");
		}
		
		
		return retval;
	}
	
	/**
	 * Get the length-weighted coverage of all marker as an approximation for this path's coverage
	 * @return average depth of this path
	 */
	public double averageCov(){
		int len=0;
		double res=0;
		for(Node n:getNodePath()){
			Sequence seq = (Sequence) n.getAttribute("seq");
			len+=(n==getRoot())?seq.length():seq.length()-BDGraph.getKmerSize();
			res+=seq.length()*n.getNumber("cov");
		}
		return res/len;
	}
	public void revert() {
		// TODO Auto-generated method stub
		
	}
	
	public BDNode getFirstNode(){return (BDNode) getRoot();} 
	public BDNode getLastNode(){return (BDNode) peekNode();}
	public boolean getFirstNodeDirection() throws Exception{
		if(getEdgeCount()<1)
			throw new Exception("Path has no edge!");
		return ((BDEdge)getEdgePath().get(0)).getDir(getFirstNode());
	}
	public boolean getLastNodeDirection() throws Exception{
		if(getEdgeCount()<1)
			throw new Exception("Path has no edge!");
		return ((BDEdge) peekEdge()).getDir(getLastNode());
	}
	
	public String getEndingID() {
		try {
			return BDEdge.createID(getFirstNode(), getLastNode(), getFirstNodeDirection(), getLastNodeDirection());
		} catch (Exception e) {
			return null;
		}
	}
	
	private double getLikelihood(double nc, double bc, long nl, double phred ){
//		Estats =  log [ Prob (X is 1-copy | coverage (X)) / Prob (X is alpha-copy | coverage (X)) ] (alpha < 1.0)
		double Estats=-Math.log(BDGraph.ALPHA)*nc*nl/BDGraph.ILLUMINA_READ_LENGTH + (BDGraph.ALPHA-1.0)*nl*bc/BDGraph.ILLUMINA_READ_LENGTH; 
		return phred*Estats/10;
	}
	//Get likelihood of extending this path to node, based on its current coverage, length 
	//and alignment phred score
	public double getExtendLikelihood(Node node, double phred){
		double 	nc=node.getNumber("cov"),
				bc=(uniqueBin!=null?uniqueBin.estCov:averageCov());
		
		//count number of occurrences for this node in this path and reduce its cov respectively  
		long ncount=nodes().filter(n->n.equals(node)).count();
		nc-=ncount*bc;
		long 	nl=(long) node.getNumber("len");
		
		return getLikelihood(nc, bc, nl, phred);
	}
	//if the alignment not available. just want to have score of single node contained in a path to increase
	public double getExtendLikelihood(Node node){
		return getExtendLikelihood(node, Alignment.MIN_QUAL);//phred doesn't matter, just pick a positive one!
	}
	
	//get likelihood of the whole path: unpractical, not in use!
	public double getPathLikelihood(){
		double retval=Double.MIN_VALUE;
		HashMap<Node,Double> occur = new HashMap<>();
		double bc=(uniqueBin!=null?uniqueBin.estCov:averageCov());
		for(Node n:getNodePath()){
			//not considering 2 tips
			if(n==getRoot()||n==peekNode())
				continue;
			
			double 	nc=n.getNumber("cov");
			if(occur.containsKey(n)){
				nc=occur.get(n);
			}
			occur.put(n, nc-bc);

			long nl=(long) n.getNumber("len");
			if(retval==Double.MIN_VALUE)
				retval=getLikelihood(nc, bc, nl, Alignment.MIN_QUAL);
			else
				retval+=getLikelihood(nc, bc, nl, Alignment.MIN_QUAL);	
				
		}
		return retval;
	}
}
