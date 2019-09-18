package org.rtassembly.npgraph;

import japsa.seq.Alphabet;
import japsa.seq.JapsaAnnotation;
import japsa.seq.JapsaFeature;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;

import java.util.HashMap;
import java.util.List;

import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.graphstream.graph.Path;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class BDPath extends Path{
    private static final Logger LOG = LoggerFactory.getLogger(BDPath.class);

	private int diff; //deviation from path to aligned read (now just length): the more the worse
	private double pathEstats=0.0; //sum of estats of inbetween nodes
    private long len=0;
    private PopBin uniqueBin;//the unique population bin that this path belongs to (can only be set from outside)
    @Override
    public void add(Edge edge) {
    	Node lastNode = edge.getOpposite(peekNode());
    	//score calculation is expensive for long path!
    	pathEstats+=getExtendLikelihood(lastNode);
    	len+=((long)lastNode.getNumber("len"))+ ((BDEdge)edge).getLength();
    	super.add(edge);  	
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
				super.add(e);
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
		BDNode 	curNode = (BDNode) graph.getNode(curID.substring(0,curID.length()-1)),
				nextNode;
		if(curNode==null){
			System.out.printf("Node %s already removed from graph! Stop reading path here!\n",curID.substring(0,curID.length()-1));
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
				System.out.printf("Node %s already removed from graph! Stop reading path here!\n",nextID.substring(0,nextID.length()-1));
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

		String 	retval=curNode.getId();
		boolean	curDir=getFirstNodeDirection();

		retval+=curDir?"+":"-";
		for(Edge e:getEdgePath()){
			curNode=(BDNode) e.getOpposite(curNode);
			retval+=","+curNode.getId();
			//if this is infinity loop, take previous direction
			if(((BDEdge) e).getNodeDirection(curNode)!=null)
				curDir=((BDEdge) e).getNodeDirection(curNode);
			retval+=curDir?"-":"+";			
		}

		return retval.trim();
	}
	
	public String toString(){
		return "path:(" + getId() + ")";
	}
	//Get the norm path of edges and nodes from a recursive path 
	//Do not calculate Estats here: very expensive for whole long path
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
	public Sequence spelling(JapsaAnnotation annotation){
		BDPath realPath=getPrimitivePath();
		JapsaFeature feature;
		BDNode curNode = (BDNode) realPath.getRoot();
		Sequence curSeq = (Sequence) curNode.getAttribute("seq");
		if(realPath.getEdgeCount()==0){
			if(annotation!=null){
				feature=new JapsaFeature(1, curSeq.length(),"CONTIG",curSeq.getName(),'+',"");
				feature.addDesc(curSeq.getName()+"+[1->"+curSeq.length()+"]");
				annotation.add(feature);
			}
			return curSeq;
		}
		
		SequenceBuilder seq = new SequenceBuilder(Alphabet.DNA5(), 1024*1024, toString());
		seq.setDesc(realPath.toString());
		boolean curDir=realPath.getFirstNodeDirection();
		curSeq = curDir?curSeq:Alphabet.DNA.complement(curSeq);
		//If path is circular: don't need to duplicate the closing node
		if(realPath.getFirstNode()!=realPath.getLastNode()){ //linear
			seq.append(curSeq);
			if(annotation!=null){
				feature=new JapsaFeature(1, curSeq.length(),"CONTIG",(String)curNode.getAttribute("name"),curDir?'+':'-',"");
				feature.addDesc((String)curNode.getAttribute("name")+(curDir?"+":"-")+"[1->"+curSeq.length()+"]");
				annotation.add(feature);

			}
		}
		BDNode nextNode=null;
		for(Edge e:realPath.getEdgePath()){
			nextNode=(BDNode) e.getOpposite(curNode);

			curSeq= (Sequence) nextNode.getAttribute("seq");
			if(((BDEdge) e).getNodeDirection(nextNode)!=null)
				curDir=!((BDEdge) e).getNodeDirection(nextNode);

			curSeq = curDir?curSeq:Alphabet.DNA.complement(curSeq);
			
			int overlap=((BDEdge) e).getLength();
			if(annotation!=null){
				feature=new JapsaFeature(seq.length()+1, seq.length()+overlap+curSeq.length(),"CONTIG",(String)nextNode.getAttribute("name"),curDir?'+':'-',"");
				feature.addDesc((String)nextNode.getAttribute("name")+(curDir?"+":"-")+"["+ (seq.length()+overlap+1) +"->"+(seq.length()+overlap+curSeq.length())+"]");
				feature.setScore(overlap);
				annotation.add(feature);

			}
			//if length of edge > 0: should add NNNN...NN to seq (in case there are gaps in NGS assembly graph)
			if(overlap < 0){
				seq.append(curSeq.subSequence(-overlap, curSeq.length())); 
			}
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
		for(int i = 1; i<edges.size()-1; i++){
			retval.add(edges.get(i));
		}
		return retval;
	}
	 /*
	  * Add a path to the current path. The path to be added must start with the last node
	  * of the current path. Return TRUE if the joining valid and succeed.
	  */
	public BDPath join(BDPath newPath) {
		
		if(newPath==null || newPath.empty()){
			return new BDPath(this); 
		}else if(this.empty() && this.getRoot() == newPath.getRoot())
			return new BDPath(newPath);
		
		if(newPath.getRoot() != peekNode()){
			if(HybridAssembler.VERBOSE)
				LOG.error("Cannot join path {} to path {} with disagreed first node: {} != {}", newPath.getId(), this.getId(), newPath.getRoot().getId() ,peekNode().getId());

			return null;
		}

		if(newPath.getNodeCount() > 1 && this.getNodeCount() > 1 && newPath.getFirstNodeDirection()==getLastNodeDirection()){
			if(HybridAssembler.VERBOSE)
				LOG.error("Conflict direction from the first node " + newPath.getRoot().getId());

			return null;
		}
		BDPath retval=new BDPath(this);
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
		dirOfTo = dirOfFrom = ref.getFirstNodeDirection();

		BDNode curNode= (BDNode)from;
		for(Edge e:ref.getEdgePath()){
			curNode=(BDNode) e.getOpposite(curNode);
			curDistance+=((BDEdge) e).getLength();
			if(curNode==to){
				if(((BDEdge) e).getNodeDirection((BDNode) curNode) != null)
					dirOfTo=!((BDEdge) e).getNodeDirection((BDNode) curNode);
				if((dirOfFrom == dirOfTo) == direction){
					if(Math.abs(curDistance-distance) < BDGraph.A_TOL || GraphUtil.approxCompare(curDistance, distance)==0){
						if(HybridAssembler.VERBOSE)
							LOG.info("|-> agree distance between node {}, node {}: {} and given distance {}",
								from.getId(), to.getId(), curDistance, distance);
						if(Math.abs(retval) > Math.abs(curDistance-distance))
							retval=curDistance-distance;
					}

				}else if(HybridAssembler.VERBOSE)
					LOG.info("!-> inconsistence direction between node {} ({}), node {} ({}) and given orientation: {}",
										from.getId(), dirOfFrom?"+":"-", to.getId(), dirOfTo?"+":"-", direction?"same":"opposite" );
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
	public boolean getFirstNodeDirection(){
		Boolean retval=null;
		List<Edge> edges=getEdgePath();
		for(int i=0;i<edges.size();i++){
			retval=((BDEdge) edges.get(i)).getNodeDirection(getFirstNode());
			if( retval!=null){
				break;
			}
		}

		return retval==null?true:retval;
	}
	public boolean getLastNodeDirection(){
		Boolean retval=null;
		List<Edge> edges=getEdgePath();
		for(int i=edges.size()-1;i>=0;i--){
			retval=((BDEdge) edges.get(i)).getNodeDirection(getLastNode());
			if( retval!=null){
				break;
			}
		}
		return retval==null?false:retval;
	}
	
	public String getEndingID() {
		return BDEdge.createID(getFirstNode(), getLastNode(), getFirstNodeDirection(), getLastNodeDirection());
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

}
