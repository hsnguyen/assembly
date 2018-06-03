package org.rtassembly.npgraph;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;


import java.util.List;

import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.graphstream.graph.Path;
import org.graphstream.graph.implementations.AbstractNode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class BidirectedPath extends Path{
	int deviation; //how this path differ to long read data (todo: by multiple-alignment??)
    private long len=0;
	private static final Logger LOG = LoggerFactory.getLogger(BidirectedPath.class);
    private PopBin uniqueBin;//the unique population bin that this path belongs to
    @Override
    public void add(Edge edge) {
    	super.add(edge);
    	Node lastNode = peekNode();
    	
    	len+=((long)lastNode.getNumber("len"))-BidirectedGraph.getKmerSize();
    	  	
    }
    @Override
    public void setRoot(Node root) {
    	super.setRoot(root);
    	len=(long) root.getNumber("len");
    }
    
	public BidirectedPath(){
		super();
	}
	public BidirectedPath(BidirectedEdge e) {
		super();
		setRoot(e.getNode0());
		add(e);
	}
	public BidirectedPath(BidirectedPath p){
		super();
		if(p!=null && !p.empty()){
			setRoot(p.getRoot());
			for(Edge e:p.getEdgePath())
				add(e);
		}
		deviation=p.deviation;
	}
	
	//This constructor is used e.g. to load in contigs.path from SPAdes
	//So no recursive path here (path contains all primitive edges)
	public BidirectedPath(BidirectedGraph graph, String paths){
		super();
		paths=paths.replace(";", ","); //optimized it!
		String[] comps = paths.split(",");
		if(comps.length<1)
			return;
		String curID = comps[0], nextID;
		boolean curDir = curID.contains("+")?true:false,
				nextDir;
		BidirectedNode curNode = (BidirectedNode) graph.getNode(curID.substring(0,curID.length()-1)),
						nextNode;
		setRoot(curNode);
		for(int i=1; i<comps.length; i++){
			nextID = comps[i];
			nextDir = nextID.contains("+")?true:false;
			nextNode = (BidirectedNode) graph.getNode(nextID.substring(0,nextID.length()-1));		
			
			BidirectedEdge curEdge=(BidirectedEdge) graph.getEdge(BidirectedEdge.createID(curNode, nextNode, curDir, !nextDir));

			add(curEdge);
			curDir=nextDir;
			curNode=nextNode;
		}
	}
	public BidirectedPath getReversedComplemented(){
		BidirectedPath rcPath = new BidirectedPath();
		rcPath.setRoot(this.peekNode());
		List<Edge> edges = this.getEdgePath();
		for(int i = edges.size()-1; i>=0; i--)
			rcPath.add(edges.get(i));
		return rcPath;
	}
	//It is not really ID because Path doesn't need an ID
	public String getId(){
		//need to make the Id unique for both sense and antisense spelling???
		BidirectedNode curNode = (BidirectedNode) getRoot();
		if(getEdgeCount()<1)
			return curNode.getId();

		String 	retval=curNode.getId(),
				curDir=((BidirectedEdge) getEdgePath().get(0)).getDir(curNode)?"+":"-";
		retval+=curDir;
		for(Edge e:getEdgePath()){
			curNode=(BidirectedNode) e.getOpposite(curNode);
			retval+=","+curNode.getId();
			curDir=((BidirectedEdge) e).getDir(curNode)?"-":"+"; //note that curNode is target node
			retval+=curDir;
		}

		return retval.trim();
	}
	
	public String toString(){
		return "(" + getId() + ")";
	}
	 
	public Sequence spelling(){
	
		BidirectedNode curNode = (BidirectedNode) getRoot();
		Sequence curSeq = (Sequence) curNode.getAttribute("seq");
	
		if(getEdgeCount()<1)
			return curSeq;
		
		SequenceBuilder seq = new SequenceBuilder(Alphabet.DNA16(), 1024*1024,  this.toString());
		boolean curDir=((BidirectedEdge) getEdgePath().get(0)).getDir(curNode);
		curSeq = curDir?curSeq:Alphabet.DNA.complement(curSeq);
	
		seq.append(curSeq.subSequence(0, curSeq.length()-BidirectedGraph.getKmerSize()));
		for(Edge e:getEdgePath()){
			curNode=(BidirectedNode) e.getOpposite(curNode);
			curSeq= (Sequence) curNode.getAttribute("seq");
			curDir=!((BidirectedEdge) e).getDir(curNode);
			curSeq = curDir?curSeq:Alphabet.DNA.complement(curSeq);
	
			seq.append(curSeq.subSequence(0, curSeq.length()-(curNode==peekNode()?
					0:BidirectedGraph.getKmerSize()))); //TODO: consider cases of overlap != kmer
			
		}
	 return seq.toSequence();
	}
	 /*
	  * Add a path to the current path. The path to be added must start with the last node
	  * of the current path. Return TRUE if the joining valid and succeed.
	  */
	public boolean join(BidirectedPath newPath) {
		if(newPath==null || newPath.size() <=1)
			return true;//old:false
		
		if(newPath.getRoot() != peekNode()){
			LOG.error("Cannot join path with disagreed first node: " + newPath.getRoot().getId() + " != " + peekNode().getId());
			return false;
		}
		if(((BidirectedEdge) newPath.getEdgePath().get(0)).getDir((AbstractNode) newPath.getRoot())
			== ((BidirectedEdge) peekEdge()).getDir((AbstractNode) peekNode())){
			LOG.error("Conflict direction from the first node " + newPath.getRoot().getId());
			return false;
		}
		//TODO: need a way to check coverage consistent

			
		for(Edge e:newPath.getEdgePath()){
			add(e);
		}
		
		return true;
	}
	
	public int getDeviation(){
		return this.deviation;
	}
	public void setDeviation(int deviation){
		this.deviation=deviation;
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
//	/**
//	 * Get the length-weighted coverage of all marker as an approximation for this path's coverage
//	 * @return average depth of this path
//	 */
//	public double averageCov(){
//		int len=0;
//		double res=0;
//		for(Node n:getNodePath()){
//			if(BidirectedGraph.isMarker(n)){
//				Sequence seq = (Sequence) n.getAttribute("seq");
//				len+=(n==getRoot())?seq.length():seq.length()-BidirectedGraph.getKmerSize();
//				res+=seq.length()*n.getNumber("cov");
//			}
//		}
//		return res/len;
//	}
	
}
