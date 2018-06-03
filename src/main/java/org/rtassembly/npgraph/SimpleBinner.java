package org.rtassembly.npgraph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

import org.apache.commons.math3.ml.clustering.Cluster;
import org.apache.commons.math3.ml.clustering.DBSCANClusterer;
import org.apache.commons.math3.ml.clustering.DoublePoint;
import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.collect.Sets;

import japsa.seq.Sequence;


public class SimpleBinner {
    private static final Logger LOG = LoggerFactory.getLogger(SimpleBinner.class);
	public static volatile int SIGNIFICANT_CTG_LEN=10000;
	
	BidirectedGraph graph;
	ArrayList<PopBin> binList;
	HashMap<Edge, HashMap<PopBin,Integer>> edge2BinMap;
	HashMap<Node, HashMap<PopBin, Integer>> node2BinMap;
	ArrayList<Edge> unresolvedEdges;
	SimpleBinner(BidirectedGraph graph){
		this.graph = graph;
		binList = new ArrayList<PopBin>();
		edge2BinMap = new HashMap<Edge, HashMap<PopBin,Integer>>();
		node2BinMap = new HashMap<Node, HashMap<PopBin, Integer>>();
		unresolvedEdges = new ArrayList<Edge>(graph.edges().collect(Collectors.toList()));
	}

	/*
	 * When a node bin is set, traverse the graph to assign edges if possible
	 */
	private void exploringFromNode(Node node){
		HashMap<PopBin, Integer> nbins = node2BinMap.get(node);
		if(nbins == null || nbins.isEmpty())
			return;
		
		ArrayList<Edge> unknownLeavingEdges = new ArrayList<>(),
						unknownEnteringEdges = new ArrayList<>();
		HashMap<PopBin,Integer> 	leavingEdgeBinCount=new HashMap<PopBin,Integer>(),
								enteringEdgeBinCount=new HashMap<PopBin,Integer>();
		// Leaving edges
		node.leavingEdges().forEach(e -> 
		{
			HashMap<PopBin, Integer> ebins = edge2BinMap.get(e); 
			if(ebins == null || ebins.size()==0)
				unknownLeavingEdges.add(e);
			else{
				for(PopBin b:ebins.keySet()) {
					leavingEdgeBinCount.put(b, ebins.get(b)+ (leavingEdgeBinCount.get(b)==null?0:leavingEdgeBinCount.get(b)));
				}
			}
		});
		if(unknownLeavingEdges.size()==1){
			Edge e = unknownLeavingEdges.get(0);
			edge2BinMap.put(e, substract(nbins, leavingEdgeBinCount));
			System.out.println("From node " + node.getId() + " firing edge " + e.getId());
			exploringFromEdge(e);
		}

		// Entering edges
		node.enteringEdges().forEach(e -> 
		{
			HashMap<PopBin, Integer> ebins = edge2BinMap.get(e); 
			if(ebins == null || ebins.size()==0)
				unknownEnteringEdges.add(e);
			else{
				for(PopBin b:ebins.keySet()) {
					enteringEdgeBinCount.put(b, ebins.get(b)+ (enteringEdgeBinCount.get(b)==null?0:enteringEdgeBinCount.get(b)));
				}
			}
		});
		if(unknownEnteringEdges.size()==1){
			Edge e = unknownEnteringEdges.get(0);
			edge2BinMap.put(e, substract(nbins, enteringEdgeBinCount));
			System.out.println("From node " + node.getId() + " firing edge " + e.getId());
			exploringFromEdge(e);
		}
		
	}
	/*
	 * When an edge bin is set, traverse the graph to assign nodes if possible
	 */
	private void exploringFromEdge(Edge edge){
		if(!edge2BinMap.containsKey(edge))
			return;
		unresolvedEdges.remove(edge);
		Node 	n0 = edge.getNode0(),
				n1 = edge.getNode1();
	
		if(!node2BinMap.containsKey(n0)){
			boolean dir0 = ((BidirectedEdge)edge).getDir0();
			Stream<Edge> edgeSet0 = dir0?n0.leavingEdges():n0.enteringEdges();		
			
			ArrayList<Edge> edgeList0 = (ArrayList<Edge>) edgeSet0.collect(Collectors.toList());
			HashMap<PopBin,Integer> binCounts0 = new HashMap<PopBin,Integer>();
			boolean fully = true;
			for(Edge e:edgeList0){
				if(edge2BinMap.containsKey(e)){
					binCounts0=sum(binCounts0, edge2BinMap.get(e));
				}else{
//					System.out.println("...node " + n0.getId() + " not yet fully solved at: " + e.getId());
					fully = false;
					break;
				}
			}
			if(fully){
				node2BinMap.put(n0, binCounts0);
//				System.out.println("From edge " + edge.getId() + " updating node " + n0.getId());
				exploringFromNode(n0);
			}
		}else{
//			System.out.println("From edge " + edge.getId() + " try firing node " + n0.getId());
			exploringFromNode(n0);
		}
		
		
		if(!node2BinMap.containsKey(n1)){
			boolean dir1 = ((BidirectedEdge)edge).getDir1();
			Stream<Edge> edgeSet1 = dir1?n1.leavingEdges():n1.enteringEdges();
			
			ArrayList<Edge> edgeList1 = (ArrayList<Edge>) edgeSet1.collect(Collectors.toList());
			HashMap<PopBin,Integer> binCounts1 = new HashMap<PopBin,Integer>();
			boolean fully = true;
			for(Edge e:edgeList1){
				if(edge2BinMap.containsKey(e)){
					binCounts1=sum(binCounts1, edge2BinMap.get(e));
				}else{
					fully = false;
//					System.out.println("...node " + n1.getId() + " not yet fully solved at: " + e.getId());
					break;
				}
			}
			if(fully){
				node2BinMap.put(n1, binCounts1);
//				System.out.println("From edge " + edge.getId() + " updating node " + n1.getId());
				exploringFromNode(n1);
			}
			
		}else{
//			System.out.println("From edge " + edge.getId() + " try firing node " + n1.getId());
			exploringFromNode(n1);
		}
	}
	
	private static HashMap<PopBin, Integer> substract(HashMap<PopBin, Integer> minuend, HashMap<PopBin, Integer> subtrahend){
		HashMap<PopBin, Integer> difference = new HashMap<PopBin, Integer>();
		for(PopBin b:minuend.keySet()){
			int _minuend = minuend.get(b),
				_diff = _minuend - (subtrahend.containsKey(b)?subtrahend.get(b):0);	
			if(_diff > 0)
				difference.put(b, _diff);
		}
		
		return difference;
	}
	
	private static HashMap<PopBin, Integer> sum(HashMap<PopBin, Integer> augend, HashMap<PopBin, Integer> addend){
		if(augend==null||addend==null)
			return null;
		
		HashMap<PopBin, Integer> retval = new HashMap<PopBin, Integer>();
		for(PopBin b:Sets.union(augend.keySet(), addend.keySet()))
			retval.put(b, (augend.containsKey(b)?augend.get(b):0) + (addend.containsKey(b)?addend.get(b):0));
		
		return retval;
	}
	
	private PopBin scanAndGuess(double cov) {
		PopBin target=null;
		double tmp=100;
		for(PopBin b:binList) {
			double distance=GraphUtil.metric(cov, b.estCov);
			if(distance<tmp) {
				tmp=distance;
				target=b;
			}
		}
		if(tmp>GraphUtil.DISTANCE_THRES) {
			target = null;
		}
		
		return target;
	}

	//A simple clustering with DBSCAN. used internal
	@SuppressWarnings({ "rawtypes", "unchecked" })
	private void nodesClustering() {
		binList = new ArrayList<PopBin>();
		List<DoublePoint> points = new ArrayList<DoublePoint>();
		
		for(Node n:graph) {
			if(n.getNumber("len") >= SIGNIFICANT_CTG_LEN) {
				points.add(new DoublePoint(new double[]{n.getNumber("cov"), new Double(n.getId())}));
			}
		}
		// FIXME: tricky epsilon. need to loop to find best value??? 
		DBSCANClusterer dbscan = new DBSCANClusterer(GraphUtil.DISTANCE_THRES, 0, (a,b)->GraphUtil.metric(a[0], b[0]));

		List<Cluster<DoublePoint>> cluster = dbscan.cluster(points);
		for(Cluster<DoublePoint> c:cluster) {
			PopBin bin=new PopBin();
			for(DoublePoint p:c.getPoints()) {
				Node tmp = graph.getNode(((int)p.getPoint()[1])+"");
				if(tmp.getDegree() <= 2){
					bin.addCoreNode(tmp);
					HashMap<PopBin, Integer> entry = new HashMap<PopBin, Integer>();
					entry.put(bin, 1);
					node2BinMap.put(tmp, entry);
				}
			}
			if(bin.getNodesList().size() > 0){
				binList.add(bin);
			}
			
		}

	}
	
	
	//Now traverse and clustering the edges (+assign cov)
	public void estimatePathsByCoverage() {
		nodesClustering();
		GraphUtil.gradientDescent(graph);
		//1.First round of assigning unit cov: from binned significant nodes
		
		for(PopBin b:binList) {
			for(Node n:b.getNodesList()) {
				exploringFromNode(n);
			}			
		}
		HashMap<PopBin, ArrayList<Edge>> highlyPossibleEdges = new HashMap<PopBin, ArrayList<Edge>>();

		for(Edge e:unresolvedEdges){
			if(e.getNode0().getDegree() <=2 || e.getNode1().getDegree() <=2){
				PopBin tmp = scanAndGuess(e.getNumber("cov"));
				if(tmp!=null){
					if(!highlyPossibleEdges.containsKey(tmp))
						highlyPossibleEdges.put(tmp, new ArrayList<Edge>());
					
					highlyPossibleEdges.get(tmp).add(e);

				}
					
			}
		}
//		2. Second round of thorough assignment: suck it deep!

		while(!unresolvedEdges.isEmpty()) {
			LOG.info("Starting assigning " + unresolvedEdges.size() + " unresolved edges");
			//sort the unresolved edges based on confidence and guess until all gone...
			if(!highlyPossibleEdges.keySet().isEmpty()){
				for(PopBin b:binList){
					if(!highlyPossibleEdges.containsKey(b))
						continue;
					else if(highlyPossibleEdges.get(b).isEmpty()){
						highlyPossibleEdges.remove(b);
						continue;
					}
					else{
						Edge guess = highlyPossibleEdges.get(b).remove(0);
						HashMap<PopBin, Integer> bc = new HashMap<>();
						bc.put(b, 1);
						edge2BinMap.put(guess, bc);
						exploringFromEdge(guess);
					}
				}

			}
			else{
				//we can go further with random guessing but let's stop here for now
				LOG.info("GUESS NO MORE!!!");
				break;
			}
		}
		
	}

//	public boolean isMarker(Node node){
//		boolean retval = false;
//		
//		if(node2BinMap.containsKey(node)){
//			HashMap<SimpleBin, Integer> bc = node2BinMap.get(node);
//			ArrayList<Integer> counts = new ArrayList<Integer>(bc.values());
//			if(counts.size()==1 && counts.get(0)==1){ //and should check for any conflict???
//				if(node.getNumber("len") >= 1000)
//					retval=true;
//			}
//		}
//				
//		return retval;
//	}

	synchronized public PopBin getUniqueBin(Node node){
		PopBin retval = null;
		
		if(node2BinMap.containsKey(node)){
			HashMap<PopBin, Integer> bc = node2BinMap.get(node);
			ArrayList<PopBin> counts = new ArrayList<PopBin>(bc.keySet());
			if(counts.size()==1 && bc.get(counts.get(0))==1){ //and should check for any conflict???
				if(node.getNumber("len") >= 1000)
					retval=counts.get(0);
			}
		}
				
		return retval;
	}
	//Traversal along a unique path (unique ends) and return list of unique edges
	//Must only be called from BidirectedGraph.reduce()
	synchronized public ArrayList<BidirectedEdge> reducedUniquePath(BidirectedPath path) {
		Node  	curNode = path.getRoot(), nextNode;
		PopBin 	uniqueBin=path.getConsensusUniqueBinOfPath(),
				startBin=getUniqueBin(path.getRoot()),
				endBin=getUniqueBin(path.peekNode());;
		if(uniqueBin==null||startBin==null||endBin==null){
			LOG.info("Ignored: population bin of the path (either at ending node or global) is not known!");
			return null;
		}else if(uniqueBin!=startBin && uniqueBin!=endBin){
			LOG.info("Ignored: consensus bin must be one of the endings bin: Ignored!");
			//clean from bin map here...
			return null;
		}else if(!uniqueBin.isCloseTo(startBin)){
			LOG.info("Ignored: consensus bin doesn't agree with one of the endings bin at node {}", path.getRoot());
			node2BinMap.remove(path.getRoot());
			//clean from bin map here...
			return null;
		}else if(!uniqueBin.isCloseTo(endBin)){
			LOG.info("Ignored: consensus bin doesn't agree with one of the endings bin at node {}", path.peekNode());
			node2BinMap.remove(path.peekNode());
			//clean from bin map here...
			return null;
		}
		HashMap<PopBin, Integer> oneBin=new HashMap<>();
		oneBin.put(uniqueBin, 1);

		ArrayList<BidirectedEdge> retval = new ArrayList<BidirectedEdge>();	
		double aveCov=uniqueBin.estCov;
	
		
		for(Edge ep:path.getEdgePath()){
			nextNode=ep.getOpposite(curNode);
			HashMap<PopBin, Integer> edgeBinsCount, bcMinusOne,  
										nodeBinsCount;	

			if(edge2BinMap.containsKey(ep)) {
				edgeBinsCount=edge2BinMap.get(ep);
				if(edgeBinsCount.containsKey(uniqueBin)) {
					bcMinusOne=substract(edgeBinsCount, oneBin);
					if(!bcMinusOne.isEmpty()) {
						edge2BinMap.replace(ep, bcMinusOne);
						
					}
					else {
						//delete here???
						edge2BinMap.remove(ep);
						retval.add((BidirectedEdge) ep);
					}
				}else {//TODO: different due to mis-binning: need to rectify here... (same to node)
				//E.g. b2 vs b1 =>  b2==b1	
					System.out.println("On path: " + path.getId());
					System.out.println("...conflict binning on " + ep.getId() + ": " +edgeBinsCount+ " vs " + oneBin + "("+uniqueBin.estCov+")");
						//...
					
				}
		
				
			}
			
//			LOG.info("--edge {} coverage:{} to {}",ep.getId(),ep.getNumber("cov"),ep.getNumber("cov") - aveCov);
			ep.setAttribute("cov", ep.getNumber("cov") - aveCov);	


//			if(ep.getNumber("cov")/aveCov < .5 && (isMarker(ep.getSourceNode()) || isMarker(ep.getTargetNode())) ) //plasmid coverage is different!!!
//				retval.add((BidirectedEdge) ep);
			
			if(getUniqueBin(curNode)==null) {
				if(node2BinMap.containsKey(curNode) && node2BinMap.get(curNode).containsKey(uniqueBin)) {
					nodeBinsCount=node2BinMap.get(curNode);
					bcMinusOne=substract(nodeBinsCount, oneBin);
					if(!bcMinusOne.isEmpty()) {
						node2BinMap.replace(curNode, bcMinusOne);
					}
					else {
						node2BinMap.remove(curNode);
					}
					
				}				
				
				curNode.setAttribute("cov", curNode.getNumber("cov")>aveCov?curNode.getNumber("cov")-aveCov:0);
				
			}
			
			curNode=nextNode;
		}
		
		return retval;
	}
	
	
	
	public static void main(String[] args) throws IOException {
		HybridAssembler hbAss = new HybridAssembler();
		hbAss.setShortReadsInput(GraphExplore.spadesFolder+"EcK12S-careful/assembly_graph.fastg");
		hbAss.setShortReadsInputFormat("fastg");
		hbAss.prepareShortReadsProcess();
		
		BidirectedGraph graph = hbAss.simGraph;		
		SimpleBinner binner = new SimpleBinner(graph);
		binner.estimatePathsByCoverage();
		System.out.println("=> number of bin = " + binner.binList.size());
		for(PopBin b:binner.binList) {
			System.out.println("Bin " + b.binID + " estCov=" + b.estCov + " totLen=" + b.totLen);
			for(Node n:b.getNodesList())
				System.out.println(n.getAttribute("name"));
		}
			
			
		for(Node n:graph) {
			Iterator<Edge> ite = n.edges().iterator();
			while(ite.hasNext()) {
				Edge e = ite.next();
				System.out.println("Edge "+e.getId() + " cov=" + e.getNumber("cov") );
				HashMap<PopBin,Integer> tmp = binner.edge2BinMap.get(e);
				if(tmp==null) {
					System.out.println("...has not yet assigned!");
					continue;
				}
				for(PopBin b:tmp.keySet()) {
					System.out.printf(" bin %d : %d, ", b.binID, tmp.get(b));
				}
				System.out.println();
			}
		}
	}
}
