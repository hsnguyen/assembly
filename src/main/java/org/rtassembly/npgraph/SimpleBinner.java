package org.rtassembly.npgraph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
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


public class SimpleBinner {
    private static final Logger LOG = LoggerFactory.getLogger(SimpleBinner.class);
	public static volatile int SIGNIFICANT_CTG_LEN=10000;
	
	BidirectedGraph graph;
	ArrayList<SimpleBin> binList;
	HashMap<Edge, HashMap<SimpleBin,Integer>> edge2BinMap;
	HashMap<Node, HashMap<SimpleBin, Integer>> node2BinMap;
	ArrayList<Edge> unresolvedEdges;
	SimpleBinner(BidirectedGraph graph){
		this.graph = graph;
		binList = new ArrayList<SimpleBin>();
		edge2BinMap = new HashMap<Edge, HashMap<SimpleBin,Integer>>();
		node2BinMap = new HashMap<Node, HashMap<SimpleBin, Integer>>();
		unresolvedEdges = new ArrayList<Edge>(graph.edges().collect(Collectors.toList()));
	}

	/*
	 * When a node bin is set, traverse the graph to assign edges if possible
	 */
	private void exploringFromNode(Node node){
		HashMap<SimpleBin, Integer> nbins = node2BinMap.get(node);
		if(nbins == null || nbins.isEmpty())
			return;
		
		ArrayList<Edge> unknownLeavingEdges = new ArrayList<>(),
						unknownEnteringEdges = new ArrayList<>();
		HashMap<SimpleBin,Integer> 	leavingEdgeBinCount=new HashMap<SimpleBin,Integer>(),
								enteringEdgeBinCount=new HashMap<SimpleBin,Integer>();
		// Leaving edges
		node.leavingEdges().forEach(e -> 
		{
			HashMap<SimpleBin, Integer> ebins = edge2BinMap.get(e); 
			if(ebins == null || ebins.size()==0)
				unknownLeavingEdges.add(e);
			else{
				for(SimpleBin b:ebins.keySet()) {
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
			HashMap<SimpleBin, Integer> ebins = edge2BinMap.get(e); 
			if(ebins == null || ebins.size()==0)
				unknownEnteringEdges.add(e);
			else{
				for(SimpleBin b:ebins.keySet()) {
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
			HashMap<SimpleBin,Integer> binCounts0 = new HashMap<SimpleBin,Integer>();
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
			HashMap<SimpleBin,Integer> binCounts1 = new HashMap<SimpleBin,Integer>();
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
	
	private static HashMap<SimpleBin, Integer> substract(HashMap<SimpleBin, Integer> minuend, HashMap<SimpleBin, Integer> subtrahend){
		HashMap<SimpleBin, Integer> difference = new HashMap<SimpleBin, Integer>();
		for(SimpleBin b:minuend.keySet()){
			int _minuend = minuend.get(b),
				_diff = _minuend - (subtrahend.containsKey(b)?subtrahend.get(b):0);	
			if(_diff > 0)
				difference.put(b, _diff);
		}
		
		return difference;
	}
	
	private static HashMap<SimpleBin, Integer> sum(HashMap<SimpleBin, Integer> augend, HashMap<SimpleBin, Integer> addend){
		if(augend==null||addend==null)
			return null;
		
		HashMap<SimpleBin, Integer> retval = new HashMap<SimpleBin, Integer>();
		for(SimpleBin b:Sets.union(augend.keySet(), addend.keySet()))
			retval.put(b, (augend.containsKey(b)?augend.get(b):0) + (addend.containsKey(b)?addend.get(b):0));
		
		return retval;
	}
	
	private SimpleBin scanAndGuess(double cov) {
		SimpleBin target=null;
		double tmp=100;
		for(SimpleBin b:binList) {
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
		binList = new ArrayList<SimpleBin>();
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
			SimpleBin bin=new SimpleBin();
			for(DoublePoint p:c.getPoints()) {
				Node tmp = graph.getNode(((int)p.getPoint()[1])+"");
				if(tmp.getDegree() <= 2){
					bin.addCoreNode(tmp);
					HashMap<SimpleBin, Integer> entry = new HashMap<SimpleBin, Integer>();
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
		
		for(SimpleBin b:binList) {
			for(Node n:b.getNodesList()) {
				exploringFromNode(n);
			}			
		}
		HashMap<SimpleBin, ArrayList<Edge>> highlyPossibleEdges = new HashMap<SimpleBin, ArrayList<Edge>>();

		for(Edge e:unresolvedEdges){
			if(e.getNode0().getDegree() <=2 || e.getNode1().getDegree() <=2){
				SimpleBin tmp = scanAndGuess(e.getNumber("cov"));
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
				for(SimpleBin b:binList){
					if(!highlyPossibleEdges.containsKey(b))
						continue;
					else if(highlyPossibleEdges.get(b).isEmpty()){
						highlyPossibleEdges.remove(b);
						continue;
					}
					else{
						Edge guess = highlyPossibleEdges.get(b).remove(0);
						HashMap<SimpleBin, Integer> bc = new HashMap<>();
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

	public boolean isMarker(Node node){
		boolean retval = false;
		
		if(node2BinMap.containsKey(node)){
			HashMap<SimpleBin, Integer> bc = node2BinMap.get(node);
			ArrayList<Integer> counts = new ArrayList<Integer>(bc.values());
			if(counts.size()==1 && counts.get(0)==1){ //and should check for any conflict???
				if(node.getNumber("len") >= 500)
					retval=true;
			}
		}
				
		return retval;
	}

	public static void main(String[] args) throws IOException {
		HybridAssembler hbAss = new HybridAssembler();
		hbAss.setShortReadsInput(GraphExplore.spadesFolder+"W303-careful/assembly_graph.fastg");
		hbAss.setShortReadsInputFormat("fastg");
		hbAss.prepareShortReadsProcess();
		
		BidirectedGraph graph = hbAss.simGraph;		
		SimpleBinner binner = new SimpleBinner(graph);
		binner.estimatePathsByCoverage();
		System.out.println("=> number of bin = " + binner.binList.size());
		for(SimpleBin b:binner.binList) {
			System.out.println("Bin " + b.binID + " estCov=" + b.estCov + " totLen=" + b.totLen);
			for(Node n:b.getNodesList())
				System.out.println(n.getAttribute("name"));
		}
			
			
		for(Node n:graph) {
			Iterator<Edge> ite = n.edges().iterator();
			while(ite.hasNext()) {
				Edge e = ite.next();
				System.out.println("Edge "+e.getId() + " cov=" + e.getNumber("cov") );
				HashMap<SimpleBin,Integer> tmp = binner.edge2BinMap.get(e);
				if(tmp==null) {
					System.out.println("...has not yet assigned!");
					continue;
				}
				for(SimpleBin b:tmp.keySet()) {
					System.out.printf(" bin %d : %d, ", b.binID, tmp.get(b));
				}
				System.out.println();
			}
		}
	}
}

class SimpleBin{
	static int lastID=1;
	int binID;
	double estCov; //or range???
	long totLen;
	ArrayList<Node> coreNodes;
	public SimpleBin() {
		binID=lastID++;
		coreNodes = new ArrayList<Node>();
	
	}
	
	public int getId() {
		return binID;
	}
	public void addCoreNode(Node node) {
		if(!coreNodes.isEmpty() && coreNodes.contains(node))
			return;
		
		coreNodes.add(node);
		estCov=(estCov*totLen+node.getNumber("cov")*node.getNumber("len"))/(node.getNumber("len")+totLen);
		totLen+=node.getNumber("len");
		
	}
	public void removeCoreNode(Node node) {
		if(coreNodes.remove(node)) {
			estCov=(estCov*totLen-node.getNumber("cov")*node.getNumber("len"))/(-node.getNumber("len")+totLen);
			totLen-=node.getNumber("len");
		}else {
			System.err.println("Node "+ node.getId() + " not found to remove!");
		}
	}
	public ArrayList<Node> getNodesList(){
		return coreNodes;
		
	}
	/*
	 * Return true if this bin also covers cov value
	 */
	public boolean isCloseTo(SimpleBin b) {
		return GraphUtil.metric(this.estCov, b.estCov) < GraphUtil.DISTANCE_THRES;
	} 
}