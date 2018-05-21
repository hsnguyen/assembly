package org.rtassembly.npgraph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
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
	ArrayList<Bin> binList;
	HashMap<Edge, HashMap<Bin,Integer>> edge2BinMap;
	HashMap<Node, HashMap<Bin, Integer>> node2BinMap;
	ArrayList<Edge> unresolvedEdges;
	SimpleBinner(BidirectedGraph graph){
		this.graph = graph;
		binList = new ArrayList<Bin>();
		edge2BinMap = new HashMap<Edge, HashMap<Bin,Integer>>();
		node2BinMap = new HashMap<Node, HashMap<Bin, Integer>>();
		unresolvedEdges = new ArrayList<Edge>(graph.edges().collect(Collectors.toList()));
	}
	private Bin getLeastAbundanceBin() {//bin that have least coverage
		Bin retval=null;
		double cov=Double.MAX_VALUE;
		for(Bin b:binList)
			if(b.estCov < cov) {
				cov=b.estCov;
				retval=b;
			}
		return retval;
	}
	private Bin getMostSignificantBin() {//bin that have cover longest fraction in the genome
		Bin retval=null;
		long maxLen=0;
		for(Bin b:binList)
			if(b.totLen > maxLen) {
				maxLen=b.totLen;
				retval=b;
			}
		return retval;
	}
//	private void assignEdgeBin(Edge e, Bin b, int multiplicity) {
//		//Add to map: edge->occurences in bin
//		HashMap<Bin, Integer> ebEntry=edge2BinMap.get(e);
//		if(ebEntry==null) {
//			ebEntry=new HashMap<Bin,Integer>();
//			ebEntry.put(b, new Integer(multiplicity));
//			edge2BinMap.put(e, ebEntry);
//		}else if(ebEntry.get(b) == null){
//			ebEntry.put(b, multiplicity);
//		}else {
//			if(ebEntry.get(b)!=multiplicity) {
//				System.out.println("WARNING: conflict edge multiplicity!");
//				//ebEntry.put(b, ebEntry.get(b)+multiplicity);
//			}
//		}
//
//	}
//	private void assignNodeBin(Node n, Bin b, int multiplicity) {
//		//Add to map: edge->occurences in bin
//		HashMap<Bin, Integer> nbEntry=node2BinMap.get(n);
//		if(nbEntry==null) {
//			nbEntry=new HashMap<Bin,Integer>();
//			nbEntry.put(b, new Integer(multiplicity));
//			node2BinMap.put(n, nbEntry);
//		}else if(nbEntry.get(b) == null){
//			nbEntry.put(b, multiplicity);
//		}else {
//			if(nbEntry.get(b)!=multiplicity) {
//				System.out.println("WARNING: conflict node multiplicity!");
//				//nbEntry.put(b, ebEntry.get(b)+multiplicity);
//			}
//		}
//
//	}
	/*
	 * When a node bin is set, traverse the graph to assign edges if possible
	 */
	private void exploringFromNode(Node node){
		HashMap<Bin, Integer> nbins = node2BinMap.get(node);
		if(nbins == null || nbins.isEmpty())
			return;
		
		ArrayList<Edge> unknownLeavingEdges = new ArrayList<>(),
						unknownEnteringEdges = new ArrayList<>();
		HashMap<Bin,Integer> 	leavingEdgeBinCount=new HashMap<Bin,Integer>(),
								enteringEdgeBinCount=new HashMap<Bin,Integer>();
		// Leaving edges
		node.leavingEdges().forEach(e -> 
		{
			HashMap<Bin, Integer> ebins = edge2BinMap.get(e); 
			if(ebins == null || ebins.size()==0)
				unknownLeavingEdges.add(e);
			else{
				for(Bin b:ebins.keySet()) {
					leavingEdgeBinCount.put(b, ebins.get(b)+ (leavingEdgeBinCount.get(b)==null?0:leavingEdgeBinCount.get(b)));
				}
			}
		});
		if(unknownLeavingEdges.size()==1){
			Edge e = unknownLeavingEdges.get(0);
			edge2BinMap.put(e, substract(nbins, leavingEdgeBinCount));
			exploringFromEdge(e);
		}

		// Entering edges
		node.enteringEdges().forEach(e -> 
		{
			HashMap<Bin, Integer> ebins = edge2BinMap.get(e); 
			if(ebins == null || ebins.size()==0)
				unknownEnteringEdges.add(e);
			else{
				for(Bin b:ebins.keySet()) {
					enteringEdgeBinCount.put(b, ebins.get(b)+ (enteringEdgeBinCount.get(b)==null?0:enteringEdgeBinCount.get(b)));
				}
			}
		});
		if(unknownEnteringEdges.size()==1){
			Edge e = unknownEnteringEdges.get(0);
			edge2BinMap.put(e, substract(nbins, enteringEdgeBinCount));
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
			
			edgeSet0.map(e -> edge2BinMap.get(e))
						.reduce(SimpleBinner::sum)
						.ifPresent(s -> {
							node2BinMap.put(n0, s);
							exploringFromNode(n0);
						}
			);
		}
		
		
		if(!node2BinMap.containsKey(n1)){
			boolean dir1 = ((BidirectedEdge)edge).getDir0();
			Stream<Edge> edgeSet1 = dir1?n1.leavingEdges():n1.enteringEdges();
			
			edgeSet1.map(e -> edge2BinMap.get(e))
						.reduce(SimpleBinner::sum)
						.ifPresent(s -> {
							node2BinMap.put(n1, s);
							exploringFromNode(n1);
						}
			);
		}
	}
	
	private static HashMap<Bin, Integer> substract(HashMap<Bin, Integer> minuend, HashMap<Bin, Integer> subtrahend){
		HashMap<Bin, Integer> difference = new HashMap<Bin, Integer>();
		for(Bin b:minuend.keySet()){
			int _minuend = minuend.get(b),
				_diff = _minuend - (subtrahend.containsKey(b)?subtrahend.get(b):0);	
			if(_diff > 0)
				difference.put(b, _diff);
		}
		
		return difference;
	}
	
	private static HashMap<Bin, Integer> sum(HashMap<Bin, Integer> opt1, HashMap<Bin, Integer> opt2){
		if(opt1==null||opt2==null)
			return null;
		
		HashMap<Bin, Integer> retval = new HashMap<Bin, Integer>();
		for(Bin b:Sets.union(opt1.keySet(), opt2.keySet()))
			retval.put(b, (opt1.containsKey(b)?opt1.get(b):0) + (opt2.containsKey(b)?opt2.get(b):0));
		
		return retval;
	}
	
	private Bin scanAndGuess(double cov) {
		Bin target=null, pop = getLeastAbundanceBin();
		if(cov < pop.estCov)
			return pop;
		
		double tmp=100;
		for(Bin b:binList) {
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
		binList = new ArrayList<Bin>();
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
			Bin bin=new Bin();
			for(DoublePoint p:c.getPoints()) {
				Node tmp = graph.getNode(((int)p.getPoint()[1])+"");
				if(tmp.getDegree() <= 2)
					bin.addNode(tmp);
			}
			if(bin.getNodesList().size() > 0)
				binList.add(bin);
			
		}
		
		Bin mainBin = getMostSignificantBin();
		for(Node n:mainBin.nodesIncludedList)
			n.setAttribute("marked");
	}
	
	//Now traverse and clustering the edges (+assign cov)
	public void estimatePathsByCoverage() {
		nodesClustering();
		GraphUtil.gradientDescent(graph);
		//1.First round of assigning unit cov: from binned significant nodes
		
		for(Bin b:binList) {
			for(Node n:b.getNodesList()) {
				exploringFromNode(n);
			}			
		}
//		ArrayList<Edge> highlyPossibleEdges = new ArrayList<>();
//		for(Edge e:unresolvedEdges){
//			if(scanAndGuess(e.getNumber("cov"))==bb && (e.getNode0().getDegree() <=2 || e.getNode1().getDegree() <=2)){
//				highlyPossibleEdges.add(e);
//			}
//		}
		//2. Second round of thorough assignment: suck it deep!
		while(!unresolvedEdges.isEmpty()) {
			LOG.info("Starting assigning " + unresolvedEdges.size() + " unresolved edges");
			//sort the unresolved edges based on confidence and guess until all gone...
			
		}
		
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
		for(Bin b:binner.binList) {
			System.out.println("Bin " + b.binID + " estCov=" + b.estCov + " totLen=" + b.totLen);
			for(Node n:b.nodesIncludedList)
				System.out.println(n.getAttribute("name"));
		}
			
			
		for(Node n:graph) {
			Iterator<Edge> ite = n.edges().iterator();
			while(ite.hasNext()) {
				Edge e = ite.next();
				System.out.println("Edge "+e.getId() + " cov=" + e.getNumber("cov") );
				HashMap<Bin,Integer> tmp = binner.edge2BinMap.get(e);
				if(tmp==null) {
					System.out.println("...has not yet assigned!");
					continue;
				}
				for(Bin b:tmp.keySet()) {
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
	public boolean isCloseTo(Bin b) {
		return GraphUtil.metric(this.estCov, b.estCov) < GraphUtil.DISTANCE_THRES;
	} 
}