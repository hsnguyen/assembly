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
	
	SimpleBinner(BidirectedGraph graph){
		this.graph = graph;
		binList = new ArrayList<Bin>();
		edge2BinMap = new HashMap<Edge, HashMap<Bin,Integer>>();
		node2BinMap = new HashMap<Node, HashMap<Bin, Integer>>();
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
	private void assignEdgeBin(Edge e, Bin b, int multiplicity) {
		//Add to map: edge->occurences in bin
		HashMap<Bin, Integer> ebEntry=edge2BinMap.get(e);
		if(ebEntry==null) {
			ebEntry=new HashMap<Bin,Integer>();
			ebEntry.put(b, new Integer(multiplicity));
			edge2BinMap.put(e, ebEntry);
		}else if(ebEntry.get(b) == null){
			ebEntry.put(b, multiplicity);
		}else {
			if(ebEntry.get(b)!=multiplicity) {
				System.out.println("WARNING: conflict edge multiplicity!");
				//ebEntry.put(b, ebEntry.get(b)+multiplicity);
			}
		}

	}
	private void assignNodeBin(Node n, Bin b, int multiplicity) {
		//Add to map: edge->occurences in bin
		HashMap<Bin, Integer> nbEntry=node2BinMap.get(n);
		if(nbEntry==null) {
			nbEntry=new HashMap<Bin,Integer>();
			nbEntry.put(b, new Integer(multiplicity));
			node2BinMap.put(n, nbEntry);
		}else if(nbEntry.get(b) == null){
			nbEntry.put(b, multiplicity);
		}else {
			if(nbEntry.get(b)!=multiplicity) {
				System.out.println("WARNING: conflict node multiplicity!");
				//nbEntry.put(b, ebEntry.get(b)+multiplicity);
			}
		}

	}
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
		if(unknownEnteringEdges.size()==1){
			Edge e = unknownEnteringEdges.get(0);
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
//		nodesClustering();
//		GraphUtil.coverageOptimizer(graph);
		//Store the unresolved nodes(edges)..
		
		Bin bb=getMostSignificantBin();
		List<Edge> 	unresolvedEdges = 
				graph.edges().sorted((o1,o2)->(int)(o1.getNumber("cov")-o2.getNumber("cov"))).collect(Collectors.toList());
		
		List<Edge>	nextSetOfUnresolvedEdges=null;
		Collections.sort(unresolvedEdges, (o1,o2)->(int)(o1.getNumber("cov")-o2.getNumber("cov")));
		//now do the task here:
		Collections.sort(binList, (o1,o2)->(int)(o1.estCov-o2.estCov));
		
		//1.First round of assigning unit cov: from binned significant nodes
		
		for(Bin b:binList) {
			for(Node n:b.getNodesList()) {
				Iterator<Edge> ite=n.edges().iterator();
				while(ite.hasNext()) {
					Edge e=ite.next();
					assignEdgeBin(e, b, 1);
					unresolvedEdges.remove(e);
				}
			}			
		}
		ArrayList<Edge> highlyPossibleEdges = new ArrayList<>();
		for(Edge e:unresolvedEdges){
			if(scanAndGuess(e.getNumber("cov"))==bb && (e.getNode0().getDegree() <=2 || e.getNode1().getDegree() <=2)){
				highlyPossibleEdges.add(e);
			}
		}
		//2. Second round of thorough assignment: suck it deep!
		int numberOfResolvedEdges=0;
		while(true) {
			LOG.info("Starting assigning " + unresolvedEdges.size() + " unresolved edges");
			nextSetOfUnresolvedEdges = new ArrayList<>();
			for(Edge e:unresolvedEdges) {
				//fuk its hard here!!  
				//TODO: should move to update(Edge justUpdatedEdge) to update graph everytime an edge is classified
				if(!assignEdgeBins(e)) 
					nextSetOfUnresolvedEdges.add(e);
				else
					System.out.println("Solved edge " + e + " cov=" + e.getNumber("cov") + " :" + edge2BinMap.get(e).toString());
			}
			
			numberOfResolvedEdges=-nextSetOfUnresolvedEdges.size()+unresolvedEdges.size();
			LOG.info(numberOfResolvedEdges + " edges has been solved! ");
			// now just guess
//			if(numberOfResolvedEdges==0){
//				ArrayList<Edge> tmp = new ArrayList<>();
//				for(Edge e:nextSetOfUnresolvedEdges) {
//					Bin b = scanAndGuess(e.getNumber("cov"));
//					if(b!=null) {
//						addToMaps(e, b, 1);
//						tmp.add(e);
//					}	
//				}
//				for(Edge e:tmp)
//					nextSetOfUnresolvedEdges.remove(e);
//			}
			
			if(nextSetOfUnresolvedEdges.size()==0)
				break;
			else if(nextSetOfUnresolvedEdges.size()==unresolvedEdges.size()){
//				Bin pop = getLeastAbundanceBin();
//				for(Edge e:nextSetOfUnresolvedEdges) {
//					addToMaps(e, pop, (int)((e.getNumber("cov")/pop.estCov)));
//				}
				if(highlyPossibleEdges.size()>0){
					Edge e = highlyPossibleEdges.remove(0);
					assignEdgeBin(e,bb,1);
					nextSetOfUnresolvedEdges.remove(e);
					continue;
				}else{
					LOG.info("IMPOTENT!!!");
					break;
				}

			}
			unresolvedEdges=nextSetOfUnresolvedEdges;
			
		}
		
	}
	private boolean assignEdgeBins(Edge e) {
		//1. considering information from neighbors
		//refer to the old balance() function!
		BidirectedNode n0 = (BidirectedNode) e.getNode0(), n1=(BidirectedNode) e.getNode1();
		boolean dir0 = ((BidirectedEdge) e).getDir0(), dir1 = ((BidirectedEdge) e).getDir1();
		Iterator<Edge> 	in0 = n0.enteringEdges().iterator(), out0 = n0.leavingEdges().iterator(),
						in1 = n1.enteringEdges().iterator(), out1 = n1.leavingEdges().iterator();
		int unknwIn0 = 0, unknwOut0 = 0, unknwIn1  = 0, unknwOut1 = 0;
		HashMap<Bin, Integer> 	inBins0 = new HashMap<Bin, Integer>(),
		 						outBins0 = new HashMap<Bin, Integer>(), 
	 							inBins1 = new HashMap<Bin, Integer>(), 
	 							outBins1 = new HashMap<Bin, Integer>();
		while(in0.hasNext()) {
			HashMap<Bin,Integer> tmp = edge2BinMap.get(in0.next());
			if(tmp!=null) {
				for(Bin key:tmp.keySet()) {
					if(inBins0.get(key)==null)
						inBins0.put(key, tmp.get(key));
					else
						inBins0.replace(key, inBins0.get(key)+tmp.get(key));
				}
			}else {
				unknwIn0++;
			}
		}
		while(out0.hasNext()) {
			HashMap<Bin,Integer> tmp = edge2BinMap.get(out0.next());
			if(tmp!=null) {
				for(Bin key:tmp.keySet()) {
					if(outBins0.get(key)==null)
						outBins0.put(key, tmp.get(key));
					else
						outBins0.replace(key, outBins0.get(key)+tmp.get(key));
				}
			}else
				unknwOut0++;
		}
		
		while(in1.hasNext()) {
			HashMap<Bin,Integer> tmp = edge2BinMap.get(in1.next());
			if(tmp!=null) {
				for(Bin key:tmp.keySet()) {
					if(inBins1.get(key)==null)
						inBins1.put(key, tmp.get(key));
					else
						inBins1.replace(key, inBins1.get(key)+tmp.get(key));
				}
			}else
				unknwIn1++;
		}
		while(out1.hasNext()) {
			HashMap<Bin,Integer> tmp = edge2BinMap.get(out1.next());
			if(tmp!=null) {
				for(Bin key:tmp.keySet()) {
					if(outBins1.get(key)==null)
						outBins1.put(key, tmp.get(key));
					else
						outBins1.replace(key, outBins1.get(key)+tmp.get(key));
				}
			}else
				unknwOut1++;
		}
		//map containing assigned bins and corresponding multiplicities
		HashMap<Bin,Integer> 	results0=new HashMap<Bin,Integer>(),
								results1=new HashMap<Bin,Integer>();
		if(dir0) {
			if(unknwOut0 == 1) {
				if(unknwIn0==0) {
					for(Bin b:inBins0.keySet()) {
						if(inBins0.get(b) != outBins0.get(b))
							results0.put(b, inBins0.get(b)-(outBins0.get(b)==null?0:outBins0.get(b)));
					}
				}
			}
		}else {
			if(unknwIn0 == 1) {
				if(unknwOut0==0) {
					for(Bin b:outBins0.keySet()) {
						if(outBins0.get(b) != inBins0.get(b))
							results0.put(b, outBins0.get(b)-(inBins0.get(b)==null?0:inBins0.get(b)));
					}
				}
			}
		}
		
		if(dir1) {
			if(unknwOut1 == 1) {
				if(unknwIn1==0) {
					for(Bin b:inBins1.keySet()) {
						if(inBins1.get(b) != outBins1.get(b))
							results1.put(b, inBins1.get(b)-(outBins1.get(b)==null?0:outBins1.get(b)));
					}
				}
			}
		}else {
			if(unknwIn1 == 1) {
				if(unknwOut1==0) {
					for(Bin b:outBins1.keySet()) {
						if(outBins1.get(b) != inBins1.get(b))
							results1.put(b, outBins1.get(b)-(inBins1.get(b)==null?0:inBins1.get(b)));
					}
				}
				
			}
		}
		
		if(results0.isEmpty() && results1.isEmpty()) {
			LOG.info("Could not resolve edge " + e.getId() + " with estimated cov=" + e.getNumber("cov"));
			return false;
		}else if(!results0.isEmpty() && results1.isEmpty()) {
			for(Bin b:results0.keySet()) 
				assignEdgeBin(e, b, results0.get(b));
			
		}else if(results0.isEmpty() && !results1.isEmpty()) {
			for(Bin b:results1.keySet()) 
				assignEdgeBin(e, b, results1.get(b));
		}else {
			for(Bin b:results0.keySet()) {
				if(results0.get(b)!=results1.get(b)) {
					LOG.info("There is conflict prediction on edge " + e.getId()  + " for bin " + b.binID + ": " + results0.get(b) + " vs " + results1.get(b));
					return false;
				}
				assignEdgeBin(e, b, results0.get(b));
			}
		}

		return true;
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