package org.rtassembly.npgraph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.apache.commons.math3.ml.clustering.Cluster;
import org.apache.commons.math3.ml.clustering.DBSCANClusterer;
import org.apache.commons.math3.ml.clustering.DoublePoint;
import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;

import japsa.seq.Sequence;

public class CoverageBinner {
	BidirectedGraph graph;
	ArrayList<Bin> binList;
	HashMap<Edge, HashMap<Bin,Integer>> edge2BinMap;
	HashMap<Bin,ArrayList<Edge>> bin2EdgeMap;
	
	CoverageBinner(BidirectedGraph graph){
		this.graph = graph;
		binList = new ArrayList<Bin>();
		edge2BinMap = new HashMap<Edge, HashMap<Bin,Integer>>();
		bin2EdgeMap = new HashMap<Bin,ArrayList<Edge>>();
	}
	
	private void addToMaps(Edge e, Bin b, int multiplicity) {
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
		//Add to the map: bin -> edges 
		ArrayList<Edge> beEntry=bin2EdgeMap.get(b);
		if(beEntry==null) {
			beEntry=new ArrayList<Edge>();
		}
		if(!beEntry.contains(e))
			beEntry.add(e);
	}

	private Bin scanAdd(double cov) {
		Bin target=null;
		double tmp=100;
		for(Bin b:binList) {
			double distance=GraphUtil.metric(cov, b.estCov);
			if(distance<tmp) {
				tmp=distance;
				target=b;
			}
		}
		if(tmp>GraphUtil.DISTANCE_THRES) {
			target=new Bin();
		}
		for(Node node:graph) {
			if(node.getInDegree()<=1 && node.getOutDegree()<=1 && GraphUtil.metric(node.getNumber("cov"), cov) < GraphUtil.DISTANCE_THRES)
				target.addNode(node);
		}
		if(target.getNodesList().size()>0)
			binList.add(target);
		
		return target;
	}

	//A simple clustering with DBSCAN. used internal
	@SuppressWarnings({ "rawtypes", "unchecked" })
	private void nodesClustering() {
		int minLen=10000;
		List<DoublePoint> points = new ArrayList<DoublePoint>();
		
		for(Node n:graph) {
			if(n.getNumber("len") > minLen) {
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
				bin.addNode(tmp);
			}
			binList.add(bin);
			
		}
		
	}
	
	//Now traverse and clustering the edges (+assign cov)
	public void estimatePathsByCoverage() {
		GraphUtil.gradientDescent(graph);
//		GraphUtil.coverageOptimizer(graph);
		nodesClustering();
		//Store the unresolved nodes(edges)..
		List<Edge> 	unresolvedEdges = 
				graph.edges().sorted((o1,o2)->(int)(o1.getNumber("cov")-o2.getNumber("cov"))).collect(Collectors.toList());
		
		List<Edge>	nextSetOfUnresolvedEdges=null;
		Collections.sort(unresolvedEdges, (o1,o2)->(int)(o1.getNumber("cov")-o2.getNumber("cov")));
		//now do the task here:
		Collections.sort(binList, (o1,o2)->(int)(o1.estCov-o2.estCov));
		//1.First round of assigning unit cov
		
		for(Bin b:binList) {
			ArrayList<Node> tobeRemoved = new ArrayList<Node>();
			for(Node n:b.getNodesList()) {
				if(n.getInDegree() > 1 || n.getOutDegree() > 1) {
					System.err.println("...node " + n.getId() + " excluded from bin " + b.getId());
					tobeRemoved.add(n);
					//create/split to 2 new bins here???
				}else {
					Iterator<Edge> ite=n.edges().iterator();
					while(ite.hasNext()) {
						Edge e=ite.next();
						addToMaps(e, b, 1);
						unresolvedEdges.remove(e);
					}
				}
			}
			for(Node n:tobeRemoved)
				b.removeNode(n);
		}

		
		//2. Second round of thorough assignment: suck it deep!
		int numberOfResolvedEdges=0;
		while(true) {
			System.out.println("=====================================================================");
			System.out.println("Starting assigning " + unresolvedEdges.size() + " unresolved edges");
			nextSetOfUnresolvedEdges = new ArrayList<>();
			for(Edge e:unresolvedEdges) {
				//fuk its hard here!!    		
				//not sure about the order of step
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
	    			nextSetOfUnresolvedEdges.add(e);
	    			System.out.println("Could not resolve edge " + e.getId() + " with estimated cov=" + e.getNumber("cov"));
	    			continue;
	    		}else if(!results0.isEmpty() && results1.isEmpty()) {
	    			for(Bin b:results0.keySet()) 
	    				addToMaps(e, b, results0.get(b));
	    			
	    		}else if(results0.isEmpty() && !results1.isEmpty()) {
	    			for(Bin b:results1.keySet()) 
	    				addToMaps(e, b, results1.get(b));
	    		}else {
	    			boolean conflict=false;
	    			for(Bin b:results0.keySet()) {
	    				if(results0.get(b)!=results1.get(b)) {
	    					conflict=true;
	    					break;
	    				}
	    			}
	    			if(conflict) {
	    				System.out.println("There is conflict prediction on edge " + e.getId());
	    				nextSetOfUnresolvedEdges.add(e);
	    				continue;
	    			}
	    		}
	    		
//				//2. scan and guess...
//				Bin b = scanAdd(e.getNumber("cov"));
//				if(b!=null) {
////					wait...
////					addToMaps(e, b);
//				}else {					
//		    		
//					//If still cannot resolve it...
//					nextSetOfUnresolvedEdges.add(e);
//				}
			}
			
			numberOfResolvedEdges=nextSetOfUnresolvedEdges.size()-unresolvedEdges.size();
			unresolvedEdges=nextSetOfUnresolvedEdges;
			if(numberOfResolvedEdges==0)
				break;
		}
		
	}
	
	
	public static void main(String[] args) throws IOException {
		HybridAssembler hbAss = new HybridAssembler(GraphExplore.spadesFolder+"EcK12S-careful/assembly_graph.fastg");
		BidirectedGraph graph = hbAss.simGraph;		
		CoverageBinner binner = new CoverageBinner(graph);
		binner.estimatePathsByCoverage();
		for(Node n:graph) {
			Iterator<Edge> ite = n.edges().iterator();
			while(ite.hasNext()) {
				Edge e = ite.next();
				System.out.println("Edge "+e.getId() + ":");
				HashMap<Bin,Integer> tmp = binner.edge2BinMap.get(e);
				for(Bin b:tmp.keySet()) {
					System.out.printf(" bin %d : %d, ", b.binID, tmp.get(b));
				}
				System.out.println();
			}
		}
	}
}

class Bin{
	static int lastID=1;
	int binID;
	double estCov; //or range???
	long totLen;
	ArrayList<Node> nodesIncludedList;
	public Bin() {
		binID=lastID++;
		nodesIncludedList = new ArrayList<Node>();
	
	}
	
	public int getId() {
		return binID;
	}
	public void addNode(Node node) {
		if(!nodesIncludedList.isEmpty() && nodesIncludedList.contains(node))
			return;
		
		nodesIncludedList.add(node);
		estCov=(estCov*totLen+node.getNumber("cov")*node.getNumber("len"))/(node.getNumber("len")+totLen);
		totLen+=node.getNumber("len");
		
	}
	public void removeNode(Node node) {
		if(nodesIncludedList.remove(node)) {
			estCov=(estCov*totLen-node.getNumber("cov")*node.getNumber("len"))/(-node.getNumber("len")+totLen);
			totLen-=node.getNumber("len");
		}else {
			System.err.println("Node "+ node.getId() + " not found to remove!");
		}
	}
	public ArrayList<Node> getNodesList(){
		return nodesIncludedList;
		
	}
	/*
	 * Return true if this bin also covers cov value
	 */
	public boolean isCloseTo(Bin b) {
		return GraphUtil.metric(this.estCov, b.estCov) < GraphUtil.DISTANCE_THRES;
	} 
}