package org.rtassembly.npgraph;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
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

import com.google.common.collect.Iterables;
import com.google.common.collect.Sets;

public class SimpleBinner {
    private static final Logger LOG = LoggerFactory.getLogger(SimpleBinner.class);
	public static volatile int 	UNIQUE_CTG_LEN=10000,
								ANCHOR_CTG_LEN=1000, //surely length for the anchors; shorter anchors wouldn't be used until double-checked
								TRANSFORMED_ANCHOR_CTG_LEN=2000;
	
	BDGraph graph;
	ArrayList<PopBin> binList;
	PopBin leastBin;
	HashMap<Edge, Multiplicity> edge2BinMap;
	static HashMap<Node, Multiplicity> node2BinMap;
	ArrayList<Edge> unresolvedEdges;
	public SimpleBinner(BDGraph graph){
		this.graph = graph;
		binList = new ArrayList<PopBin>();
		leastBin = null;
		edge2BinMap = new HashMap<Edge, Multiplicity>();
		node2BinMap = new HashMap<Node, Multiplicity>();
		unresolvedEdges = new ArrayList<Edge>(graph.edges().collect(Collectors.toList()));
	}

	public SimpleBinner(BDGraph graph, String binFileName) {
		this(graph);
		if(binFileName!=null && !binFileName.isEmpty())
			try {
				readBinningFromFile(binFileName);
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
	}


	/*
	 * When a node bin is set, traverse the graph to assign edges if possible
	 */
	private void exploringFromNode(Node node){
		Multiplicity nbins = node2BinMap.get(node);
		if(nbins == null || nbins.isEmpty())
			return;
		
		ArrayList<Edge> unknownLeavingEdges = new ArrayList<>(),
						unknownEnteringEdges = new ArrayList<>();
		Multiplicity 	leavingEdgeBinCount=new Multiplicity(),
						enteringEdgeBinCount=new Multiplicity(),
						induceEdgeBin=new Multiplicity();
		// Leaving edges
		node.leavingEdges().forEach(e -> 
		{
			Multiplicity ebins = edge2BinMap.get(e); 
			if(ebins == null || ebins.isEmpty())
				unknownLeavingEdges.add(e);
			else
				leavingEdgeBinCount.add(ebins);
			
		});
		if(unknownLeavingEdges.size()==1){
			Edge e = unknownLeavingEdges.get(0);
			induceEdgeBin=nbins.substract(leavingEdgeBinCount);
			if(induceEdgeBin!=null && induceEdgeBin.getSum() > 0){
				edge2BinMap.put(e, induceEdgeBin);
				System.out.printf("From node %s%s firing edge %s%s\n",node.getId(),getBinsOfNode(node), e.getId(), getBinsOfEdge(e));
				exploringFromEdge(e);
			}
		}

		// Entering edges
		node.enteringEdges().forEach(e -> 
		{
			Multiplicity ebins = edge2BinMap.get(e); 
			if(ebins == null || ebins.isEmpty())
				unknownEnteringEdges.add(e);
			else
				enteringEdgeBinCount.add(ebins);

		});
		if(unknownEnteringEdges.size()==1){
			Edge e = unknownEnteringEdges.get(0);
			induceEdgeBin=nbins.substract(enteringEdgeBinCount);
			if(induceEdgeBin!=null && induceEdgeBin.getSum() > 0){
				edge2BinMap.put(e, induceEdgeBin);
				System.out.printf("From node %s%s firing edge %s%s\n",node.getId(),getBinsOfNode(node), e.getId(), getBinsOfEdge(e));
				exploringFromEdge(e);
			}
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
		System.out.printf("From edge %s%s: ", edge.getId(), getBinsOfEdge(edge));
		if(!node2BinMap.containsKey(n0)){
			boolean dir0 = ((BDEdge)edge).getDir0();
			Stream<Edge> edgeSet0 = dir0?n0.leavingEdges():n0.enteringEdges();		
			
			ArrayList<Edge> edgeList0 = (ArrayList<Edge>) edgeSet0.collect(Collectors.toList());
			Multiplicity binCounts0 = new Multiplicity();
			boolean fully = true;
			for(Edge e:edgeList0){
				if(edge2BinMap.containsKey(e)){
					binCounts0.add(edge2BinMap.get(e));
				}else{
					fully = false;
					break;
				}
			}
			//check if the total cov of inferred pop bins agree with its original coverage or not
			double covSum=binCounts0.getCoverage();

			if(fully && GraphUtil.approxCompare(covSum, n0.getNumber("cov"))==0){
				node2BinMap.put(n0, binCounts0);
//				System.out.println("From edge " + edge.getId() + " updating node " + n0.getId());
				System.out.printf("completing node %s%s\n", n0.getId(),getBinsOfNode(n0));
				exploringFromNode(n0);
			}else
				System.out.printf("skip node %s%s\n", n0.getId(),getBinsOfNode(n0));
		}else{
			//TODO: check consistent here			
			System.out.printf("firing node %s%s\n", n0.getId(),getBinsOfNode(n0));
			exploringFromNode(n0);
		}
		
		if(!node2BinMap.containsKey(n1)){
			boolean dir1 = ((BDEdge)edge).getDir1();
			Stream<Edge> edgeSet1 = dir1?n1.leavingEdges():n1.enteringEdges();
			
			ArrayList<Edge> edgeList1 = (ArrayList<Edge>) edgeSet1.collect(Collectors.toList());
			Multiplicity binCounts1 = new Multiplicity();
			boolean fully = true;
			for(Edge e:edgeList1){
				if(edge2BinMap.containsKey(e)){
					binCounts1.add(edge2BinMap.get(e));
				}else{
					fully = false;
//					System.out.println("...node " + n1.getId() + " not yet fully solved at: " + e.getId());
					break;
				}
			}
			//check if the total cov of inferred pop bins agree with its original coverage or not
			double covSum=binCounts1.getCoverage();
			
			if(fully && GraphUtil.approxCompare(covSum, n1.getNumber("cov"))==0){				
				node2BinMap.put(n1, binCounts1);
				System.out.printf("completing node %s%s\n", n1.getId(),getBinsOfNode(n1));
				exploringFromNode(n1);
			}else
				System.out.printf("skip node %s%s\n", n1.getId(),getBinsOfNode(n1));
			
		}else{ 
			//TODO: check consistent here also
			System.out.printf("firing node %s%s\n", n1.getId(),getBinsOfNode(n1));
			exploringFromNode(n1);
		}
	}
	
	private PopBin scanAndGuess(double cov) {
		assert leastBin!=null:"Populations not decided yet!";
		PopBin target=null;
		double 	tmp=Double.MAX_VALUE;
		for(PopBin b:binList) {
			double distance=GraphUtil.metric(cov, b.estCov);
			if(distance<tmp) {
				tmp=distance;
				target=b;
			}
			
		}
		//TODO: 1.5 is important! need to find a way for more robust guess (+self-correct)!!!	
		if(cov > 1.5*leastBin.estCov && tmp>GraphUtil.DISTANCE_THRES) {
			System.out.println("...reason: " + tmp + " > " + GraphUtil.DISTANCE_THRES);
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
			if(	n.getNumber("len") >= UNIQUE_CTG_LEN 
				&& Math.max(n.getInDegree(), n.getOutDegree()) <= 1
				){
				points.add(new DoublePoint(new double[]{n.getNumber("cov"), new Double(n.getId())}));
			}
		}
		DBSCANClusterer dbscan = new DBSCANClusterer(GraphUtil.DISTANCE_THRES, 0, (a,b)->GraphUtil.metric(a[0], b[0]));

		List<Cluster<DoublePoint>> cluster = dbscan.cluster(points);
		for(Cluster<DoublePoint> c:cluster) {
			PopBin bin=new PopBin();
			for(DoublePoint p:c.getPoints()) {
				Node tmp = graph.getNode(((int)p.getPoint()[1])+"");
				bin.addCoreNode(tmp);
				Multiplicity entry = new Multiplicity(bin,1);
				node2BinMap.put(tmp, entry);
				
			}
			if(bin.getCoreNodes().size() > 0){
				binList.add(bin);
			}
			
			if(leastBin==null || bin.estCov<leastBin.estCov)
				leastBin=bin;

		}

	}
	
	//read binning file from metaBAT
	private void readBinningFromFile(String binFileName) throws Exception {

		BufferedReader binReader = new BufferedReader(new FileReader(binFileName));
		HashMap<String, PopBin> bins = new HashMap<>();
		String line="", nodeID="", binID="";
		Node node=null;
		PopBin bin=null;
		//Read contigs from contigs.paths
		while((line=binReader.readLine()) != null){			
			String[] comps= line.split("\\s+");
			binID=comps[1];
			//FIXME: contigs binned "0" from metaBAT is unspecified? unique?
			if(binID.equals("0"))
				continue;
			
			nodeID=GraphUtil.getIDFromName(comps[0]);
			
			if(bins.containsKey(binID))
				bin=bins.get(binID);
			else {
				bin=new PopBin(Integer.parseInt(binID));
				bins.put(binID, bin);
				binList.add(bin);
			}
			
			node = (BDNode) graph.getNode(nodeID);
			bin.addCoreNode(node);
			
			Multiplicity entry = new Multiplicity(bin,1);
			node2BinMap.put(node, entry);
		}
		
		binList.stream().forEach(b->{if(leastBin==null || leastBin.estCov>b.estCov) leastBin=b;});
		binReader.close();
		
	}
	
	//Now traverse and clustering the edges (+assign cov)
	public void estimatePathsByCoverage() {
		//1.first clustering the biggest nodes into populations
		if(binList.isEmpty())
			nodesClustering();
		
		//2. assign edges and nodes coverage based on original nodes coverage and graph topo
		GraphUtil.gradientDescent(graph);
		
		for(PopBin bin:binList) {
			System.out.println("Bin " + bin.binID + ": " + bin.estCov);		
			for(Node n:bin.getCoreNodes())
				System.out.println("...core node " + n.getAttribute("name"));
		}
		graph.edges().forEach(e->System.out.println("Edge " + e.getId() + " cov=" + e.getNumber("cov")));		

		//3.1.First round of assigning unit cov: from binned significant nodes
		for(PopBin b:binList) {
			for(Node n:b.getCoreNodes()) {
				exploringFromNode(n);
			}			
		}

		//3.2. Second round of thorough assignment: suck it deep!
		//if big node have more than 1 edge going in/out and the edges have significant less coverage (.5+.5)
		//then it is due to sequencing error and this node should be unique!!!
		HashMap<PopBin, ArrayList<Edge>> highlyPossibleEdges = new HashMap<PopBin, ArrayList<Edge>>();
		unresolvedEdges.sort((a,b)->Double.compare(a.getNumber("cov"),b.getNumber("cov")));

		for(Edge e:unresolvedEdges){
				System.out.print("...scanning edge " + e.getId() + "cov=" + e.getNumber("cov"));
				PopBin tmp = scanAndGuess(e.getNumber("cov"));
				if(tmp!=null){
					
					//FIXME: need more robust binning
					//should build bridge as stronger end contains weaker end!
					
					if(tmp==leastBin && GraphUtil.approxCompare(e.getNumber("cov"), leastBin.estCov) < 0){
						Node n0=e.getNode0(), n1=e.getNode1();
						Multiplicity unitBinMap = new Multiplicity(leastBin,1);
						
						if(n0.getNumber("len") > UNIQUE_CTG_LEN && getBinIfUnique(n0)==null){
							n0.setAttribute("unique", leastBin);
							node2BinMap.put(n0, unitBinMap);
							System.out.printf(" :node %s is unique but have degree=%d, length=%d\n", n0.getId(), n0.getDegree(), (int)n0.getNumber("len"));
							continue;
						}
						
						if(n1.getNumber("len") > UNIQUE_CTG_LEN && getBinIfUnique(n1)==null){
							n1.setAttribute("unique", leastBin);
							node2BinMap.put(n1, unitBinMap);
							System.out.printf(" :node %s is unique but have degree=%d, length=%d\n", n1.getId(), n1.getDegree(), (int)n1.getNumber("len"));
							continue;
						}
						
					}
					
					if(!highlyPossibleEdges.containsKey(tmp))
						highlyPossibleEdges.put(tmp, new ArrayList<Edge>());
					
					highlyPossibleEdges.get(tmp).add(e);
					System.out.printf(": highly possible in bin %d!", tmp.getId());
				}else
					System.out.print(": none!");
				System.out.println();
		}

		while(!unresolvedEdges.isEmpty()) {
			System.out.println("Starting assigning " + unresolvedEdges.size() + " unresolved edges");
			//sort the unresolved edges based on abundance and guess until all gone...
			if(!highlyPossibleEdges.keySet().isEmpty()){
				for(PopBin b:binList){
					if(highlyPossibleEdges.containsKey(b)){
						while(!highlyPossibleEdges.get(b).isEmpty()) {
							Edge guess = highlyPossibleEdges.get(b).remove(0);
							System.out.print("...assigning " + guess.getId());
							if(unresolvedEdges.contains(guess)){
								Multiplicity bc = new Multiplicity(b,1);
								edge2BinMap.put(guess, bc);
								System.out.println(": start explore");
								exploringFromEdge(guess);
							}else{
								//TODO: check consistent here or just take it?
								System.out.println(": already in bin " + getBinsOfEdge(guess));
							}
						
						}
						highlyPossibleEdges.remove(b);
					}
				}

			}
			else{
				//we can go further with random guessing but let's stop here for now
				LOG.info("GUESS NO MORE!!!");
				break;
			}
		}
		
		//3.3 Assign unique nodes here: need more tricks
		for(Node node:graph) {
			if(	node2BinMap.containsKey(node) && node.getNumber("len") > ANCHOR_CTG_LEN ){ 
				Multiplicity bc = node2BinMap.get(node);
				if(Math.max(node.getInDegree(), node.getOutDegree()) <= 1){ //not true if e.g. sequencing errors inside unique contig
					Set<PopBin> counts = bc.getBinsSet();
	//				if(counts.size()==1 && bc.get(counts.get(0))==1){ //and should check for any conflict???
					int totOcc = bc.getSum();
					if(totOcc==1) //and should check for any conflict???
						node.setAttribute("unique", counts.toArray()[0]);
				} else{ // REMOVE NODES WITH MULTIPLICITY > 3 SINCE THEY'RE NOT SO CONFIDENT 
					if(bc.getSum() > 3)
						node2BinMap.remove(node);
				}
			}else if(node.getNumber("len") > TRANSFORMED_ANCHOR_CTG_LEN //big unbinned nodes
					&& Math.max(node.getInDegree(), node.getOutDegree()) <= 1 //to investigate unique nodes only
					)
				graph.addUnknownNodes((BDNode) node);
		}
		
		
	}
	/******************************************************
	 ************** Utility functions *********************
	 *****************************************************/
	public boolean checkIfBinContainingNode(PopBin bin, Node node){
		if(node2BinMap.containsKey(node)){
			return node2BinMap.get(node).getBinCount(bin) > 0;
		}else
			return GraphUtil.approxCompare(bin.estCov, node.getNumber("cov"))>=0;
	}
	
	//unique node from the beginning
	static public PopBin getBinIfUnique(Node node){
		return (PopBin)node.getAttribute("unique");
	}
	//also take into account nodes that transformed to unique after reduced
	static public PopBin getBinIfUniqueNow(Node node){
		PopBin retval=getBinIfUnique(node);
		if(retval==null && node.getNumber("len") > TRANSFORMED_ANCHOR_CTG_LEN){
			if(node2BinMap.containsKey(node)){
				Multiplicity bc = node2BinMap.get(node);
				if(bc.getSum() == 1){
					retval=Iterables.getOnlyElement(bc.getBinsSet());
				}
			}else //for unknown nodes only
				retval=BDGraph.getUniqueBinFromLongReads((BDNode) node);
		}
		return retval;
		  
	}
	
	static public boolean isPotentialAnchorNode(BDNode node){
		return getBinIfUnique(node)!=null||(BDGraph.isSuspectedNode(node)&&node.getInDegree()<=1&&node.getOutDegree()<=1);
	}
	public boolean checkRemovableNode(Node node) {
//		LOG.info("Checking node {}", node.getAttribute("name"));
		if(node.getInDegree()*node.getOutDegree()!=0 || SimpleBinner.getBinIfUnique(node)!=null) {
//			LOG.info("1 not removable!");
			return false;
		}
		else if(node2BinMap.containsKey(node)) {
			if(node2BinMap.get(node).getSum() != 0){
//				LOG.info("2 not removable!");
				return false;
			}
		}
		else if(GraphUtil.approxCompare(node.getNumber("cov"),leastBin.estCov)>=0) {
//			LOG.info("3 not removable!");
			return false;
		}
//		LOG.info(" removable!");
		return true;
	}
	//Traversal along a unique path (unique ends) and return list of unique edges
	//Must only be called from BDGraph.reduce()
	synchronized public Set<Edge> walkAlongUniquePath(BDPath path) {
		System.out.printf("Path %s being processed based on binning info...\n", path.getId());
		Node  	curNode = path.getRoot(), nextNode;
		PopBin 	consensusBin=path.getConsensusUniqueBinOfPath(),
				startBin=getBinIfUniqueNow(path.getRoot()),
				endBin=getBinIfUniqueNow(path.peekNode());
		
		
		if(consensusBin==null||startBin==null||endBin==null){
			System.err.println("Ignored: population bin of the path (either at ending node or global) is not known!");
			return null;
		}
//		else if(uniqueBin!=startBin && uniqueBin!=endBin){//at least one end must agree with the wholepath bin
//			System.err.println("Ignored: consensus bin must be one of the endings bin: " + uniqueBin.getId() + " != " + startBin.getId() + " or " + endBin.getId() );
//			//clean from bin map here...
//			node2BinMap.remove(path.getRoot());
//			node2BinMap.remove(path.peekNode());
//			return null;
//		}
		else if(!PopBin.isCloseBins(consensusBin,startBin)){
			System.err.printf("Ignored: consensus bin %s doesn't agree with one of the endings bin %s at node %s\n", consensusBin, startBin, path.getRoot());
			node2BinMap.remove(path.getRoot());
			//clean from bin map here...
			return null;
		}else if(!PopBin.isCloseBins(consensusBin,endBin)){
			System.err.printf("Ignored: consensus bin %s doesn't agree with one of the endings bin %s at node %s\n", consensusBin, endBin, path.peekNode());
			node2BinMap.remove(path.peekNode());
			//clean from bin map here...
			return null;
		}
		PopBin other=(consensusBin==startBin?endBin:startBin);
		Multiplicity 	oneBin=new Multiplicity(consensusBin,1),
						otherBin=new Multiplicity(other,1);
		
		Set<Edge> retval = new HashSet<Edge>();	
		double aveCov=consensusBin.estCov;
		Boolean dir=null;
		for(Edge ep:path.getEdgePath()){
			nextNode=ep.getOpposite(curNode);
			Multiplicity edgeBinsCount, bcMinusOne=null,  
							nodeBinsCount;	
			//remove corresponding edges from the unique nodes
			if(getBinIfUniqueNow(curNode)!=null){
					dir=((BDEdge)ep).getNodeDirection((BDNode)curNode);
					if(dir!=null){
						Stream<Edge> streamEdges=dir?curNode.leavingEdges():curNode.enteringEdges();
						streamEdges.forEach(retval::add);
					}

			}
			if(getBinIfUniqueNow(nextNode)!=null){
					dir=((BDEdge)ep).getNodeDirection((BDNode)nextNode);
					if(dir!=null){
						Stream<Edge> streamEdges=dir?nextNode.leavingEdges():nextNode.enteringEdges();
						streamEdges.forEach(retval::add);
					}
			}
			
			if(edge2BinMap.containsKey(ep)) {
				edgeBinsCount=edge2BinMap.get(ep);
				if(edgeBinsCount.getBinCount(consensusBin) > 0) {
					bcMinusOne=edgeBinsCount.substract(oneBin);
					edge2BinMap.replace(ep, bcMinusOne);

				}else if(edgeBinsCount.getBinCount(other) > 0){
				//E.g. b2 vs b1 =>  b2==b1							//...
					bcMinusOne=edgeBinsCount.substract(otherBin);
					edge2BinMap.replace(ep, bcMinusOne);				

				}else {
					System.err.printf("WARNING: not found appropriate binning information of edge %s: %s...remove conflict bin!\n", ep.getId(), getBinsOfEdge(ep));
					edge2BinMap.remove(ep);
				}
		
				
			}
			
//			LOG.info("--edge {} coverage:{} to {}",ep.getId(),ep.getNumber("cov"),ep.getNumber("cov") - aveCov);
			ep.setAttribute("cov", ep.getNumber("cov")>aveCov?ep.getNumber("cov")-aveCov:0);	

//			//Heuristic attempts:
//			if(bcMinusOne!=null && bcMinusOne.values().stream().mapToInt(Integer::intValue).sum() == 0) {
//				retval.add((BDEdge) ep);
//			}
//			
////			if(ep.getNumber("cov") < 0  && !unresolvedEdges.contains(ep)) //plasmid coverage is different!!!
//			if(ep.getNumber("cov") <= BPOP*0.5  && !edge2BinMap.containsKey(ep)) //plasmid coverage is different!!!
//				retval.add((BDEdge) ep);
			
			bcMinusOne=null;
			if(curNode!=path.getRoot() && curNode!=path.peekNode()) {
				if(node2BinMap.containsKey(curNode)) {
					if(node2BinMap.get(curNode).getBinCount(consensusBin) > 0) {
						nodeBinsCount=node2BinMap.get(curNode);
						bcMinusOne=nodeBinsCount.substract(oneBin);
						node2BinMap.replace(curNode, bcMinusOne);

					}else if(node2BinMap.get(curNode).getBinCount(other) > 0) {
						nodeBinsCount=node2BinMap.get(curNode);
						bcMinusOne=nodeBinsCount.substract(otherBin);
						node2BinMap.replace(curNode, bcMinusOne);
//						if(!bcMinusOne.isEmpty()) {
//							node2BinMap.replace(curNode, bcMinusOne);
//						}
//						else {
//							node2BinMap.remove(curNode);
//						}
					}
					
				}				
				
				curNode.setAttribute("cov", curNode.getNumber("cov")>aveCov?curNode.getNumber("cov")-aveCov:0);
//				if(	(bcMinusOne!=null && bcMinusOne.values().stream().mapToInt(Integer::intValue).sum() == 0)
//					|| curNode.getNumber("cov") < .1*aveCov) 
//					curNode.edges().forEach(e->retval.add((BDEdge) e));
				
			}
			
			curNode=nextNode;
		}
		
		//check again
//		retval.removeIf(e->	getUniqueBin(e.getNode0())==null && getUniqueBin(e.getNode0())==null &&
//							!(checkEdgeSafeToRemove(e)));
		return retval;
	}
	
//	private boolean checkEdgeSafeToRemove(BDEdge edge) {
//		boolean retval=true;
//		BDNode 	n0=(BDNode) edge.getNode0(),
//						n1=(BDNode) edge.getNode1();
//		boolean dir0=edge.getDir0(),
//				dir1=edge.getDir1();
//		double remainCov=0.0;
//		Optional<Double> tmp;
//		System.out.printf("Checking edge " + edge.getId() + ": remainCov=");
//		boolean onlyFlag=false;
//		if(getUniqueBin(n0)==null && getUniqueBin(n1)==null){
//			if((dir0?n0.getOutDegree():n0.getInDegree()) == 1){//the only in/out edge for n0
//				tmp=(dir0?n0.enteringEdges():n0.leavingEdges()).map(e->(e.getNumber("cov"))).reduce(Double::sum);
//				if(tmp.isPresent())
//					remainCov+=tmp.get();
//				onlyFlag=true;
//			}
//			if((dir1?n1.getOutDegree():n1.getInDegree()) == 1){//the only in/out edge for n1
//				tmp=(dir1?n1.enteringEdges():n1.leavingEdges()).map(e->(e.getNumber("cov"))).reduce(Double::sum);
//				if(tmp.isPresent())
//					remainCov+=tmp.get();			
//				onlyFlag=true;
//			}
//			if(onlyFlag && remainCov > edge.getNumber("cov"))
//				retval=false;
//		}
//		System.out.println(remainCov + " => " + retval);
//		return retval;
//	}
	
	public String getBinsOfNode(Node node) {
		String retval="[";
		Multiplicity binCount=node2BinMap.get(node);
		if(binCount==null)
			retval+="unknown";
		else
			retval+=binCount.toString();
		
		retval+="]";
		return retval;
	}
	public String getBinsOfEdge(Edge edge) {
		String retval="[";
		Multiplicity binCount=edge2BinMap.get(edge);
		if(binCount==null||binCount.isEmpty())
			retval+="unknown";
		else
			retval+=binCount.toString();
		
		retval+="]";
		return retval;
	}
	public static void main(String[] args) throws IOException {
		HybridAssembler hbAss = new HybridAssembler();
		hbAss.setShortReadsInput("/home/sonhoanghguyen/Projects/scaffolding/data/porecamp/metaSPAdes/assembly_graph.fastg");
		hbAss.setBinReadsInput("/home/sonhoanghguyen/Projects/scaffolding/data/porecamp/metabat/bin");
		hbAss.setShortReadsInputFormat("fastg");
		hbAss.prepareShortReadsProcess();
		
		SimpleBinner binner = hbAss.simGraph.binner;
		//binner.estimatePathsByCoverage();
		System.out.println("=> number of bin = " + binner.binList.size());
		for(PopBin b:binner.binList) {
			System.out.println("Bin " + b.binID + " estCov=" + b.estCov + " totLen=" + b.estLen);
			for(Node n:b.getCoreNodes())
				System.out.println(n.getAttribute("name"));
		}
			
			
		for(Node n:hbAss.simGraph) {
			Iterator<Edge> ite = n.edges().iterator();
			while(ite.hasNext()) {
				Edge e = ite.next();
				System.out.println("Edge "+e.getId() + " cov=" + e.getNumber("cov") );
				Multiplicity tmp = binner.edge2BinMap.get(e);
				if(tmp==null) {
					System.out.println("...has not yet assigned!");
					continue;
				}
				
				System.out.println(tmp);
			}
		}
	}
}
