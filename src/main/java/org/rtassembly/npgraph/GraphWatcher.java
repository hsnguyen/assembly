package org.rtassembly.npgraph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.stream.Collectors;

import org.graphstream.algorithm.ConnectedComponents;
import org.graphstream.algorithm.ConnectedComponents.ConnectedComponent;
import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;

import com.google.common.util.concurrent.AtomicDouble;

import japsa.seq.Sequence;

public class GraphWatcher {
	BDGraph inputGraph, outputGraph;
	ConnectedComponents rtComponents;
	HashSet<BDEdge> cutEdges;
	int numberOfComponents=0;
	Set<Integer> rubbish; //save id of insignificant components
	
	public GraphWatcher(BDGraph graph) {
		this.inputGraph=graph;
		rtComponents = new ConnectedComponents();
		rtComponents.init(graph);
		rtComponents.setCutAttribute("cut");
		rubbish=new TreeSet<Integer>();	
		numberOfComponents=rtComponents.getConnectedComponentsCount();
	}

	//Remove nodes with degree <=1 and length || cov low
	//TODO: should just hide them
	private void cleanInsignificantNodes(){
		List<Node> badNodes = inputGraph.nodes()
						.filter(n->(inputGraph.binner.checkRemovableNode(n)))
						.collect(Collectors.toList());
		while(!badNodes.isEmpty()) {
			Node node = badNodes.remove(0);
			List<Node> neighbors = node.neighborNodes().collect(Collectors.toList());	

			inputGraph.removeNode(node);
			neighbors.stream()
				.filter(n->(inputGraph.binner.checkRemovableNode(n)))
				.forEach(n->{if(!badNodes.contains(n)) badNodes.add(n);});
		}
	}

	synchronized boolean checkGoodComponent(ConnectedComponent comp) {
//		LOG.info("==========================================================================");
//		LOG.info("\nTotal number of components: {} \ncomponents containing more than 1: {} \nsize of biggest component: {}", 
//					rtComponents.getConnectedComponentsCount(),rtComponents.getConnectedComponentsCount(2),rtComponents.getGiantComponent().size());					    		
//		LOG.info("==========================================================================");    		

		//Hide components with no markers! Optimize it to work dynamically
		boolean retval=true;
		ArrayList<Node> cleanup = new ArrayList<>();

		AtomicDouble lengthWeightedCov = new AtomicDouble(0.0);
		ArrayList<Node> tmp = new ArrayList<>();
		comp.nodes().forEach(n->{
			lengthWeightedCov.getAndAdd(n.getNumber("cov")*(n.getNumber("len")-BDGraph.getKmerSize()));
			tmp.add(n);
		});
		if(lengthWeightedCov.get() < 10000*inputGraph.binner.leastBin.estCov){
			cleanup.addAll(tmp);
			retval=false;
		}		
		return retval;
	}
	
	
	/*
	 * TODO: replace linearComponentsDecomposition() with this + real-time + threads...
	 * Update the outputGraph to show statistics and current output
	 * Should merge with updating the GUI (colors, labels...)???
	 */
	synchronized void update() {
		//reset
		cutEdges = new HashSet<BDEdge>();
		inputGraph.edges().filter(e->e.hasAttribute("cut")).forEach(e->{e.removeAttribute("cut"); e.removeAttribute("ui.hide");});;

		//the set the cut edges
		inputGraph.nodes()
		.forEach(n->{
			if(n.getInDegree()>=2)
				n.enteringEdges().forEach(e->{e.setAttribute("ui.hide");e.setAttribute("cut");cutEdges.add((BDEdge) e);});
			if(n.getOutDegree()>=2)
				n.leavingEdges().forEach(e->{e.setAttribute("ui.hide");e.setAttribute("cut");cutEdges.add((BDEdge) e);});

		});
		
		outputGraph=new BDGraph();
		BDPath repPath=null; //representative path of a component
		System.out.println("Another round of updating connected components: " + rtComponents.getConnectedComponentsCount());

		for (Iterator<ConnectedComponent> compIter = rtComponents.iterator(); compIter.hasNext(); ) {
			ConnectedComponent comp = compIter.next();
//			System.out.printf("... id=%s edges=%d nodes=%d \n", comp.id, comp.getEdgeCount(), comp.getNodeCount());
			
			if(rubbish.contains(comp.id))
				continue;
			if(!checkGoodComponent(comp)){
				//hide all nodes and edges from this comp
				comp.nodes().forEach(n->n.setAttribute("ui.hide"));
				comp.edges().forEach(e->e.setAttribute("ui.hide"));
				rubbish.add(comp.id);
				continue;
			}
			
			//Start analyzing significant components from here
			//check comp: should be linear paths, should start with node+
			 Node node = comp.nodes().toArray(Node[]::new)[0];
			 repPath = new BDPath(node);
			 boolean isCircular=false;
			 
			 if(comp.getEdgeCount()>=1){
				 //extend to
				 Node curNode=node;
				 boolean curDir=true;
				 List<Edge> ways = (curDir?curNode.leavingEdges():curNode.enteringEdges()).filter(e->!e.hasAttribute("cut")).collect(Collectors.toList());
				 while(ways.size()==1){
					 Edge edge = ways.get(0);
					 repPath.add(edge);
					 curNode=edge.getOpposite(curNode);
					 
					 if(curNode==node){//circular
						 isCircular=true;
						 break;
					 }
					 
					 curDir=!((BDEdge) edge).getDir((BDNode)curNode);
					 ways = (curDir?curNode.leavingEdges():curNode.enteringEdges()).filter(e->!e.hasAttribute("cut")).collect(Collectors.toList());

				 }
				 
				 //if linear: reverse
				 if(!isCircular){
					 repPath=repPath.reverse();
					 //extend in opposite direction
					 curNode=node;
					 curDir=false;
					 ways = (curDir?curNode.leavingEdges():curNode.enteringEdges()).filter(e->!e.hasAttribute("cut")).collect(Collectors.toList());

					 while(ways.size()==1){
						 Edge edge = ways.get(0);
						 repPath.add(edge);
						 curNode=edge.getOpposite(curNode);
						 curDir=!((BDEdge) edge).getDir((BDNode)curNode);
						 ways = (curDir?curNode.leavingEdges():curNode.enteringEdges()).filter(e->!e.hasAttribute("cut")).collect(Collectors.toList());

					 }
				 }
				 
			 }
			 //now we have repPath
			 Sequence seq=repPath.spelling();
			 double cov=GraphUtil.getRealCoverage(repPath.averageCov());
			 Node n=outputGraph.addNode(Integer.toString(comp.id));
			 seq.setName("Contig_"+comp.id+"_"+(isCircular?"circular":"linear")+"_length_"+seq.length()+"_cov_"+cov);
			 n.setAttribute("seq", seq);
			 n.setAttribute("len", seq.length());
			 n.setAttribute("cov",cov);
			 n.setAttribute("path", repPath);
			 if(isCircular)
				 n.setAttribute("circular");
//			 System.out.println("\n" + seq.getName() + ":" + repPath.getId() + "\n=> "+ repPath.getPrimitivePath().getId());

		}
		//now set the edge of outputGraph based on the cut edges
		for(Edge e:cutEdges) {
			Node n0=e.getNode0(), n1=e.getNode1();
			//get corresponding grouped nodes in outputGraph
			Node 	nn0=outputGraph.getNode(Integer.toString(rtComponents.getConnectedComponentOf(n0).id)),
					nn1=outputGraph.getNode(Integer.toString(rtComponents.getConnectedComponentOf(n1).id));
			if(nn0!=null && nn1!=null) {
				boolean dir0=((BDEdge)e).getDir((BDNode)n0),
						dir1=((BDEdge)e).getDir((BDNode)n1);
				if(((BDPath)nn0.getAttribute("path")).getNodeCount()>1) 
					dir0=(n0==((BDPath)nn0.getAttribute("path")).peekNode())?true:false;
				
				if(((BDPath)nn1.getAttribute("path")).getNodeCount()>1) 
					dir1=(n1==((BDPath)nn1.getAttribute("path")).getRoot())?false:true;	
				
				outputGraph.addEdge((BDNode)nn0, (BDNode)nn1 , dir0, dir1);
//				Edge newEdge=outputGraph.addEdge((BDNode)nn0, (BDNode)nn1 , dir0, dir1);
//				System.out.printf("Cut edge %s == New edge %s (%s=%s and %s=%s)\n", e.getId(), newEdge.getId(), 
//							nn0.getId(),((BDPath)nn0.getAttribute("path")).getId(),  
//							nn1.getId(),((BDPath)nn1.getAttribute("path")).getId());
			}
		}
		
		System.out.printf("Output stats: %d sequences (%d circular) N50=%.2f\n", getNumberOfSequences(), getNumberOfCircularSequences(), getN50());
	}
	synchronized int getNumberOfSequences() {
		return outputGraph.getNodeCount();
	}
	synchronized double getN50() {
		return outputGraph.getN50();
	}
	synchronized int getNumberOfCircularSequences() {
		return (int) outputGraph.nodes().filter(n->n.hasAttribute("circular")).count();
	}
	synchronized void outputGFA() {
		//TODO: output GFA of outputGraph or inputGraph?
	}
	synchronized void outputFASTA(String fileName) throws IOException {
		outputGraph.outputFASTA(fileName);
	}
}
