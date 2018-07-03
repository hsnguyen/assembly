package org.rtassembly.npgraph;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

import org.graphstream.algorithm.ConnectedComponents;
import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;

import com.google.common.util.concurrent.AtomicDouble;

public class GraphObserver {
	BidirectedGraph inputGraph, outputGraph;
	ConnectedComponents rtComponents;
	HashSet<BidirectedEdge> cutEdges;
	int numberOfComponents=0;
	
	public GraphObserver(BidirectedGraph graph) {
		this.inputGraph=graph;
		rtComponents = new ConnectedComponents();
		rtComponents.init(graph);
		rtComponents.setCutAttribute("cut");

		numberOfComponents=rtComponents.getConnectedComponentsCount();
	}
	synchronized void forFunUpdate() {
//		This is just for fun
//		LOG.info("==========================================================================");
//		LOG.info("\nTotal number of components: {} \ncomponents containing more than 1: {} \nsize of biggest component: {}", 
//					rtComponents.getConnectedComponentsCount(),rtComponents.getConnectedComponentsCount(2),rtComponents.getGiantComponent().size());					    		
//		LOG.info("==========================================================================");    		
		if(numberOfComponents != rtComponents.getConnectedComponentsCount()) {
			numberOfComponents = rtComponents.getConnectedComponentsCount();

    		//Hide components with no markers! Optimize it to work dynamically
    		ArrayList<Node> cleanup = new ArrayList<>();
    		for (Iterator<ConnectedComponents.ConnectedComponent> compIter = rtComponents.iterator(); compIter.hasNext(); ) {
    			ConnectedComponents.ConnectedComponent comp = compIter.next();

    			AtomicDouble lengthWeightedCov = new AtomicDouble(0.0);
    			ArrayList<Node> tmp = new ArrayList<>();
    			comp.nodes().forEach(n->{
    				lengthWeightedCov.getAndAdd(n.getNumber("cov")*(n.getNumber("len")-BidirectedGraph.getKmerSize()));
    				tmp.add(n);
    			});
    			if(lengthWeightedCov.get() < 10000*BidirectedGraph.RCOV)
    				cleanup.addAll(tmp);
    		}
    		for(Node n:cleanup) {
				n.setAttribute("ui.hide");
				n.edges().forEach(e->e.setAttribute("ui.hide"));
////			simGraph.removeNode(n); //this faster but careful here!!!
    		}
		}
	}
	
	synchronized void scanAndUpdate() {
		//1. clean it first
		
		//2. then decompose it (using cut attribute instead of removing edges)
		//reset
		cutEdges = new HashSet<BidirectedEdge>();
		inputGraph.edges().filter(e->e.hasAttribute("cut")).forEach(e->{e.removeAttribute("cut"); e.removeAttribute("ui.hide");});;

		inputGraph.nodes()
//		.filter(n->(n.clean))
		.forEach(n->{
//			if(n.getInDegree()>=2)
//				n.enteringEdges().forEach(e->{e.setAttribute("ui.class","cut");e.setAttribute("cut");cutEdges.add((BidirectedEdge) e);});
//			if(n.getOutDegree()>=2)
//				n.leavingEdges().forEach(e->{e.setAttribute("ui.class","cut");e.setAttribute("cut");cutEdges.add((BidirectedEdge) e);});
		
			if(n.getInDegree()>=2)
				n.enteringEdges().forEach(e->{e.setAttribute("ui.hide");e.setAttribute("cut");cutEdges.add((BidirectedEdge) e);});
			if(n.getOutDegree()>=2)
				n.leavingEdges().forEach(e->{e.setAttribute("ui.hide");e.setAttribute("cut");cutEdges.add((BidirectedEdge) e);});

		});
		
		outputGraph=new BidirectedGraph();
		BidirectedPath repPath=null; //representative path of a component
		for (Iterator<ConnectedComponents.ConnectedComponent> compIter = rtComponents.iterator(); compIter.hasNext(); ) {
			ConnectedComponents.ConnectedComponent comp = compIter.next();
			//check comp: should be linear paths, should start with node+
			 repPath = new BidirectedPath();
			 Node node = comp.nodes().toArray(Node[]::new)[0];
			 repPath.setRoot(node);
			 if(comp.getEdgeCount()>1){
				 //extend to
				 Node curNode=node;
				 boolean curDir=true, isCircular=false;
				 while(curDir?curNode.getOutDegree()==1:curNode.getInDegree()==1){
					 Edge e = curDir?curNode.leavingEdges().toArray(Edge[]::new)[0]:curNode.enteringEdges().toArray(Edge[]::new)[0];
					 repPath.add(e);
					 curNode=e.getOpposite(curNode);
					 curDir=!((BidirectedEdge) e).getDir((BidirectedNode)curNode);
					 if(curNode==node){//circular
						 isCircular=true;
						 break;
					 }
				 }
				 
				 //if linear: reverse
				 if(!isCircular){
					 repPath=repPath.getReversedComplemented();
					 //extend in opposite direction
					 curNode=node;
					 curDir=false;
					 while(curDir?curNode.getOutDegree()==1:curNode.getInDegree()==1){
						 Edge e = curDir?curNode.leavingEdges().toArray(Edge[]::new)[0]:curNode.enteringEdges().toArray(Edge[]::new)[0];
						 repPath.add(e);
						 curNode=e.getOpposite(curNode);
						 curDir=!((BidirectedEdge) e).getDir((BidirectedNode)curNode);
					 }
				 }
				 
			 }
			 //now we have repPath
			 System.out.println(repPath.getId());
		}
		
	}


	synchronized int getNumberOfSequences() {
		
		return 0;
	}
	synchronized double getN50() {
		
		return 0.0;
	}
	synchronized int getNumberOfCircularSequences() {
		
		return 0;
	}
	synchronized void outputGFA() {
		
	}
	synchronized void outputFASTA() {
		
	}
}
