package org.rtassembly.npgraph;

import java.util.ArrayList;
import java.util.Iterator;

import org.graphstream.algorithm.ConnectedComponents;
import org.graphstream.graph.Node;

import com.google.common.util.concurrent.AtomicDouble;

public class GraphObserver {
	BidirectedGraph graph;
	ConnectedComponents rtComponents;
	int numberOfComponents=0;
	
	public GraphObserver(BidirectedGraph graph) {
		this.graph=graph;
		rtComponents = new ConnectedComponents();
		rtComponents.init(graph);
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
		
	}
	synchronized void componentTraversal(ConnectedComponents.ConnectedComponent component) {
		//1. clean it first
		
		//2. then decompose it 
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
