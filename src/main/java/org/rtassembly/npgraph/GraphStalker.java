package org.rtassembly.npgraph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.stream.Collectors;

import org.graphstream.algorithm.ConnectedComponents;
import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;

import com.google.common.util.concurrent.AtomicDouble;

import japsa.seq.Sequence;

public class GraphStalker {
	BidirectedGraph inputGraph, outputGraph;
	ConnectedComponents rtComponents;
	HashSet<BidirectedEdge> cutEdges;
	int numberOfComponents=0;
	
	public GraphStalker(BidirectedGraph graph) {
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
	
	synchronized void linearComponentsDecomposition() {
		//1. clean it first
		cleanGraphByTraversal();
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
					 
					 curDir=!((BidirectedEdge) edge).getDir((BidirectedNode)curNode);
					 ways = (curDir?curNode.leavingEdges():curNode.enteringEdges()).filter(e->!e.hasAttribute("cut")).collect(Collectors.toList());

				 }
				 
				 //if linear: reverse
				 if(!isCircular){
					 repPath=repPath.getReversedComplemented();
					 //extend in opposite direction
					 curNode=node;
					 curDir=false;
					 ways = (curDir?curNode.leavingEdges():curNode.enteringEdges()).filter(e->!e.hasAttribute("cut")).collect(Collectors.toList());

					 while(ways.size()==1){
						 Edge edge = ways.get(0);
						 repPath.add(edge);
						 curNode=edge.getOpposite(curNode);
						 curDir=!((BidirectedEdge) edge).getDir((BidirectedNode)curNode);
						 ways = (curDir?curNode.leavingEdges():curNode.enteringEdges()).filter(e->!e.hasAttribute("cut")).collect(Collectors.toList());

					 }
				 }
				 
			 }
			 //now we have repPath
			 System.out.println(repPath.getId() + " => "+ repPath.getPrimitivePath().getId());
			 Sequence seq=repPath.spelling();
			 double cov=repPath.averageCov();
			 Node n=outputGraph.addNode(Integer.toString(comp.id));
			 seq.setName("Contig_"+comp.id+"_"+(isCircular?"cicurlar":"linear")+"_length_"+seq.length()+"_cov_"+cov);
			 n.setAttribute("seq", seq);
			 n.setAttribute("len", seq.length());
			 n.setAttribute("cov",cov);
			 n.setAttribute("path", repPath);
		}
		//now set the edge of outputGraph based on the cut edges
		for(Edge e:cutEdges) {
			Node n0=e.getNode0(), n1=e.getNode1();
			//get corresponding grouped nodes in outputGraph
			Node 	nn0=outputGraph.getNode(Integer.toString(rtComponents.getConnectedComponentOf(n0).id)),
					nn1=outputGraph.getNode(Integer.toString(rtComponents.getConnectedComponentOf(n1).id));
			if(nn0!=null && nn1!=null) {
				boolean dir0=((BidirectedEdge)e).getDir((BidirectedNode)n0),
						dir1=((BidirectedEdge)e).getDir((BidirectedNode)n1);
				if(((BidirectedPath)nn0.getAttribute("path")).getNodeCount()>1) 
					dir0=(n0==((BidirectedPath)nn0.getAttribute("path")).peekNode())?true:false;
				
				if(((BidirectedPath)nn1.getAttribute("path")).getNodeCount()>1) 
					dir1=(n1==((BidirectedPath)nn1.getAttribute("path")).getRoot())?false:true;	
				
				outputGraph.addEdge((BidirectedNode)nn0, (BidirectedNode)nn1 , dir0, dir1);
//				Edge newEdge=outputGraph.addEdge((BidirectedNode)nn0, (BidirectedNode)nn1 , dir0, dir1);
//				System.out.printf("Cut edge %s == New edge %s (%s=%s and %s=%s)\n", e.getId(), newEdge.getId(), 
//							nn0.getId(),((BidirectedPath)nn0.getAttribute("path")).getId(),  
//							nn1.getId(),((BidirectedPath)nn1.getAttribute("path")).getId());
			}
		}
		
	}

	//Traverse graph to find longest linear path possible, remove redundant paths on the way...
	private void cleanGraphByTraversal(){
		
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
	synchronized void outputFASTA(String fileName) throws IOException {
		outputGraph.outputFASTA(fileName);
	}
}
