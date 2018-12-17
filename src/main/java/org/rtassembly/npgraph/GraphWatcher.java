package org.rtassembly.npgraph;

import java.io.IOException;
import java.time.LocalTime;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.concurrent.atomic.AtomicInteger;
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
	
	public GraphWatcher(BDGraph graph) {
		this.inputGraph=graph;
		rtComponents = new ConnectedComponents();
		rtComponents.init(graph);
		rtComponents.setCutAttribute("cut");
		numberOfComponents=rtComponents.getConnectedComponentsCount();
	}

	//Remove nodes with degree <=1 and length || cov low
	//TODO: use for posprocess only???
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

	synchronized private void removeBadComponents() {
		List<Node> 	removeNodes=new ArrayList<Node>();
		
		for (Iterator<ConnectedComponent> compIter = rtComponents.iterator(); compIter.hasNext(); ) {
			ConnectedComponent comp = compIter.next();
			AtomicDouble lengthWeightedCov = new AtomicDouble(0.0);
			AtomicInteger length = new AtomicInteger(0);
			comp.nodes().forEach(n->{
				int len = (int) (n.getNumber("len")-BDGraph.getKmerSize());
				length.getAndAdd(len);
				lengthWeightedCov.getAndAdd(n.getNumber("cov")*len);
			});
			double aveCov=lengthWeightedCov.get()/length.get();
			if(GraphUtil.approxCompare(aveCov, inputGraph.binner.leastBin.estCov) < 0 || length.get() < SimpleBinner.ANCHOR_CTG_LEN)
				comp.nodes().forEach(n->removeNodes.add(n));
				
		}
		//Remove abundant components here
		removeNodes.stream().forEach(n->inputGraph.removeNode(n));
	}
	
	
	/*
	 * TODO: replace linearComponentsDecomposition() with this + real-time + threads...
	 * Update the outputGraph to show statistics and current output
	 * Should merge with updating the GUI (colors, labels...)???
	 */
	synchronized void update(boolean lastTime) {
		//cleaning...
		removeBadComponents();
//		if(last)
//			cleanInsignificantNodes();
		
		cutEdges = new HashSet<BDEdge>();
		//then set the cut edges: just for outputGraph stats (will reset after)
		inputGraph.nodes()
		.forEach(n->{
			if(n.getInDegree()>=2)
				n.enteringEdges().forEach(e->{if(lastTime) e.setAttribute("ui.hide");e.setAttribute("cut");cutEdges.add((BDEdge) e);});
			if(n.getOutDegree()>=2)
				n.leavingEdges().forEach(e->{if(lastTime) e.setAttribute("ui.hide");e.setAttribute("cut");cutEdges.add((BDEdge) e);});

		});
		
		if(lastTime){
			//TODO: only remove low cov edges+nodes
			removeBadComponents();
		}
		
		outputGraph=new BDGraph();
		BDPath repPath=null; //representative path of a component
		System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
		System.out.println("Current time: " + LocalTime.now());

		for (Iterator<ConnectedComponent> compIter = rtComponents.iterator(); compIter.hasNext(); ) {
			ConnectedComponent comp = compIter.next();
//			System.out.printf("... id=%s edges=%d nodes=%d \n", comp.id, comp.getEdgeCount(), comp.getNodeCount());
					
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
			ConnectedComponent 	comp0=rtComponents.getConnectedComponentOf(n0),
								comp1=rtComponents.getConnectedComponentOf(n1);
			if(comp0==null || comp1==null)
				continue;
			Node 	nn0=outputGraph.getNode(Integer.toString(comp0.id)),
					nn1=outputGraph.getNode(Integer.toString(comp1.id));
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
		outputGraph.updateStats();
		System.out.printf("Output stats: %d sequences (%d circular) N50=%d N75=%d Max=%d\n", getNumberOfSequences(), getNumberOfCircularSequences(), getN50(), getN75(), getLongestContig());
		if(lastTime)
			System.out.println("FINISH!");
		//reset the cutting attributes
//		inputGraph.edges().filter(e->e.hasAttribute("cut")).forEach(e->{e.removeAttribute("cut"); e.removeAttribute("ui.hide");});
		inputGraph.edges().filter(e->e.hasAttribute("cut")).forEach(e->{e.removeAttribute("cut");});

	}
	
	synchronized public int getN50() {
		return outputGraph==null?inputGraph.n50:outputGraph.n50;
	}
	synchronized public int getN75() {
		return outputGraph==null?inputGraph.n75:outputGraph.n75;
	}
	synchronized public int getLongestContig() {
		return outputGraph==null?inputGraph.maxl:outputGraph.maxl;
	}
	synchronized public int getNumberOfSequences() {
		return outputGraph==null?inputGraph.numOfCtgs:outputGraph.numOfCtgs;
	}
	synchronized public int getNumberOfCircularSequences() {
		return outputGraph==null?inputGraph.numOfCircularCtgs:outputGraph.numOfCircularCtgs;
	}
	synchronized public void outputGFA(String fileName) throws IOException{
		inputGraph.outputGFA(fileName);
	}
	synchronized public void outputFASTA(String fileName) throws IOException {
		outputGraph.outputFASTA(fileName);
	}
}
