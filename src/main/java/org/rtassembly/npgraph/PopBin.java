package org.rtassembly.npgraph;

import java.util.ArrayList;

import org.graphstream.graph.Node;

public class PopBin{
	static int lastID=1;
	int binID;
	double estCov; //or range???
	long totLen;
	ArrayList<Node> coreNodes;
	public PopBin() {
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
	 * A-stats here?
	 */
	public boolean isCloseTo(PopBin b) {
		return GraphUtil.approxCompare(this.estCov, b.estCov)==0;
	} 
}
