package org.rtassembly.npgraph;

import java.text.DecimalFormat;
import java.util.ArrayList;

import org.graphstream.graph.Node;

public class PopBin{
	static int lastID=1;
	int binID;
	double estCov; //or range???
	long estLen; //bad
	ArrayList<Node> coreNodes;
	private static DecimalFormat df2 = new DecimalFormat(".##");
	
	public PopBin() {
		binID=lastID++;
		coreNodes = new ArrayList<Node>();
	
	}
	//this constructor for reading from file with assigned binID
	public PopBin(int binID) {
		this.binID=binID;
		coreNodes = new ArrayList<Node>();

	}
	public int getId() {
		return binID;
	}
	public void addCoreNode(Node node) {
		if(!coreNodes.isEmpty() && coreNodes.contains(node))
			return;
		
		coreNodes.add(node);
		estCov=(estCov*estLen+node.getNumber("cov")*node.getNumber("len"))/(node.getNumber("len")+estLen);
		estLen+=node.getNumber("len");
		
	}
	public void removeCoreNode(Node node) {
		if(coreNodes.remove(node)) {
			estCov=(estCov*estLen-node.getNumber("cov")*node.getNumber("len"))/(-node.getNumber("len")+estLen);
			estLen-=node.getNumber("len");
		}else {
			System.err.println("Node "+ node.getId() + " not found to remove!");
		}
	}
	public ArrayList<Node> getCoreNodes(){
		return coreNodes;
		
	}
	public String toString(){
		return "B-" + binID + "(cov=" + df2.format(estCov)  + " totLen=" + estLen +")";
	}
	/*
	 * A-stats here?
	 */
	public static boolean isCloseBins(PopBin a, PopBin b) {
		if(a==null||b==null)
			return false;
		else
			return BDGraph.isMetagenomics||GraphUtil.approxCompare(a.estCov, b.estCov)==0;
	} 

	public static PopBin getDominateBin(PopBin a, PopBin b){
		if(a==null)
			return b;
		else if(b==null)
			return a;
		else
			return a.estLen>b.estLen?a:b;
	}
}
