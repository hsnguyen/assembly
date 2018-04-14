package org.rtassembly.npgraph;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.math3.ml.clustering.Cluster;
import org.apache.commons.math3.ml.clustering.DBSCANClusterer;
import org.apache.commons.math3.ml.clustering.DoublePoint;
import org.graphstream.graph.Edge;
import org.graphstream.graph.Graph;
import org.graphstream.graph.Node;

public class CoverageBinner {
	Graph graph;
	HashMap<HashMap<Edge,Bin>, Integer> edge2BinMap;
	HashMap<Bin,ArrayList<Edge>> bin2EdgeMap;
	
	CoverageBinner(Graph graph){
		this.graph = graph;
		edge2BinMap = new HashMap<HashMap<Edge,Bin>, Integer>();
		bin2EdgeMap = new HashMap<Bin,ArrayList<Edge>>();
	}
	
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public void clustering() {
		int minLen=10000;
		List<DoublePoint> points = new ArrayList<DoublePoint>();
		
		for(Node n:graph) {
			if(n.getNumber("len") > minLen) {
				points.add(new DoublePoint(new double[]{n.getNumber("cov"), new Double(n.getId())}));
			}
		}
		// FIXME: tricky epsilon. need to loop to find best value??? 
		DBSCANClusterer dbscan = new DBSCANClusterer(1.0, 0, (a,b)->(.5*(a[0]-b[0])*(Math.log(a[0]) -Math.log(b[0]))));

		List<Cluster<DoublePoint>> cluster = dbscan.cluster(points);
		int count=0;
		for(Cluster<DoublePoint> c:cluster) {
			System.out.println("Cluster " + count++ + ":");
			for(DoublePoint p:c.getPoints())
				System.out.printf(" Node %d:%.2f,", (int)p.getPoint()[1], p.getPoint()[0]);
			System.out.println();
		}
	}
}

class Bin{
	int binID;
	double estCov; //or range???
	
}