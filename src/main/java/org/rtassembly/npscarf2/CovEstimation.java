package org.rtassembly.npscarf2;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;


public class CovEstimation {
	void initialGuess(BidirectedGraph graph) {
		for(Edge e:graph.getEdgeSet()) {
    		BidirectedNode n0 = e.getNode0(), n1=e.getNode1();
    		boolean dir0 = ((BidirectedEdge) e).getDir0(), dir1 = ((BidirectedEdge) e).getDir1();
    		n0.getInDegree();
    		double guess = 	(n0.getNumber("len")*n0.getNumber("cov")/(dir0?n0.getOutDegree():n0.getInDegree())
    						+n1.getNumber("len")*n1.getNumber("cov")/(dir1?n1.getOutDegree():n1.getInDegree()))
    						/(n0.getNumber("len")+n1.getNumber("len"));
    		
    		e.setAttribute("cov", guess);

    		
		}
	}
	void newton(BidirectedGraph graph) {
		
		int maxIterations=500;
		double epsilon=.01;
		for(int i=0;i<maxIterations;i++) {
//			System.out.println("======================= Iteration " + i + " ==========================");
			HashMap<String,Double> stepMap = new HashMap<>();
			for(Edge e:graph.getEdgeSet()) {
	    		BidirectedNode n0 = e.getNode0(), n1=e.getNode1();
	    		boolean dir0 = ((BidirectedEdge) e).getDir0(), dir1 = ((BidirectedEdge) e).getDir1();
	    		Iterator<Edge> 	ite0 = dir0?n0.getLeavingEdgeIterator():n0.getEnteringEdgeIterator(),
	    						ite1 = dir1?n1.getLeavingEdgeIterator():n1.getEnteringEdgeIterator();
	    		double sum0=0, sum1=0, tmp;
	    		while(ite0.hasNext()) {
	    			tmp=ite0.next().getNumber("cov");
	    			sum0+=Double.isNaN(tmp)?1.0:tmp;
	    		}
	    		while(ite1.hasNext())	    		 {
	    			tmp=ite1.next().getNumber("cov");
	    			sum1+=Double.isNaN(tmp)?1.0:tmp;
	    		}	    		
	    		double value=.5*(n0.getNumber("len")*(sum0-n0.getNumber("cov"))+n1.getNumber("len")*(sum1-n1.getNumber("cov")))/(n0.getNumber("len")+n1.getNumber("len"));
	    			
	    		stepMap.put(e.getId(), value);
			}
			boolean isConverged=true;
			for(Edge e:graph.getEdgeSet()) {
				double delta=stepMap.get(e.getId()),
						curCov=Double.isNaN(e.getNumber("cov"))?1.0:e.getNumber("cov");
				if(Math.abs(delta/curCov) > epsilon) {
					isConverged=false;
				}
				e.setAttribute("cov", curCov-delta);
			}
			if(isConverged) {
				System.out.println("Estimation CONVERGED at iteration " + i + "th");
				break;
			}
		}
	}
	
	

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		HybridAssembler hbAss = new HybridAssembler(GraphExplore.spadesFolder+"Kp13883-careful/assembly_graph.fastg");
		BidirectedGraph graph = hbAss.simGraph;
		CovEstimation est = new CovEstimation();
		
//		est.initialGuess(graph);
//		System.out.println("=========================BEFORE============================");
//		for(Node node:graph) {
//			Iterator<Edge> 	inIte = node.getEnteringEdgeIterator(),
//							outIte = node.getLeavingEdgeIterator();
//			double inCov=0, outCov=0;
//			while(inIte.hasNext()) 
//				inCov+=inIte.next().getNumber("cov");			
//			while(outIte.hasNext()) 
//				outCov+=outIte.next().getNumber("cov");	
//			
//			System.out.println(node.getAttribute("name") + " cov= " + node.getNumber("cov") + " inCov=" + inCov + " outCov=" + outCov);
//
//		}
		
		est.newton(graph);
		
		System.out.println("=========================AFTER============================");
		for(Node node:graph) {
			Iterator<Edge> 	inIte = node.getEnteringEdgeIterator(),
							outIte = node.getLeavingEdgeIterator();
			double inCov=0, outCov=0;
			while(inIte.hasNext()) 
				inCov+=inIte.next().getNumber("cov");			
			while(outIte.hasNext()) 
				outCov+=outIte.next().getNumber("cov");	
			
			System.out.println(node.getAttribute("name") + " len=" + node.getNumber("len") + " cov=" + node.getNumber("cov") + " inCov=" + inCov + " outCov=" + outCov);

		}
		System.out.println("=============================================================");
		
		HashMap<Double,Double> lengthWeightedCoverageDistribution = new HashMap<>();
		for(Node n:graph) {
			if(n.getNumber("len") > 10000)
				lengthWeightedCoverageDistribution.put(n.getNumber("cov"), n.getNumber("cov")*n.getNumber("len"));
		}
		
		List<Map.Entry<Double, Double>> covDistribution = new ArrayList<>(lengthWeightedCoverageDistribution.entrySet());
		Collections.sort(covDistribution, (o1,o2)->o2.getValue().compareTo(o1.getValue()));
		
//		System.out.println(covDistribution);
		ArrayList<Ranger> confidentIntervals=new ArrayList<>();
		double alpha=.1;
		for(Map.Entry<Double,Double> entry:covDistribution) {
//			System.out.println(entry.getKey() + ":" + entry.getValue());
			int k = entry.getKey().intValue();
			ChiSquaredDistribution 	x1dist=new ChiSquaredDistribution(2*k),
									x2dist=new ChiSquaredDistribution(2*k+2);
			
			double  lower=.5*x1dist.inverseCumulativeProbability(alpha/2),
					upper=.5*x2dist.inverseCumulativeProbability(1-.5*alpha);
			Ranger tmp=new Ranger(lower,upper);
			if(confidentIntervals.size()==0)
				confidentIntervals.add(tmp);
			else {
				boolean flag=false;
				for(Ranger r:confidentIntervals) {
					if(r.intersect(tmp)!=null) {
						r=tmp;
						flag=true;
						break;
					}
						
				}
				if(!flag)
					confidentIntervals.add(tmp);
			}
			System.out.println(entry + " :" + lower + " to " + upper);
		}
		
		System.out.println("Clustering result: ");
		System.out.println(confidentIntervals);
	}

}
class Ranger{
	double left, right;
	Ranger(){
		left=right=0.0;
	}
	Ranger(double x, double y){
		left=x; right=y;
	}
	Ranger intersect(Ranger r) {
		Ranger retval=null;
		if(left < r.right && right > r.left) 
			return new Ranger(Math.min(left, r.left), Math.max(right, r.right));
		
		return retval;
	}
	public String toString() {
		return left+"->"+right;
	}
}
