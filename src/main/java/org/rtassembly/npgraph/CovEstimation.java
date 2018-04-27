package org.rtassembly.npgraph;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math3.ml.clustering.*;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.LUDecomposition;
import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.jfree.util.Log;


public class CovEstimation {
	static double alpha=.05; //confident level ~95%
	
	public static void initialGuess(BidirectedGraph graph) {
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
	Ranger getConfidentInterval(double k) {
		ChiSquaredDistribution 	x1dist=new ChiSquaredDistribution(2*k),
				x2dist=new ChiSquaredDistribution(2*k+2);

		double  lower=.5*x1dist.inverseCumulativeProbability(alpha/2),
		upper=.5*x2dist.inverseCumulativeProbability(1-.5*alpha);
		return new Ranger(lower,upper);
	}
	
	void gradientDescent(BidirectedGraph graph) {
		//function: \sum{i}{len_i*deg_i*((\sum{edges_in} - cov_i)^2 + (\sum{edges_out} - cov_i)^2)/4}
		int 	maxIterations=500, 
				eIteCount=0, nIteCount=0;
		double epsilon=.1;
		while(true) {
			nIteCount++;
			eIteCount=0;
//			System.out.println("======================= Iteration " + i + " ==========================");
			//1. Updating edges' coverage			
			while(true) {
				eIteCount++;
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
		    		//gamma_ij=1/(len_i+len_j) -> small enough!
		    		double 	w0=n0.getNumber("len")*n0.getDegree()/2,
		    				w1=n1.getNumber("len")*n1.getDegree()/2;
		    		
		    		double value=.5*(w0*(sum0-n0.getNumber("cov"))+w1*(sum1-n1.getNumber("cov")))/(w0+w1);
	
		    		stepMap.put(e.getId(), value);
				}
				boolean isConverged=true;
				ArrayList<Edge> negative = new ArrayList<Edge>();
				for(Edge e:graph.getEdgeSet()) {
					double delta=stepMap.get(e.getId()),
							curCov=Double.isNaN(e.getNumber("cov"))?1.0:e.getNumber("cov");
					if(Math.abs(delta/curCov) > epsilon) {
						isConverged=false;
					}
					e.setAttribute("cov", curCov-delta);
					if(curCov<=delta)
						negative.add(e);
//					e.setAttribute("cov", curCov>delta? curCov-delta:0);
				}
				if(isConverged || eIteCount >= maxIterations) {
					System.out.println("======== edges coverage CONVERGED at iteration " + eIteCount + "th =========");
					for(Edge e:negative) {
						System.out.println("-negative edge " + e.getId() + " cov=" + e.getNumber("cov"));
					}
					break;
				}
			}
			//2. Updating nodes' coverage
			boolean isConverged=true;
			for(Node n:graph) {
				Iterator<Edge> 	in=n.getEnteringEdgeIterator(),
								out=n.getLeavingEdgeIterator();
				long inWeight=0, outWeight=0;
				double inCov=0, outCov=0;
				while(in.hasNext()) {
					Edge tmp=in.next();
					inWeight+=tmp.getOpposite(n).getNumber("len");
					inCov+=tmp.getNumber("cov");
				}
				while(out.hasNext()) {
					Edge tmp=out.next();
					outWeight+=tmp.getOpposite(n).getNumber("len");
					outCov+=tmp.getNumber("cov");
				}
				double newCovEst=(inCov*inWeight+outCov*outWeight)/(inWeight+outWeight);
				if(Math.abs(newCovEst/n.getNumber("cov")) > epsilon)
					isConverged=false;
				n.setAttribute("cov", newCovEst);
			}
			if(isConverged || nIteCount >= maxIterations) {
//			if(isConverged) {
				System.out.println("Node coverage CONVERGED at iteration " + nIteCount + "th");
				System.out.println("======================================================");
				break;
			}
		}
	}
	
	

	public static void main(String[] args) throws IOException {
		HybridAssembler hbAss = new HybridAssembler(GraphExplore.spadesFolder+"W303-careful/assembly_graph.fastg");
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
		
		/*
		 * Estimate edge coverage by using gradient descent method
		 */
		est.gradientDescent(graph);
		
		System.out.println("=========================AFTER============================");
		int nc=0, ec=0;
		for(Node node:graph) {
			Iterator<Edge> 	inIte = node.getEnteringEdgeIterator(),
							outIte = node.getLeavingEdgeIterator();
			double inCov=0, outCov=0;
			while(inIte.hasNext()) 
				inCov+=inIte.next().getNumber("cov");			
			while(outIte.hasNext()) 
				outCov+=outIte.next().getNumber("cov");	
			if(node.getDegree()>0)
				nc++;
			System.out.println(node.getAttribute("name") + " len=" + node.getNumber("len") + " cov=" + node.getNumber("cov") + " inCov=" + inCov + " outCov=" + outCov);

		}
		System.out.println("Nodes: " + nc + " Edges: " + graph.getEdgeCount());
		/*
		 * Estimate dominant peaks by searching for union of confident intervals of Poisson distributions
		 */
		System.out.println("=============================================================");
		
		HashMap<Double,Double> lengthWeightedCoverageDistribution = new HashMap<>();
		int minLen=10000;
		List<DoublePoint> points = new ArrayList<DoublePoint>();
		
		for(Node n:graph) {
			if(n.getNumber("len") > minLen) {
//				System.out.printf("\t%.2f", n.getNumber("cov"));
				points.add(new DoublePoint(new double[]{n.getNumber("cov"), new Double(n.getId())}));
			}
		}
		
//		System.out.println();
//		for(Node n1:graph) {
//			if(n1.getNumber("len") > minLen) {
//				lengthWeightedCoverageDistribution.put(n1.getNumber("cov"), n1.getNumber("cov")*n1.getNumber("len"));
//				double cov1=n1.getNumber("cov");
//				System.out.printf("%.2f\t", cov1);
//				for(Node n2:graph) {
//					double cov2=n2.getNumber("cov");
//					if(n2.getNumber("len") > minLen)
//						System.out.printf("%.2f\t", .5*(cov1*Math.log(cov1/cov2) + cov2*Math.log(cov2/cov1)));
//				}
//				System.out.println();
//				
//			}
//		}
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
		
		
		/*
		 * Trying to work out coverage components for edges for their fuking multiplicities
		 */
		System.out.println("=============================================================");
		
		ArrayList<Node> unresolvedNodes = new ArrayList<>();
		ArrayList<Node> allNodes = new ArrayList<Node>(graph.getNodeSet());
		Collections.sort(allNodes, (o1,o2)->(int)(((Node)o2).getNumber("len")-((Node)o1).getNumber("len")));
		for(Node node:allNodes) {
//			System.out.println(n.getId() + ":" + n.getNumber("len"));
			Iterator<Edge> 	inIte = node.getEnteringEdgeIterator(),
							outIte = node.getLeavingEdgeIterator();
			int m = node.getInDegree(), n = node.getOutDegree();
			
			if(m*n==0)
				continue;
			ArrayList<Double> 	inComponents = new ArrayList<>(),
								outComponents = new ArrayList<>();
			while(inIte.hasNext()) {
				Edge in = inIte.next();
				
				if(in.getAttribute("covdis")==null) {
					in.setAttribute("covdis", new CoverageDistribution(in.getNumber("cov")));
				}
				CoverageDistribution inCovDist = in.getAttribute("covdis");
				inComponents.addAll(inCovDist.getAllComponents());
				
			}
			while(outIte.hasNext()) {
				Edge out = outIte.next();
				CoverageDistribution outCovDist = out.getAttribute("covdis");
				if(outCovDist!=null)
					outComponents.addAll(outCovDist.getAllComponents());
			}
			
			
			System.out.println("Node " + node.getAttribute("name"));
			System.out.println("IN components: " + inComponents);
			System.out.println("OUT components: " + outComponents);

			
			
			
//			/******************************************************
//			 * Another Newton
//			 * Aij-=(\sum{k}{Akj} + \sum{k}{Aik} - Ii - Oj)
//			 ******************************************************/
//			ArrayList<Edge> 	inComponents = new ArrayList<>(),
//								outComponents = new ArrayList<>();
//			int maxIterations=500, count=0;
//			double epsilon=1;
////			double[][] 	comps = new double[m][n];
////			
//			while(inIte.hasNext()) {
//				Edge in = inIte.next();
//				inComponents.add(in);			
//			}
//			while(outIte.hasNext()) {
//				Edge out = outIte.next();
//				outComponents.add(out);	
//			}
//			boolean isContinue=true;
////			while(isContinue) {
////				double[][] delta = new double[m][n];
////				count++;
////				if(count>maxIterations){
////					break;
////				}
////				System.out.println("Iteration " + count + "th:");
////				for(int i=0;i<m;i++) {
////					for(int j=0;j<n;j++)
////						System.out.printf("%.2f ",comps[i][j]);
////					System.out.println();
////				}
////				
////				isContinue=false;
////				for(int i=0;i<m;i++) {
////					for(int j=0;j<n;j++) {
////						for(int k=0;k<m;k++)
////							delta[i][j]+=comps[k][j];
////						for(int k=0;k<n;k++)
////							delta[i][j]+=comps[i][k];
////						delta[i][j]-=inComponents.get(i).getNumber("cov")+outComponents.get(j).getNumber("cov");
////						delta[i][j]/=2;
////						if(Math.abs(delta[i][j])>epsilon)
////							isContinue=true;
////					}
////				}
////				
////				for(int i=0;i<m;i++) 
////					for(int j=0;j<n;j++)
////						comps[i][j]-=delta[i][j];
////				
////			}
////			
////			System.out.println("Node " + node.getAttribute("name") + " comps estimate after " + count + " iterations!");
////			for(int i=0;i<m;i++) {
////				Edge e = inComponents.get(i);
////				System.out.print("...edge " + e.getId() + " cov=" + e.getNumber("cov") + " components: ");
////				for(int j=0;j<n;j++)
////					System.out.print(comps[i][j] + "; ");
////				System.out.println();
////			}
////			for(int j=0;j<n;j++) {
////				Edge e = outComponents.get(j);
////				System.out.print("...edge " + e.getId() + " cov=" + e.getNumber("cov") + " components: ");
////				for(int i=0;i<m;i++)
////					System.out.print(comps[i][j] + "; ");
////				System.out.println();
////			}
//			
//			/******************************************************************************
//			 * Real Newton from wiki: doesn't work with this function!!!
//			 ******************************************************************************/
//			
//			//First build Hessian matrix from submatrices A and B
//			//		A 	B ... B
//			//		B	A ... B
//			// H = 	...
//			//		B 	B ... A
//			
//			double[][] 	A=new double[n][n],
//						B=new double[n][n];
//			for(int i=0;i<n;i++) {
//				for(int j=0;j<n;j++) {
//					if(i==j) {
//						A[i][j]=4;
//						B[i][j]=1;
//					}else {
//						A[i][j]=1;
//						B[i][j]=0;
//					}
//				}
//			}
//			
//			RealMatrix H = MatrixUtils.createRealMatrix(m*n,m*n);
//			for(int i=0;i<m;i++) {
//				for(int j=0;j<m;j++) {
//					if(i==j) {
//						H.setSubMatrix(A, i*n, j*n);
//					}else {
//						H.setSubMatrix(B, i*n, j*n);
//					}
//				}
//			}
//			// Now the gradient
//			double[] 	grad = new double[m*n],
//						xn = new double[m*n];
//
//			while(isContinue) {
//				count++;
//				if(count > maxIterations)
//					break;
//				//set the gradient at xn
//				for(int k=0;k<m*n;k++) {
//					int i=Math.floorDiv(k, n), 	//i=0...m-1
//						j=k%n;					//j=0...n-1
//					for(int p=0;p<n;p++) {
//						grad[k]+=xn[i*n+p];
//					}
//					grad[k]-=inComponents.get(i).getNumber("cov");
//					
//					for(int q=0;q<m;q++) {
//						grad[k]+=xn[j+q*n];
//					}
//					grad[k]-=outComponents.get(j).getNumber("cov");
//					
//					grad[k]*=2;
//				}
//				
//				//calculate the step
//				RealMatrix 	Grad = MatrixUtils.createRealMatrix(m*n,1),
//							step = MatrixUtils.createRealMatrix(m*n,1);
//				Grad.setColumn(0, grad);
//				
//				step=new LUDecomposition(H).getSolver().getInverse().multiply(Grad);
//				
//				isContinue=false;
//				for(int index=0;index<m*n;index++) {
//					xn[index]-=step.getEntry(index,0);
//					if(Math.abs(step.getEntry(index,0))>epsilon)
//						isContinue=true;
//				}
//				
//				System.out.println("Node " + node.getAttribute("name") + " comps estimate after " + count + " iterations!");
//				for(int i=0;i<m;i++) {
//					Edge e = inComponents.get(i);
//					System.out.print("...edge " + e.getId() + " cov=" + e.getNumber("cov") + " components: ");
//					for(int j=0;j<n;j++)
//						System.out.print(xn[i*n+j] + "; ");
//					System.out.println();
//				}
//				for(int j=0;j<n;j++) {
//					Edge e = outComponents.get(j);
//					System.out.print("...edge " + e.getId() + " cov=" + e.getNumber("cov") + " components: ");
//					for(int i=0;i<m;i++)
//						System.out.print(xn[i*n+j] + "; ");
//					System.out.println();
//				}
//				
//			}
//
		}
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
	Ranger union(Ranger r) {
		Ranger retval=null;
		if(left < r.right && right > r.left) 
			return new Ranger(Math.min(left, r.left), Math.max(right, r.right));

		return retval;
	}
	Ranger intersect(Ranger r) {
		Ranger retval=null;
		if(left < r.right && right > r.left) 
			return new Ranger(Math.max(left, r.left), Math.min(right, r.right));

		return retval;
	}
	public void setCooridnates(double left, double right) {
		this.left=left;
		this.right=right;
	}
	public double getLength() {
		return right-left;
	}
	public String toString() {
		return left+"->"+right;
	}
}
//This is actually a tree
class CoverageDistribution{
	double sum;
	ArrayList<CoverageDistribution> components;
	
	public CoverageDistribution(double sum) {
		this.sum=sum;
		components=null;
	}
	
	public boolean isLeaf() {
		return components==null;
	}
	public void setComponents (ArrayList<CoverageDistribution> components) {
		this.components=components;
	}
	public Set<Double> getAllComponents(){
		Set<Double> leafNodes = new HashSet<Double>();
		if(this.isLeaf())
			leafNodes.add(sum);
		else
			for(CoverageDistribution comp:components)
				leafNodes.addAll(comp.getAllComponents());
		
		return leafNodes;
	}
	public String toString() {
		String retval=sum+"=";
		for(Double cov:getAllComponents())
			retval+=cov+",";
		return retval;
	}

}

class PopGroup{
	int id;
	Ranger popRange; 
	ArrayList<Edge> edgesList;
	PopGroup(){}
	
	void addEdge(Edge e) {
		//...
	}
	
	void removeEdge(Edge e) {
		//...
	}
}
