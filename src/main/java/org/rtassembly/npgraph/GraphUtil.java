package org.rtassembly.npgraph;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;
import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.graphstream.graph.implementations.AbstractNode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.joptimizer.functions.ConvexMultivariateRealFunction;
import com.joptimizer.functions.LinearMultivariateRealFunction;
import com.joptimizer.functions.PDQuadraticMultivariateRealFunction;
import com.joptimizer.optimizers.JOptimizer;
import com.joptimizer.optimizers.NewtonUnconstrained;
import com.joptimizer.optimizers.OptimizationRequest;

import htsjdk.samtools.util.CigarUtil;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;

public class GraphUtil {

	private static final Logger LOG = LoggerFactory.getLogger(GraphUtil.class);
	/**********************************************************************************
	 * ****************************Algorithms go from here*****************************
	 */
    //TODO: read from ABySS assembly graph (graph of final contigs, not like SPAdes)
	public static volatile double DISTANCE_THRES=1.0;
    
    public static void loadFromFASTG(String graphFileName, BidirectedGraph graph, boolean spadesBridging) throws IOException{
        graph.setAutoCreate(true);
        graph.setStrict(false);
		/*
		 * 1. next iterate over again to read the connections
		 */
		SequenceReader reader = new FastaReader(graphFileName);
		Sequence seq;
		int shortestLen = 10000;
		ArrayList<EdgeComponents> potentialEdgeSet = new ArrayList<EdgeComponents>();
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			if(seq.length()<shortestLen)
				shortestLen=seq.length();
			
			String[] adjList = seq.getName().split(":");
			String name = adjList[0];
			
			
			boolean dir0=name.contains("'")?false:true;
			
			name=name.replaceAll("[^a-zA-Z0-9_.]", "").trim(); //EDGE_X_length_Y_cov_Z
			
			String nodeID = name.split("_")[1];
			AbstractNode node = (AbstractNode) graph.addNode(nodeID); //or get the existing node prototype created below (to set the attributes)
			node.setAttribute("name", name);
			
			if(dir0){
				seq.setName(name);
				node.setAttribute("seq", seq);
				node.setAttribute("len", seq.length());

				double cov = Double.parseDouble(name.split("_")[5]);

				//Convert from kmer coverage (read kmer per contig kmer)
				//to read coverage (read bp per contig bp)
				//Cx=Ck*L/(L-k+1)
//				cov=cov*BidirectedGraph.ILLUMINA_READ_LENGTH/(BidirectedGraph.ILLUMINA_READ_LENGTH-BidirectedGraph.KMER+1);
				
				node.setAttribute("cov", cov);

//				totReadsLen += cov*seq.length();
//				totContigsLen += seq.length();
			}
			if (adjList.length > 1){
				String[] nbList = adjList[1].split(",");
				for(int i=0; i < nbList.length; i++){
					String neighbor = nbList[i];
					// note that the direction is read reversely in the dest node
					boolean dir1=neighbor.contains("'")?true:false;
					neighbor=neighbor.replaceAll("[^a-zA-Z0-9_.]", "").trim();
					
					String neighborID = neighbor.split("_")[1];
					AbstractNode nbr = (AbstractNode) graph.addNode(neighborID); //just need a prototype, attributes can be set later...

					potentialEdgeSet.add(new EdgeComponents(node,nbr,dir0,dir1)); //edges' prototype
					
				}
			}
			
		}

		reader.close();
		
		for(EdgeComponents ec:potentialEdgeSet) {
			BidirectedEdge e = graph.addEdge(ec.n1, ec.n2, ec.dir1, ec.dir2);
//			System.out.println("...adding edge " + e.getId() + " -> total no of edges = " + graph.getEdgeCount());
//			graph.updateGraphMap(e, new BidirectedPath(e));
		}
		//rough estimation of kmer used
		if((shortestLen-1) != BidirectedGraph.getKmerSize()){
			BidirectedGraph.setKmerSize(shortestLen-1);
//			for(Edge e:graph.getEdgeSet()){
//				((BidirectedEdge)e).changeKmerSize(BidirectedGraph.KMER);
//			}
			graph.edges().forEach(e -> ((BidirectedEdge)e).changeKmerSize(BidirectedGraph.KMER));
		}
		
		double totReadsLen=0, totContigsLen=0;

		for(Node node:graph) {
			if(node.getDegree()==0 && node.getNumber("len")<2*BidirectedGraph.ILLUMINA_READ_LENGTH) {
//				removeNode(node);
			}
			else {
				totReadsLen += node.getNumber("cov")*node.getNumber("len");
				totContigsLen += node.getNumber("len");
			}
		}
		BidirectedGraph.RCOV = totReadsLen/totContigsLen;
		
		
		/**
		 * A-statistics: A(delta,r)=log(e)*delta*n/G -r*log(2)
		 * delta: contig length, r: number of reads comprise this contig
		 * n: total number of reads, G: genome size
		 * Recalculated by Cx, contig_len, read_len, RCOV (average read coverage over the genome)
		 */
		for (Node n:graph) {
			Sequence nseq = (Sequence) n.getAttribute("seq");
			double astats = nseq.length()*BidirectedGraph.RCOV/BidirectedGraph.ILLUMINA_READ_LENGTH
							-Math.log(2)*n.getNumber("cov")*nseq.length()/BidirectedGraph.ILLUMINA_READ_LENGTH;
			astats*=Math.log10(Math.E);
			n.setAttribute("astats", astats);
			LOG.info("{} Read coverage = {} Length = {} A-stats = {}", n.getAttribute("name"), n.getNumber("cov"), n.getNumber("len"), astats );
		}
		LOG.info("No of nodes= {} No of edges = {} Estimated avg. read coverage = {} Total contigs length = {}", graph.getNodeCount(), graph.getEdgeCount(), BidirectedGraph.RCOV, totContigsLen );
		
		
		/*
		 * 2. Use a binner to estimate graph multiplicity
		 */
		graph.binning();
		/*
		 * 3. Now scan for the contigs.path file in SPAdes folder for the paths if specified
		 */
		if(spadesBridging){
			File graphFile = new File(graphFileName);
			String pathsFile = FilenameUtils.getFullPathNoEndSeparator(graphFile.getAbsolutePath()) + "/contigs.paths";
			if(! new File(pathsFile).exists()){
				LOG.warn("Path file {} couldn't be found in SPAdes output! Skipped!", pathsFile);
				return;
			}
			BufferedReader pathReader = new BufferedReader(new FileReader(pathsFile));
			
			String s="", curpath="";
			//Read contigs from contigs.paths
			boolean flag=false, changed=false;
			while((s=pathReader.readLine()) != null){
				if(s.contains("NODE")){
					if(flag){
						BidirectedPath path=new BidirectedPath(graph, curpath);
				    	if(graph.reduce(path))
				    		changed=true;
					}
					flag=s.contains("'")?false:true;
					curpath=new String();
					continue;
				}else if(flag){
					curpath+=s;
				}	
					

			}
			pathReader.close();
			if(changed)
				GraphExplore.redrawGraphComponents(graph);
		}
    }
    
    
    
    public static void loadFromGFA(String graphFile, BidirectedGraph graph, boolean spadesBridging) throws IOException{
        graph.setAutoCreate(true);
        graph.setStrict(false);
		/*
		 * 1. next iterate over again to read the connections
		 */
        BufferedReader reader = new BufferedReader(new FileReader(graphFile));
        String line=null;
		Sequence seq;
		int shortestLen = 10000;
		
		ArrayList<BidirectedPath> spadesPaths = new ArrayList<>();
		while ((line=reader.readLine()) != null){
			String[] gfaFields = line.split("\\s");
			String type = gfaFields[0];
			switch (type.toUpperCase().trim()) {
			case "#"://comments
				break;
			case "H"://header
				break;
			case "C"://containment
				break;
			case "S"://segment
				String 	nodeID=gfaFields[1].trim(),
						optField=gfaFields[3].trim();
				String[] toks=optField.split(":");
				assert toks[0].equals("KC")&&toks[1].equals("i"):"Invalid k-mer count field!";
				
				seq = new Sequence(Alphabet.DNA5(), gfaFields[2], nodeID);
				AbstractNode node = (AbstractNode) graph.addNode(nodeID); //or get the existing node prototype created below (to set the attributes)
				node.setAttribute("name", "Node "+nodeID);
				node.setAttribute("seq", seq);
				node.setAttribute("len", seq.length());
				//node.setAttribute("kc", Integer.parseInt(toks[2]));
				//TODO: replace read cov with kc calculation?
				//1.kmer count to kmer cov
				double cov=Double.parseDouble(toks[2])/(seq.length()-BidirectedGraph.getKmerSize()); 
				//2.kmer cov to read cov
				cov*=BidirectedGraph.ILLUMINA_READ_LENGTH/(BidirectedGraph.ILLUMINA_READ_LENGTH-BidirectedGraph.getKmerSize());
				node.setAttribute("cov", cov);		
				
				break;
			case "L"://links
				BidirectedNode 	n0=(BidirectedNode) graph.getNode(gfaFields[1]),
								n1=(BidirectedNode) graph.getNode(gfaFields[3]);
				boolean dir0=gfaFields[2].equals("+")?true:false,
						dir1=gfaFields[4].equals("+")?false:true;
				BidirectedEdge e = graph.addEdge(n0, n1, dir0, dir1);
				
				//just do it simple for now when the last field of Links line is xxM (kmer=xx)
				String cigar=gfaFields[5];
				for(int i=0; i < cigar.length();i++) {
					char c = cigar.charAt(i);
					if(c >= '0' && c <= '9')
						continue;
					else {
						int len = Integer.parseInt(cigar.substring(0, i));
						if(shortestLen > len)
							shortestLen=len;
						break;
					}
				}
									
				break;
			case "P"://path
				if(spadesBridging){
					if(gfaFields.length>3 && gfaFields[2].contains(",")) {
						spadesPaths.add(new BidirectedPath(graph, gfaFields[2]));
					}
				}
				break;

				default:throw new IllegalStateException("Unrecognized GFA field: " + type.toUpperCase().trim());
			}
			
		
		}

		reader.close();

		//rough estimation of kmer used
		if((shortestLen-1) != BidirectedGraph.getKmerSize()){
			BidirectedGraph.setKmerSize(shortestLen-1);
//			for(Edge e:graph.getEdgeSet()){
//				((BidirectedEdge)e).changeKmerSize(BidirectedGraph.KMER);
//			}
			graph.edges().forEach(e -> ((BidirectedEdge)e).changeKmerSize(BidirectedGraph.KMER));

		}
		
		double totReadsLen=0, totContigsLen=0;

		for(Node node:graph) {
			if(node.getDegree()==0 && node.getNumber("len")<2*BidirectedGraph.ILLUMINA_READ_LENGTH) {
//				removeNode(node);
			}
			else {
				totReadsLen += node.getNumber("cov")*node.getNumber("len");
				totContigsLen += node.getNumber("len");
			}
		}
		BidirectedGraph.RCOV = totReadsLen/totContigsLen;
		
		
		/**
		 * A-statistics: A(delta,r)=log(e)*delta*n/G -r*log(2)
		 * delta: contig length, r: number of reads comprise this contig
		 * n: total number of reads, G: genome size
		 * Recalculated by Cx, contig_len, read_len, RCOV (average read coverage over the genome)
		 */
		for (Node n:graph) {
			Sequence nseq = (Sequence) n.getAttribute("seq");
			double astats = nseq.length()*BidirectedGraph.RCOV/BidirectedGraph.ILLUMINA_READ_LENGTH
							-Math.log(2)*n.getNumber("cov")*nseq.length()/BidirectedGraph.ILLUMINA_READ_LENGTH;
			astats*=Math.log10(Math.E);
			n.setAttribute("astats", astats);
			LOG.info("{} Read coverage = {} Length = {} A-stats = {}", n.getAttribute("name"), n.getNumber("cov"), n.getNumber("len"), astats );
		}
		LOG.info("No of nodes= {} No of edges = {} Estimated avg. read coverage = {} Total contigs length = {}", graph.getNodeCount(), graph.getEdgeCount(), BidirectedGraph.RCOV, totContigsLen );
		
		/*
		 * 2. Binning the graph
		 */
		graph.binning();
		
		/*
		 * 3. Reduce the SPAdes path if specified
		 */
		boolean changed=false;
		if(spadesBridging){
			for(BidirectedPath p:spadesPaths)
				if(graph.reduce(p))
					changed=true;
			if(changed)
				GraphExplore.redrawGraphComponents(graph);
		}
    }
    
	public static void gradientDescent(BidirectedGraph graph) {
		int 	maxIterations=100, 
				eIteCount=0, nIteCount=0;
		double epsilon=.01;
		while(true) {
			nIteCount++;
			eIteCount=0;
			//1. Updating edges' coverage			
			while(true) {
				eIteCount++;
				HashMap<String,Double> stepMap = new HashMap<>();
				Iterator<Edge> edgesIterator=graph.edges().iterator();
				while(edgesIterator.hasNext()) {
					Edge e = edgesIterator.next();
//					System.out.println("Working on edge "+e.getId());
		    		BidirectedNode n0 = (BidirectedNode) e.getNode0(), n1=(BidirectedNode) e.getNode1();
		    		boolean dir0 = ((BidirectedEdge) e).getDir0(), dir1 = ((BidirectedEdge) e).getDir1();
		    		Iterator<Edge> 	ite0 = dir0?n0.leavingEdges().iterator():n0.enteringEdges().iterator(),
		    						ite1 = dir1?n1.leavingEdges().iterator():n1.enteringEdges().iterator();
		    		double sum0=0, sum1=0, tmp;
		    		int deg0=0, deg1=0;
		    		while(ite0.hasNext()) {
		    			Edge e0=ite0.next();
		    			tmp=e0.getNumber("cov");
		    			sum0+=Double.isNaN(tmp)?1.0:tmp;
		    			deg0++;
//		    			System.out.printf("\tedge %s cov=%.2f sum0=%.2f\n",e0.getId(),tmp,sum0);
		    		}
		    		while(ite1.hasNext()) {
		    			Edge e1=ite1.next();
		    			tmp=e1.getNumber("cov");		    			
		    			sum1+=Double.isNaN(tmp)?1.0:tmp;
		    			deg1++;
//		    			System.out.printf("\tedge %s cov=%.2f sum1=%.2f\n",e1.getId(),tmp,sum1);

		    		}	    		
		    		//gamma_ij=1/*(len_i+len_j) -> failed!
		    		//gamma_ij=1/2*(len_i+len_j) -> small enough! (explanation???)
		    		double value=.5*(n0.getNumber("len")*(sum0-n0.getNumber("cov"))/(deg0*deg0) + n1.getNumber("len")*(sum1-n1.getNumber("cov"))/(deg1*deg1))/(n0.getNumber("len")+n1.getNumber("len"));
//		    		System.out.println("=> step="+value);
		    		stepMap.put(e.getId(), value);
				}
				boolean isConverged=true, isZero=false;
				
				edgesIterator=graph.edges().iterator();
				while(edgesIterator.hasNext()) {
					Edge e = edgesIterator.next();
					double delta=stepMap.get(e.getId()),
							curCov=Double.isNaN(e.getNumber("cov"))?1.0:e.getNumber("cov");
					if(Math.abs(delta/curCov) > epsilon) {
						isConverged=false;
					}
					
					if(curCov<=delta) {
						LOG.warn("Edge " + e.getId() + " coverage is not positive : curCov=" + curCov + ", delta=" + delta);
//						isZero=true;
					}else
						e.setAttribute("cov", curCov-delta);
				}
//				if(isZero)
//					return;
				if(isConverged || eIteCount >= maxIterations) {
//					System.out.println("...edges coverage CONVERGED at iteration " + eIteCount + "th");
					break;
				}
			}
			//2. Updating nodes' coverage: keep significant node info intact!
			boolean isConverged=true;
			for(Node n:graph) {
//				if(n.getAttribute("marked") != null)//do smt???
//					continue;
				
				Iterator<Edge> 	in=n.enteringEdges().iterator(),
								out=n.leavingEdges().iterator();
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
				if(Math.abs(newCovEst-n.getNumber("cov"))/n.getNumber("cov") > epsilon)
					isConverged=false;
				n.setAttribute("cov", newCovEst);
			}
			if(isConverged || nIteCount >= maxIterations) {
//				System.out.println("Node coverage CONVERGED at iteration " + nIteCount + "th");
//				System.out.println("======================================================");
				break;
			}
			
		}
	}
	
	public static void coverageOptimizer(BidirectedGraph graph) {
		org.apache.log4j.BasicConfigurator.configure();
		int nIteCount=0;
		while(true) {
			nIteCount++;
			//1. Updating edges' coverage
			ArrayList<Edge> edges = new ArrayList<Edge>(graph.edges().collect(Collectors.toList()));
			int edgesNumber=edges.size();
			HashMap<Edge,Integer> idToIndex = new HashMap<Edge,Integer>();
			double[] init = new double[edgesNumber];
			for(int i=0;i<edgesNumber;i++) {
				idToIndex.put(edges.get(i), i);
				init[i]=Double.isNaN(edges.get(i).getNumber("cov"))?0:edges.get(i).getNumber("cov");
			}
			
			
			/**	A minimization problem in the form of:
			***	  minimize_{x} (1/2)x^{T}Px+q^{T}x+r  s.t. 
			***	    Gx ≤ h 
			***	    Ax = b,  
			***	where P ∈ S+^{n} (symmetric nxn), G ∈ R^{mxn} and A ∈ R^{pxn}
			***
			***	is called a quadratic program (QP). In a quadratic program we minimize a (convex) quadratic objective function with affine constraint functions. Quadratic programs include linear program as a special case by taking P=0.
			***/
			//Objective function
			double[][] PMatrix = new double[edgesNumber][edgesNumber];
			double[] qVector = new double[edgesNumber];
			double r = 0;
			
			//equalities
			double[][] A = new double[graph.getNodeCount()][edgesNumber];
			double[]  b = new double[graph.getNodeCount()];
			int p=0;
			
		
			for(Node node:graph) {
				if(node.getDegree()==0)
					continue;
				double cov=node.getNumber("cov");
				int len=(int)node.getNumber("len");
				Iterator<Edge> 	in=node.enteringEdges().iterator(),
								out=node.leavingEdges().iterator();
				ArrayList<Integer> 	inIndices = new ArrayList<Integer>(), 
									outIndices = new ArrayList<Integer>();
				while(in.hasNext()) {
					inIndices.add(idToIndex.get(in.next()));
				}
				while(out.hasNext()) {
					outIndices.add(idToIndex.get(out.next()));
				}
				for(int i:inIndices) {
					qVector[i]-=cov*len;
					A[p][i]=1; 
					for(int j:inIndices)
						PMatrix[i][j]+=len;
				}
				for(int i:outIndices) {
					qVector[i]-=cov*len;
					A[p][i]=-1;
					for(int j:outIndices)
						PMatrix[i][j]+=len;
				}
				b[p]=0;
				p++;
				r+=len*cov*cov;
			}
			
			PDQuadraticMultivariateRealFunction objectiveFunction = new PDQuadraticMultivariateRealFunction(PMatrix, qVector, r);
			
//			ConvexMultivariateRealFunction[] inequalities = new ConvexMultivariateRealFunction[p];
//			for(int i=0;i<p;i++) {
//				double[] tmp=new double[edgesNumber];
//				tmp[i]=-1;
//				inequalities[i] = new LinearMultivariateRealFunction(tmp, 0);
//			}
			
			OptimizationRequest or = new OptimizationRequest();
			or.setF0(objectiveFunction);
			or.setInitialPoint(init);
//			or.setFi(inequalities);
//			or.setA(A);
//			or.setB(b);
//			or.setToleranceFeas(1.E-8);
			or.setTolerance(1.E-8);
			
//			JOptimizer opt = new JOptimizer();
			NewtonUnconstrained opt = new NewtonUnconstrained();
			opt.setOptimizationRequest(or);
			try {
				opt.optimize();
				double[] sol = opt.getOptimizationResponse().getSolution();
				for(int i=0;i<edgesNumber;i++) {
					edges.get(i).setAttribute("cov", sol[i]);
				}
				
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			
			//2. Updating nodes' coverage
			boolean isConverged=true;
			for(Node n:graph) {
				Iterator<Edge> 	in=n.enteringEdges().iterator(),
								out=n.leavingEdges().iterator();
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
				if(Math.abs(newCovEst/n.getNumber("cov")) > 1.E-2)
					isConverged=false;
				n.setAttribute("cov", newCovEst);
			}
			if(isConverged || nIteCount >= 10) {
				System.out.println("STOP at iteration " + nIteCount + "th");
				System.out.println("======================================================");
				break;
			}
			
		}
	}
	
	//The distance between 2 Poisson distribution based on Kullback-Leibner divergence
	//distance=0.5*(KL(a,b) + KL(b,a))
	//TODO: take density, not just mean of distribution into consideration
	public static double metric(double a, double b) {
		return (.5*(a-b)*(Math.log(a) -Math.log(b)));
	}
    
	
    /*
     * Compare 2 double values x,y
     * Return 0 if x~=y, 1 if x>>y, -1 if x<<y
     */
    public static int approxCompare(double x, double y) {
    	int retval=0;
    	double ratio=Math.abs(x-y)/(Math.max(Math.abs(x), Math.abs(y)));
    	if(ratio > .33)
    		retval=x>y?1:-1;
    	
    	return retval;
    }
    
    /*
     * Get string representative of a path starting node. 
     * Note that there are 2 of them (template start + reversed-complementary end)
     * E.g: 103-82- would return 103- and 82+
     */
    public static String[] getHeadOfPathString(String brg){
    	String[] retval = null;
    	Pattern pattern = Pattern.compile("(\\d+)([\\+\\-]),(\\d+)([\\+\\-])");
		Matcher matcher =pattern.matcher(brg);
		if (matcher.find()){
			retval = new String[2];
			retval[0]=matcher.group(1)+matcher.group(2);
			retval[1]=matcher.group(3)+(matcher.group(4).trim().equals("+")?"-":"+");
		}
		
		return retval;
    }
    
}
class EdgeComponents{
	AbstractNode n1,n2;
	boolean dir1,dir2;

	EdgeComponents(AbstractNode n1, AbstractNode n2, boolean dir1, boolean dir2){
		this.n1=n1;
		this.n2=n2;
		this.dir1=dir1;
		this.dir2=dir2;
	}
}
