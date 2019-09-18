package org.rtassembly.npgraph;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
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

import com.joptimizer.functions.PDQuadraticMultivariateRealFunction;
import com.joptimizer.optimizers.NewtonUnconstrained;
import com.joptimizer.optimizers.OptimizationRequest;

import htsjdk.samtools.SAMRecord;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceReader;

public class GraphUtil {

	private static final Logger LOG = LoggerFactory.getLogger(GraphUtil.class);
	/**********************************************************************************
	 * ****************************Algorithms go from here*****************************
	 */
    //TODO: read from ABySS assembly graph (graph of final contigs, not like SPAdes)
	public static volatile double DISTANCE_THRES=3.0; //i like number 3
	public static HashMap<Node,Double> originalCoverageValues = new HashMap<>();
    
    public static void loadFromFASTG(String graphFileName, String binFileName, BDGraph graph, boolean spadesBridging) throws IOException{
        graph.setAutoCreate(true);
        graph.setStrict(false);
		/*
		 * 1. next iterate over again to read the connections
		 */
		SequenceReader reader = new FastaReader(graphFileName);
		Sequence seq;
		int shortestLen = 10000;
		ArrayList<BDEdgePrototype> potentialEdgeSet = new ArrayList<BDEdgePrototype>();
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			if(seq.length()<shortestLen)
				shortestLen=seq.length();
			
			String[] adjList = seq.getName().split(":");
			String name = adjList[0];
			
			
			boolean dir0=name.contains("'")?false:true;
			
			name=name.replaceAll("[^a-zA-Z0-9_.]", "").trim(); //EDGE_X_length_Y_cov_Z
			
			String nodeID = name.split("_")[1];
			BDNode node = (BDNode) graph.addNode(nodeID); //or get the existing node prototype created below (to set the attributes)
			node.setAttribute("name", name);
			
			if(dir0){
				seq.setName(name);
				node.setAttribute("seq", seq);
				node.setAttribute("len", seq.length());

				double cov = Double.parseDouble(name.split("_")[5]);
				
				node.setAttribute("cov", cov);

			}
			if (adjList.length > 1){
				String[] nbList = adjList[1].split(",");
				for(int i=0; i < nbList.length; i++){
					String neighbor = nbList[i];
					// note that the direction is read reversely in the dest node
					boolean dir1=neighbor.contains("'")?true:false;
					neighbor=neighbor.replaceAll("[^a-zA-Z0-9_.]", "").trim();
					
					String neighborID = neighbor.split("_")[1];
					BDNode nbr = (BDNode) graph.addNode(neighborID); //just need a prototype, attributes can be set later...

					potentialEdgeSet.add(new BDEdgePrototype(node,nbr,dir0,dir1)); //edges' prototype
					
				}
			}
			
		}

		reader.close();
		
		for(BDEdgePrototype ec:potentialEdgeSet) {
			try {
				graph.addEdge(ec.getNode0(), ec.getNode1(), ec.getDir0(), ec.getDir1());
			} catch (Exception e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				System.exit(1);
			}
		}
		//rough estimation of kmer used
		if((shortestLen-1) < BDGraph.getKmerSize()){
			BDGraph.setKmerSize(shortestLen-1);
		}
		
		double totReadsLen=0, totContigsLen=0;

		for(Node node:graph) {
			//Convert from kmer coverage (read kmer per contig kmer)
			//to read coverage (read bp per contig bp)
			//Cx=Ck*L/(L-k+1)
//			cov=cov*BDGraph.ILLUMINA_READ_LENGTH/(BDGraph.ILLUMINA_READ_LENGTH-BDGraph.KMER+1);
			double cov=node.getNumber("cov")*BDGraph.ILLUMINA_READ_LENGTH/(BDGraph.ILLUMINA_READ_LENGTH-BDGraph.getKmerSize()+1);
			node.setAttribute("cov", cov);	
			
			//ignore noises
			if(node.getDegree()==0 && node.getNumber("len")<2*BDGraph.ILLUMINA_READ_LENGTH) {
				continue;
			}
			else {
				totReadsLen += node.getNumber("cov")*node.getNumber("len");
				totContigsLen += node.getNumber("len");
			}
		}
		BDGraph.RCOV = totReadsLen/totContigsLen;
		
		
		/**
		 * A-statistics: A(delta,r)=log(e)*delta*n/G -r*log(2)
		 * delta: contig length, r: number of reads comprise this contig
		 * n: total number of reads, G: genome size
		 * Recalculated by Cx, contig_len, read_len, RCOV (average read coverage over the genome)
		 */
		
		//Try print sequence of Astat
		//A(delta,r,k)=log(e)*delta*n/G +r*log(n/n+1)
		//Celera Astats = A(delta,r,1)
		
		for (Node node:graph) {
			Sequence nseq = (Sequence) node.getAttribute("seq");
			double astats=-1;

//			int estcov=(int) Math.round(node.getNumber("cov")/BDGraph.RCOV);
//			//get the first positive astats, if > 10 then assign its multiplicity (not work with different pops) 
//			int i=0;
//			while(astats<0){
//				i++;
//				astats = nseq.length()*BDGraph.RCOV/BDGraph.ILLUMINA_READ_LENGTH
//				+Math.log((i+1)*1.0/(i+2))*node.getNumber("cov")*nseq.length()/BDGraph.ILLUMINA_READ_LENGTH;
//				astats*=Math.log10(Math.E);		
//				System.out.printf("...A-stat at iteration %d: %.2f\n", i, astats);
//			}
//			normalizedCoverage(node);			
//			if(astats>10)
//				LOG.info("{} Normalized coverage={} Length={} Coverage={} A-stats={}", node.getAttribute("name"), node.getNumber("cov"), node.getNumber("len"), estcov, astats );
//			else
//				LOG.info("{} Normalized coverage={} Length={}", node.getAttribute("name"), node.getNumber("cov"), node.getNumber("len") );

			//calculate 1 astat as normal
			astats = nseq.length()*BDGraph.RCOV/BDGraph.ILLUMINA_READ_LENGTH
							-Math.log(2)*node.getNumber("cov")*nseq.length()/BDGraph.ILLUMINA_READ_LENGTH;
			astats*=Math.log10(Math.E);
			node.setAttribute("astats", astats);
			normalizedCoverage(node);
			if(HybridAssembler.VERBOSE)
				LOG.info("{} Normalized coverage = {} Length = {} \nA-stats = {}", node.getAttribute("name"), node.getNumber("cov"), node.getNumber("len"), astats);
		}
		
		LOG.info("No of nodes= {} No of edges = {} Estimated avg. read coverage = {} (normalized to 100.0) Total contigs length = {}", graph.getNodeCount(), graph.getEdgeCount(), BDGraph.RCOV, totContigsLen );
		
		
		/*
		 * 2. Use a binner to estimate graph multiplicity
		 */
//		graph.nodes().filter(n->n.getNumber("cov") < .2*BDGraph.RCOV).forEach(n->{n.edges().forEach(e->graph.removeEdge(e));});
		
		graph.fixDeadEnds();
		graph.binning(binFileName);

		/*
		 * 3. Now scan for the contigs.path file in SPAdes folder for the paths if specified
		 */
		if(spadesBridging){
			File graphFile = new File(graphFileName);
			String pathsFile = FilenameUtils.getFullPathNoEndSeparator(graphFile.getAbsolutePath()) + File.separator + "contigs.paths";
			if(! new File(pathsFile).exists()){
				LOG.warn("Path file {} couldn't be found in SPAdes output! Skipped!", pathsFile);
				return;
			}
			BufferedReader pathReader = new BufferedReader(new FileReader(pathsFile));
			
			String s="", curpath="";
			//Read contigs from contigs.paths
			boolean flag=false;
			while((s=pathReader.readLine()) != null){
				if(s.contains("NODE")){
					if(flag){
						String[] consecutivePaths = curpath.split(";");
						BDPath path=null;
						for(int i=0;i<consecutivePaths.length;i++){
							//TODO: make use of the information about gapped paths (separated by ";")
							path=new BDPath(graph, consecutivePaths[i]); //only after binning
					    	graph.reduceFromSPAdesPath(path);
						}
					}
					flag=s.contains("'")?false:true;
					curpath=new String();
					continue;
				}else if(flag){
					curpath+=s;
				}	
					

			}
			pathReader.close();
		}
    }
    
    
    
    public static void loadFromGFA(String graphFile, String binFileName, BDGraph graph, boolean spadesBridging) throws IOException{
        graph.setAutoCreate(true);
        graph.setStrict(false);
		/*
		 * 1. next iterate over again to read the connections
		 */
        BufferedReader reader = new BufferedReader(new FileReader(graphFile));
        String line=null;
		Sequence seq;
		int overlapLen = 10000;
		
		ArrayList<String> spadesPaths = new ArrayList<>();
		boolean covFlag=true; //true if GFA from SPAdes contains KC:i:xx; false if GFA from Unicycler dp:f:xx
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
					String 	nodeID=gfaFields[1].trim();
					assert gfaFields.length>3 && gfaFields[2].trim().length()>1:"Invalid GFA v1 file!";
					seq = new Sequence(Alphabet.DNA5(), gfaFields[2], nodeID);
					AbstractNode node = (AbstractNode) graph.addNode(nodeID); //or get the existing node prototype created below (to set the attributes)
					node.setAttribute("name", "Contig_"+nodeID);
					seq.setName("Contig_"+nodeID);
					node.setAttribute("seq", seq);
					node.setAttribute("len", seq.length());
					
					for(int i=3; i<gfaFields.length; i++){
						String	optField=gfaFields[i].trim();
						String[] toks=optField.split(":");
						if(toks[0].equals("KC") && toks[1].equals("i")){
							//here is the kmer coverage
							node.setAttribute("cov", Integer.parseInt(toks[2]));
						}else if(toks[0].equals("df") && toks[1].equals("f")){
							covFlag=false;
							node.setAttribute("cov", Double.parseDouble((toks[2])));
						}
					}
					
					break;
				case "L"://links
					BDNode 	n0=(BDNode) graph.getNode(gfaFields[1]),
							n1=(BDNode) graph.getNode(gfaFields[3]);
					boolean dir0=gfaFields[2].equals("+")?true:false,
							dir1=gfaFields[4].equals("+")?false:true;
					graph.addEdge(n0, n1, dir0, dir1);
					
					//just do it simple for now when the last field of Links line is xxM (kmer=xx)
					String cigar=gfaFields[5];
					for(int i=0; i < cigar.length();i++) {
						char c = cigar.charAt(i);
						if(c >= '0' && c <= '9')
							continue;
						else {
							int len = Integer.parseInt(cigar.substring(0, i));
							if(overlapLen > len)
								overlapLen=len;
							break;
						}
					}
										
					break;
				case "P"://path
					if(spadesBridging){
						if(gfaFields.length>3 && gfaFields[2].contains(",")) {
							spadesPaths.add(gfaFields[2]);
						}
					}
					break;
	
				default:
					LOG.warn("Unrecognized GFA field: " + type.toUpperCase().trim());
			}
			
		
		}

		reader.close();

		//rough estimation of kmer used
		if((overlapLen) < BDGraph.getKmerSize()){
			BDGraph.setKmerSize(overlapLen);

		}

		
		double totReadsLen=0, totContigsLen=0;

		for(Node node:graph) {
			
			if(covFlag){//only do this for SPAdes' GFA to convert kmer count to read cov
				//1.kmer count to kmer cov
				double cov=node.getNumber("cov")/(node.getNumber("len")-BDGraph.getKmerSize()); 
				//2.kmer cov to read cov
				cov*=BDGraph.ILLUMINA_READ_LENGTH/(BDGraph.ILLUMINA_READ_LENGTH-BDGraph.getKmerSize()+1);
				node.setAttribute("cov", cov);	
			}
			
			//ignore noises
			if(node.getDegree()==0 && node.getNumber("len")<2*BDGraph.ILLUMINA_READ_LENGTH) {
				continue;
			}
			else {
				totReadsLen += node.getNumber("cov")*node.getNumber("len");
				totContigsLen += node.getNumber("len");
			}
		}
		BDGraph.RCOV = totReadsLen/totContigsLen;
		
		
		/**
		 * A-statistics: A(delta,r)=log(e)*delta*n/G -r*log(2)
		 * delta: contig length, r: number of reads comprise this contig
		 * n: total number of reads, G: genome size
		 * Recalculated by Cx, contig_len, read_len, RCOV (average read coverage over the genome)
		 */
		for (Node node:graph) {
			Sequence nseq = (Sequence) node.getAttribute("seq");
			double astats = nseq.length()*BDGraph.RCOV/BDGraph.ILLUMINA_READ_LENGTH
							-Math.log(2)*node.getNumber("cov")*nseq.length()/BDGraph.ILLUMINA_READ_LENGTH;
			astats*=Math.log10(Math.E);
			node.setAttribute("astats", astats);
			normalizedCoverage(node);
			if(HybridAssembler.VERBOSE)		
				LOG.info("{} Normalized coverage = {} Length = {} A-stats = {}", node.getAttribute("name"), node.getNumber("cov"), node.getNumber("len"), astats );
		}
		
		LOG.info("No of nodes= {} No of edges = {} Estimated avg. read coverage = {} (normalized to 100.0) Total contigs length = {}", graph.getNodeCount(), graph.getEdgeCount(), BDGraph.RCOV, totContigsLen );
		
		/*
		 * 2. Binning the graph
		 */

		graph.fixDeadEnds();
		graph.binning(binFileName);
		
		/*
		 * 3. Reduce the SPAdes path if specified
		 */
		if(spadesBridging){
			for(String pID:spadesPaths)
				graph.reduceFromSPAdesPath(new BDPath(graph,pID));
		}
    }
    
   
    public static int overlap(Sequence s0, Sequence s1){
    	int retval;
    	for(retval=BDGraph.getKmerSize() ; retval>0; retval--){
    		int match=0;
    		while(match < retval && s0.getBase(s0.length()-retval+match)==s1.getBase(match))
    			match++;
    		if(match==retval)
    			return match;
    	}
    	return 0;
    }
    
    public static String getIDFromName(String readName){
    	String retval=readName; 
		if(readName.contains("_"))
			retval=readName.split("_")[1];	//SPAdes style
   	
    	return retval;
    }
    
    //Normalized the average cov to 100
    public static void normalizedCoverage(Node node){
    	assert BDGraph.RCOV!=0 && node.hasAttribute("cov"):"Cannot normalize coverage of node " + node.getId();
    	node.setAttribute("cov", node.getNumber("cov")*100.0/BDGraph.RCOV);
    	originalCoverageValues.put(node, node.getNumber("cov"));
    }
    
    //Get real value for output
    public static double getRealCoverage(double cov){
    	return cov*BDGraph.RCOV/100.0;
    }
    
    //Check if a node has been likely used up its coverage. Return true iff its current coverage still suggest its residence somewhere.
    public static boolean isLikelyStillPresented(Node node){
    	return node.getNumber("cov") > originalCoverageValues.get(node)/3;
    }
    
    
	public static void gradientDescent(BDGraph graph) {
		int 	maxIterations=21, 
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
		    		BDNode n0 = (BDNode) e.getNode0(), n1=(BDNode) e.getNode1();
		    		boolean dir0 = ((BDEdge) e).getDir0(), dir1 = ((BDEdge) e).getDir1();
		    		Iterator<Edge> 	ite0 = dir0?n0.leavingEdges().iterator():n0.enteringEdges().iterator(),
		    						ite1 = dir1?n1.leavingEdges().iterator():n1.enteringEdges().iterator();
		    		double sum0=0, sum1=0, tmp;
		    		int deg0=0, deg1=0;
		    		while(ite0.hasNext()) {
		    			Edge e0=ite0.next();
		    			tmp=e0.getNumber("cov");
		    			sum0+=Double.isNaN(tmp)?1.0:tmp;
		    			deg0++;
		    		}
		    		while(ite1.hasNext()) {
		    			Edge e1=ite1.next();
		    			tmp=e1.getNumber("cov");		    			
		    			sum1+=Double.isNaN(tmp)?1.0:tmp;
		    			deg1++;
		    		}	    		
		    		//gamma_ij=1/*(len_i+len_j) -> failed!
		    		//gamma_ij=1/2*(len_i+len_j) -> small enough! (explanation???)
		    		double value=.5*(n0.getNumber("len")*(sum0-n0.getNumber("cov"))/(deg0*deg0) + n1.getNumber("len")*(sum1-n1.getNumber("cov"))/(deg1*deg1))/(n0.getNumber("len")+n1.getNumber("len"));
		    		stepMap.put(e.getId(), value);
				}
				boolean isConverged=true;
				
				edgesIterator=graph.edges().iterator();
				while(edgesIterator.hasNext()) {
					Edge e = edgesIterator.next();
					double delta=stepMap.get(e.getId()),
							curCov=Double.isNaN(e.getNumber("cov"))?1.0:e.getNumber("cov");
					if(Math.abs(delta/curCov) > epsilon) {
						isConverged=false;
					}
					
					if(curCov<=delta && HybridAssembler.VERBOSE)
						LOG.warn("Edge " + e.getId() + " coverage is not positive : curCov=" + curCov + ", delta=" + delta);
					else
						e.setAttribute("cov", curCov-delta);
				}

				if(isConverged || eIteCount >= maxIterations)
					break;
			}
			//2. Updating nodes' coverage: keep significant node info intact!
			boolean isConverged=true;
			for(Node n:graph) {
				if(n.getInDegree()==0 && n.getOutDegree()==0)//do smt???
					continue;
				
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
				break;
			}
			
		}
	}
	
	public static void coverageOptimizer(BDGraph graph) {
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
				if(HybridAssembler.VERBOSE) 				
					LOG.info("STOP at iteration " + nIteCount + "th");
				
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
     * Compare 2 double values x,y >=0
     * Return 0 if x~=y, 1 if x>>y, -1 if x<<y
     */
    public static int approxCompare(double x, double y) {
    	int retval=0;
    	if(x==0 && y==0)
    		return 0;
    	double ratio=Math.abs(x-y)/(Math.max(Math.abs(x), Math.abs(y)));
//    	TODO: check this!
//    	double ratio=Math.abs(x-y)/(Math.min(Math.abs(x), Math.abs(y))); <>1.5?
    	
    	if(ratio > BDGraph.R_TOL)
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
			retval[0]=matcher.group(1)+(matcher.group(2).trim().equals("+")?"o":"i");
			retval[1]=matcher.group(3)+(matcher.group(4).trim().equals("+")?"i":"o");
		}
		
		return retval;
    }
   
    

    public static void sleep() {
        try { Thread.sleep(1000); } catch (Exception e) {}
    }

    	
	public static String styleSheet =				// 1
			"node { size: 7px; fill-color: rgb(150,150,150); }" +
			"edge { fill-color: rgb(255,50,50); size: 2px; }" +
			"edge.cut { fill-color: rgba(200,200,200,128); }";
    
	/* 
	 * Adapt from japsa without need to output fasta input file
	 */
	public static Sequence consensusSequence(String prefix, int distance) throws IOException, InterruptedException{
		Sequence consensus = null;	
		String msa = BDGraph.MSA;
		String 	faoFile = AlignedRead.tmpFolder+File.separator+prefix+"_"+msa+".fasta";
		File consFile = new File(faoFile);
		if(!consFile.exists()){
			String faiFile=AlignedRead.tmpFolder+File.separator+prefix+".fasta";
			//1.0 Check faiFile?
			FastaReader faiReader =  new FastaReader(faiFile);
			int count=0, minGap=Integer.MAX_VALUE;
			Sequence seq=null;
			while((seq=faiReader.nextSequence(Alphabet.DNA5()))!=null){
				//if there is no MSA detected, output the sequence with the least non-Illumina bases
				if(msa.isEmpty() || msa.startsWith("none")){
					String[] toks=seq.getDesc().split("=");
					try{
						int curGap=Integer.parseInt(toks[1]);
						if( curGap < minGap){
							minGap=curGap;
							consensus=seq;
						}		
					}catch(NumberFormatException e){
						if(HybridAssembler.VERBOSE)
							LOG.info("Pattern %s=%d not found in the header of input FASTA file {}", faiFile);
						if(consensus==null || consensus.length() > seq.length())
							consensus=seq;
					}
				}
				count++;
			}
			
			faiReader.close();
			if(count<BDGraph.MIN_SUPPORT) //at least 3 of the good read count is required
				return null;
			
			//2.0 Run multiple alignment. 
			// For now, only kalign were chosen due to its speedy operation
			String cmd  = "";
			if (msa.startsWith("kalign")){
				cmd = "kalign -gpo 60 -gpe 10 -tgpe 0 -bonus 0 -q -i " + faiFile	+ " -o " + faoFile;
//			}else if (msa.startsWith("poa")){
//				String poaDir="/home/sonhoanghguyen/sw/poaV2/";//test
//				cmd = poaDir+"poa -read_fasta " + faiFile + " -clustal " + faoFile + " -hb " + poaDir+"blosum80.mat";
//			}else if (msa.startsWith("muscle")){
//				cmd = "muscle -in " + faiFile + " -out " + faoFile + " -maxiters 5 -quiet";				
//			}else if (msa.startsWith("clustal")) {
//				cmd = "clustalo --force -i " + faiFile + " -o " + faoFile;
//			}else if (msa.startsWith("msaprobs")){
//				cmd = "msaprobs -o " + faoFile + " " + faiFile;
//			}else if (msa.startsWith("mafft")){
//				cmd = "mafft_wrapper.sh  " + faiFile + " " + faoFile;
			}else if(msa.isEmpty() || msa.startsWith("none")){
				consensus.setName(prefix+"_consensus");
				return consensus;
			}else{
				if(HybridAssembler.VERBOSE)
					LOG.info("Unknown msa function " + msa);
				return null;
			}
			if(HybridAssembler.VERBOSE)
				LOG.info("Running " + cmd);
			Process process = Runtime.getRuntime().exec(cmd);
			process.waitFor();
			if(HybridAssembler.VERBOSE)
				LOG.info("Done " + cmd);
		}
				
//		if ("poa".equals(msa)){
//			SequenceBuilder sb = new SequenceBuilder(Alphabet.DNA(), (int) ((1+BDGraph.R_TOL)*distance)+2*BDGraph.getKmerSize());
//			BufferedReader bf =  FastaReader.openFile(faoFile);
//			String line = bf.readLine();
//			while ( (line = bf.readLine()) != null){
//				if (line.startsWith("CONSENS0")){
//					for (int i = 10;i < line.length();i++){
//						char c = line.charAt(i);
//						int base = DNA.DNA().char2int(c);
//						if (base >= 0 && base < 4)
//							sb.append((byte)base);
//					}//for							
//				}//if
//			}//while
//			sb.setName(prefix+"_consensus");
//			LOG.info(sb.getName() + "  " + sb.length());
//			return sb.toSequence();
//		}

		//3.0 Read in multiple alignment
		ArrayList<Sequence> seqList = new ArrayList<Sequence>();
		{
			SequenceReader msaReader = FastaReader.getReader(faoFile);
			Sequence nSeq = null;
			while ((nSeq = msaReader.nextSequence(Alphabet.DNA())) != null) {
				seqList.add(nSeq);						
			}
			msaReader.close();
		}

		//4.0 get consensus			
		{
			int [] coef = new int[seqList.size()];			
			for (int y = 0; y < seqList.size(); y++){
				coef[y] = 1;				
			}
			//TODO: combine error profiles?? 
			int [] counts = new int[6];

			SequenceBuilder sb = new SequenceBuilder(Alphabet.DNA(), seqList.get(0).length());
			for (int x = 0; x < seqList.get(0).length();x++){		
				Arrays.fill(counts, 0);
				for (int y = 0; y < seqList.size(); y++){					
					byte base = seqList.get(y).getBase(x);
					if (base >= 6) 
						counts[4] += coef[y];//N
					else
						counts[base] += coef[y];
				}//for y
				int maxIdx = 0;
				for (int y = 1; y < counts.length; y++){
					if (counts[y] > counts[maxIdx])
						maxIdx = y;					
				}//for y
				if (maxIdx < Alphabet.DNA.GAP){//not a gap
					sb.append((byte)maxIdx);
				}//if
			}//for x
			sb.setName(prefix+"_consensus");
			if(HybridAssembler.VERBOSE)
				LOG.info(sb.getName() + "  " + sb.length());
			consensus = sb.toSequence();
		}

		return consensus;
	}
	
	
	//no check: use with care
	public static void saveReadToDisk(AlignedRead read){
		if(BDGraph.getReadsNumOfBrg(read.getEndingsID())>=BDGraph.MAX_LISTING)
			return;
		
		if(read.getEFlag()>=3)
			read.saveCorrectedSequenceInBetween();
		else
			read.splitAtPotentialAnchors().forEach(r->r.saveCorrectedSequenceInBetween());
			
	}
	
	/*
	 * Note that read sequence from SAMRecord could be reverse-complemented
	 */
	public static Sequence getQueryReadFromSAMRecord(SAMRecord sam){
		Sequence retval = new Sequence(Alphabet.DNA5(), sam.getReadString(), sam.getReadName());
		if(sam.getReadNegativeStrandFlag())
			retval=Alphabet.DNA.complement(retval);
		
		return retval;
	}
}
