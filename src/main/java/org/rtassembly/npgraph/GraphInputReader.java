package org.rtassembly.npgraph;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.graphstream.graph.implementations.AbstractNode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.util.CigarUtil;
import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;

public class GraphInputReader {

	private static final Logger LOG = LoggerFactory.getLogger(GraphInputReader.class);
	/**********************************************************************************
	 * ****************************Algorithms go from here*****************************
	 */
    //TODO: read from ABySS assembly graph (graph of final contigs, not like SPAdes)
    
    public static void loadFromFASTG(String graphFile, BidirectedGraph graph) throws IOException{
        graph.setAutoCreate(true);
        graph.setStrict(false);
		//1. next iterate over again to read the connections
		SequenceReader reader = new FastaReader(graphFile);
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
			AbstractNode node = graph.addNode(nodeID); //or get the existing node prototype created below (to set the attributes)
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
					AbstractNode nbr = graph.addNode(neighborID); //just need a prototype, attributes can be set later...

					potentialEdgeSet.add(new EdgeComponents(node,nbr,dir0,dir1)); //edges' prototype
					
				}
			}
			
		}

		reader.close();
		
		for(EdgeComponents ec:potentialEdgeSet) {
			BidirectedEdge e = graph.addEdge(ec.n1, ec.n2, ec.dir1, ec.dir2);
			graph.updateGraphMap(e, new BidirectedPath(e));
		}
		//rough estimation of kmer used
		if((shortestLen-1) != BidirectedGraph.getKmerSize()){
			BidirectedGraph.setKmerSize(shortestLen-1);
			for(Edge e:graph.getEdgeSet()){
				((BidirectedEdge)e).changeKmerSize(BidirectedGraph.KMER);
			}
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
			Sequence nseq = n.getAttribute("seq");
			double astats = nseq.length()*BidirectedGraph.RCOV/BidirectedGraph.ILLUMINA_READ_LENGTH
							-Math.log(2)*n.getNumber("cov")*nseq.length()/BidirectedGraph.ILLUMINA_READ_LENGTH;
			astats*=Math.log10(Math.E);
			n.setAttribute("astats", astats);
			LOG.info("{} Read coverage = {} Length = {} A-stats = {}", n.getAttribute("name"), n.getNumber("cov"), n.getNumber("len"), astats );
		}
		LOG.info("Estimated read coverage = {} Total contigs length = {}", BidirectedGraph.RCOV,totContigsLen );
		
    }
    
    
    
    public static void loadFromGFA(String graphFile, BidirectedGraph graph) throws IOException{
        graph.setAutoCreate(true);
        graph.setStrict(false);
		//1. next iterate over again to read the connections
        BufferedReader reader = new BufferedReader(new FileReader(graphFile));
        String line=null;
		Sequence seq;
		int shortestLen = 10000;
		
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
				AbstractNode node = graph.addNode(nodeID); //or get the existing node prototype created below (to set the attributes)
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
				BidirectedNode 	n0=graph.getNode(gfaFields[1]),
								n1=graph.getNode(gfaFields[3]);
				boolean dir0=gfaFields[2].equals("+")?true:false,
						dir1=gfaFields[4].equals("+")?false:true;
				BidirectedEdge e = graph.addEdge(n0, n1, dir0, dir1);
				graph.updateGraphMap(e, new BidirectedPath(e));
				
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
//				if(gfaFields.length>3 && gfaFields[2].contains(",")) {
//					BidirectedPath path=new BidirectedPath(graph, gfaFields[2]);
//			    	if(graph.reduce(path))
//			    		GraphExplore.redrawGraphComponents(graph);
//				}
				break;

				default:throw new IllegalStateException("Unrecognized GFA field: " + type.toUpperCase().trim());
			}
			
		
		}

		reader.close();

		//rough estimation of kmer used
		if((shortestLen-1) != BidirectedGraph.getKmerSize()){
			BidirectedGraph.setKmerSize(shortestLen-1);
			for(Edge e:graph.getEdgeSet()){
				((BidirectedEdge)e).changeKmerSize(BidirectedGraph.KMER);
			}
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
			Sequence nseq = n.getAttribute("seq");
			double astats = nseq.length()*BidirectedGraph.RCOV/BidirectedGraph.ILLUMINA_READ_LENGTH
							-Math.log(2)*n.getNumber("cov")*nseq.length()/BidirectedGraph.ILLUMINA_READ_LENGTH;
			astats*=Math.log10(Math.E);
			n.setAttribute("astats", astats);
			LOG.info("{} Read coverage = {} Length = {} A-stats = {}", n.getAttribute("name"), n.getNumber("cov"), n.getNumber("len"), astats );
		}
		LOG.info("Estimated read coverage = {} Total contigs length = {}", BidirectedGraph.RCOV,totContigsLen );
		
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
