package org.rtassembly.npgraph;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.util.ArrayList;
import java.util.Date;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;

import org.graphstream.algorithm.ConnectedComponents;
import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;


public class HybridAssembler {
    private static final Logger LOG = LoggerFactory.getLogger(HybridAssembler.class);
	
//	final BidirectedGraph origGraph;
	public BidirectedGraph simGraph; //original and simplified graph should be separated, no???
	public ConnectedComponents rtComponents;
	public HybridAssembler(){
//		origGraph=new BidirectedGraph("batch");
		simGraph=new BidirectedGraph("real");
		rtComponents = new ConnectedComponents();
	}
	
	
	public HybridAssembler(String graphInputFile) throws IOException{
		this();

		if(graphInputFile.toLowerCase().endsWith(".gfa")) 
			GraphUtil.loadFromGFA(graphInputFile, simGraph);
		else if(graphInputFile.toLowerCase().endsWith(".fastg"))
			GraphUtil.loadFromFASTG(graphInputFile, simGraph);
		else {
			System.err.println("Assembly graph file must have .gfa or .fastg extension!");
			System.exit(1);
		}
		
		rtComponents.init(simGraph);

	}
	/**
	 * SHN modified the default aligner to minimap2
	 * @param bamFile
	 * @param readNumber
	 * @param timeNumber
	 * @param minCov
	 * @param qual
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void assembly(String inFile, String inFormat, String mm2Path, String mm2Preset, int mm2Threads, String mm2Index) 
			throws IOException, InterruptedException{

		LOG.info("Scaffolding ready at {}", new Date());

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = null;

		Process mm2Process = null;

		if (inFormat.endsWith("am")){//bam or sam
			if ("-".equals(inFile))
				reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
			else
				reader = SamReaderFactory.makeDefault().open(new File(inFile));	
		}else{
			LOG.info("Starting alignment by minimap2 at {}", new Date());
			ProcessBuilder pb = null;
			if ("-".equals(inFile)){
				pb = new ProcessBuilder(mm2Path+"minimap2", 
						"-t",
						"" + mm2Threads,
						"-ax",
						mm2Preset,
//						"-I",
//						"40g",
						"-K",
						"20000",
						mm2Index,
						"-"
						).
						redirectInput(Redirect.INHERIT);
			}else{
				pb = new ProcessBuilder(mm2Path+"minimap2", 
						"-t",
						"" + mm2Threads,
						"-ax",
						mm2Preset,
//						"-I",
//						"40g",
						"-K",
						"20000",
						mm2Index,
						inFile
						);
			}

			mm2Process  = pb.redirectError(ProcessBuilder.Redirect.to(new File("/dev/null"))).start();
//			mm2Process  = pb.redirectError(ProcessBuilder.Redirect.to(new File("mm2.err"))).start();

			LOG.info("minimap2 started!");			

			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(mm2Process.getInputStream()));

		}
		SAMRecordIterator iter = reader.iterator();

		String readID = "";
		int currentNumOfComponents=rtComponents.getConnectedComponentsCount();

		ArrayList<Alignment> samList =  new ArrayList<Alignment>();;// alignment record of the same read;	
		
		while (iter.hasNext()) {
			SAMRecord rec = iter.next();
			if (rec.getReadUnmappedFlag())
				continue;
			
			if (rec.getMappingQuality() < Alignment.MIN_QUAL)
				continue;
			
			String refName = rec.getReferenceName();
			String refID = refName.split("_").length > 1 ? refName.split("_")[1]:refName;
			
			if (simGraph.getNode(refID)==null)
				continue;
			Alignment myRec = new Alignment(rec, simGraph.getNode(refID)); //FIXME: optimize

			//////////////////////////////////////////////////////////////////
			// make list of alignments of the same (Nanopore) read. 

			//not the first occurrance				
			if (!readID.equals("") && !readID.equals(myRec.readID)) {	
				synchronized(simGraph) {
					List<BidirectedPath> paths=simGraph.uniqueBridgesFinding(samList);
					if(paths!=null)						
						for(BidirectedPath p:paths) 
						{
					    	if(simGraph.reduce(p)) {

					    		GraphExplore.redrawGraphComponents(simGraph);

////					    		LOG.info("==========================================================================");
////					    		LOG.info("\nTotal number of components: {} \ncomponents containing more than 1: {} \nsize of biggest component: {}", 
////					    					rtComponents.getConnectedComponentsCount(),rtComponents.getConnectedComponentsCount(2),rtComponents.getGiantComponent().size());
////					    		
////					    		LOG.info("==========================================================================");    		
//					    		if(currentNumOfComponents != rtComponents.getConnectedComponentsCount()) {
//						    		currentNumOfComponents = rtComponents.getConnectedComponentsCount();
//
//						    		//Hide components with no markers! Optimize it to work dynamically
//						    		ArrayList<Node> cleanup = new ArrayList<>();
//						    		for (Iterator<ConnectedComponents.ConnectedComponent> compIter = rtComponents.iterator(); compIter.hasNext(); ) {
//						    			ConnectedComponents.ConnectedComponent comp = compIter.next();
//	
//						    			int numOfMarker=0;
//						    			ArrayList<Node> tmp = new ArrayList<>();
//						    			for(Node n:comp.getEachNode()) {
//						    				if(BidirectedGraph.isMarker(n)) 
//						    					numOfMarker++;
//						    				tmp.add(n);
//						    			}
//						    			if(numOfMarker==0)
//						    				cleanup.addAll(tmp);
//						    		}
//						    		for(Node n:cleanup) {
//					    				n.addAttribute("ui.hide");
//						    			for(Edge e:n.getEachEdge())
//						    				e.addAttribute("ui.hide");
//	//					    			simGraph.removeNode(n); //this faster but careful here!!!
//						    		}
//					    		}
////					    		promptEnterKey();
					    	}
						}
				}
				samList = new ArrayList<Alignment>();
			}	
			readID = myRec.readID;
			samList.add(myRec); 
		}// while
		iter.close();
		reader.close();

		if (mm2Process != null){
			mm2Process.destroy();
		}	

	}
	
	@Deprecated
	public void assembly(String bamFile) throws IOException{
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);

		SamReader reader;
		if ("-".equals(bamFile))
			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			reader = SamReaderFactory.makeDefault().open(new File(bamFile));	

		SAMRecordIterator iter = reader.iterator();

		String readID = "";
		int currentNumOfComponents=rtComponents.getConnectedComponentsCount();
		//ReadFilling readFilling = null;
		ArrayList<Alignment> samList =  new ArrayList<Alignment>();;// alignment record of the same read;	
		while (iter.hasNext()) {
			SAMRecord rec = iter.next();
			if (rec.getReadUnmappedFlag())
				continue;
//			if (rec.getMappingQuality() < qual)
//				continue;
			
			String refID = rec.getReferenceName().split("_")[1];
			if (simGraph.getNode(refID)==null)
				continue;
			Alignment myRec = new Alignment(rec, simGraph.getNode(refID)); //FIXME: optimize

			//////////////////////////////////////////////////////////////////
			// make list of alignments of the same (Nanopore) read. 

			//not the first occurrance				
			if (!readID.equals("") && !readID.equals(myRec.readID)) {	
				synchronized(simGraph) {
					//p=origGraph.pathFinding(samList);
//					List<BidirectedPath> paths=simGraph.pathFinding(samList);// the graph MUST be the same as from new Alignment(...)
					List<BidirectedPath> paths=simGraph.uniqueBridgesFinding(samList);
					if(paths!=null)						
						for(BidirectedPath p:paths) 
						{
					    	if(simGraph.reduce(p)) {

					    		GraphExplore.redrawGraphComponents(simGraph);

////					    		LOG.info("==========================================================================");
////					    		LOG.info("\nTotal number of components: {} \ncomponents containing more than 1: {} \nsize of biggest component: {}", 
////					    					rtComponents.getConnectedComponentsCount(),rtComponents.getConnectedComponentsCount(2),rtComponents.getGiantComponent().size());
////					    		
////					    		LOG.info("==========================================================================");    		
//					    		if(currentNumOfComponents != rtComponents.getConnectedComponentsCount()) {
//						    		currentNumOfComponents = rtComponents.getConnectedComponentsCount();
//
//						    		//Hide components with no markers! Optimize it to work dynamically
//						    		ArrayList<Node> cleanup = new ArrayList<>();
//						    		for (Iterator<ConnectedComponents.ConnectedComponent> compIter = rtComponents.iterator(); compIter.hasNext(); ) {
//						    			ConnectedComponents.ConnectedComponent comp = compIter.next();
//	
//						    			int numOfMarker=0;
//						    			ArrayList<Node> tmp = new ArrayList<>();
//						    			for(Node n:comp.getEachNode()) {
//						    				if(BidirectedGraph.isMarker(n)) 
//						    					numOfMarker++;
//						    				tmp.add(n);
//						    			}
//						    			if(numOfMarker==0)
//						    				cleanup.addAll(tmp);
//						    		}
//						    		for(Node n:cleanup) {
//					    				n.addAttribute("ui.hide");
//						    			for(Edge e:n.getEachEdge())
//						    				e.addAttribute("ui.hide");
//	//					    			simGraph.removeNode(n); //this faster but careful here!!!
//						    		}
//					    		}
////					    		promptEnterKey();
					    	}
						}
				}


//				reduce2(p);
				samList = new ArrayList<Alignment>();
				//readID = myRec.readID;	
			}	
			readID = myRec.readID;
			samList.add(myRec); 
		}// while
		iter.close();

		//outOS.close();
		reader.close();		
	
	}
	
	
    /*
     * Read paths from contigs.path and reduce the graph
     */
    public void reduceFromSPAdesPaths(String paths) throws IOException{

		BufferedReader pathReader = new BufferedReader(new FileReader(paths));
		
		String s="", curpath="";
		//Read contigs from contigs.paths
		boolean flag=false;
		while((s=pathReader.readLine()) != null){
			if(s.contains("NODE")){
				if(flag){
					BidirectedPath path=new BidirectedPath(simGraph, curpath);
			    	if(simGraph.reduce(path))
			    		GraphExplore.redrawGraphComponents(simGraph);
//			    	reduce2(path);
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
	


    
    @SuppressWarnings("resource")
	public static void promptEnterKey(){
    	   System.out.println("Press \"ENTER\" to continue...");
    	   Scanner scanner = new Scanner(System.in);
    	   scanner.nextLine();
    	}
    
    protected void sleep() {
        try { Thread.sleep(1000); } catch (Exception e) {}
    }

	public static void main(String[] argv) throws IOException{
		HybridAssembler hbAss = new HybridAssembler(GraphExplore.spadesFolder+"EcK12S-careful/assembly_graph.fastg");
		//For SAM file, run bwa first on the edited assembly_graph.fastg by running:
		//awk -F '[:;]' -v q=\' 'BEGIN{flag=0;}/^>/{if(index($1,q)!=0) flag=0; else flag=1;}{if(flag==1) print $1;}' ../EcK12S-careful/assembly_graph.fastg > Eck12-careful.fasta
		//TODO: need to make this easier

		hbAss.assembly(GraphExplore.spadesFolder+"EcK12S-careful/assembly_graph.sam");

	}
	
}
