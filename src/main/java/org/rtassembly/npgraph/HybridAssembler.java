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
import java.util.concurrent.atomic.AtomicInteger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.graphstream.algorithm.ConnectedComponents;
import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import com.google.common.util.concurrent.AtomicDouble;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;


public class HybridAssembler {
    private static final Logger LOG = LoggerFactory.getLogger(HybridAssembler.class);
	//setting parameter for the GUI
    private boolean ready=false, overwrite=false;
    private String prefix = "/tmp/";
	private String 	mm2Path="", 
					mm2Opt="-t4 -x map-ont";

	private String shortReadsInput="", longReadsInput="";
	private String shortReadsInputFormat="", longReadsInputFormat="";
	Process mm2Process = null;
	private boolean stop=false;
	//Getters and Setters
	//==============================================================================================//
	public void setReady(boolean isReady) {ready=isReady;}
	public boolean getReady() {return ready;}
	
	public void setOverwrite(boolean overwrite) {this.overwrite=overwrite;}
	public boolean getOverwrite() {return overwrite;}
	
	public void setPrefix(String prefix) {this.prefix=prefix;}
	public String getPrefix() {return prefix;}
	
	public void setMinimapPath(String path) {mm2Path=path;}
	public String getMinimapPath() {return mm2Path;}
	
	public void setMinimapOpts(String setting) {mm2Opt=setting;}
	public String getMinimapOpts() {return mm2Opt;}
	
	
	public void setShortReadsInput(String srInput) {
		shortReadsInput=srInput;
//		String fn = srInput.toLowerCase();
//		if(	fn.endsWith(".fastg") || fn.endsWith(".gfa") ) 
//			setShortReadsInputFormat("fastg");
//		else if(fn.endsWith(".sam") || fn.endsWith(".bam")) 
//			setShortReadsInputFormat("gfa");
//		else
//			setShortReadsInputFormat("");
	}
	public String getShortReadsInput() {return shortReadsInput;}
	
	public void setLongReadsInput(String lrInput) {
		longReadsInput=lrInput;
//		String fn = lrInput.toLowerCase();
//		if(	fn.endsWith(".fasta") || fn.endsWith(".fa") || fn.endsWith("fna")
//			|| fn.endsWith(".fastq") || fn.endsWith(".fq") 
//			) 
//			setLongReadsInputFormat("fasta/fastq");
//		else if(fn.endsWith(".sam") || fn.endsWith(".bam")) 
//			setLongReadsInputFormat("sam/bam");
//		else
//			setLongReadsInputFormat("");
	}
	public String getLongReadsInput() {return longReadsInput;}
	
	public void setShortReadsInputFormat(String srInputFormat) {shortReadsInputFormat=srInputFormat;}
	public String getShortReadsInputFormat() {return shortReadsInputFormat;}
	
	public void setLongReadsInputFormat(String lrInputFormat) {longReadsInputFormat=lrInputFormat;}
	public String getLongReadsInputFormat() {return longReadsInputFormat;}
	
	public synchronized void setStopSignal(boolean stop) {this.stop=stop;}
	public synchronized boolean getStopSignal() {return stop;}
	//===============================================================================================//
	
	//Operational variables
//	final BidirectedGraph origGraph;
	public BidirectedGraph simGraph; //original and simplified graph should be separated, no???
	public ConnectedComponents rtComponents;
	
	public HybridAssembler(){
//		origGraph=new BidirectedGraph("batch");
		simGraph=new BidirectedGraph("real");
		rtComponents = new ConnectedComponents();
		
		simGraph.setAttribute("ui.quality");
		simGraph.setAttribute("ui.antialias");
		
	}
		
	
	//Indexing reference, prepare for alignment...
	public boolean prepareLongReadsProcess(){
		//if long reads data not given in SAM/BAM, need to invoke minimap2
        if(longReadsInputFormat.toLowerCase().startsWith("fast")) {
			File indexFile=new File(prefix+"/assembly_graph.mmi");
			if(overwrite || !indexFile.exists()) {						
				try{
					simGraph.printNodeSequencesToFile(prefix+"/assembly_graph.fasta");
//					if(!checkMinimap2()) {
//							LOG.error("Dependancy check failed! Please config to the right version of minimap2!");
//							return false;
//					}
					ProcessBuilder pb = new ProcessBuilder(mm2Path+"/minimap2", mm2Opt,"-d", prefix+"/assembly_graph.mmi",prefix+"/assembly_graph.fasta");
					Process indexProcess =  pb.start();
					indexProcess.waitFor();
					
				}catch (IOException | InterruptedException e){
					System.err.println("Issue when indexing with minimap2: \n" + e.getMessage());
					return false;
				}
			}
        }
        return true;
	}
	//Loading the graph, doing preprocessing
	//binning, ...
	public boolean prepareShortReadsProcess(boolean useSPAdesPaths) {
		try {
			if(shortReadsInputFormat.toLowerCase().equals("gfa")) 
				GraphUtil.loadFromGFA(shortReadsInput, simGraph, useSPAdesPaths);
			else if(shortReadsInputFormat.toLowerCase().equals("fastg"))
				GraphUtil.loadFromFASTG(shortReadsInput, simGraph, useSPAdesPaths);
			else 				
				throw new IOException("assembly graph file must have .gfa or .fastg extension!");
			
		}catch(IOException e) {
			System.err.println("Issue when loading pre-assembly: \n" + e.getMessage());
			return false;
		}
		rtComponents.init(simGraph);
		return true;
	}

	/**
	 * SHN modified the default aligner to minimap2
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void assembly() 
			throws IOException, InterruptedException{

		LOG.info("Scaffolding ready at {}", new Date());

		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader = null;

		if (longReadsInputFormat.endsWith("am")){//bam or sam
			if ("-".equals(longReadsInput))
				reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
			else
				reader = SamReaderFactory.makeDefault().open(new File(longReadsInput));	
		}else{
			LOG.info("Starting alignment by minimap2 at {}", new Date());
			ProcessBuilder pb = null;
			if ("-".equals(longReadsInput)){
				pb = new ProcessBuilder(mm2Path+"/minimap2", 
						"-a",
						mm2Opt,
						"-K",
						"20000",
						prefix+"./assembly_graph.mmi",
						"-"
						).
						redirectInput(Redirect.INHERIT);
			}else{
				pb = new ProcessBuilder(mm2Path+"/minimap2", 
						"-a",
						mm2Opt,
						"-K",
						"20000",
						prefix+"./assembly_graph.mmi",
						longReadsInput
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

		ArrayList<Alignment> samList =  new ArrayList<Alignment>();// alignment record of the same read;	
		
		while (iter.hasNext()) {
			if(getStopSignal())
				break;
			
			SAMRecord rec = iter.next();
			if (rec.getReadUnmappedFlag())
				continue;
			
			if (rec.getMappingQuality() < Alignment.MIN_QUAL)
				continue;
			
			String refName = rec.getReferenceName();
			String refID = refName.split("_").length > 1 ? refName.split("_")[1]:refName;
			
			if (simGraph.getNode(refID)==null)
				continue;
			Alignment myRec = new Alignment(rec, (BidirectedNode) simGraph.getNode(refID)); 

			//////////////////////////////////////////////////////////////////
			//FIXME: optimize
			// make list of alignments of the same (Nanopore) read. 
			//not the first occurrance				
			if (!readID.equals("") && !readID.equals(myRec.readID)) {	
				synchronized(simGraph) {
					List<BidirectedPath> paths=simGraph.uniqueBridgesFinding(samList);
					if(paths!=null)						
						for(BidirectedPath path:paths) 
						{
							//path here is already unique! (2 unique ending nodes)
					    	if(simGraph.reduceUniquePath(path)) {

					    		GraphExplore.redrawGraphComponents(simGraph);
					    		
//								This is just for fun
//					    		LOG.info("==========================================================================");
//					    		LOG.info("\nTotal number of components: {} \ncomponents containing more than 1: {} \nsize of biggest component: {}", 
//					    					rtComponents.getConnectedComponentsCount(),rtComponents.getConnectedComponentsCount(2),rtComponents.getGiantComponent().size());					    		
//					    		LOG.info("==========================================================================");    		
					    		if(currentNumOfComponents != rtComponents.getConnectedComponentsCount()) {
						    		currentNumOfComponents = rtComponents.getConnectedComponentsCount();

						    		//Hide components with no markers! Optimize it to work dynamically
						    		ArrayList<Node> cleanup = new ArrayList<>();
						    		for (Iterator<ConnectedComponents.ConnectedComponent> compIter = rtComponents.iterator(); compIter.hasNext(); ) {
						    			ConnectedComponents.ConnectedComponent comp = compIter.next();
	
						    			AtomicDouble lengthWeightedCov = new AtomicDouble(0.0);
						    			ArrayList<Node> tmp = new ArrayList<>();
						    			comp.nodes().forEach(n->{
						    				lengthWeightedCov.getAndAdd(n.getNumber("cov")*(n.getNumber("len")-BidirectedGraph.getKmerSize()));
						    				tmp.add(n);
						    			});
						    			if(lengthWeightedCov.get() < 10000*BidirectedGraph.RCOV)
						    				cleanup.addAll(tmp);
						    		}
						    		for(Node n:cleanup) {
					    				n.setAttribute("ui.hide");
					    				n.edges().forEach(e->e.setAttribute("ui.hide"));
////					    			simGraph.removeNode(n); //this faster but careful here!!!
						    		}
					    		}
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
	
	public void lastAttempt(){
		//TODO: traverse for the last time,remove redundant edges, infer the path...
		//may want to run consensus to determine the final path
		for(BidirectedBridge brg:simGraph.getUnsolvedBridges()){
			System.out.println("Last attempt: " + brg.getBridgeString());
			if(brg.getBridgeStatus()==0)
				simGraph.reduceUniquePath(brg.fullPaths.get(0));
			else
				System.out.println("bridge contain no path! ignored");
		}
		GraphExplore.redrawGraphComponents(simGraph);
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
			    	if(simGraph.reduceUniquePath(path))
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
    
    public boolean checkMinimap2() {

		ProcessBuilder pb = new ProcessBuilder(mm2Path+"/minimap2","-V").redirectErrorStream(true);
		Process process;
		try {
			process = pb.start();

			//Allen changes: BWA process doesn't produce gzip-compressed output
			BufferedReader bf = SequenceReader.openInputStream(process.getInputStream());
	
	
			String line;
			String version = "";
			Pattern versionPattern = Pattern.compile("^(\\d+\\.\\d+).*");
			Matcher matcher=versionPattern.matcher("");
			
			while ((line = bf.readLine())!=null){				
				matcher.reset(line);
				if (matcher.find()){
				    version = matcher.group(1);
				    break;//while
				}
				
								
			}	
			bf.close();
			
			if (version.length() == 0){
				System.err.println("ERROR: minimap2 command not found. Please install minimap2 and set the appropriate PATH variable;\n"
									+ "	or run the alignment yourself and provide the SAM file instead of FASTA/Q file.");
				return false;
			}else{
				System.out.println("minimap version: " + version);
				if (version.compareTo("2.0") < 0){
					System.err.println(" ERROR: require minimap version 2 or above!");
					return false;
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			System.err.println("Error running: " + mm2Path + "/minimap2 \n" + e.getMessage());
			e.printStackTrace();			
			return false;
		}
		
		return true;
			
    }

	public static void main(String[] argv) throws IOException, InterruptedException{
		HybridAssembler hbAss = new HybridAssembler();
		
		hbAss.setShortReadsInput(GraphExplore.spadesFolder+"EcK12S-careful/assembly_graph.fastg");
		hbAss.setShortReadsInputFormat("fastg");
		hbAss.prepareShortReadsProcess(false);
		hbAss.setLongReadsInput(GraphExplore.spadesFolder+"EcK12S-careful/assembly_graph.sam");
		hbAss.setLongReadsInputFormat("sam/bam");
		hbAss.prepareLongReadsProcess();
		hbAss.assembly();

	}
	
}
