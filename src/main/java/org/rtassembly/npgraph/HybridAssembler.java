package org.rtassembly.npgraph;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.lang.ProcessBuilder.Redirect;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Scanner;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

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
import japsa.seq.SequenceReader;


public class HybridAssembler {
    private static final Logger LOG = LoggerFactory.getLogger(HybridAssembler.class);
	//setting parameter for the GUI
    private boolean ready=false, overwrite=false;
    private String prefix = "/tmp/";
	private String 	mm2Path="", 
					mm2Opt="-t4 -x map-ont -k15 -w5";

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
		File shortReadsInputFile = new File(srInput);
		if(shortReadsInputFile.isFile()) {
			try {
				shortReadsInput=shortReadsInputFile.getCanonicalPath();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return;
			}
			String fn = srInput.toLowerCase();
			if(	fn.endsWith(".fastg")) 
				setShortReadsInputFormat("fastg");
			else if(fn.endsWith(".gfa"))
				setShortReadsInputFormat("gfa");
			else
				setShortReadsInputFormat("");
		}
	}
	public String getShortReadsInput() {return shortReadsInput;}
	
	public void setLongReadsInput(String lrInput) {
		File longReadsInputFile = new File(lrInput);
		if(longReadsInputFile.isFile()) {
			try {
				longReadsInput=longReadsInputFile.getCanonicalPath();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return;
			}
			String fn = lrInput.toLowerCase();
			if(	fn.endsWith(".fasta") || fn.endsWith(".fa") || fn.endsWith("fna")
				|| fn.endsWith(".fastq") || fn.endsWith(".fq")
				|| fn.endsWith(".fasta.gz") || fn.endsWith(".fa.gz") || fn.endsWith("fna.gz")
				|| fn.endsWith(".fastq.gz") || fn.endsWith(".fq.gz") 
				) 
				setLongReadsInputFormat("fasta/fastq");
			else if(fn.endsWith(".sam") || fn.endsWith(".bam")) 
				setLongReadsInputFormat("sam/bam");
			else
				setLongReadsInputFormat("");
		}
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
	public GraphWatcher observer;
	
	public HybridAssembler(){
//		origGraph=new BidirectedGraph("batch");
		simGraph=new BidirectedGraph("real");
//		rtComponents = new ConnectedComponents();
		
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
					simGraph.outputFASTA(prefix+"/assembly_graph.fasta");
					if(!checkMinimap2()) {
							LOG.error("Dependancy check failed! Please config to the right version of minimap2!");
							return false;
					}
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
		
		observer = new GraphWatcher(simGraph);
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
						prefix+"/assembly_graph.mmi",
						"-"
						).
						redirectInput(Redirect.INHERIT);
			}else{
				pb = new ProcessBuilder(mm2Path+"/minimap2", 
						"-a",
						mm2Opt,
						"-K",
						"20000",
						prefix+"/assembly_graph.mmi",
						longReadsInput
						);
			}

//			mm2Process  = pb.redirectError(ProcessBuilder.Redirect.to(new File("/dev/null"))).start();
			mm2Process  = pb.redirectError(ProcessBuilder.Redirect.to(new File(prefix+"/mm2.log"))).start();

			LOG.info("minimap2 started!");			

			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(mm2Process.getInputStream()));

		}
		SAMRecordIterator iter = reader.iterator();

		String readID = "";
		Sequence nnpRead = null;
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
			
			if (!readID.equals("") && !readID.equals(myRec.readID)) {	
				synchronized(simGraph) {
					List<BidirectedPath> paths=simGraph.uniqueBridgesFinding(nnpRead, samList);
					if(paths!=null){	
						for(BidirectedPath path:paths) 
						{
							//path here is already unique! (2 unique ending nodes)
					    	if(simGraph.reduceUniquePath(path)) {
					    		//TODO: replace with proper observer.update()
//					    		observer.forFunUpdate();					    		
					    		GraphUtil.redrawGraphComponents(simGraph);
					    	}
						}
					}
				}
				samList = new ArrayList<Alignment>();
				nnpRead = new Sequence(Alphabet.DNA5(), rec.getReadString(), "R" + readID);
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
	
	public void postProcessGraph() throws IOException{
		//Take the current best path among the candidate of a bridge and connect the bridge(greedy)
		for(GoInBetweenBridge brg:simGraph.getUnsolvedBridges()){
			System.out.printf("Last attempt on incomplete bridge %s : anchors=%d \n %s \n", brg.getEndingsID(), brg.getNumberOfAnchors(), brg.getAllPossiblePaths());
			if(brg.getCompletionLevel()<=2) {//examine bridge with completion level = 2 that unable to connected
				brg.steps.connectBridgeSteps(true);
			}
			
			if(brg.getCompletionLevel()>=3) {
				for(BidirectedPath p:brg.getBestPath().chopPathAtAnchors())
					simGraph.reduceUniquePath(p);
			}
			else
				System.out.printf("Last attempt failed \n");


		}
		
		//TODO: traverse for the last time,remove redundant edges
		//may want to run consensus to determine the final path
		
		
		//Finally redraw the graph and output
		GraphUtil.redrawGraphComponents(simGraph);
		
        observer.linearComponentsDecomposition();
		observer.outputFASTA(getPrefix()+"npgraph_assembly.fasta");

	}
	


    
    @SuppressWarnings("resource")
	public static void promptEnterKey(){
    	   System.out.println("Press \"ENTER\" to continue...");
    	   Scanner scanner = new Scanner(System.in);
    	   scanner.nextLine();
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
		
		hbAss.setShortReadsInput("/home/sonhoanghguyen/Projects/scaffolding/data/spades_3.7/EcK12S-careful/assembly_graph.fastg");
		hbAss.setShortReadsInputFormat("fastg");
		hbAss.prepareShortReadsProcess(false);
		hbAss.setLongReadsInput("/home/sonhoanghguyen/Projects/scaffolding/data/spades_3.7/EcK12S-careful/assembly_graph.sam");
		hbAss.setLongReadsInputFormat("sam/bam");
		hbAss.prepareLongReadsProcess();
		hbAss.assembly();

	}
	
}
