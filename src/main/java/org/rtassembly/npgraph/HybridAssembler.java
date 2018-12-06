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
import javafx.beans.property.StringProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.SimpleBooleanProperty;
@SuppressWarnings("restriction")
public class HybridAssembler {
    private static final Logger LOG = LoggerFactory.getLogger(HybridAssembler.class);
	//setting parameter for the GUI
    private boolean ready=false;
    private BooleanProperty overwrite;
    private StringProperty 	prefix,
    						aligner,
    						alignerPath, 
							alignerOpt,
							errorLog,
							shortReadsInput, 
							binReadsInput, 
							longReadsInput,
							shortReadsInputFormat, 
							longReadsInputFormat;
	
	Process alignmentProcess = null;
	private boolean stop=false;
	//Getters and Setters
	//==============================================================================================//
	public void setReady(boolean isReady) {ready=isReady;}
	public boolean getReady() {return ready;}
	
	public final void setOverwrite(boolean owr) {overwrite.set(owr);}
	public final boolean getOverwrite() {return overwrite.get();}
	public BooleanProperty overwriteProperty() {return overwrite;}
	
	public final void setPrefix(String output) {prefix.set(output);}
	public final String getPrefix() {return prefix.get();}
	public StringProperty prefixProperty() {return prefix;}
	
	public final void setAligner(String tool) {	aligner.set(tool);}
	public final String getAligner() {return aligner.get();}
	public StringProperty alignerProperty(){return aligner;}
	
	public final void setAlignerPath(String path) {alignerPath.set(path);}
	public final String getAlignerPath() {return alignerPath.get();}
	public StringProperty alignerPathProperty(){return alignerPath;}
	
	public final void setAlignerOpts(String setting) {alignerOpt.set(setting);}
	public final String getAlignerOpts() {return alignerOpt.get();}
	public StringProperty alignerOptProperty(){return alignerOpt;}
	
	public final void setErrorLog(String log) {errorLog.set(log);}
	public final String getErrorLog() {return errorLog.get();}
	public StringProperty errorLogProperty(){return errorLog;}
	
	public final String getFullPathOfAligner() {return getAlignerPath()+"/"+getAligner();}
	
	public final void setBinReadsInput(String brInput) {binReadsInput.set(brInput);}
	public final String getBinReadsInput() {return binReadsInput.get();}
	public StringProperty binReadsInputProperty() {return binReadsInput;}
	
	public final void setShortReadsInput(String srInput) {
		shortReadsInput.set(srInput);
	}
	public final String getShortReadsInput() {return shortReadsInput.get();}
	public StringProperty shortReadsInputProperty() {return shortReadsInput;}
	
	public final void setLongReadsInput(String lrInput) {
		longReadsInput.set(lrInput);
	}
	public final String getLongReadsInput() {return longReadsInput.get();}
	public StringProperty longReadsInputProperty() {return longReadsInput;}
	
	public final void setShortReadsInputFormat(String srInputFormat) {shortReadsInputFormat.set(srInputFormat);}
	public final String getShortReadsInputFormat() {return shortReadsInputFormat.get();}
	public StringProperty shortReadsInputFormatProperty() {return shortReadsInputFormat;}
	
	public final void setLongReadsInputFormat(String lrInputFormat) {longReadsInputFormat.set(lrInputFormat);}
	public final String getLongReadsInputFormat() {return longReadsInputFormat.get();}
	public StringProperty longReadsInputFormatProperty() {return longReadsInputFormat;}
	
	public synchronized void setStopSignal(boolean stop) {this.stop=stop;}
	public synchronized boolean getStopSignal() {return stop;}
	//===============================================================================================//
	
	//Operational variables
//	final BDGraph origGraph;
	public BDGraph simGraph; //original and simplified graph should be separated, no???
	public GraphWatcher observer;
	
	public HybridAssembler(){
//		origGraph=new BDGraph("batch");
		simGraph=new BDGraph("real");
//		rtComponents = new ConnectedComponents();
		simGraph.setAttribute("ui.quality");
		simGraph.setAttribute("ui.antialias");
		
		overwrite = new SimpleBooleanProperty(true);
	    prefix = new SimpleStringProperty("/tmp/");
	    aligner = new SimpleStringProperty("");
		alignerPath = new SimpleStringProperty(); 
		alignerOpt = new SimpleStringProperty();
		errorLog = new SimpleStringProperty();	
		
		shortReadsInput = new SimpleStringProperty(""); 
		binReadsInput = new SimpleStringProperty(""); 
		longReadsInput = new SimpleStringProperty("");
		shortReadsInputFormat = new SimpleStringProperty(); 
		longReadsInputFormat = new SimpleStringProperty();
		
		//set all binding options here...
        shortReadsInput.addListener((observable, oldValue, newValue) -> 
			{
				String fn = ((String)observable.getValue()).toLowerCase();
				if(	fn.endsWith(".fastg")) 
					setShortReadsInputFormat("fastg");
				else if(fn.endsWith(".gfa"))
					setShortReadsInputFormat("gfa");
			}	 

        );
        
        shortReadsInputFormat.addListener((observable, oldValue, newValue) -> 
			{
				if(!getShortReadsInput().toLowerCase().endsWith(newValue))
					setShortReadsInput("");
			}	 

        );
		
        longReadsInput.addListener( (observable, oldValue, newValue) -> 
    		{
				String fn = ((String)observable.getValue()).toLowerCase();
				if(	fn.endsWith(".fasta") || fn.endsWith(".fa") || fn.endsWith("fna")
					|| fn.endsWith(".fastq") || fn.endsWith(".fq")
					|| fn.endsWith(".fasta.gz") || fn.endsWith(".fa.gz") || fn.endsWith("fna.gz")
					|| fn.endsWith(".fastq.gz") || fn.endsWith(".fq.gz") 
					) 
					setLongReadsInputFormat("fasta/fastq");
				else if(fn.endsWith(".sam") || fn.endsWith(".bam")) 
					setLongReadsInputFormat("sam/bam");		
    		}	 
        );
        
        longReadsInputFormat.addListener((observable, oldValue, newValue) -> 
			{
				String oldFile=getLongReadsInput().toLowerCase();
				if(	newValue.equals("fasta/fastq") && !oldFile.endsWith(".fasta") && !oldFile.endsWith(".fa") && !oldFile.endsWith("fna")
							 && !oldFile.endsWith(".fastq") && !oldFile.endsWith(".fq")
							 && !oldFile.endsWith(".fasta.gz") && !oldFile.endsWith(".fa.gz") && !oldFile.endsWith("fna.gz")
							 && !oldFile.endsWith(".fastq.gz") && !oldFile.endsWith(".fq.gz") 
							) 
					setLongReadsInput("");
						
				if(newValue.equals("sam/bam") && !oldFile.endsWith(".sam") && !oldFile.endsWith(".bam"))
					setLongReadsInput("");
			}	 

        );
        
        aligner.addListener( (observable, oldValue, newValue) ->
        	{
				String aligner=(String)observable.getValue();
				if(aligner.toLowerCase().equals("minimap2"))
					setAlignerOpts("-t4 -x map-ont -k15 -w5");
				else if (aligner.toLowerCase().equals("bwa"))
					setAlignerOpts("-t4 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y");			
			}	 

        );
        

	}
		
	
	//Indexing reference, prepare for alignment...
	public boolean prepareLongReadsProcess(){
		//TODO: check all the fields filled up correctly...
		
		
		try{
			System.setProperty("usr.dir", getPrefix());
		}
		catch(NullPointerException | IllegalArgumentException | SecurityException exception ){
			LOG.error("Fail to set working directory usr.dir to {}", getPrefix());
			return false;
		}
		
		//if long reads data not given in SAM/BAM, need to invoke minimap2
        if(getLongReadsInputFormat().toLowerCase().startsWith("fast")) {
        	if(getAligner().equals("minimap2")) { 	
				File indexFile=new File(getPrefix()+"/assembly_graph.mmi");
				if(getOverwrite() || !indexFile.exists()) {						
					try{
						simGraph.outputFASTA(getPrefix()+"/assembly_graph.fasta");
						if(!checkMinimap2()) {
								LOG.error("Dependancy check failed! Please config to the right version of minimap2!");
								return false;
						}
						ProcessBuilder pb = new ProcessBuilder(getFullPathOfAligner(), getAlignerOpts(),"-d", getPrefix()+"/assembly_graph.mmi",prefix+"/assembly_graph.fasta");
						Process indexProcess =  pb.start();
						indexProcess.waitFor();
						
					}catch (IOException | InterruptedException e){
						LOG.error("Issue when indexing with minimap2: \n" + e.getMessage());
						return false;
					}
				}
        	}else if(getAligner().equals("bwa")) {
				File indexFile=new File(getPrefix()+"/assembly_graph.fasta.bwt");
				if(getOverwrite() || !indexFile.exists()) {						
					try{
						simGraph.outputFASTA(getPrefix()+"/assembly_graph.fasta");
						if(!checkBWA()) {
								LOG.error("Dependancy check failed! Please config to the right version of bwa!");
								return false;
						}
						ProcessBuilder pb = new ProcessBuilder(getFullPathOfAligner(),"index", getPrefix()+"/assembly_graph.fasta");
						Process indexProcess =  pb.start();
						indexProcess.waitFor();
						
					}catch (IOException | InterruptedException e){
						LOG.error("Issue when indexing with bwa: \n" + e.getMessage());
						return false;
					}
				}
        	}else {
        		LOG.error("Unknown aligner!");
        		return false;
        	}
			
			
        }
        return true;
	}
	//Loading the graph, doing preprocessing
	//binning, ...
	public boolean prepareShortReadsProcess(boolean useSPAdesPaths) {
		//TODO: check all the fields filled up correctly...

		//try to read input file
		try {
			if(getShortReadsInputFormat().toLowerCase().equals("gfa")) 
				GraphUtil.loadFromGFA(getShortReadsInput(), getBinReadsInput(), simGraph, useSPAdesPaths);
			else if(getShortReadsInputFormat().toLowerCase().equals("fastg"))
				GraphUtil.loadFromFASTG(getShortReadsInput(), getBinReadsInput(), simGraph, useSPAdesPaths);
			else 				
				throw new IOException("assembly graph file must have .gfa or .fastg extension!");
			
		}catch(IOException e) {
			System.err.println("Issue when loading pre-assembly: \n" + e.getMessage());
			return false;
		}
		
		
		simGraph.updateStats();
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

		if (getLongReadsInputFormat().endsWith("am")){//bam or sam
			if ("-".equals(getLongReadsInput()))
				reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
			else
				reader = SamReaderFactory.makeDefault().open(new File(getLongReadsInput()));	
		}else{
			LOG.info("Starting alignment by {} at {}", getAligner(), new Date());
			ProcessBuilder pb = null;
			if ("-".equals(getLongReadsInput())){
				if(getAligner().equals("minimap2"))
					pb = new ProcessBuilder(getFullPathOfAligner(), 
							"-a",
							getAlignerOpts(),
							"-K",
							"20000",
							getPrefix()+"/assembly_graph.mmi",
							"-"
							).
							redirectInput(Redirect.INHERIT);
				else if(getAligner().equals("bwa"))
					pb = new ProcessBuilder(getFullPathOfAligner(), 
							"mem",
							getAlignerOpts(),
							"-K",
							"20000",
							getPrefix()+"/assembly_graph.fasta",
							"-"
							).
							redirectInput(Redirect.INHERIT);
			}else{
				if(getAligner().equals("minimap2"))
					pb = new ProcessBuilder(getFullPathOfAligner(), 
							"-a",
							getAlignerOpts(),
							"-K",
							"20000",
							getPrefix()+"/assembly_graph.mmi",
							getLongReadsInput()
							);
				else if(getAligner().equals("bwa"))
					pb = new ProcessBuilder(getFullPathOfAligner(), 
							"mem",
							getAlignerOpts(),
							"-K",
							"20000",
							getPrefix()+"/assembly_graph.fasta",
							getLongReadsInput()
							);
			}

//			alignmentProcess  = pb.redirectError(ProcessBuilder.Redirect.to(new File("/dev/null"))).start();
			alignmentProcess  = pb.redirectError(ProcessBuilder.Redirect.to(new File(getPrefix()+"/alignment.log"))).start();

			LOG.info("{} started!", getAligner());			

			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(alignmentProcess.getInputStream()));

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
			
			if (simGraph.getNode(refID)==null) {
				LOG.error("Node {} not found from the graph!", refID);
				continue;
			}
			Alignment myRec = new Alignment(rec, (BDNode) simGraph.getNode(refID)); 

			//////////////////////////////////////////////////////////////////
			
			if (!readID.equals("") && !readID.equals(myRec.readID)) {	
				synchronized(simGraph) {
					List<BDPath> paths=simGraph.uniqueBridgesFinding(nnpRead, samList);
					if(paths!=null){	
						for(BDPath path:paths) 
						{
							//path here is already unique! (2 unique ending nodes)
					    	if(simGraph.reduceUniquePath(path)) {
					    		observer.update(false);					    		
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

		if (alignmentProcess != null){
			alignmentProcess.destroy();
		}	

	}
	
	public void postProcessGraph() throws IOException{
		//Take the current best path among the candidate of a bridge and connect the bridge(greedy)
		for(GoInBetweenBridge brg:simGraph.getUnsolvedBridges()){
			System.out.printf("Last attempt on incomplete bridge %s : anchors=%d \n %s \n", brg.getEndingsID(), brg.getNumberOfAnchors(), brg.getAllPossiblePaths());
//			if(brg.getCompletionLevel()<=2) {//examine bridge with completion level = 2 that unable to connected
//				brg.steps.connectBridgeSteps(true);
//			}
			
			if(brg.getCompletionLevel()>=3) 
				simGraph.chopPathAtAnchors(brg.getBestPath(brg.pBridge.getNode0(),brg.pBridge.getNode1())).stream().forEach(p->simGraph.reduceUniquePath(p));
			else{
				brg.scanForAnEnd(true);
				//selective connecting
				brg.steps.connectBridgeSteps(true);
				//return appropriate path
				if(brg.segments!=null)
					simGraph.chopPathAtAnchors(brg.getBestPath(brg.steps.start.getNode(),brg.steps.end.getNode())).stream().forEach(p->simGraph.reduceUniquePath(p));
				else
					System.out.printf("Last attempt failed \n");
			}


		}
		
        //update for the last time
        observer.update(true);
		observer.outputFASTA(getPrefix()+"npgraph_assembly.fasta");

	}
	


    
    @SuppressWarnings("resource")
	public static void promptEnterKey(){
    	   System.out.println("Press \"ENTER\" to continue...");
    	   Scanner scanner = new Scanner(System.in);
    	   scanner.nextLine();
    	}
    
    
    public boolean checkMinimap2() {    		
		ProcessBuilder pb = new ProcessBuilder(getFullPathOfAligner(),"-V").redirectErrorStream(true);
		Process process;
		try {
			process = pb.start();
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
				LOG.error(getFullPathOfAligner() + " is not the right path to minimap2!");
				return false;
			}else{
				System.out.println("minimap version: " + version);
				if (version.compareTo("2.0") < 0){
					LOG.error("Require minimap version 2 or above!");
					return false;
				}
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			LOG.error("Error running: " + getFullPathOfAligner() + "\n" + e.getMessage());
			return false;
		}
		
		return true;
			
    }
    
    public boolean checkBWA() {    		
		try{
			ProcessBuilder pb = new ProcessBuilder(getFullPathOfAligner()).redirectErrorStream(true);
			Process process =  pb.start();
			BufferedReader bf = SequenceReader.openInputStream(process.getInputStream());


			String line;
			String version = "";
			Pattern versionPattern = Pattern.compile("^Version:\\s(\\d+\\.\\d+\\.\\d+).*");
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
				LOG.error(getFullPathOfAligner() + " is not the right path to BWA!");
				return false;
			}else{
				LOG.info("bwa version: " + version);
				if (version.compareTo("0.7.11") < 0){
					LOG.error(" Require bwa of 0.7.11 or above");
					return false;
				}
			}

		}catch (IOException e){
			LOG.error("Error running: " + getFullPathOfAligner() + "\n" + e.getMessage());
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
