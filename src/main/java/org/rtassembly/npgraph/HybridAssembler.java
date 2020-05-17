package org.rtassembly.npgraph;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.lang.ProcessBuilder.Redirect;
import java.lang.invoke.MethodHandles;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.PAFRecord;
import japsa.seq.Sequence;
import javafx.beans.property.StringProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.BooleanProperty;
import javafx.beans.property.SimpleBooleanProperty;
public class HybridAssembler {
    private static final Logger logger = LogManager.getLogger(MethodHandles.lookup().lookupClass());
	//setting parameter for the GUI
    private boolean ready=false;
    private BooleanProperty overwrite;
    private StringProperty 	prefix;
	public InputData input;
	Process alignmentProcess = null;
	private boolean stop=false;
	public int currentReadCount = 0;
	public long currentBaseCount = 0;	
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
	
	public synchronized void setStopSignal(boolean stop) {this.stop=stop;}
	public synchronized boolean getStopSignal() {return stop;}
	//===============================================================================================//
	
	//Operational variables
	public volatile BDGraph simGraph; //original and simplified graph should be separated, no???
	public RealtimeGraphWatcher observer;
	
	public HybridAssembler(){
		simGraph=new BDGraph("real");
		simGraph.setAttribute("ui.quality");
		simGraph.setAttribute("ui.antialias");
		
		overwrite = new SimpleBooleanProperty(true);
	    prefix = new SimpleStringProperty(System.getProperty("java.io.tmpdir"));
	    
	    input = new InputData();       
	}	
	
	//Indexing reference, prepare for alignment...
	public boolean prepareLongReadsProcess(){
		if(!input.validateLongReadInput())
			return false;
		
		if(!GraphUtil.checkFolder(getPrefix()))
			return false;
		try{
			System.setProperty("usr.dir", getPrefix());			
		}
		catch(NullPointerException | IllegalArgumentException | SecurityException exception ){
			logger.error("Fail to set working directory usr.dir to " + getPrefix());
			return false;
		}
		
		//create temporary folder to store bridging reads
		File tmpFolder=new File(getPrefix()+File.separator+"npGraph_tmp");
		if(tmpFolder.exists()){
		  try {
			Files.walk(tmpFolder.toPath())
			    .sorted(Comparator.reverseOrder())
			    .map(Path::toFile)
			    .forEach(File::delete);
			} catch (IOException e) {
				logger.error("Cannot remove existed temporary folder: {}", e);
			}
  		}
		
		if(!tmpFolder.mkdir()){
			logger.error("Cannot set temporary folder " + tmpFolder.getAbsolutePath());
			return false;
		}else{
			AlignedRead.tmpFolder=tmpFolder.getAbsolutePath();
		}
		
		//if long reads data not given in SAM/BAM, need to invoke minimap2
        if(input.getLongReadsInputFormat().contains("fast")) {
        	File indexFile=null;
        	ArrayList<String> idxCmd = new ArrayList<>();
        	idxCmd.add(input.getAligner());
        	if(input.getAligner().endsWith("minimap2")) { 	
				indexFile=new File(getPrefix()+File.separator+"assembly_graph.mmi");												
				if(!GraphUtil.checkMinimap2(input.getAligner())) 
						return false;
				idxCmd.addAll(Arrays.asList(input.getAlignerOpts().split("\\s")));
				idxCmd.add("-d");
				idxCmd.add(getPrefix()+File.separator+"assembly_graph.mmi");
														
        	}else if(input.getAligner().endsWith("bwa")) {
				indexFile=new File(getPrefix()+File.separator+"assembly_graph.fasta.bwt");
				if(!GraphUtil.checkBWA(input.getAligner())) 
						return false;
				idxCmd.add("index");

        	}else {
        		logger.error("Invalid aligner! Set to BWA or minimap2 please!");
        		return false;
        	}
			idxCmd.add(getPrefix()+File.separator+"assembly_graph.fasta");

			if(getOverwrite() || !indexFile.exists()) {						
				try{
					simGraph.outputFASTA(getPrefix()+File.separator+"assembly_graph.fasta");
					
					ProcessBuilder pb = new ProcessBuilder(idxCmd);
					Process indexProcess =  pb.start();
					indexProcess.waitFor();
					
				}catch (IOException | InterruptedException e){
					logger.error("Issue when indexing the pre-assemblies:\n{}", e);
					return false;
				}
			}
			
        }
        
        //check consensus tool
    	if(!GraphUtil.checkMSA(input.getMSA())){
    		logger.warn("WARNING: MSA tool \"{}\" not found! Set to none", input.getMSA());
    		input.setMSA("none");
    	}else
    		logger.info("MSA for consensus calling is set to {}", input.getMSA());
    	
		simGraph.consensus.setConsensusMSA(input.getMSA());
    	
        return true;
	}
	//Loading the graph, doing preprocessing
	//binning, ...
	public boolean prepareShortReadsProcess() {
	
		//try to read input file
		try {
			if(input.getShortReadsInputFormat().toLowerCase().equals("gfa")) 
				GraphUtil.loadFromGFA(input.getShortReadsInput(), input.getBinReadsInput(), simGraph, input.getUseSPAdesPath());
			else if(input.getShortReadsInputFormat().toLowerCase().equals("fastg"))
				GraphUtil.loadFromFASTG(input.getShortReadsInput(), input.getBinReadsInput(), simGraph, input.getUseSPAdesPath());
			else 				
				throw new IOException("Assembly graph file must have .gfa or .fastg extension!");
			
		}catch(IOException e) {
			logger.error("Issue when loading pre-assembly:\n{}", e);
			return false;
		}
		
		simGraph.updateStats();
		observer = new RealtimeGraphWatcher(this);
		
		observer.setReadPeriod(RealtimeGraphWatcher.R_INTERVAL);
		observer.setTimePeriod(RealtimeGraphWatcher.T_INTERVAL * 1000);
		//re-estimate timely report based on graph complexity
//		int timeInterval=(int) (Math.round(Math.log10(simGraph.getNodeCount()))-1); //estimated interval time based on graph complexity
//		timeInterval=(timeInterval>1?timeInterval:1)*10;
		logger.info("The results will be reported every {} seconds and {} reads pass", RealtimeGraphWatcher.T_INTERVAL, RealtimeGraphWatcher.R_INTERVAL);
		return true;
	}

	/**
	 * Alignment (if needed) running to update the graph in parallel with an observer thread 
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void assembly() 
			throws IOException, InterruptedException{
		
		if(input.getLongReadsInput().isEmpty()) {
//			LOG.info("Scaffolding is ignored due to lack of long-read input!");
			return;
		}else
			logger.info("Scaffolding ready at ", new Date());

		Thread thread = new Thread(observer);
		thread.start();	
		/*PAF*/
		if (input.getLongReadsInputFormat().endsWith("paf")){//bam or sam
			InputStreamReader inputStream=null;
			if ("-".equals(input.getLongReadsInput()))
				inputStream = new InputStreamReader(System.in);
			else
				inputStream = new FileReader(input.getLongReadsInput());	
			
			try(BufferedReader reader=new BufferedReader(inputStream)){
				String readID = "";
				Sequence read = null;
				ArrayList<Alignment> hits =  new ArrayList<Alignment>();// alignment record of the same read;	
				PAFRecord curRecord=null;
				
				String line;
				while ((line=reader.readLine()) != null) {
					if(getStopSignal())
						break;
					
					try {
						curRecord = new PAFRecord(line);
					}catch(Exception e) {
						logger.error("Error reading PAF record: \n {}", e);
		//				continue;
						break;
					}
					
					if (curRecord.qual < Alignment.MIN_QUAL){		
						logger.debug("Ignore low-quality map record!");
						if (!readID.equals(curRecord.qname)){
							update(read, hits);
							hits = new ArrayList<Alignment>();
							readID = curRecord.qname;
							read = GraphUtil.getNSequence(curRecord.qname, curRecord.qlen);//there is no read data from PAF, so just fake one!

						}
						continue;		
					}
					
					String refName = curRecord.tname;
					String refID = refName.split("_").length > 1 ? refName.split("_")[1]:refName;
					
					//check if this node still in. FIXME: do not remove nodes for metagenomics' graph?
					if (simGraph.getNode(refID)==null) {
						logger.debug("Ignore record with reference {} not found (removed) from the graph!", refID);
						if (!readID.equals(curRecord.qname)){
							update(read, hits);
							hits = new ArrayList<Alignment>();
							readID = curRecord.qname;
							read = GraphUtil.getNSequence(curRecord.qname, curRecord.qlen);//there is no read data from PAF, so just fake one!

						}
						continue;
					}
					Alignment curAlignment = new Alignment(curRecord, (BDNode) simGraph.getNode(refID)); 
		
					//////////////////////////////////////////////////////////////////
					
					if (!readID.equals("") && !readID.equals(curRecord.qname)) {	
						update(read, hits);
						hits = new ArrayList<Alignment>();
						readID = curRecord.qname;
						read = GraphUtil.getNSequence(curRecord.qname, curRecord.qlen);//there is no read data from PAF, so just fake one!

					}	
					hits.add(curAlignment); 
				}// while
	
			}
		/*SAM/BAM/FASTA/FASTQ*/
		}else {
			SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
			SamReader reader = null;
	
			if (input.getLongReadsInputFormat().endsWith("am")){//bam or sam
				if ("-".equals(input.getLongReadsInput()))
					reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
				else
					reader = SamReaderFactory.makeDefault().open(new File(input.getLongReadsInput()));	
			}else{
				logger.info("Starting alignment by {} at {}", input.getAligner(), new Date());
				ProcessBuilder pb = null;
				List<String> command = new ArrayList<>();
				command.add(input.getAligner());
				if(input.getAligner().endsWith("minimap2")) {
					command.add("-a");
					command.addAll(Arrays.asList(input.getAlignerOpts().split("\\s")));
					command.add("-K20000");
					command.add(getPrefix()+File.separator+"assembly_graph.mmi");
					command.add(input.getLongReadsInput());
				}
				else if(input.getAligner().endsWith("bwa")) {
					command.add("mem");
					command.addAll(Arrays.asList(input.getAlignerOpts().split("\\s")));
					command.add("-K20000");
					command.add(getPrefix()+File.separator+"assembly_graph.fasta");
					command.add(input.getLongReadsInput());
				}
				
				if ("-".equals(input.getLongReadsInput())){
					pb = new ProcessBuilder(command).redirectInput(Redirect.INHERIT);
				}else{
					pb = new ProcessBuilder(command);
				}
	
				alignmentProcess  = pb.redirectError(ProcessBuilder.Redirect.to(new File(getPrefix()+File.separator+"alignment.log"))).start();
	
				logger.info("{} started!", input.getAligner());			
	
				reader = SamReaderFactory.makeDefault().open(SamInputResource.of(alignmentProcess.getInputStream()));
	
			}
			SAMRecordIterator iter = reader.iterator();
	
			String readID = "";
			Sequence read = null;
			ArrayList<Alignment> hits =  new ArrayList<Alignment>();// alignment record of the same read;	
			SAMRecord curRecord=null;
						
			while (iter.hasNext()) {
				if(getStopSignal())
					break;
				
				try {
					curRecord = iter.next();
				}catch(Exception e) {
					logger.error("Error SAM record: \n {}", e);
					break;
				}
				
				if (curRecord.getReadUnmappedFlag() || curRecord.getMappingQuality() < Alignment.MIN_QUAL){		
					logger.debug("Ignore one unmapped or low-quality map record!");
					if (!readID.equals(curRecord.getReadName())){
						update(read, hits);
						hits = new ArrayList<Alignment>();
						read = GraphUtil.getQueryReadFromSAMRecord(curRecord);
						readID = curRecord.getReadName();
					}
					continue;		
				}
				
				String refName = curRecord.getReferenceName();
				String refID = refName.split("_").length > 1 ? refName.split("_")[1]:refName;
				
				//check if this node still in. FIXME: do not remove nodes for metagenomics' graph?
				if (simGraph.getNode(refID)==null) {
					logger.debug("Ignore record with reference {} not found (removed) from the graph!", refID);
					if (!readID.equals(curRecord.getReadName())){
						update(read, hits);
						hits = new ArrayList<Alignment>();
						read = GraphUtil.getQueryReadFromSAMRecord(curRecord);
						readID = curRecord.getReadName();
					}
					continue;
				}
				Alignment curAlignment = new Alignment(curRecord, (BDNode) simGraph.getNode(refID)); 
	
				//////////////////////////////////////////////////////////////////
				
				if (!readID.equals("") && !readID.equals(curRecord.getReadName())) {	
					update(read, hits);
					hits = new ArrayList<Alignment>();
					read = GraphUtil.getQueryReadFromSAMRecord(curRecord);
					readID = curRecord.getReadName();
	
				}	
				hits.add(curAlignment); 
			}// while
			iter.close();
			reader.close();
		}
		
		
		observer.stopWaiting();
		thread.join();
		terminateAlignmentProcess();
	}
	
	/**
	 * Version 2 using PAF format instead of SAM/BAM
	 * @throws IOException
	 * @throws InterruptedException
	 */
	public void assemblyPAF() 
			throws IOException, InterruptedException{
	
		if(input.getLongReadsInput().isEmpty()) {
			logger.info("Scaffolding is ignored due to lack of long-read input!");
			return;
		}else
			logger.info("Scaffolding ready at {}", new Date());


		InputStreamReader inputStream = null;

		if (input.getLongReadsInputFormat().endsWith("paf")){//bam or sam
			if ("-".equals(input.getLongReadsInput()))
				inputStream = new InputStreamReader(System.in);
			else
				inputStream = new FileReader(input.getLongReadsInput());	
		}else if(input.getLongReadsInputFormat().contains("fast")){
			logger.info("Alignment command: {} {} start at {}", input.getAligner(), input.getAlignerOpts(), new Date());
			ProcessBuilder pb = null;
			List<String> command = new ArrayList<>();
			command.add(input.getAligner());
//			command.add("-a");
			command.addAll(Arrays.asList(input.getAlignerOpts().split("\\s")));
			command.add("-K20000");
			command.add(getPrefix()+File.separator+"assembly_graph.mmi");
			command.add(input.getLongReadsInput());
			

			
			if ("-".equals(input.getLongReadsInput())){
				pb = new ProcessBuilder(command).redirectInput(Redirect.INHERIT);
			}else{
				pb = new ProcessBuilder(command);
			}

			alignmentProcess  = pb.redirectError(ProcessBuilder.Redirect.to(new File(getPrefix()+File.separator+"alignment.log"))).start();

			logger.info("{} started!", input.getAligner());			

			inputStream = new InputStreamReader(alignmentProcess.getInputStream());

		}
		
		try(BufferedReader reader=new BufferedReader(inputStream)){
			String readID = "";
			Sequence read = null;
			ArrayList<Alignment> hits =  new ArrayList<Alignment>();// alignment record of the same read;	
			PAFRecord curRecord=null;
			
			Thread thread = new Thread(observer);
			thread.start();	
			String line;
			while ((line=reader.readLine()) != null) {
				if(getStopSignal())
					break;
				
				try {
					curRecord = new PAFRecord(line);
				}catch(Exception e) {
					logger.error("Error reading PAF record: \n {}", e);
	//				continue;
					break;
				}
				
				if (curRecord.qual < Alignment.MIN_QUAL){		
					logger.debug("Ignore low-quality map record!");
					if (!readID.equals(curRecord.qname)){
						update(read, hits);
						hits = new ArrayList<Alignment>();
						readID = curRecord.qname;
						read = GraphUtil.getNSequence(curRecord.qname, curRecord.qlen);//there is no read data from PAF, so just fake one!

					}
					continue;		
				}
				
				String refName = curRecord.tname;
				String refID = refName.split("_").length > 1 ? refName.split("_")[1]:refName;
				
				//check if this node still in. FIXME: do not remove nodes for metagenomics' graph?
				if (simGraph.getNode(refID)==null) {
					logger.debug("Ignore record with reference {} not found (removed) from the graph!", refID);
					if (!readID.equals(curRecord.qname)){
						update(read, hits);
						hits = new ArrayList<Alignment>();
						readID = curRecord.qname;
						read = GraphUtil.getNSequence(curRecord.qname, curRecord.qlen);//there is no read data from PAF, so just fake one!

					}
					continue;
				}
				Alignment curAlignment = new Alignment(curRecord, (BDNode) simGraph.getNode(refID)); 
	
				//////////////////////////////////////////////////////////////////
				
				if (!readID.equals("") && !readID.equals(curRecord.qname)) {	
					update(read, hits);
					hits = new ArrayList<Alignment>();
					readID = curRecord.qname;
					read = GraphUtil.getNSequence(curRecord.qname, curRecord.qlen);//there is no read data from PAF, so just fake one!

				}	
				hits.add(curAlignment); 
			}// while
			
			observer.stopWaiting();
			thread.join();
			terminateAlignmentProcess();	
		}

	}
	
	//update when more read alignments coming in
	synchronized void update(Sequence nnpRead, ArrayList<Alignment> alignments){
		if(alignments.isEmpty() || nnpRead==null)
			return;
		
		currentReadCount ++;
		currentBaseCount += nnpRead.length();

		List<BDPath> paths=simGraph.uniqueBridgesFinding(nnpRead, alignments);
		if(paths!=null)
		    paths.stream().forEach(p->simGraph.reduceUniquePath(p));
	} 		
	
	public void terminateAlignmentProcess() {
 		if (alignmentProcess != null){
 			alignmentProcess.destroy();
 		}		
 	}
	
	
	@Deprecated
	//last attempt to connect bridges greedily. Now move to RealtimeGraphWatcher
	public void postProcessGraph() throws IOException{
		System.out.printf("Post-processing the graph by greedy path-finding algorithm. Please wait...\n");
		HashSet<GoInBetweenBridge> 		unsolved=simGraph.getUnsolvedBridges(),
										solved=new HashSet<>();
		while(true){
			boolean changed=false;
			for(GoInBetweenBridge brg:unsolved){
				logger.debug("Last attempt on incomplete bridge {} : anchors={}\n{}",
								brg.getEndingsID(), brg.getNumberOfAnchors(), brg.getAllPossiblePaths());
				//Take the current best path among the candidate of a bridge and connect the bridge(greedy)
				if(brg.getCompletionLevel()>=3){ 
					simGraph.getNewSubPathsToReduce(brg.getBestPath(brg.pBridge.getNode0(),brg.pBridge.getNode1())).stream().forEach(p->simGraph.reduceUniquePath(p));
					solved.add(brg);
					changed=true;
				}else{
					brg.scanForAnEnd(true);	
					changed=brg.steps.connectBridgeSteps(true);
					
					//return appropriate path
					if(changed){
						simGraph.getNewSubPathsToReduce(brg.getBestPath(brg.steps.start.getNode(),brg.steps.end.getNode())).stream().forEach(p->simGraph.reduceUniquePath(p));
						solved.add(brg);
					}
					else
						logger.debug("Last attempt failed!");
				}
	
			}
			if(solved.isEmpty()&&!changed)
				break;
			else{
				unsolved.removeAll(solved);
				solved.clear();
			}
				
		}
        //update for the last time
        observer.update(true);
		System.out.printf("Input stats: read count=%d base count=%d\n", currentReadCount, currentBaseCount);
		
		observer.outputFASTA(getPrefix()+File.separator+"npgraph_assembly.fasta");
		observer.outputJAPSA(getPrefix()+File.separator+"npgraph_assembly.japsa");
		observer.outputAssGFA(getPrefix()+File.separator+"npgraph_assembly.gfa");
		observer.outputOrigGFA(getPrefix()+File.separator+"npgraph_components.gfa");

		
		//delete temporary files
//		File tmpFolder=new File(getPrefix()+File.separator+"npGraph_tmp");
//		if(tmpFolder.exists()){
//		  try {
//			Files.walk(tmpFolder.toPath())
//			    .sorted(Comparator.reverseOrder())
//			    .map(Path::toFile)
//			    .forEach(File::delete);
//			} catch (IOException e) {
//				LOG.info("Cannot remove existed temporary folder {}!", tmpFolder.getPath());
//				e.printStackTrace();
//			}
//  		}
		
	}

    
//	public static void main(String[] argv) throws IOException, InterruptedException{
//		HybridAssembler hbAss = new HybridAssembler();
//		
//		hbAss.input.setShortReadsInput("/home/sonhoanghguyen/Projects/scaffolding/data/spades_3.7/EcK12S-careful/assembly_graph.fastg");
//		hbAss.input.setShortReadsInputFormat("fastg");
//		hbAss.prepareShortReadsProcess();
//		hbAss.input.setLongReadsInput("/home/sonhoanghguyen/Projects/scaffolding/data/spades_3.7/EcK12S-careful/assembly_graph.sam");
//		hbAss.input.setLongReadsInputFormat("sam/bam");
//		hbAss.prepareLongReadsProcess();
//		hbAss.assembly();
//
//	}
	
}
