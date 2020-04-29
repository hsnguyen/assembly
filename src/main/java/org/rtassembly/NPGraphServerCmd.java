package org.rtassembly;
import java.io.File;
import java.io.IOException;
import java.lang.invoke.MethodHandles;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.rtassembly.gui.NPGraphFX;
import org.rtassembly.npgraph.Alignment;
import org.rtassembly.npgraph.BDGraph;
import org.rtassembly.npgraph.RealtimeGraphWatcher;
import org.rtassembly.npgraph.HybridAssembler;
import org.rtassembly.npgraph.SimpleBinner;
import org.rtassembly.npgraph.grpc.AssemblyGuideServer;

import japsa.util.CommandLine;
import javafx.application.Application;


public class NPGraphServerCmd extends CommandLine{
    private static final Logger logger = LogManager.getLogger(MethodHandles.lookup().lookupClass());
	public NPGraphServerCmd(){
		super();
		//Server setttings
		addInt("port",2105, "Port number for which the server will listen to.");
		
		//Input settings
		addString("si", "", "Name of the short-read assembly file.");
		addString("sf", "", "Format of the assembly input file. Accepted format are FASTG, GFA");
		addString("output", "/tmp/", "Output folder for temporary files and the final assembly npgraph_assembly.fasta");
				
		addString("sb", "", "Name of the metaBAT file for binning information (experimental).");

		addBoolean("overwrite", true, "Whether to overwrite or reuse the intermediate file");
		addBoolean("sp", false, "Whether to use SPAdes contigs.paths for bridging.");
		
		//Algorithm-wise
		addInt("qual", 10, "Minimum quality of alignment to considered");
		addInt("mincov", 3, "Minimum number of reads spanning a confident bridge");
		addInt("maxcov", 20, "Cut-off number of reads spanning a confident bridge");
		addInt("depth", 300, "Limit depth for searching path between 2 neighbors");
		addInt("anchor", 1000, "Lowerbound of an anchor contig's length.");
		addInt("unique", 10000, "Lowerbound of a unique contig's length.");
		addInt("time", 10, "Time interval to considered for real-time assembly.");
		addInt("read", 50, "Read interval to considered for real-time assembly.");
		
		addBoolean("gui", false, "Whether using GUI or not.");
		addBoolean("keep", false, "Whether to keep extremely-low-coveraged contigs.");
//		addBoolean("verbose", false, "For debugging.");
		addStdHelp();
	}
	
	public static void main(String[] args) throws IOException{
		CommandLine cmdLine = new NPGraphServerCmd();		
		args = cmdLine.stdParseLine(args);

		/***********************************************************************/
		int port = cmdLine.getIntVal("port");
		String 	shortReadsInput = cmdLine.getStringVal("si"),
				shortReadsInputFormat = cmdLine.getStringVal("sf"),
				outputDir = cmdLine.getStringVal("output"),
				shortReadsBinInput = cmdLine.getStringVal("sb");
		boolean overwrite = cmdLine.getBooleanVal("overwrite"),
				spaths = cmdLine.getBooleanVal("sp"),
				gui = cmdLine.getBooleanVal("gui");
			
		Alignment.MIN_QUAL = cmdLine.getIntVal("qual");
		SimpleBinner.ANCHOR_CTG_LEN=cmdLine.getIntVal("anchor");
		SimpleBinner.UNIQUE_CTG_LEN=cmdLine.getIntVal("unique");
		RealtimeGraphWatcher.KEEP=cmdLine.getBooleanVal("keep");
		BDGraph.MIN_SUPPORT=cmdLine.getIntVal("min");
		BDGraph.GOOD_SUPPORT=cmdLine.getIntVal("max");
		BDGraph.S_LIMIT=cmdLine.getIntVal("depth");
		RealtimeGraphWatcher.R_INTERVAL=cmdLine.getIntVal("read");
		RealtimeGraphWatcher.T_INTERVAL=cmdLine.getIntVal("time");

		//Default output dir 
		if(outputDir == null) {
			outputDir = new File(shortReadsInput).getAbsoluteFile().getParent();
		}
		File outDir = new File(outputDir);
		if(!outDir.exists())
			outDir.mkdirs();
			
		//1. Create an assembler object with appropriate file loader
		HybridAssembler hbAss = new HybridAssembler();
		if(shortReadsInput!=null && !shortReadsInput.isEmpty())
			hbAss.setShortReadsInput(shortReadsInput);
		if(shortReadsInputFormat!=null && !shortReadsInputFormat.isEmpty())
			hbAss.setShortReadsInputFormat(shortReadsInputFormat);
		
		hbAss.setLongReadsInputFormat("paf"); //just to pass prepareLongReadProcess()
		
		hbAss.setPrefix(outputDir);
		if(shortReadsBinInput!=null && !shortReadsBinInput.isEmpty())
			hbAss.setBinReadsInput(shortReadsBinInput);
		
	
		hbAss.setOverwrite(overwrite);
		hbAss.setUseSPAdesPath(spaths);
		        
		//4. Call the assembly function or invoke GUI to do so
        if(gui) {
//			NPGraphFX.setAssembler(hbAss);
//			Application.launch(NPGraphFX.class,args);
        	System.out.println("Sorry it's easy but not implemented yet!");
        	System.exit(1);
        }else if(shortReadsInput.isEmpty()) {
			System.out.println(cmdLine.usageString());			
			System.exit(-1);
        }else {
			try {
				if(hbAss.prepareShortReadsProcess() &&	hbAss.prepareLongReadsProcess()) {
					AssemblyGuideServer assServer = new AssemblyGuideServer(port, hbAss);
					assServer.start();
				}
				else{
					logger.error("Error with pre-processing step: \n{}", hbAss.getCheckLog());
					System.exit(1);
				}
					
			} catch (Exception e) {
				logger.error("Error {}", e);
				System.exit(1);
			}
        }
		
	}
	
}
