package org.rtassembly;
import java.io.File;
import java.io.IOException;
import org.rtassembly.gui.NPGraphFX;
import org.rtassembly.npgraph.Alignment;
import org.rtassembly.npgraph.BDGraph;
import org.rtassembly.npgraph.HybridAssembler;

import japsa.util.CommandLine;
import javafx.application.Application;



public class NPGraphCmd extends CommandLine{
	public NPGraphCmd(){
		super();

		addString("si", null, "Name of the short-read assembly file.",true);
		addString("sf", "gfa", "Format of the assembly input file. Accepted format are FASTG, GFA", true);
		addString("li", "-", "Name of the long-read data input file, - for stdin.", true);
		addString("lf", "fastq", "Format of the long-read data input file. This may be FASTQ/FASTA (MinION reads) or SAM/BAM (aligned with the assembly graph already)", true);
		addString("output", "npgraph", "Name of the output folder.", true);
				
		addString("alg","","Absolute path to the folder containing binary minimap2",true);

		addString("algPath","","Absolute path to the binary aligner file");
		addString("algOpt", "", "Settings used by aligner to align long reads to the contigs");
		
		addBoolean("overwrite", false, "Whether to overwrite or reuse the intermediate file");
		addInt("qual", 1, "Minimum quality of alignment to considered");
		addInt("dfs", 15, "Number of DFS steps to search");

		addBoolean("gui", false, "Whether using GUI or not.");
		
		addStdHelp();
	}
	
	@SuppressWarnings("restriction")
	public static void main(String[] args) throws IOException{
		CommandLine cmdLine = new NPGraphCmd();		
		args = cmdLine.stdParseLine(args);

		/***********************************************************************/
		String 	shortReadsInput = cmdLine.getStringVal("si"),
				shortReadsInputFormat = cmdLine.getStringVal("sf"),
				longReadsInput = cmdLine.getStringVal("li"),
				longReadsInputFormat = cmdLine.getStringVal("lf"),
				outputDir = cmdLine.getStringVal("output"),
				alg=cmdLine.getStringVal("alg"),
				algPath = cmdLine.getStringVal("algPath"),
				algOpt = cmdLine.getStringVal("algOpt");
		boolean overwrite = cmdLine.getBooleanVal("overwrite"),
				gui = cmdLine.getBooleanVal("gui");
			
		Alignment.MIN_QUAL = cmdLine.getIntVal("qual");
		BDGraph.S_LIMIT=cmdLine.getIntVal("dfs");
		//Default output dir 
		if(outputDir == null) {
			outputDir = new File(shortReadsInput).getAbsoluteFile().getParent();
		}
		File outDir = new File(outputDir);
		if(!outDir.exists())
			outDir.mkdirs();
			
		//1. Create an assembler object with appropriate file loader
		HybridAssembler hbAss = new HybridAssembler();
		hbAss.setShortReadsInput(shortReadsInput);
		hbAss.setShortReadsInputFormat(shortReadsInputFormat);
		hbAss.setLongReadsInput(longReadsInput);
		hbAss.setLongReadsInputFormat(longReadsInputFormat);
		
		hbAss.setPrefix(outputDir);
		
		hbAss.setAlignerPath(algPath);
		hbAss.setAlignerOpts(algOpt);
		hbAss.setAligner(alg);
		hbAss.setOverwrite(overwrite);
		        
		//4. Call the assembly function or invoke GUI to do so
        if(gui) {
			NPGraphFX.setAssembler(hbAss);
			Application.launch(NPGraphFX.class,args);
        }else {
	        
			try {
				if(hbAss.prepareShortReadsProcess(false) &&	hbAss.prepareLongReadsProcess()) {
					hbAss.assembly();
					hbAss.postProcessGraph();
				}
				else{
					System.err.println("Error with pre-processing step: \n" + hbAss.getErrorLog());
					System.exit(1);
				}
					
			} catch (InterruptedException|IOException e) {
				System.err.println("Issue when assembly: \n" + e.getMessage());
				e.printStackTrace();
				System.exit(1);
			}
        }
		
	}
	
}
