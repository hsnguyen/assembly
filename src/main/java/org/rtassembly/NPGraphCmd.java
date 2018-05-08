package org.rtassembly;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.graphstream.ui.view.Viewer;
import org.rtassembly.gui.NPGraphFX;
import org.rtassembly.npgraph.Alignment;
import org.rtassembly.npgraph.BidirectedGraph;
import org.rtassembly.npgraph.GraphExplore;
import org.rtassembly.npgraph.HybridAssembler;

import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import javafx.application.Application;



public class NPGraphCmd extends CommandLine{
	public NPGraphCmd(){
		super();

		addString("si", null, "Name of the short-read assembly file.",true);
		addString("sf", "gfa", "Format of the assembly input file. Accepted format are FASTG, GFA", true);
		addString("li", "-", "Name of the long-read data input file, - for stdin.", true);
		addString("lf", "fastq", "Format of the long-read data input file. This may be FASTQ/FASTA (MinION reads) or SAM/BAM (aligned with the assembly graph already)", true);
		addString("output", null, "Name of the output folder.", true);
		
		addBoolean("overwrite", false, "Whether to overwrite or reuse the intermediate file.");
		
		addString("mm2Path","","Path to the folder containing binary minimap2 if you don't have it included in $PATH");
		addString("mm2Preset", "map-ont", "Preset used by minimap2 to align long reads to the contigs");
		addInt("mm2Threads", 4, "Theads used by minimap2");
		addInt("qual", 1, "Minimum quality of alignment to considered");
 
		addBoolean("gui", false, "Whether using GUI or not.");
		
		addStdHelp();
	}
	
	public static void main(String[] args) throws IOException{
		CommandLine cmdLine = new NPGraphCmd();		
		args = cmdLine.stdParseLine(args);

		/***********************************************************************/
		String 	shortReadsAssembly = cmdLine.getStringVal("si"),
				shortReadsAssemblyFormat = cmdLine.getStringVal("sf"),
				longReadsInput = cmdLine.getStringVal("li"),
				longReadsInputFormat = cmdLine.getStringVal("lf"),
				outputDir = cmdLine.getStringVal("output"),
				mm2Path = cmdLine.getStringVal("mm2Path"),
				mm2Preset = cmdLine.getStringVal("mm2Preset");
		boolean overwrite = cmdLine.getBooleanVal("overwrite"),
				gui = cmdLine.getBooleanVal("gui");
		int 	mm2Threads = cmdLine.getIntVal("mm2Threads");
				
		Alignment.MIN_QUAL = cmdLine.getIntVal("qual");
		
		//Default output dir 
		if(outputDir == null) {
			outputDir = new File(shortReadsAssembly).getAbsoluteFile().getParent();
		}
		File outDir = new File(outputDir);
		if(!outDir.exists())
			outDir.mkdirs();
		
		mm2Path+=mm2Path.isEmpty()?"":"/";
		
		//1. Create an assembler object with appropriate file loader
		HybridAssembler hbAss = new HybridAssembler();
		hbAss.setShortReadsInput(shortReadsAssembly);
		hbAss.setShortReadsInputFormat(shortReadsAssemblyFormat);
		hbAss.setLongReadsInput(longReadsInput);
		hbAss.setLongReadsInputFormat(longReadsInputFormat);
		
		hbAss.setPrefix(outputDir);
		
		hbAss.setMinimapPath(mm2Path);
		hbAss.setMinimapPreset(mm2Preset);
		hbAss.setMinimapThreads(mm2Threads);
		
		hbAss.setOverwrite(overwrite);
		
		
		BidirectedGraph graph = hbAss.simGraph;      
        
		//4. Call the assembly function or invoke GUI to do so
        if(gui) {
        	//settings...
//        	GraphExplore.redrawGraphComponents(graph);
//            graph.display();
			NPGraphFX.setAssembler(hbAss);
			Application.launch(NPGraphFX.class,args);
        }else {
	        
			try {
				hbAss.prepareShortReadsProcess();
				hbAss.prepareLongReadsProcess();

				hbAss.assembly();
			} catch (InterruptedException e) {
				// TODO Auto-generated catch block
				System.err.println("Issue when assembly: \n" + e.getMessage());
				System.exit(1);
			}
        }
		
	}
	
	
//	public static void main(String[] args) throws IOException{
//		CommandLine cmdLine = new NPGraphCmd ();
//		args = cmdLine.stdParseLine(args);
//
//		Alignment.MIN_QUAL = cmdLine.getIntVal("qual");
//		String fastgFile = cmdLine.getStringVal("fastg");
//		String samFile = cmdLine.getStringVal("sam");
//		String pathFile = cmdLine.getStringVal("path");
//		String name = cmdLine.getStringVal("title");
//		
//		String styleSheet =
//			        "node {" +
//			        "	fill-color: black; z-index: 0;" +
//			        "}" +
//			        "edge {" +
//			        "	text-alignment: along;" +
//			        "}" +
//			        "node.marked {" +
//			        "	fill-color: red;" +
//			        "}" +
//			        "edge.marked {" +
//			        "	fill-color: red;" +
//			        "}";
//		System.setProperty("java.awt.headless", "false");
//		HybridAssembler hbAss = new HybridAssembler(fastgFile);
//		//For SAM file, run bwa first on the edited assembly_graph.fastg by running:
//		//awk -F '[:;]' -v q=\' 'BEGIN{flag=0;}/^>/{if(index($1,q)!=0) flag=0; else flag=1;}{if(flag==1) print $1;}' ../EcK12S-careful/assembly_graph.fastg > Eck12-careful.fasta
//		//TODO: need to make this easier
//		BidirectedGraph graph= hbAss.simGraph;
//		
//        //graph.addAttribute("ui.quality");
//        //graph.addAttribute("ui.antialias");
//        graph.addAttribute("ui.stylesheet", styleSheet);
//        graph.addAttribute("ui.default.title", name);
//
//        Viewer viewer = graph.display();
//        // Let the layout work ...
//        
//        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
//
//        
//        for (Node node : graph) {
//            node.addAttribute("ui.label", node.getId());
//            node.setAttribute("ui.style", "text-offset: -10;"); 
//            node.addAttribute("layout.weight", 10); 
//
//            if(BidirectedGraph.isMarker(node))
//            	node.setAttribute("ui.class", "marked");
//        }
//
//
//        try {
//        	if(pathFile!=null)
//        		hbAss.reduceFromSPAdesPaths(pathFile);
//        	hbAss.assembly(samFile);
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//        
//        
//        
//        HybridAssembler.promptEnterKey();
//        viewer.disableAutoLayout();
//        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
//	}
}
