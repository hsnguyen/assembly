package org.rtassembly;
import java.io.File;
import java.io.IOException;
import org.rtassembly.gui.NPGraphFX;
import org.rtassembly.npgraph.Alignment;
import org.rtassembly.npgraph.BidirectedGraph;
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
		
		addBoolean("overwrite", false, "Whether to overwrite or reuse the intermediate file");
		
		addString("mm2Path","","Absolute path to the folder containing binary minimap2");
		addString("mm2Opt", "-t4 -x map-ont", "Preset used by minimap2 to align long reads to the contigs");
		
		addInt("qual", 1, "Minimum quality of alignment to considered");
		addInt("dfs", 15, "Number of DFS steps to search");

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
				mm2Opt = cmdLine.getStringVal("mm2Opt");
		boolean overwrite = cmdLine.getBooleanVal("overwrite"),
				gui = cmdLine.getBooleanVal("gui");
			
		Alignment.MIN_QUAL = cmdLine.getIntVal("qual");
		BidirectedGraph.S_LIMIT=cmdLine.getIntVal("dfs");
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
		hbAss.setMinimapOpts(mm2Opt);
		
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
				if(hbAss.prepareShortReadsProcess(true) &&	hbAss.prepareLongReadsProcess()) {
					hbAss.assembly();
					hbAss.postProcessGraph();
				}
				else{
					System.err.println("Error with pre-processing step! Config and try again!");
					System.exit(1);
				}
					
			} catch (InterruptedException|IOException e) {
				// TODO Auto-generated catch block
				System.err.println("Issue when assembly: \n" + e.getMessage());
				e.printStackTrace();
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
