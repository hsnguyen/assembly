package org.rtassembly;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.graphstream.ui.view.Viewer;
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

		addString("graph", null, "Name of the assembly graph file. Accepted format are FASTG, GFA",true);
		addString("input", "-", "Name of the input file, - for stdin.", true);
		addString("format", "fastq", "Format of the input file. This may be FASTQ/FASTA (MinION reads) or SAM/BAM (aligned with the assembly graph).", true);
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
		String 	graphFileName = cmdLine.getStringVal("graph"),
				dataInputFileName = cmdLine.getStringVal("input"),
				inputFormat = cmdLine.getStringVal("format"),
				outputDir = cmdLine.getStringVal("output"),
				mm2Path = cmdLine.getStringVal("mm2Path"),
				mm2Setting = cmdLine.getStringVal("mm2Preset");
		boolean overwrite = cmdLine.getBooleanVal("overwrite"),
				gui = cmdLine.getBooleanVal("gui");
		int 	mm2Threads = cmdLine.getIntVal("mm2Threads");
				
		Alignment.MIN_QUAL = cmdLine.getIntVal("qual");
		
		//Default output dir 
		if(outputDir == null) {
			outputDir = new File(dataInputFileName).getAbsoluteFile().getParent();
		}
		File outDir = new File(outputDir);
		if(!outDir.exists())
			outDir.mkdirs();
		mm2Path+=mm2Path.isEmpty()?"":"/";
		
		//1. Create an assembler object with appropriate file loader
		HybridAssembler hbAss = new HybridAssembler(graphFileName);
		BidirectedGraph graph = hbAss.simGraph;
		
    	GraphExplore.redrawGraphComponents(graph);
        graph.display();
        
        //2. Create minimap2 index file for alignment if needed
        if(inputFormat.toLowerCase().startsWith("fast")) {
			File indexFile=new File(outputDir+"/assembly_graph.mmi");
			if(overwrite || !indexFile.exists()) {
				
				graph.printNodeSequencesToFile(outputDir+"/assembly_graph.fasta");
		
				
				try{
					ProcessBuilder pb = new ProcessBuilder(mm2Path+"minimap2","-V").redirectErrorStream(true);
					Process process =  pb.start();
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
						System.exit(1);
					}else{
						System.out.println("minimap version: " + version);
						if (version.compareTo("2.0") < 0){
							System.err.println(" ERROR: require minimap version 2 or above!");
							System.exit(1);
						}
					}
					
					
					ProcessBuilder pb2 = new ProcessBuilder(mm2Path+"minimap2", "-t", Integer.toString(mm2Threads), "-x", mm2Setting,"-d", outputDir+"/assembly_graph.mmi",outputDir+"/assembly_graph.fasta");
					Process indexProcess =  pb2.start();
					indexProcess.waitFor();
					
				}catch (IOException | InterruptedException e){
					System.err.println("Issue when indexing with minimap2: \n" + e.getMessage());
					e.printStackTrace();
					System.exit(1);
				}
			}
        }
		//4. Call the assembly function or invoke GUI to do so
        if(gui) {
        	//settings...
        	
//			NPGraphFX.setAssembler(hbAss);
//			Application.launch(NPGraphFX.class,args);
        }else {
	        
			try {
				hbAss.assembly(dataInputFileName,inputFormat, mm2Path, mm2Setting, mm2Threads, outputDir+"/assembly_graph.mmi");
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
