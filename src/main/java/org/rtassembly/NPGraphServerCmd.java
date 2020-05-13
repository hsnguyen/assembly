package org.rtassembly;
import java.io.File;
import java.io.IOException;
import java.lang.invoke.MethodHandles;

import org.apache.logging.log4j.Logger;
import org.graphstream.stream.thread.ThreadProxyPipe;
import org.graphstream.ui.fx_viewer.FxDefaultView;
import org.graphstream.ui.fx_viewer.FxViewer;
import org.graphstream.ui.javafx.FxGraphRenderer;
import org.graphstream.ui.view.Viewer;
import org.graphstream.ui.view.Viewer.CloseFramePolicy;
import org.apache.logging.log4j.LogManager;
import org.rtassembly.npgraph.Alignment;
import org.rtassembly.npgraph.BDGraph;
import org.rtassembly.npgraph.GraphUtil;
import org.rtassembly.npgraph.RealtimeGraphWatcher;
import org.rtassembly.npgraph.HybridAssembler;
import org.rtassembly.npgraph.SimpleBinner;
import org.rtassembly.npgraph.grpc.AssemblyGuideServer;

import japsa.util.CommandLine;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.ColumnConstraints;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.RowConstraints;
import javafx.stage.Stage;


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
			hbAss.input.setShortReadsInput(shortReadsInput);
		if(shortReadsInputFormat!=null && !shortReadsInputFormat.isEmpty())
			hbAss.input.setShortReadsInputFormat(shortReadsInputFormat);
		
		hbAss.input.setLongReadsInputFormat("paf"); //to pass prepareLongReadProcess()
		
		hbAss.setPrefix(outputDir);
		if(shortReadsBinInput!=null && !shortReadsBinInput.isEmpty())
			hbAss.input.setBinReadsInput(shortReadsBinInput);
		
	
		hbAss.setOverwrite(overwrite);
		hbAss.input.setUseSPAdesPath(spaths);
		        
		//4. Call the assembly function or invoke GUI to do so
		if(shortReadsInput.isEmpty() && !gui) {
			System.out.println(cmdLine.usageString());			
			System.exit(-1);
        }else {
			try {
				if(hbAss.prepareShortReadsProcess() &&	hbAss.prepareLongReadsProcess()) {
					if(gui) {
						GraphView.hbAss=hbAss;
						Application.launch(GraphView.class, args);
					}
					AssemblyGuideServer assServer = new AssemblyGuideServer(port, hbAss);
					assServer.start();
				}
				else{
					logger.error("Error(s) happened! Check the log above for more info.");
					System.exit(1);
				}
					
			} catch (Exception e) {
				logger.error("Error {}", e);
				System.exit(1);
			}
        }
		
	}
}


///////////////////////////////////////////////////////////////////////////////////////////
class GraphView extends Application{

	static HybridAssembler hbAss;	
    private GridPane createAutoresizeGridPane(int ncols, int nrows){
        GridPane gridpane = new GridPane();
        for (int i = 0; i < ncols; i++) {
            ColumnConstraints column = new ColumnConstraints();
            column.setPercentWidth(100/ncols);
            gridpane.getColumnConstraints().add(column);
        }
        for (int i = 0; i < nrows; i++) {
            RowConstraints row = new RowConstraints();
            row.setPercentHeight(100/nrows);
            gridpane.getRowConstraints().add(row);
        }
        gridpane.setPadding(new Insets(5, 5, 5, 5));
        gridpane.setVgap(5);
        gridpane.setHgap(5);
        return gridpane;
    }
    private GridPane addGraphResolverPane(){  		
    	GridPane mainGrid = createAutoresizeGridPane(1,1);
		mainGrid.setStyle("-fx-background-color: #C0C0C0;");
		
		
		ThreadProxyPipe pipe = new ThreadProxyPipe() ;
		pipe.init(hbAss.simGraph);
		Viewer graphViewer = new FxViewer(pipe);
		System.setProperty("org.graphstream.ui", "javafx");

		FxDefaultView view = new FxDefaultView(graphViewer, "npGraph", new FxGraphRenderer());
		graphViewer.addView(view);
		graphViewer.enableAutoLayout();
		graphViewer.setCloseFramePolicy(CloseFramePolicy.CLOSE_VIEWER);
		
		mainGrid.getChildren().add(view);
		
		return mainGrid;
    }
	@Override
	public void start(Stage primaryStage) throws Exception {
        BorderPane gBorder = new BorderPane();
        gBorder.setCenter(addGraphResolverPane());
        Scene gscene = new Scene(gBorder);
        primaryStage.setScene(gscene);
        primaryStage.setTitle("Assembly graph");
        primaryStage.setOnCloseRequest(e -> {
        	hbAss.terminateAlignmentProcess();
            Platform.exit();
            System.exit(0);
        });
        
    	BDGraph graph= hbAss.simGraph;
		
    	graph.setAttribute("ui.style", GraphUtil.styleSheet);
    	graph.setAttribute("layout.force", .5);
    	graph.setAttribute("layout.weight", .1);
    	

//        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
//        for (Node node : graph) {
//        	node.setAttribute("ui.label", node.getId());
//        }        

        primaryStage.show();

        /*
         * Testing reduce function
         */
        /*
         * A thread to run the algorithm
         */
	    new Thread(new Runnable(){

			@Override
			public void run() {

				try{
//					hbAss.assemblyPAF();
					hbAss.assembly();
				}catch (Exception e){
					System.err.println(e.getMessage());
					e.printStackTrace();
				}

			}
			
		}).start();
      
	}   

}