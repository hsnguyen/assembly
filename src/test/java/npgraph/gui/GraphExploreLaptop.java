package npgraph.gui;

import org.graphstream.graph.Node;
import org.graphstream.stream.thread.ThreadProxyPipe;
import org.graphstream.ui.fx_viewer.FxDefaultView;
import org.graphstream.ui.fx_viewer.FxViewer;
import org.graphstream.ui.javafx.FxGraphRenderer;
import org.graphstream.ui.view.Viewer;
import org.graphstream.ui.view.Viewer.CloseFramePolicy;
import org.rtassembly.npgraph.BDGraph;
import org.rtassembly.npgraph.GraphUtil;
import org.rtassembly.npgraph.HybridAssembler;

import javafx.application.Application;
import javafx.application.Platform;
import javafx.geometry.Insets;
import javafx.scene.Scene;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.ColumnConstraints;
import javafx.scene.layout.GridPane;
import javafx.scene.layout.RowConstraints;
import javafx.stage.Stage;

public class GraphExploreLaptop extends Application{

	static HybridAssembler hbAss = new HybridAssembler();	
	//	public static String dataFolder="/home/hoangnguyen/workspace/data/spades/"; //sony

	public static void main(String args[]) {
    	System.setProperty("org.graphstream.ui", "javafx");
    	//System.setProperty("org.graphstream.ui.renderer", "org.graphstream.ui.j2dviewer.J2DGraphRenderer"); 
    	String sInput="", lInput="", output="", binFile="";
    	boolean useSPAdesPath=false;
    	String dataFolder="/home/sonnguyen/Projects/npGraph/test/"; //dell FASTG
    	
		sInput=dataFolder+"spades/assembly_graph.fastg";
    	lInput=dataFolder+"alignment.paf";
//		lInput=dataFolder+"EcK12S_ONT.fastq.gz";
    	
    	/*******************************************************************************
    	 ****************************** Share code *************************************
    	 *******************************************************************************/
//		hbAss.setUseSPAdesPath(useSPAdesPath);
		
		if(!sInput.isEmpty())
			hbAss.setShortReadsInput(sInput);
		if(!output.isEmpty())
			hbAss.setPrefix(output);
		if(!binFile.isEmpty())
			hbAss.setBinReadsInput(binFile);
		
		hbAss.setAligner("/home/sonnguyen/sw/minimap2/minimap2");
		hbAss.setMSA("kalign"); //spoa: out of memory with reads > 18kbp
		hbAss.prepareShortReadsProcess();
      	hbAss.setLongReadsInput(lInput);
//      	hbAss.setLongReadsInputFormat("fastq");
      	hbAss.prepareLongReadsProcess();
      	
		Application.launch(GraphExploreLaptop.class,args);
		
    }
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
    	

        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
        for (Node node : graph) {
        	node.setAttribute("ui.label", node.getId());
        }        

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
					hbAss.assembly2();
//					hbAss.postProcessGraph();
				}catch (Exception e){
					System.err.println(e.getMessage());
					e.printStackTrace();
				}

			}
			
		}).start();
      
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
	}   

}