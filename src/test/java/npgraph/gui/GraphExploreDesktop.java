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

/*
 * Follow this to run tests on Eclipse:
 * https://stackoverflow.com/questions/52144931/how-to-add-javafx-runtime-to-eclipse-in-java-11
 */
public class GraphExploreDesktop extends Application{

	static HybridAssembler hbAss = new HybridAssembler();	
	//	public static String dataFolder="/home/hoangnguyen/workspace/data/spades/"; //sony

	public static void main(String args[]) {
    	System.setProperty("org.graphstream.ui", "javafx");
    	//System.setProperty("org.graphstream.ui.renderer", "org.graphstream.ui.j2dviewer.J2DGraphRenderer"); 
    	String sInput="", lInput="", output="", binFile="";
    	boolean useSPAdesPath=false;
    	/*
    	 * npScarf data set
    	 */
//    	String dataFolder="/home/sonhoanghguyen/Projects/scaffolding/data/spades_3.7/";
//    	
////    	String sample="EcK12S-careful";
//    	String sample="Kp2146-careful";
////    	String sample="Kp13883-careful";
// 
////    	String sample="W303-careful";
////    	String sample="meta-careful";
////    	String sample="cp_S5";
//    	sInput=dataFolder+sample+"/assembly_graph.fastg";
//		output=dataFolder+sample+"/";
//		lInput=dataFolder+sample+"/assembly_graph.sam";
    	    	
//    	/*
//    	 * unicycler data set
//    	 */
//    	String dataFolder="/home/sonhoanghguyen/Projects/scaffolding/data/unicycler/synthetic/";
////    	String sample="Acinetobacter_A1/";
////    	String sample="Acinetobacter_AB30/"; //SimpleBinner.ANCHOR_CTG_LEN=2000 || Alignment.MIN_QUAL=10
////    	String sample="E_coli_K-12_MG1655/";
////    	String sample="E_coli_O25b_H4-ST131/";
//    	String sample="Klebsiella_30660_NJST258_1/";
////    	String sample="Klebsiella_MGH_78578/"; //SimpleBinner.ANCHOR_CTG_LEN=500
////    	String sample="Klebsiella_NTUH-K2044/";
////    	String sample="Mycobacterium_tuberculosis_H37Rv/";
////    	String sample="random_sequences_many_repeats/";
////    	String sample="random_sequences_no_repeats/";
////    	String sample="random_sequences_some_repeats/";
////    	String sample="Saccharomyces_cerevisiae_S288c/";
////    	String sample="Shigella_dysenteriae_Sd197/"; //SimpleBinner.ANCHOR_CTG_LEN=2000
////    	String sample="Shigella_sonnei_53G/";
////    	String sample="Streptococcus_suis_BM407/";
////    	
////    	
////    	String quality="bad/";
////    	String quality="medium/";
//    	String quality="good/";
//
//    	sInput=dataFolder+sample+quality+"spades/assembly_graph.fastg";
//    	output=dataFolder+sample+quality+"/";
//    	lInput=dataFolder+sample+quality+"mm2.sam";	
    
    	
    	/*
    	 * To test real-time assembly
    	 */
    	String dataFolder="/home/sonhoanghguyen/Projects/scaffolding/npgraph/realtime/data/";
    	String sample="B09_Dap8";
    	sInput=dataFolder+sample+"/assembly_graph.fastg";
    	output=dataFolder+sample;
    	lInput=dataFolder+sample+"/assembly_graph.bam";	
    	BDGraph.GOOD_SUPPORT=10;
    	/*
    	 * Metagenomics data:
    	 */

//    	//1. porecamp
//    	String sass="metaSPAdes";
////    	String sass="megaHIT";
//    	sInput="/home/sonhoanghguyen/Projects/scaffolding/data/porecamp/"+sass+"/assembly_graph.fastg";
//    	output="/home/sonhoanghguyen/Projects/scaffolding/data/porecamp/";
//    	lInput="/home/sonhoanghguyen/Projects/scaffolding/data/porecamp/"+sass+"/assembly_graph.sam";
//		binFile="/home/sonhoanghguyen/Projects/scaffolding/data/porecamp/metabat/"+sass+"_contigs.bin";	
    	
//		//2.zymo
//    	sInput="/home/sonhoanghguyen/Projects/scaffolding/zymo/metaSPAdes/assembly_graph.fastg";
//    	output="/home/sonhoanghguyen/Projects/scaffolding/zymo/metaSPAdes/npgraph";
//    	lInput="/home/sonhoanghguyen/Projects/scaffolding/zymo/metaSPAdes/assembly_graph_G.bam";
////    	lInput="/media/sonhoanghguyen/Seagate Backup Plus Drive/Data/zymo/assembly_graph_P.bam";
//		binFile="/home/sonhoanghguyen/Projects/scaffolding/zymo/metaSPAdes/bin";	
//	
//		useSPAdesPath=true;
    	
    	/*
    	 * MRSA day 0
    	 */
//    	String 	sInput="/home/sonhoanghguyen/Projects/scaffolding/data/spades_v3.10/S.aureus_day0/spades/assembly_graph.fastg",
//    			output="/home/sonhoanghguyen/Projects/scaffolding/data/spades_v3.10/S.aureus_day0",
//    			lInput="/home/sonhoanghguyen/Projects/scaffolding/data/spades_v3.10/S.aureus_day0/MRSA_Rapid_230916.fastq";  
    	/*******************************************************************************
    	 ****************************** Share code *************************************
    	 *******************************************************************************/
    	HybridAssembler.VERBOSE=true;
//		hbAss.setUseSPAdesPath(useSPAdesPath);
		
		if(!sInput.isEmpty())
			hbAss.setShortReadsInput(sInput);
		if(!output.isEmpty())
			hbAss.setPrefix(output);
		if(!binFile.isEmpty())
			hbAss.setBinReadsInput(binFile);
		
		hbAss.setAligner("minimap2");
		hbAss.setMSA("kalign"); //spoa: out of memory with reads > 18kbp
		hbAss.prepareShortReadsProcess();
      	hbAss.setLongReadsInput(lInput);
      	hbAss.prepareLongReadsProcess();
      	
		Application.launch(GraphExploreDesktop.class,args);

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
					hbAss.assembly();
//					myass.postProcessGraph();
				}catch (Exception e){
					System.err.println(e.getMessage());
					e.printStackTrace();
				}

			}
			
		}).start();
      
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
	}
    

}