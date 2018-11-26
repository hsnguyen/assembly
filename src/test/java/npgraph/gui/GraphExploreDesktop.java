package npgraph.gui;

import java.io.IOException;

import org.graphstream.ui.view.Viewer;
import org.rtassembly.npgraph.BDGraph;
import org.rtassembly.npgraph.GraphUtil;
import org.rtassembly.npgraph.HybridAssembler;

public class GraphExploreDesktop {


	//	public static String dataFolder="/home/hoangnguyen/workspace/data/spades/"; //sony

	public static void main(String args[]) {
    	try {
			new GraphExploreDesktop();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


    }

    public GraphExploreDesktop() throws IOException{
    	System.setProperty("org.graphstream.ui", "javafx");
    	//System.setProperty("org.graphstream.ui.renderer", "org.graphstream.ui.j2dviewer.J2DGraphRenderer"); 
    	
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
//    	String sInput=dataFolder+sample+"/assembly_graph.fastg",
//    			output=dataFolder+sample+"/",
//    			lInput=dataFolder+sample+"/assembly_graph.sam";
    	
    	
    	
    	
    	/*
    	 * unicycler data set
    	 */
    	String dataFolder="/home/sonhoanghguyen/Projects/scaffolding/data/unicycler/";
//    	String sample="Acinetobacter_A1/";
    	String sample="Acinetobacter_AB30/"; //SimpleBinner.ANCHOR_CTG_LEN=2000 || Alignment.MIN_QUAL=10
//    	String sample="E_coli_K-12_MG1655/";
//    	String sample="E_coli_O25b_H4-ST131/";
//    	String sample="Klebsiella_30660_NJST258_1/";
//    	String sample="Klebsiella_MGH_78578/"; //SimpleBinner.ANCHOR_CTG_LEN=500
//    	String sample="Klebsiella_NTUH-K2044/";
//    	String sample="Mycobacterium_tuberculosis_H37Rv/";
//    	String sample="random_sequences_many_repeats/";
//    	String sample="random_sequences_no_repeats/";
//    	String sample="random_sequences_some_repeats/";
//    	String sample="Saccharomyces_cerevisiae_S288c/";
//    	String sample="Shigella_dysenteriae_Sd197/";
//    	String sample="Shigella_sonnei_53G/";
//    	String sample="Streptococcus_suis_BM407/";
    	
    	
//    	String quality="bad/";
//    	String quality="medium/";
    	String quality="good/";

    	String sInput=dataFolder+sample+quality+"spades/assembly_graph.fastg",
    			output=dataFolder+sample+quality+"/",
    			lInput=dataFolder+sample+quality+"mm2.sam";	
    
    	/*
    	 * Share code
    	 */
		HybridAssembler hbAss = new HybridAssembler();
		hbAss.setShortReadsInput(sInput);
		hbAss.setPrefix(output);
		hbAss.setShortReadsInputFormat("fastg");
		hbAss.prepareShortReadsProcess(false);
		
    	BDGraph graph= hbAss.simGraph;
    	
    	GraphUtil.redrawGraphComponents(graph);
    	graph.setAttribute("ui.style", GraphUtil.styleSheet);

        Viewer viewer=graph.display();
        
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
                
        /*
         * Testing reduce function
         */
        try {
        	hbAss.setLongReadsInput(lInput);
        	hbAss.setLongReadsInputFormat("sam");
        	hbAss.prepareLongReadsProcess();
        	
        	hbAss.assembly();
			
		} catch (IOException | InterruptedException e) {
			e.printStackTrace();
		}
        
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
        
//        HybridAssembler.promptEnterKey();
        hbAss.postProcessGraph();
        
//        HybridAssembler.promptEnterKey();
//        viewer.disableAutoLayout();
        
       
    }
    

}