package npgraph.gui;

import java.io.IOException;

import org.graphstream.graph.Node;
import org.graphstream.ui.view.View;
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
////    	String sample="Klebsiella_30660_NJST258_1/";
////    	String sample="Klebsiella_MGH_78578/"; //SimpleBinner.ANCHOR_CTG_LEN=500
////    	String sample="Klebsiella_NTUH-K2044/";
////    	String sample="Mycobacterium_tuberculosis_H37Rv/";
////    	String sample="random_sequences_many_repeats/";
////    	String sample="random_sequences_no_repeats/";
////    	String sample="random_sequences_some_repeats/";
////    	String sample="Saccharomyces_cerevisiae_S288c/";
//    	String sample="Shigella_dysenteriae_Sd197/"; //SimpleBinner.ANCHOR_CTG_LEN=2000
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
    	 * To test unicycler's graph:
    	 */
//    	sInput="/home/sonhoanghguyen/Projects/scaffolding/npgraph/results_2/unicycler/citrobacter-freundii_CAV1374/003_bridges_applied.gfa";
//    	output="/home/sonhoanghguyen/Projects/scaffolding/npgraph/results_2/spades/citrobacter-freundii_CAV1374/";
//    	lInput="";	
    	
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
    	
		//2.zymo
    	sInput="/home/sonhoanghguyen/Projects/scaffolding/zymo/metaSPAdes/assembly_graph.fastg";
    	output="/home/sonhoanghguyen/Projects/scaffolding/zymo/metaSPAdes/npgraph";
    	lInput="/home/sonhoanghguyen/Projects/scaffolding/zymo/metaSPAdes/assembly_graph_G.bam";
//    	lInput="/media/sonhoanghguyen/Seagate Backup Plus Drive/Data/zymo/assembly_graph_P.bam";
		binFile="/home/sonhoanghguyen/Projects/scaffolding/zymo/metaSPAdes/bin";	
		
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
		HybridAssembler hbAss = new HybridAssembler();
		hbAss.setUseSPAdesPath(useSPAdesPath);
		
		if(!sInput.isEmpty())
			hbAss.setShortReadsInput(sInput);
		if(!output.isEmpty())
			hbAss.setPrefix(output);
		if(!binFile.isEmpty())
			hbAss.setBinReadsInput(binFile);
		
		hbAss.setAligner("minimap2");
		hbAss.setMSA("spoa");
		hbAss.prepareShortReadsProcess();
		
    	BDGraph graph= hbAss.simGraph;
    	
    	graph.setAttribute("ui.style", GraphUtil.styleSheet);
    	graph.setAttribute("layout.force", .5);
    	graph.setAttribute("layout.weight", .1);
        Viewer viewer=graph.display();
//        HybridAssembler.promptEnterKey();
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
//        for (Node node : graph) {
//            node.setAttribute("ui.label", node.getId());
//        }        
//        HybridAssembler.promptEnterKey();

        /*
         * Testing reduce function
         */
        try {
        	hbAss.setLongReadsInput(lInput);
        	hbAss.prepareLongReadsProcess();
        	
        	hbAss.assembly();
			
		} catch (IOException | InterruptedException e) {
			e.printStackTrace();
		}
        
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
                
//        HybridAssembler.promptEnterKey();
        hbAss.postProcessGraph();
        HybridAssembler.promptEnterKey();

        for (Node node : graph) {
            node.setAttribute("ui.label", node.getId());
        }
        viewer.disableAutoLayout();
       
    }
    

}