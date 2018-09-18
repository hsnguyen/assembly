package npgraph.gui;

import java.io.IOException;

import org.graphstream.ui.view.Viewer;
import org.rtassembly.npgraph.BidirectedGraph;
import org.rtassembly.npgraph.GraphUtil;
import org.rtassembly.npgraph.HybridAssembler;

public class GraphExploreDesktop {
	//imb desktop
	//npscarf
//	public static String dataFolder="/home/sonhoanghguyen/Projects/scaffolding/data/spades_3.7/"; 	
	//unicycler
	public static String dataFolder="/home/sonhoanghguyen/Projects/scaffolding/data/unicycler/"; 

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
    	
//    	npscarf
//    	String sample="EcK12S-careful";
//    	String sample="Kp2146-careful";
//    	String sample="Kp13883-careful";
 
//    	String sample="W303-careful";
//    	String sample="meta-careful";
//    	String sample="cp_S5";

    	//unicycler
//    	String sample="Acinetobacter_A1/";
//    	String sample="Acinetobacter_AB30/";
    	String sample="E_coli_K-12_MG1655/";
//    	String sample="E_coli_O25b_H4-ST131/";
//    	String sample="Klebsiella_30660_NJST258_1/";
//    	String sample="Klebsiella_MGH_78578/";
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

    	
		HybridAssembler hbAss = new HybridAssembler();
		
		//npscarf
//		hbAss.setShortReadsInput(dataFolder+sample+"/assembly_graph.fastg");
//		hbAss.setPrefix(dataFolder+sample+"/");
//		hbAss.setShortReadsInputFormat("fastg");
//		hbAss.prepareShortReadsProcess(false);//change true/false to use/not use SPAdes path
		
		
		//unicycler
		hbAss.setShortReadsInput(dataFolder+sample+quality+"spades/assembly_graph.fastg");
		hbAss.setPrefix(dataFolder+sample+quality);
		hbAss.setShortReadsInputFormat("fastg");
		hbAss.prepareShortReadsProcess(true);
		
    	BidirectedGraph graph= hbAss.simGraph;
    	
    	GraphUtil.redrawGraphComponents(graph);
    	graph.setAttribute("ui.style", GraphUtil.styleSheet);

        Viewer viewer=graph.display();
        
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
                
        /*
         * Testing reduce function
         */
        try {
        	//npscarf
//        	hbAss.setLongReadsInput(dataFolder+sample+"/assembly_graph.sam");
        	//unicycler
        	hbAss.setLongReadsInput(dataFolder+sample+quality+"bwa.sam");
        	hbAss.setLongReadsInputFormat("sam");
        	hbAss.prepareLongReadsProcess();
        	
        	hbAss.assembly();
			
		} catch (IOException | InterruptedException e) {
			e.printStackTrace();
		}
        
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
        
        HybridAssembler.promptEnterKey();
        hbAss.postProcessGraph();

//        viewer.disableAutoLayout();
        /*
         * Testing BidirectedEdge id pattern
         */
//    	String pattern = "^\\[([0-9\\+\\-]*)\\]([oi])\\[([0-9\\+\\-]*)\\]([oi])$";
//        // Create a Pattern object
//        Pattern r = Pattern.compile(pattern);
//        // Now create matcher object.
//        String id="[3]i[4+8+]o";
//        Matcher m = r.matcher(id);
//         	
//        if(m.find()){
//        	System.out.println(m.group(1)+"|"+m.group(2)+"|"+m.group(3)+"|"+m.group(4));
//        } else
//        	System.out.println("Fuck");
    }
    

}