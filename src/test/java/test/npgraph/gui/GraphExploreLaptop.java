package test.npgraph.gui;

import java.io.IOException;
import java.util.Iterator;
import org.graphstream.graph.*;
import org.graphstream.ui.view.Viewer;
import org.rtassembly.npgraph.BidirectedGraph;
import org.rtassembly.npgraph.GraphUtil;
import org.rtassembly.npgraph.HybridAssembler;

import japsa.seq.Sequence;


public class GraphExploreLaptop {
	
	public static String dataFolder="/home/s_hoangnguyen/Projects/scaffolding/test-graph/spades/"; //dell FASTG
//	public static String dataFolder="/home/s_hoangnguyen/Projects/scaffolding/test-graph/spades_v3.10/"; //dell GFA


	public static void main(String args[]) {
    	try {
			new GraphExploreLaptop();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


    }

    public GraphExploreLaptop() throws IOException{
    	System.setProperty("org.graphstream.ui", "javafx");
    	//System.setProperty("org.graphstream.ui.renderer", "org.graphstream.ui.j2dviewer.J2DGraphRenderer"); 
    	
    	//npscarf
//    	String sample="EcK12S-careful";
//    	String sample="Kp2146-careful";
    	String sample="Kp13883-careful";
 
//    	String sample="W303-careful";
//    	String sample="meta-careful";
//    	String sample="cp_S5";


    	
		HybridAssembler hbAss = new HybridAssembler();
		
		//npscarf
		hbAss.setShortReadsInput(dataFolder+sample+"/assembly_graph.fastg");
		hbAss.setPrefix(dataFolder+sample+"/");
		hbAss.setShortReadsInputFormat("fastg");
		hbAss.prepareShortReadsProcess(false);//change true/false to use/not use SPAdes path
		
		
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
        	hbAss.setLongReadsInput(dataFolder+sample+"/assembly_graph.sam");
        	hbAss.setLongReadsInputFormat("sam");
        	hbAss.prepareLongReadsProcess();
        	
        	hbAss.assembly();
			
		} catch (IOException | InterruptedException e) {
			e.printStackTrace();
		}

        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
        
        HybridAssembler.promptEnterKey();
        hbAss.postProcessGraph();

        viewer.disableAutoLayout();
    }
    

}