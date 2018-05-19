package org.rtassembly.npgraph;

import java.io.IOException;
import java.util.Iterator;
import org.graphstream.graph.*;
import japsa.seq.Sequence;


public class GraphExplore {

	public static String spadesFolder="/home/sonhoanghguyen/Projects/scaffolding/data/spades_3.7/"; //imb desktop
//	public static String spadesFolder="/home/hoangnguyen/workspace/data/spades/"; //sony
//	public static String spadesFolder="/home/s_hoangnguyen/Projects/scaffolding/test-graph/spades/"; //dell FASTG
//	public static String spadesFolder="/home/s_hoangnguyen/Projects/scaffolding/test-graph/spades_v3.10/"; //dell GFA


	public static void main(String args[]) {
    	try {
			new GraphExplore();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


    }

    public GraphExplore() throws IOException{
    	System.setProperty("org.graphstream.ui", "javafx");
    	//System.setProperty("org.graphstream.ui.renderer", "org.graphstream.ui.j2dviewer.J2DGraphRenderer"); 
//    	String sample="EcK12S-careful";
    	String sample="Kp2146-careful";
//    	String sample="Kp13883-careful";
    	
    	//memory!!!
//    	String sample="W303-careful";
//    	String sample="meta-careful";
//    	String sample="cp_S5";

		HybridAssembler hbAss = new HybridAssembler();
		hbAss.setShortReadsInput(spadesFolder+sample+"/assembly_graph.fastg");
		hbAss.setShortReadsInputFormat("fastg");
		hbAss.prepareShortReadsProcess();
		
    	BidirectedGraph graph= hbAss.simGraph;

    	redrawGraphComponents(graph);

        graph.display();
        
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());
                
        /*
         * Testing reduce function
         */
        try {
//        	HybridAssembler.promptEnterKey();
//			hbAss.reduceFromSPAdesPaths(spadesFolder+sample+"/contigs.paths");
//			HybridAssembler.promptEnterKey();
        	hbAss.setLongReadsInput(spadesFolder+sample+"/assembly_graph.sam");
        	hbAss.setLongReadsInputFormat("sam");
        	hbAss.prepareLongReadsProcess();
        	
        	hbAss.assembly();
			
		} catch (IOException | InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        
        //TODO: thorough cleaning... should have flag dead for each node
//        boolean dead=true;
//        while(dead){
//        	dead=false;
//	        for (Node node : graph) {
//	            if((node.getDegree() < 2) 
////	            	|| (node.getDegree()==0 && ((Sequence)node.getAttribute("seq")).length() < 1000)
//	            	){
//	            	graph.removeNode(node);
//	            	dead=true;
//	            }
//	            	
//	        }
//        }
        System.out.println("Node: " + graph.getNodeCount() + " Edge: " + graph.getEdgeCount());

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
    
    public static void redrawGraphComponents(BidirectedGraph graph) {
//      graph.addAttribute("ui.quality");
//      graph.addAttribute("ui.antialias");
//    	graph.addAttribute("ui.default.title", "New real-time hybrid assembler");
    	
    	graph.setAttribute("ui.style", dynamicStyle);

    	for (Node node : graph) {

    		Sequence seq = (Sequence) node.getAttribute("seq");
    		double lengthScale = 1+(Math.log10(seq.length())-2)/3.5; //100->330,000
          
			if(lengthScale<1) lengthScale=1;
			else if(lengthScale>2) lengthScale=2;
	          
			int covScale = (int) Math.round(node.getNumber("cov")/BidirectedGraph.RCOV);

//          Color[] palette= {Color.GRAY,Color.BLUE,Color.YELLOW,Color.ORANGE,Color.GREEN,Color.PINK,Color.MAGENTA,Color.RED};
//          Color color=null;
	          
			String[] palette= {"grey","blue","yellow","orange","green","pink","magenta","red"};
			String color=null;
			if(covScale>=palette.length)
				color=palette[palette.length-1];
			else
				color=palette[covScale];
          
//          node.addAttribute("ui.color", color);
//          node.addAttribute("ui.size", lengthScale+"gu");
          
//          node.addAttribute("ui.label", covScale);
          
			node.setAttribute("ui.label", node.getId());
//			node.addAttribute("ui.label", (int)(node.getNumber("cov")));


			node.setAttribute("ui.style", "	size: " + lengthScale + "gu;" +
        		  						"	fill-color: "+color+";" +
        		  			            " 	stroke-mode: plain;" +
        		  			            "	stroke-color: black;" +
        		  			            "	stroke-width: 2px;");
    	}
    	
    	graph.balancing();
//    	for(Edge e : graph.getEdgeSet()) {
//    		e.addAttribute("ui.label", (int)(e.getNumber("cov")));
//    	}
    	
    
    }
    
    //TODO: traverse the graph and print FASTA/GFA out
    public void explore(Node source) {
        Iterator<? extends Node> k = source.getBreadthFirstIterator();

        while (k.hasNext()) {
            Node next = k.next();
            next.setAttribute("ui.class", "marked");
            sleep();
        }
    }

    protected void sleep() {
        try { Thread.sleep(1000); } catch (Exception e) {}
    }

    protected String styleSheet =
        "node {" +
        "	fill-color: black;" +
        "}" +
        "node.marked {" +
        "	fill-color: red;" +
        "}" +
        "edge.marked {" +
        "	fill-color: red;" +
        "}";
    	
    protected static String dynamicStyle =
            "node {" +
    		"	fill-mode: dyn-plain;" +
            "	size-mode: dyn-size;" +
            " 	stroke-mode: plain;" +
            "	stroke-color: black;" +
            "	stroke-width: 2px;" +
            "}";
}