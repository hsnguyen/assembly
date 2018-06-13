package org.rtassembly.npgraph;

import java.io.File;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.graphstream.graph.Graph;
import org.graphstream.graph.implementations.AdjacencyListGraph;

import htsjdk.samtools.*;

public class GraphTest {
    public static void main(String[] args) throws Exception {
//        System.setProperty("org.graphstream.ui.renderer",
//                "org.graphstream.ui.j2dviewer.J2DGraphRenderer");
//        Graph g = new AdjacencyListGraph("g");
//        g.addNode("A").addAttribute("xyz", new double[] { 0, 0 });
//        g.addNode("B").addAttribute("xyz", new double[] { 10, 10 });
//
//        g.addEdge("AB", "A", "B", false)
//                .addAttribute(
//                        "ui.points",
//                        (Object) new double[] { 0, 0, 0, 0, 5, 0, 5, 10, 0, 10,
//                                10, 0 });
//
//        g.addAttribute("ui.stylesheet", "edge {shape: polyline; }"); // or shape: cubic-curve
//
//        g.display(false);
    	
//    	String data="/home/hoangnguyen/workspace/data/spades/EcK12S-careful/assembly_graph.fastg"; //sony
//
//    	final SamFileValidator validator=new SamFileValidator(new PrintWriter(System.out),8000);
//        validator.setIgnoreWarnings(true);
//        validator.setVerbose(true,1000);
//        validator.setErrorsToIgnore(Collections.singletonList(SAMValidationError.Type.MISSING_READ_GROUP));
//        SamReaderFactory factory=SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT);
//        SamReader samReader=factory.open(new File(data));
//        boolean ischeck=false;
//		SAMRecordIterator iter = samReader.iterator();
//		System.out.println(iter.next());
//        
//        samReader.close();
    	
    	Pattern pattern = Pattern.compile("(\\d+)([\\+\\-])(\\d+)([\\+\\-])");
    	String line="103-82-";
		Matcher matcher =pattern.matcher(line);
		System.out.println(line);
		String[] retval = new String[2];
		if (matcher.find()){
			retval[0]=matcher.group(1)+matcher.group(2);

			retval[1]=matcher.group(3) + matcher.group(4);
		}
		System.out.println(retval[0] + " " + retval[1]);
    }
}
