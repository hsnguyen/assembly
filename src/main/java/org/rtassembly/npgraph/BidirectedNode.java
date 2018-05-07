package org.rtassembly.npgraph;

import org.graphstream.graph.implementations.AbstractEdge;
import org.graphstream.graph.implementations.AbstractGraph;
import org.graphstream.graph.implementations.MultiNode;

public class BidirectedNode extends MultiNode {

//    private static final Logger LOG = LoggerFactory.getLogger(BidirectedNode.class);
	protected BidirectedNode(AbstractGraph graph, String id) {
		super(graph, id);
	}

	// *** Helpers ***
	@Override
	protected char edgeType(AbstractEdge e) {
		BidirectedNode opposite = (BidirectedNode) e.getOpposite(this);
		
		if(((BidirectedEdge) e).getDir(this)) {
			if(this==opposite && ((BidirectedEdge) e).getDir0()!=((BidirectedEdge) e).getDir1())
				return IO_EDGE;
			else
				return O_EDGE;
		}
		else {
			if(this==opposite && ((BidirectedEdge) e).getDir0()!=((BidirectedEdge) e).getDir1())
				return IO_EDGE;
			else
				return I_EDGE;
		}
	}

	// *** Access methods ***
	
}
