package org.rtassembly.npgraph;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.graphstream.graph.implementations.AbstractEdge;
import org.graphstream.graph.implementations.AbstractGraph;
import org.graphstream.graph.implementations.AbstractNode;
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

//	// *** Overide for the case of self-loop edges ***
	@Override
	public int getDegree() {
		int retval=degree;
		for(int i=0;i<degree;i++)
			if(edges[i].getNode0()==edges[i].getNode1())
				retval++;		
		return retval;
	}

	@Override
	public int getInDegree() {
		int retval=oStart;
		for(int i=0;i<ioStart;i++)
			if(edges[i].getNode0()==edges[i].getNode1())
				retval++;		
		return retval;
	}

	@Override
	public int getOutDegree() {
		int retval = degree - ioStart;
		for(int i=oStart;i<degree;i++)
			if(edges[i].getNode0()==edges[i].getNode1())
				retval++;	
		return retval;
	}
}
