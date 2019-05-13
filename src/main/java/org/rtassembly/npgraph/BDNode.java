package org.rtassembly.npgraph;

import org.graphstream.graph.implementations.AbstractEdge;
import org.graphstream.graph.implementations.AbstractGraph;
import org.graphstream.graph.implementations.MultiNode;

public class BDNode extends MultiNode {

//    private static final Logger LOG = LoggerFactory.getLogger(BDNode.class);
	protected BDNode(AbstractGraph graph, String id) {
		super(graph, id);
	}

	// *** Helpers ***
	@Override
	protected char edgeType(AbstractEdge e) {		
		if(((BDEdge) e).getNodeDirection(this)==null)
				return IO_EDGE;
		else if(((BDEdge) e).getNodeDirection(this))
				return O_EDGE;
		else
				return I_EDGE;
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
	
	@Override
	public String toString() {
		return String.format("[Node %s: cov=%.2f len=%d]", getId(), getNumber("cov"), (int)getNumber("len"));
	}
}
