package org.rtassembly.npgraph;

import java.io.File;
import java.io.IOException;
import java.time.LocalTime;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

import org.graphstream.algorithm.ConnectedComponents;
import org.graphstream.algorithm.ConnectedComponents.ConnectedComponent;
import org.graphstream.graph.Edge;
import org.graphstream.graph.Node;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import japsa.bio.np.RealtimeAnalysis;
import japsa.seq.JapsaAnnotation;
import japsa.seq.Sequence;

public class RealtimeGraphWatcher extends RealtimeAnalysis{
    private static final Logger LOG = LoggerFactory.getLogger(RealtimeGraphWatcher.class);

	public static boolean KEEP=false; //whether or not keeping the low-coverage nodes
	BDGraph outputGraph;
	ConnectedComponents rtComponents;
	HashSet<BDEdge> cutEdges;
	int numberOfComponents=0;
	HybridAssembler hAss;
	
	public RealtimeGraphWatcher(HybridAssembler hAss) {
		this.hAss=hAss;
		rtComponents = new ConnectedComponents();
		rtComponents.init(hAss.simGraph);
		rtComponents.setCutAttribute("cut");
		numberOfComponents=rtComponents.getConnectedComponentsCount();
		//initial cleaning
		removeDeadEdges();
		
	}
	
	//Should be applied for a Collection of Edges, or subgraph instead of the whole graph???
	synchronized private void removeDeadEdges(){
		Set<Edge> cleanedEdges = new HashSet<>();
		while(true){
			cleanedEdges.stream().forEach(e->hAss.simGraph.removeEdge(e));
			cleanedEdges.clear();
			hAss.simGraph.edges().filter(e->checkDeadEdges(e)).forEach(e->cleanedEdges.add(e));
			if(cleanedEdges.isEmpty())
				break;
		}
	}
	
	//check and remove edge if it lead to an insignificant component
	synchronized private boolean checkDeadEdges(Edge e){
		BDNode 	src=(BDNode) e.getNode0(), 
				dst=(BDNode) e.getNode1();
		if(src==dst)
			return false;
		boolean srcDir=((BDEdge)e).getDir0(),
				dstDir=((BDEdge)e).getDir1();
		e.setAttribute("cut");
		if( numberOfComponents<rtComponents.getConnectedComponentsCount() 
			&&
			((srcDir?src.getOutDegree():src.getInDegree())>1 || (dstDir?dst.getOutDegree():dst.getInDegree())>1)
			&&
			(!isSignificantComponent(rtComponents.getConnectedComponentOf(dst))) || !isSignificantComponent(rtComponents.getConnectedComponentOf(src))
			){
				return true;
		}
		e.removeAttribute("cut");
		return false;
	}
	
	synchronized private void removeBadComponents() {
		List<Node> 	removeNodes=new ArrayList<Node>();
		for (Iterator<ConnectedComponent> compIter = rtComponents.iterator(); compIter.hasNext(); ) {
			ConnectedComponent comp = compIter.next();		
			if(!isSignificantComponent(comp))
				comp.nodes().forEach(n->removeNodes.add(n));			
		}
		//Remove abundant components here
		removeNodes.stream().forEach(n->hAss.simGraph.removeNode(n));
	}
	//FIXME: find most significant path and check if it cover >90%?
	synchronized private boolean isSignificantComponent(ConnectedComponent comp){
		int clen = 1, tlen=0;
		double ccov = 0;
		double threshold=KEEP?0:hAss.simGraph.binner.leastBin.estCov; //lower-bound for coverage?!
		for(Node n:comp.getNodeSet()){
			if(SimpleBinner.getBinIfUniqueNow(n)!=null)
				return true;
			
			int len = (int) (n.getNumber("len")-BDGraph.getKmerSize());
			double cov = n.getNumber("cov");
			tlen+=len;
			if(GraphUtil.isLikelyStillPresented(n)){
				ccov+=cov*len;
				clen+=len;
			}
			
		}
		double 	aveCov=ccov/tlen;
		if(GraphUtil.approxCompare(aveCov, threshold) < 0 || clen < Math.max(SimpleBinner.ANCHOR_CTG_LEN,.5*tlen) || (!KEEP&&comp.getEdgeCount() < 1))
			return false;
		else{
//			System.out.printf("Keep non-unique components with cov=%.2f, clen=%d, tlen=%d\n", aveCov, clen, tlen);
			return true;
		}
	}
	/*
	 * FIXME: output graph updating is significantly slower as the graph getting more complex
	 */
	synchronized void update(boolean lastTime) {
		//cleaning...
		if(lastTime)
			removeDeadEdges();
		removeBadComponents();
		
		cutEdges = new HashSet<BDEdge>();
		//then set the cut edges: just for outputGraph stats (will reset after). FIXME: need a Eulerian paths finding iff lastTime
		hAss.simGraph.nodes()
		.forEach(n->{
			if(n.getInDegree()>=2)
				n.enteringEdges().forEach(e->{e.setAttribute("cut");cutEdges.add((BDEdge) e);});
			if(n.getOutDegree()>=2)
				n.leavingEdges().forEach(e->{e.setAttribute("cut");cutEdges.add((BDEdge) e);});

		});
		
		
		outputGraph=new BDGraph();
		JapsaAnnotation annotation=null;
		BDPath repPath=null; //representative path of a component
		System.out.println("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
		System.out.println("Current time: " + LocalTime.now());

		HashMap<PopBin, List<Integer>> binn50 = new HashMap<>();
		PopBin b=null;
		
		for (Iterator<ConnectedComponent> compIter = rtComponents.iterator(); compIter.hasNext(); ) {
			ConnectedComponent comp = compIter.next();
//			System.out.printf("... id=%s edges=%d nodes=%d \n", comp.id, comp.getEdgeCount(), comp.getNodeCount());
			if(comp.getNodeCount()==0)
				continue;
			//Start analyzing significant components from here
			//check comp: should be linear paths, should start with node+
			 Node node = comp.nodes().toArray(Node[]::new)[0];
			 repPath = new BDPath(node);
			 boolean isCircular=false;
			 
			 if(comp.getEdgeCount()>=1){
				 //extend to
				 Node curNode=node;
				 boolean curDir=true;
				 List<Edge> ways = (curDir?curNode.leavingEdges():curNode.enteringEdges()).filter(e->!e.hasAttribute("cut")).collect(Collectors.toList());
				 while(ways.size()==1){
					 Edge edge = ways.get(0);
					 repPath.add(edge);
					 curNode=edge.getOpposite(curNode);
					 
					 if(curNode==node){//circular
						 isCircular=true;
						 break;
					 }
					 if(((BDEdge) edge).getNodeDirection((BDNode)curNode)!=null)
						 curDir=!((BDEdge) edge).getNodeDirection((BDNode)curNode);
					 ways = (curDir?curNode.leavingEdges():curNode.enteringEdges()).filter(e->!e.hasAttribute("cut")).collect(Collectors.toList());
				 }
				 
				 //if linear: reverse
				 if(!isCircular){
					 repPath=repPath.reverse();
					 //extend in opposite direction
					 curNode=node;
					 curDir=false;
					 ways = (curDir?curNode.leavingEdges():curNode.enteringEdges()).filter(e->!e.hasAttribute("cut")).collect(Collectors.toList());

					 while(ways.size()==1){
						 Edge edge = ways.get(0);
						 repPath.add(edge);
						 curNode=edge.getOpposite(curNode);
						 if(((BDEdge) edge).getNodeDirection((BDNode)curNode)!=null)
							 curDir=!((BDEdge) edge).getNodeDirection((BDNode)curNode);
						 ways = (curDir?curNode.leavingEdges():curNode.enteringEdges()).filter(e->!e.hasAttribute("cut")).collect(Collectors.toList());
					 }
				 }
				 
			 }
			 //now we have repPath
			 if(lastTime){
				 annotation = new JapsaAnnotation();
			 }
			 Sequence seq=repPath.spelling(annotation);
			 
			 double cov=GraphUtil.getRealCoverage(repPath.averageCov());
			 Node n=outputGraph.addNode(Integer.toString(comp.id));
			 seq.setName("Contig_"+comp.id+"_"+(isCircular?"circular":"linear")+"_length_"+seq.length()+"_cov_"+cov);
			 n.setAttribute("seq", seq);
			 n.setAttribute("len", seq.length());
			 n.setAttribute("cov",cov);
			 n.setAttribute("path", repPath);
			 
			 b=SimpleBinner.getBinOfPath(repPath);
			 if(b!=null) {
				 n.setAttribute("bin", b);
				 if(!binn50.containsKey(b))
					 binn50.put(b, new ArrayList<Integer>());
				 else
					 binn50.get(b).add(seq.length());
			 }
			 
			 
			 if(isCircular){
				 n.setAttribute("circular");
			 }
			 if(lastTime){
				 annotation.setSequence(seq);
				 n.setAttribute("annotation", annotation);
			 }
//			 System.out.println("\n" + seq.getName() + ":" + repPath.getId() + "\n=> "+ repPath.getPrimitivePath().getId());

		}


		//now set the edge of outputGraph based on the cut edges
		for(Edge e:cutEdges) {
			Node n0=e.getNode0(), n1=e.getNode1();
			//get corresponding grouped nodes in outputGraph
			ConnectedComponent 	comp0=rtComponents.getConnectedComponentOf(n0),
								comp1=rtComponents.getConnectedComponentOf(n1);
			if(comp0==null || comp1==null)
				continue;
			Node 	nn0=outputGraph.getNode(Integer.toString(comp0.id)),
					nn1=outputGraph.getNode(Integer.toString(comp1.id));
			if(nn0!=null && nn1!=null) {
				boolean d0=((BDEdge)e).getDir0(),
						d1=((BDEdge)e).getDir1();
				//If it consists of a path, should be the direction of the whole path, not a particular node anymore!
				if(((BDPath)nn0.getAttribute("path")).getNodeCount()>1) 
					d0=(n0==((BDPath)nn0.getAttribute("path")).peekNode())?true:false; 
				
				if(((BDPath)nn1.getAttribute("path")).getNodeCount()>1) 
					d1=(n1==((BDPath)nn1.getAttribute("path")).getRoot())?false:true;	
				
				if(e.hasAttribute("path")){
					BDPath path = (BDPath)e.getAttribute("path");
					BDPath trimedPath=path.trimEndingNodes();
					if(trimedPath!=null){
						//Add the "middle" node
						if(lastTime)
							annotation = new JapsaAnnotation();
						Sequence seq=trimedPath.spelling(annotation);
						double cov=GraphUtil.getRealCoverage(trimedPath.averageCov());
						int id=0;
						while(outputGraph.getNode(Integer.toString(++id))!=null);
						Node n=outputGraph.addNode(Integer.toString(id));
						seq.setName("Contig_"+id+"_linear_length_"+seq.length()+"_cov_"+cov);
						n.setAttribute("seq", seq);
						n.setAttribute("len", seq.length());
						n.setAttribute("cov",cov);
						n.setAttribute("path", trimedPath);
						b=SimpleBinner.getBinOfPath(repPath);
						if(b!=null) {
							n.setAttribute("bin", b);
							if(!binn50.containsKey(b))
								binn50.put(b, new ArrayList<Integer>());
							else
								binn50.get(b).add(seq.length());
						}
						 
						if(lastTime){
							annotation.setSequence(seq);
							n.setAttribute("annotation", annotation);
						}
						//Add 2 edges
						Boolean dd0=((BDEdge)path.getEdgePath().get(0)).getNodeDirection(trimedPath.getFirstNode()), 
								dd1=((BDEdge)path.peekEdge()).getNodeDirection(trimedPath.getLastNode());
						if(dd0==null)
							dd0=!path.getFirstNodeDirection();
						if(dd1==null)
							dd1=!path.getLastNodeDirection();
						 
						//if trimedPath has more than 1 nodes and has been merged into one
						if(trimedPath.getNodeCount()>1){
							if(trimedPath.getRoot()==path.getNodePath().get(1)){
								dd0=false;
								dd1=true;
							}else{
								dd0=true;
								dd1=false;
							}
						}
						outputGraph.addEdge((BDNode)nn0, (BDNode)n , d0, dd0);
						outputGraph.addEdge((BDNode)n, (BDNode)nn1 , dd1, d1);
					}
						
				}else{
					String id=BDEdge.createID((BDNode)nn0, (BDNode)nn1 , d0, d1);
					if(outputGraph.getEdge(id)==null)
						outputGraph.addEdge((BDNode)nn0, (BDNode)nn1 , d0, d1);
				}

			}
		}
		outputGraph.updateStats();
		System.out.printf("Input: read count=%d base count=%d\n", hAss.currentReadCount, hAss.currentBaseCount);
		System.out.printf("Output: %d sequences (%d circular) N50=%d N75=%d Max=%d ", getNumberOfSequences(), getNumberOfCircularSequences(), getN50(), getN75(), getLongestContig());
		for(PopBin bb:binn50.keySet()) {
			List<Integer> lengths=binn50.get(bb);
			
			if(lengths!=null && lengths.size()>0)
				System.out.printf("%s.N50=%d ", bb.getId(), GraphUtil.getNStats(.5, lengths.stream().mapToInt(i -> i).toArray()));
			
		}
		System.out.println();
		
//		for(PopBin bb:binn50.keySet()) {
//			List<Integer> lengths=binn50.get(bb);
//			System.out.printf("%d (%d): ", bb.getId(), lengths.size());
//			for(int i:lengths)
//				System.out.printf("%d, ",i);
//			System.out.println();			
//		}
			
		if(lastTime)
			System.out.println("FINISH!");
		//reset the cutting attributes
		hAss.simGraph.edges().filter(e->e.hasAttribute("cut")).forEach(e->{e.removeAttribute("cut");});

	}
	
	synchronized public int getN50() {
		return outputGraph==null?hAss.simGraph.n50:outputGraph.n50;
	}
	synchronized public int getN75() {
		return outputGraph==null?hAss.simGraph.n75:outputGraph.n75;
	}
	synchronized public int getLongestContig() {
		return outputGraph==null?hAss.simGraph.maxl:outputGraph.maxl;
	}
	synchronized public int getNumberOfSequences() {
		return outputGraph==null?hAss.simGraph.numOfCtgs:outputGraph.numOfCtgs;
	}
	synchronized public int getNumberOfCircularSequences() {
		return outputGraph==null?hAss.simGraph.numOfCircularCtgs:outputGraph.numOfCircularCtgs;
	}
	synchronized public void outputAssGFA(String fileName) throws IOException{
		outputGraph.outputGFA(fileName);
	}
	synchronized public void outputOrigGFA(String fileName) throws IOException{
		hAss.simGraph.outputGFA(fileName);
	}
	synchronized public void outputFASTA(String fileName) throws IOException {
		outputGraph.outputFASTA(fileName);
	}
	synchronized public void outputJAPSA(String fileName) throws IOException {
		outputGraph.outputJAPSA(fileName);
	}

	@Override
	protected void close() {
		// TODO Auto-generated method stub
		System.out.printf("Post-processing the graph by greedy path-finding algorithm. Please wait...\n");
		HashSet<GoInBetweenBridge> 		unsolved=hAss.simGraph.getUnsolvedBridges(),
										solved=new HashSet<>();
		while(true){
			boolean changed=false;
			for(GoInBetweenBridge brg:unsolved){
				if(HybridAssembler.VERBOSE)
					LOG.info("Last attempt on incomplete bridge {} : anchors={} \n {}", brg.getEndingsID(), brg.getNumberOfAnchors(), brg.getAllPossiblePaths());
				//Take the current best path among the candidate of a bridge and connect the bridge(greedy)
				if(brg.getCompletionLevel()>=3){ 
					hAss.simGraph.getNewSubPathsToReduce(brg.getBestPath(brg.pBridge.getNode0(),brg.pBridge.getNode1())).stream().forEach(p->hAss.simGraph.reduceUniquePath(p));
					solved.add(brg);
					changed=true;
				}else{
					brg.scanForAnEnd(true);	
					changed=brg.steps.connectBridgeSteps(true);
					
					//return appropriate path
					if(changed){
						hAss.simGraph.getNewSubPathsToReduce(brg.getBestPath(brg.steps.start.getNode(),brg.steps.end.getNode())).stream().forEach(p->hAss.simGraph.reduceUniquePath(p));
						solved.add(brg);
					}
					else if(HybridAssembler.VERBOSE)
						LOG.info("Last attempt failed \n");
				}
	
			}
			if(solved.isEmpty()&&!changed)
				break;
			else{
				unsolved.removeAll(solved);
				solved.clear();
			}
				
		}
        //update for the last time
        update(true);		
		try {
			outputFASTA(hAss.getPrefix()+File.separator+"npgraph_assembly.fasta");
			outputJAPSA(hAss.getPrefix()+File.separator+"npgraph_assembly.japsa");
			outputAssGFA(hAss.getPrefix()+File.separator+"npgraph_assembly.gfa");
			outputOrigGFA(hAss.getPrefix()+File.separator+"npgraph_components.gfa");
		}catch(IOException e) {
			LOG.error("Could not output assembly results: {}", e.getMessage());
		}
	}

	@Override
	protected void analysis() {
		// TODO Auto-generated method stub
		update(false);
	}

	@Override
	protected int getCurrentRead() {
		// TODO Auto-generated method stub
		return hAss.currentReadCount;
	}
}
