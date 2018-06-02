package org.rtassembly.npgraph;

import java.util.ArrayList;

public class BidirectedBridge {
	ArrayList<Alignment> steps;
	BidirectedPath path;
	
	BidirectedBridge(Alignment start){
		steps=new ArrayList<>();
		steps.add(start);
		path=null;
	}

	
	public boolean append(Alignment alg) {
		Alignment last=getEndAlignment();
		if(	alg==null || last.readID.compareTo(alg.readID)!=0 || 
			Math.max(last.readStart, last.readEnd) >= Math.max(alg.readStart, alg.readEnd)) {
			System.err.println("Illegal alignment being added, skip!");
			return false;
		}else {
			steps.add(alg);
			return true;
		}
	}
//	public boolean appendAll(BidirectedBridge brg){
//		if(brg.getStartAlignment()!=getEndAlignment())
//			return false;
//		else{
//			for(int i=1;i<brg.steps.size();i++){
//				append(brg.steps.get(i));
//			}
//			brg=null;
//			return true;
//		}
//	}
	public Alignment getStartAlignment() {
		return steps.get(0);
	}
	
	public Alignment getEndAlignment() {
		return steps.get(steps.size()-1);
	}
	
	public BidirectedPath getPath() {
		return path;
	}
	public String getEndingsID() {
		if(steps.size()<=1)
			return null;
		else {
			//because it really goes from start to end
			Alignment start = getStartAlignment(), end = getEndAlignment();
			return BidirectedEdge.createID(start.node, end.node, start.strand, !end.strand);
		}
	}
	

	public void bridging(BidirectedGraph graph){
		if(steps.size()<=1) {
			return;		
		}
		Alignment 	curAlignment = getStartAlignment(),
					nextAlignment;

		BidirectedPath curPath=null;
		for(int i=1; i<steps.size();i++){
			nextAlignment = steps.get(i);
			BidirectedPath curBridge=null;			
	
			int distance = nextAlignment.readAlignmentStart()-curAlignment.readAlignmentEnd();
			if(distance<BidirectedGraph.D_LIMIT){
				ArrayList<BidirectedPath> bridges = graph.getClosestPath(curAlignment, nextAlignment, distance);
				if(bridges!=null) {
					//TODO: find the best path representing this bridge here:
					curBridge = bridges.get(0);
				}
			}else 
				continue;
	
			//join all paths from previous to the new ones
			//TODO:optimize it
			System.out.println("====Before curPath: " + curPath + "; curBridge:" + curBridge);
	
			if(curPath==null)
				curPath=curBridge;
			else{
				if(!curPath.join(curBridge)) {
					return;
				}
			}
			curAlignment=nextAlignment;
			
			System.out.println("====After curPath: " + curPath + "; curBridge:" + curBridge);
		}
		
		if(curPath!=null)
			path=curPath;
		
	}
		
		//should we overwrite compareTo() to compare 2 bridges???
}
