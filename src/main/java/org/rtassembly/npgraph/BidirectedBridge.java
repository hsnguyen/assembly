package org.rtassembly.npgraph;

import java.util.ArrayList;

public class BidirectedBridge {
	ArrayList<Alignment> steps;
	ArrayList<BidirectedPath> paths;
	boolean isSolved=false;
	
	BidirectedBridge(Alignment start){
		steps=new ArrayList<>();
		steps.add(start);
		
		paths = new ArrayList<>();
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
	
	public BidirectedPath getBestPath() {
		//return best path (unique path or best voted path); or null if undetermined
		if(isSolved)
			return paths.get(0);
		else{
			//still a list of possible paths to deal with
			return null;
		}
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
	//Using another list of steps to rectify the current bridge's steps
	public void merging(BidirectedBridge brg){
		if(isSolved)
			return;
		else{
			for(BidirectedPath p:paths)
				if(brg.agreeWith(p))
					p.elected();										
		}
		//check if one path has dominant vote -> solved!
//		if()
//			isSolved=true;
	}

	public void bridging(BidirectedGraph graph){
		if(steps.size()<=1) {
			return;		
		}
		Alignment 	curAlignment = getStartAlignment(),
					nextAlignment;

		ArrayList<BidirectedPath> 	wholePaths=new ArrayList<>(),
									stepPaths=new ArrayList<>();
		for(int i=1; i<steps.size();i++){
			nextAlignment = steps.get(i);
			
			int distance = nextAlignment.readAlignmentStart()-curAlignment.readAlignmentEnd();
			if(distance<BidirectedGraph.D_LIMIT){
				stepPaths = graph.getClosestPaths(curAlignment, nextAlignment, distance);
				if(stepPaths==null || stepPaths.isEmpty())
					continue;
			}else 
				return;
	
			//join all paths from previous to the new ones
			//TODO:optimize it

			if(wholePaths.size()==0)
				wholePaths=stepPaths;
			else{
				for(BidirectedPath curPath:wholePaths)
					for(BidirectedPath stepPath:stepPaths)
						if(!curPath.join(stepPath)) {
							return;
						}
			}
			curAlignment=nextAlignment;
			
		}
		
		if(!wholePaths.isEmpty())
			paths=wholePaths;
		if(wholePaths.size()==1)
			isSolved=true;
		
	}
		
	public String getAllPossiblePathsString(){
		String retval="";
		if(paths==null || paths.isEmpty()){
			retval="possible paths not available!";
		}else{
			retval+=paths.stream().map(p->p.getId()).reduce((a,b)->a+"\n"+b);
		}
		
		return retval;
	}
	
	//check if the steps case agree with a path or not
	private boolean agreeWith(BidirectedPath path){
		boolean retval=false;
		//1. first check the unique starting point, reverse it if necessary
		
		//2. then check the distances of nodes to the starting node 
		
		
		return retval;
	}
}
