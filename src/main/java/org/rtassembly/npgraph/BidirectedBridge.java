package org.rtassembly.npgraph;

import java.util.ArrayList;

import org.jfree.util.Log;

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
	public String getBridgeString() {
		String retval="";
		if(steps.size()>1) {
			retval=steps.stream().map(a->a.node.getId().concat(a.strand?"+":"-")).reduce("", (a,b) -> a + b);;
		}
		return retval;
	}
	
	//Using another list of steps to rectify the current bridge's steps
	public void merging(BidirectedBridge brg, BidirectedNode startNode){
		if(isSolved)
			return;
		else{
			ArrayList<BidirectedPath> tobeRemoved = new ArrayList<BidirectedPath>();
			for(BidirectedPath p:paths){
				boolean consistencyChecking=brg.agreeWith(p, startNode);
				System.out.printf("Checking consistency of bridge %s to path %s: %b\n", getBridgeString(), p.getId(), consistencyChecking);
				if(consistencyChecking) {
					p.elected();	
					System.out.printf("...path %s now has %d votes!", p.getId(), p.getVote());
				}
				else
					tobeRemoved.add(p);
			}
			if(tobeRemoved.size()==paths.size()-1)
				isSolved=true;
			
			for(BidirectedPath p:tobeRemoved)
				paths.remove(p);
			
		}
	}

	public void bridging(BidirectedGraph graph){
		if(steps.size()<=1) {
			return;		
		}
		Alignment 	curAlignment = getStartAlignment(),
					nextAlignment;

		ArrayList<BidirectedPath> 	wholePaths=new ArrayList<>(),
									tmpPaths=new ArrayList<>(),
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

			if(wholePaths.isEmpty())
				wholePaths=stepPaths;
			else{
				for(BidirectedPath curPath:wholePaths){
					for(BidirectedPath stepPath:stepPaths){
						BidirectedPath newPath=curPath.join(stepPath);
						if(newPath==null)
							return;
						else{
							tmpPaths.add(newPath);
						}
					}
				}
				wholePaths=tmpPaths;
				tmpPaths=new ArrayList<>();
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
	
	//check if the steps case agree with a unique path or not
	//must start with an unique Node (startdingNode)
	private boolean agreeWith(BidirectedPath path, BidirectedNode startingNode){
		if(steps.size()<2)
			return false;
		boolean retval=true;
		//1. first check the unique starting point, reverse it if necessary
		if(path.getRoot()!=startingNode && path.peekNode() != startingNode){
			System.err.println("Couldn't find " + startingNode.getId() + " from both ends of path " + path.getId());
			return false;
		}
		
		boolean dir=true; //true if bridge starts from left, false if from right
		Alignment root = getStartAlignment();
		if(startingNode==getEndAlignment().node){
			root=getEndAlignment();
			dir=false;
		}
		else if(startingNode!=getStartAlignment().node){
			System.err.println("Couldn't find " + startingNode.getId() + " from both ends of bridge " + this.getBridgeString());
			return false;			
		}
		//2. then check the distances of nodes to the starting node: 
		//need a function to calculate the distance from the common unique node in the path!
		//TODO: check for consistency of the direction between the two also
		if(dir){
			for(int i=1; i < steps.size(); i++){
				Alignment curAlg = steps.get(i);
				int distance = curAlg.readAlignmentStart() - root.readAlignmentEnd()+1;
				boolean direction=root.strand == curAlg.strand;
				if(!path.checkDistanceConsistency(startingNode, curAlg.node, direction, distance)){
					
					return false;
				}
			}
		}else{
			for(int i=steps.size()-2; i>=0 ; i--){
				Alignment curAlg = steps.get(i);
				int distance = root.readAlignmentStart() - curAlg.readAlignmentEnd()+1;
				boolean direction=root.strand == curAlg.strand;
				if(!path.checkDistanceConsistency(startingNode, curAlg.node, direction, distance)){
					
					return false;
				}
			}
		}		
		
		return retval;
	}
}
