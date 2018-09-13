package org.rtassembly.npgraph;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;

public class BidirectedBridge {
	public static volatile int SAFE_VOTE_DISTANCE=3;
	ArrayList<Alignment> steps;
	ArrayList<BidirectedPath> 	fullPaths; //paths connect two ends (inducing by graph traversal)
	
	//these variables to store auxillary info
	ArrayList<BidirectedPath>	halfPaths; //paths from one end only (if any): from SPAdes
	ArrayList<BidirectedBridge> halfBridges; //bridges from one end only (if any)
	
	BidirectedBridge(Alignment start){
		steps=new ArrayList<>();
		steps.add(start);
		
		fullPaths = new ArrayList<>();
	}

	//This is for SPAdes path reader only
	BidirectedBridge (BidirectedPath path, boolean isFullPath){
		if(isFullPath){
			fullPaths = new ArrayList<>();
			fullPaths.add(path);
		}else{
			halfPaths=new ArrayList<>();
			halfPaths.add(path);
		}
	}
	

	public void addHalfBridge(BidirectedBridge brg){
		if(halfBridges==null)
			halfBridges=new ArrayList<>();
		halfBridges.add(brg);
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
	//Return bridge status:
	//-1: half bridge, 0: complete bridge, unsolved; 1: solved
	public int getBridgeStatus(){
		int retval=-1;
		if(fullPaths!=null && !fullPaths.isEmpty()){
			if(fullPaths.size()==1)
				retval=1;
			else
				retval=0;
		}
		return retval;
	}
	
	public Alignment getStartAlignment() {
		return steps.get(0);
	}
	
	public Alignment getEndAlignment() {
		return steps.get(steps.size()-1);
	}
	
	public BidirectedPath getBestPath() {
		//return best path (unique path or best voted path); or null if undetermined
		if(getBridgeStatus()==1)
			return fullPaths.get(0);
		else{
			//still a list of possible paths to deal with
			return null;
		}
	}
	
	public String getEndingsID() {
		if(steps!=null && steps.size() > 1){
			//because it really goes from start to end
			Alignment start = getStartAlignment(), end = getEndAlignment();
			return BidirectedEdge.createID(start.node, end.node, start.strand, !end.strand);
		}else{//SPAdes path
			BidirectedPath repPath=null;
			if(fullPaths!=null && !fullPaths.isEmpty())
				repPath=fullPaths.get(0);
			else
				repPath=halfPaths.get(0);
			BidirectedNode 	root = (BidirectedNode) repPath.getRoot(), 
							end = (BidirectedNode) repPath.peekNode();
			boolean	 	rootDir=((BidirectedEdge) repPath.getEdgePath().get(0)).getDir(root),
						endDir=((BidirectedEdge) repPath.peekEdge()).getDir(end);
			return BidirectedEdge.createID(root, end, rootDir, endDir);
		}
		
	}
	public String getBridgeString() {
		String retval="<unknown>";
		if(getBridgeStatus()==1)
			retval="solved<"+fullPaths.get(0).getId()+">";
		else if(steps!=null && steps.size()>1) {
			retval="unsolved<"+steps.stream().map(a->a.node.getId().concat(a.strand?"+":"-")).reduce("", (a,b) -> a + b)+">";
		}
		return retval;
	}
	
	//Using another list of steps to rectify the current bridge's steps
	//This bridge HAS TO BE complete(2 unique ends), 
	//the reference need not be complete, instead one unique end is enough (half bridge)
	public void referencingTo(BidirectedBridge brg){
		if(getBridgeStatus()!=0)
			return;
		else{
			System.out.printf("Check consistency using bridge %s:\n", brg.getBridgeString());
			//look for the common unique node between 2 bridges
			BidirectedNode startNode=this.getStartAlignment().node;
			if(startNode!=brg.getStartAlignment().node && startNode!=brg.getEndAlignment().node)
				startNode=this.getEndAlignment().node;
			
			for(BidirectedPath p:fullPaths){
				brg.compareAndVote(p, startNode);
				System.out.printf("...path %s now has %d votes!\n", p.getId(), p.getVote());
			}
			
			Collections.sort(fullPaths, (p1,p2)->p2.getVote()-p1.getVote());
			int currentMaxVote=fullPaths.get(0).getVote();
			//TODO: more formal
			fullPaths.removeIf(p->(currentMaxVote-p.getVote() > SAFE_VOTE_DISTANCE));
		}
	}
	//Merging with another half bridge, inheriting all info from it
	public void merging(BidirectedBridge brg){
		if(brg!=null){
			System.out.println("Merging half bridge " + brg.getEndingsID());
			if(brg.halfPaths!=null && !brg.halfPaths.isEmpty()){
				if(halfPaths==null)
					halfPaths=brg.halfPaths;
				else{
					for(BidirectedPath p:brg.halfPaths)
						if(!halfPaths.contains(p))
							halfPaths.add(p);
				}
			}
			
			halfBridges=new ArrayList<>();
			if(brg.steps!=null && steps.size()>2)
				halfBridges.add(brg);
			if(brg.halfBridges!=null){
				System.out.println("...taking half bridge");
				for(BidirectedBridge b:brg.halfBridges){
					if(b.steps==null||b.steps.size()<3)
						continue;
					halfBridges.add(b);
					System.out.println("\t" + b.getBridgeString());
				}
			}

//			brg=null;
			
			flushInfo();
		}
	}

	//Start to build bridge (i.e. identify possible paths) when 2 unique ends are determined
	public void bridging(BidirectedGraph graph, PopBin bin){
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

			stepPaths = graph.getClosestPaths(curAlignment, nextAlignment);
			if(stepPaths==null || stepPaths.isEmpty())
				return;//must follow the steps

	
			//join all paths from previous to the new ones
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
		
		if(!wholePaths.isEmpty()){
			fullPaths=wholePaths;
			//scan for unexpected unique node in the bridge (missed by alignments) and update the map of bridges
			wholePaths.forEach(p->{
									p.nodes().filter(n-> (n!=p.getRoot() && n!=p.peekNode()))
											.forEach(n->{
															if(SimpleBinner.getUniqueBin(n)!=null) 
//																graph.updateBridgesMap(n, this)
																;
															}
											);
									}
			);
		}
		fullPaths.forEach(p->p.setConsensusUniqueBinOfPath(bin));
		flushInfo();

		
	}
		
	public String getAllPossiblePathsString(){
		String retval="";
		if(fullPaths==null || fullPaths.isEmpty()){
			retval="possible paths not available!";
		}else{
			retval+=fullPaths.stream().map(p->p.getId()).reduce((a,b)->a+"\n"+b);
		}
		
		return retval;
	}
	
	private void flushInfo(){
		if(getBridgeStatus()==0){
			if(halfPaths!=null && !halfPaths.isEmpty()){
				for(BidirectedPath halfPath:halfPaths){
					System.out.println("bridging with half path " + halfPath.getId());
	
					fullPaths.removeIf(p->(	p.getId().indexOf(halfPath.getId()) < 0 
											&& p.getId().indexOf(halfPath.getReversedComplemented().getId()) < 0));
				}
			}
			
			if(halfBridges!=null && !halfBridges.isEmpty()){
				System.out.println("bridging with half bridges " + halfBridges);
				halfBridges.forEach(brg->referencingTo(brg));
				halfBridges=null;
			}
		}
	}
	
	//check if the steps case agree with a unique path or not
	//must start with an unique Node (startdingNode)
	private void compareAndVote(BidirectedPath path, BidirectedNode startingNode){
		if(steps.size()<2)
			return;
		//1. first check the unique starting point, reverse it if necessary
		if(path.getRoot()!=startingNode && path.peekNode() != startingNode){
			System.err.println("Couldn't find " + startingNode.getId() + " from both ends of path " + path.getId());
			return;
		}
		
		boolean dir=true; //true if bridge starts from left, false if from right
		Alignment root = getStartAlignment();
		if(startingNode==getEndAlignment().node){
			root=getEndAlignment();
			dir=false;
		}
		else if(startingNode!=getStartAlignment().node){
			System.err.println("Couldn't find " + startingNode.getId() + " from both ends of bridge " + this.getBridgeString());
			return;			
		}
		//2. then check the distances of nodes to the starting node: 
		//need a function to calculate the distance from the common unique node in the path!
		if(dir){
			for(int i=1; i < steps.size(); i++){
				Alignment curAlg = steps.get(i);
				int distance = curAlg.readAlignmentStart() - root.readAlignmentEnd()+1;
				boolean direction=root.strand == curAlg.strand;
				if(!path.checkDistanceConsistency(startingNode, curAlg.node, direction, distance)){
					path.downVote();
				}else
					path.upVote();
			}
		}else{
			for(int i=steps.size()-2; i>=0 ; i--){
				Alignment curAlg = steps.get(i);
				int distance = root.readAlignmentStart() - curAlg.readAlignmentEnd()+1;
				boolean direction=root.strand == curAlg.strand;
				if(!path.checkDistanceConsistency(startingNode, curAlg.node, direction, distance)){
					path.downVote();
				}else
					path.upVote();
			}
		}		
		
	}
	//check if a bridge is worth to merge to this bridge
	public boolean checkIfMerge(BidirectedBridge brg) {
		boolean retval=false;
		if(fullPaths.size() >= 2) {
			HashSet<String> pathsIntersect = 
				fullPaths.stream()
				.map(p->{return new HashSet<String>(Arrays.asList(p.getId().split("[\\+\\-,]")));})
				.reduce((a,b)->{
								HashSet<String> intersect = new HashSet<>(a);
								intersect.remove("");
								intersect.retainAll(b);
								return intersect;
								})
				.get();
			
			HashSet<String> brgStepSet = new HashSet<String>(Arrays.asList(brg.getBridgeString().split("[\\+\\-,]")));
			brgStepSet.remove("");
			int oldSize=pathsIntersect.size();
			pathsIntersect.addAll(brgStepSet);//union

			if(pathsIntersect.size()>oldSize)
				retval=true;
			else
				retval=false;
		}
		
		return retval;
	}
}
