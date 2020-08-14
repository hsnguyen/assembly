package org.rtassembly.npgraph.grpc;

import java.util.ArrayList;
import org.rtassembly.npgraph.*;

import japsa.seq.PAFRecord;

public class AssemblyGuideUtil {
//	  public static ArrayList<Alignment> getAlignmentsFromRequest(RequestAssembly request, HybridAssembler myAss){
//		  ArrayList<Alignment> retval = new ArrayList<>();
//		  String refID;
//		  for(AlignmentMsg msg:request.getHitsListList()) {
//			  String refName=msg.getTargetName();
//			  refID = refName.split("_").length > 1 ? refName.split("_")[1]:refName;
//			  BDNode node = (BDNode) myAss.simGraph.getNode(refID);
//			  if(node==null)
//				  return retval;
//			  PAFRecord record = new PAFRecord(msg.getQueryName(), msg.getQueryLength(), msg.getQueryStart(), msg.getQueryEnd(), 
//					  							msg.getStrand(), 
//					  							msg.getTargetName(), msg.getTargetLength(), msg.getTargetStart(), msg.getTargetEnd(), 
//					  							msg.getScore(), msg.getQuality());
//			  retval.add(new Alignment(record, node));
//			  
//		  }
//		  
//		  return retval;
//	  }

}
