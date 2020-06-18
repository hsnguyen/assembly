package org.rtassembly.npgraph.grpc;

import io.grpc.Server;
import io.grpc.ServerBuilder;
import io.grpc.stub.StreamObserver;
import japsa.seq.PAFRecord;
import japsa.seq.Sequence;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.rtassembly.npgraph.Alignment;
import org.rtassembly.npgraph.BDGraph;
import org.rtassembly.npgraph.BDNode;
import org.rtassembly.npgraph.BDPath;
import org.rtassembly.npgraph.GoInBetweenBridge;
import org.rtassembly.npgraph.GraphUtil;
import org.rtassembly.npgraph.HybridAssembler;
import org.rtassembly.npgraph.SimpleBinner;

public class AssemblyGuideServer {
    private static final Logger logger = LogManager.getLogger(MethodHandles.lookup().lookupClass());

	  public int port;
	  public final Server server;
	  public static HybridAssembler myAss;
	  public static int ELEN = 10000;
	  public final Thread myThread;
	  
	  public AssemblyGuideServer(int port, HybridAssembler ass) {
		  this(ServerBuilder.forPort(port), port, ass);
	  }
	  /** Create a RouteGuide server using serverBuilder as a base */
	  public AssemblyGuideServer(ServerBuilder<?> serverBuilder, int port, HybridAssembler badAss) {
		  this.port = port;
		  server = serverBuilder.addService(new AssemblyGuideService())
				  				.build();
		  myAss=badAss;
		  myThread = new Thread(badAss.observer);
	  }

	  /** Start serving requests. */
	  public void start() throws IOException {
		  myThread.start();
		  server.start();
		  
		  logger.info("Server started, listening on {}", port);
		  Runtime.getRuntime().addShutdownHook(new Thread() {
			  @Override
			  public void run() {
				  // Use stderr here since the logger may have been reset by its JVM shutdown hook.
				  System.err.println("*** shutting down gRPC server since JVM is shutting down");
				  try {
					  AssemblyGuideServer.this.stop();
				  } catch (InterruptedException e) {
					  e.printStackTrace(System.err);
				  }
				  System.err.println("*** server shut down");
			  }
		  });
	  }

	  /** Stop serving requests and shutdown resources. */
	  public void stop() throws InterruptedException {
		  if (server != null) {
			  server.shutdown().awaitTermination(30, TimeUnit.SECONDS);
		  }
		  if(myThread.isAlive()) {
				myAss.observer.stopWaiting();
				myThread.join();
		  }
	  }

	  /**
	   * Await termination on the main thread since the grpc library uses daemon threads.
	   */
	  private void blockUntilShutdown() throws InterruptedException {
		  if (server != null) {
			  server.awaitTermination();
		  }
//		  if(myThread.isAlive()) {
//				myAss.observer.stopWaiting();
//				myThread.join();
//		  }
	  }


	  /**
	   * Our implementation of RouteGuide service.
	   *
	   * <p>See route_guide.proto for details of the methods.
	   */
	  private static class AssemblyGuideService extends AssemblyGuideGrpc.AssemblyGuideImplBase {
		  private static final HashMap<String, Alignment> lastMap = new HashMap<>(); //readID to the last unique contig it mapped to
		  private static final HashMap<String, Integer> reduceRead = new HashMap<>(); //readID to the length of chunk used to reduce

		  @Override
		  public void getAssemblyContribution(RequestAssembly request,
				StreamObserver<ResponseAssembly> responseObserver) {
			  	ArrayList<Alignment> hits = getAlignmentsFromRequest(request);
			  	boolean continueing=true;
			  	if(hits.size() > 0) {//this should be checked from client side
			  		Collections.sort(hits);
			  		Alignment 	a = hits.get(hits.size()-1), //get the last alignment of current hit list
			  					b = lastMap.get(a.readID); //previous last unique alignment
			  		//1. compare the last mapped contig, proceed if the same
			  		//TODO: index repeats of long bridges that already resolved?
			  		
			  		//with read chunk shorter than 1000bp, mapping to unique contig cannot be missed 
			  		//when sequentially considering only last Alignment
			  		if(SimpleBinner.getBinIfUniqueNow(a.node)!=null) {
			  			//if the same as prevAlg, proceed
			  			boolean investigating=true;
			  			if(b==null || a.node != b.node) {
			  				lastMap.put(a.readID, a);
			  			}else{
			  				int alignedReadLen = Math.abs(a.readEnd - a.readStart) + Math.abs(b.readEnd - b.readStart),
			  					alignedRefLen = Math.abs(a.refEnd - a.refStart) + Math.abs(b.refEnd - b.refStart);
			  				double rate = 1.0 * alignedRefLen/alignedReadLen;		

			  				int alignP = (int) ((b.readStart - a.readStart) * rate);
			  				//(rough) relative position from ref_b (contig of b) to ref_a (contig of a) in the assembled genome
			  				int gP = Math.abs((alignP + (a.strand ? a.refStart:-a.refStart) - (b.strand?b.refStart:-b.refStart)));
			  				if(gP+BDGraph.A_TOL > a.node.getNumber("len")) { //circular
			  					lastMap.put(a.readID, a);
			  					investigating=continueing=false; //terminate sequencing this one
			  				}
			  			}
			  			
			  			//3. estimate distance to the end of this unique contig to calculate usefulness
			  			if(investigating) {
			  				int eLen = a.readAlignmentEnd();
			  				eLen+=(a.strand?a.node.getNumber("len")-a.refEnd:a.refStart);
			  				
			  				BDNode prevNode, unqNode;
			  				prevNode=unqNode=a.node;
			  				boolean dir=!a.strand;
			  				GoInBetweenBridge brg=BDGraph.bridgesMap.get(unqNode.getId()+(dir?"o":"i"));
			  				while(brg!=null&&brg.getCompletionLevel()==4) {
			  					
			  					if(unqNode==brg.pBridge.getNode0() && dir==brg.pBridge.getDir0()) {
			  						unqNode=brg.pBridge.getNode1();
			  						dir=!brg.pBridge.getDir1();
			  					}else {
			  						unqNode=brg.pBridge.getNode0();
			  						dir=!brg.pBridge.getDir0();
			  					}
			  					eLen+=brg.steps.getSpanVector().distance(prevNode, unqNode);
			  					eLen+=unqNode.getNumber("len");
			  					
			  					prevNode=unqNode;
			  					brg=BDGraph.bridgesMap.get(unqNode.getId()+(dir?"o":"i"));
			  				}
			  				
			  				continueing=(eLen < ELEN);//??too simple!
			  			}
			  			
//					  	//5. send control signal back to client
//					  	responseObserver.onNext(ResponseAssembly.newBuilder().setUsefulness(continueing).setReadId(request.getReadId()).build());
//					  	responseObserver.onCompleted();
					  	
				  		//4. reduce
					  	Sequence read = GraphUtil.getNSequence(a.readID, a.readLength); 	

						List<BDPath> paths=myAss.simGraph.uniqueBridgesFinding(read, hits);
						if(paths!=null) {
						    paths.stream().forEach(p->myAss.simGraph.reduceUniquePath(p));
							if(reduceRead.containsKey(a.readID)) {//already used for reduction
								myAss.currentBaseCount += a.readLength-reduceRead.get(a.readID);
								reduceRead.replace(a.readID, a.readLength);
							}else {//new
								reduceRead.put(a.readID, a.readLength);
								myAss.currentReadCount++;
								myAss.currentBaseCount+=a.readLength;
							}
			
						}	
						
//						return;
			  		}
			  			

			  	}
			  	
			  	//5. send control signal back to client
			  	responseObserver.onNext(ResponseAssembly.newBuilder().setUsefulness(continueing).setReadId(request.getReadId()).build());
			  	responseObserver.onCompleted();
		  }
		  
		  private ArrayList<Alignment> getAlignmentsFromRequest(RequestAssembly request){
			  ArrayList<Alignment> retval = new ArrayList<>();
			  String refID;
			  for(AlignmentMsg msg:request.getHitsListList()) {
				  String refName=msg.getTargetName();
				  refID = refName.split("_").length > 1 ? refName.split("_")[1]:refName;
				  BDNode node = (BDNode) myAss.simGraph.getNode(refID);
				  if(node==null)
					  return retval;
				  //Convert the hit message to a PAFRecord: [0-based inclusive; 0-based exlusive] -> [1-based inclusive; 1-based inclusive]
				  PAFRecord record = new PAFRecord(msg.getQueryName(), msg.getQueryLength(), msg.getQueryStart()+1, msg.getQueryEnd(), 
						  							msg.getStrand(), 
						  							msg.getTargetName(), msg.getTargetLength(), msg.getTargetStart()+1, msg.getTargetEnd(), 
						  							msg.getScore(), msg.getQuality());
				  retval.add(new Alignment(record, node));
				  
			  }
			  
			  return retval;
		  }
	  }	  
}
