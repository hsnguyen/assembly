package org.rtassembly.npgraph.grpc;

import io.grpc.Server;
import io.grpc.ServerBuilder;
import io.grpc.stub.StreamObserver;
import japsa.seq.PAFRecord;
import japsa.seq.Sequence;

import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.LogManager;
import org.rtassembly.npgraph.Alignment;
import org.rtassembly.npgraph.BDNode;
import org.rtassembly.npgraph.BDPath;
import org.rtassembly.npgraph.GraphUtil;
import org.rtassembly.npgraph.HybridAssembler;

public class AssemblyGuideServer {
    private static final Logger logger = LogManager.getLogger(MethodHandles.lookup().lookupClass());

	  private final int port;
	  private final Server server;
	  private final HybridAssembler myAss;
	  private final Thread myThread;
	  
	  public AssemblyGuideServer(int port, HybridAssembler ass) {
		  this(ServerBuilder.forPort(port), port, ass);
	  }
	  /** Create a RouteGuide server using serverBuilder as a base */
	  public AssemblyGuideServer(ServerBuilder<?> serverBuilder, int port, HybridAssembler badAss) {
		  this.port = port;
		  server = serverBuilder.addService(new AssemblyGuideService(badAss))
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

		  private final HybridAssembler myAss; //redundant?
		  AssemblyGuideService(HybridAssembler badAss){ myAss=badAss;}
		  @Override
		  public void getAssemblyContribution(RequestAssembly request,
				StreamObserver<ResponseAssembly> responseObserver) {
			  	ArrayList<Alignment> hits = getAlignmentsFromRequest(request);
			  	boolean usefulness=false;
			  	if(!hits.isEmpty()) {
			  		Alignment alg = hits.get(0);
				  	Sequence read = GraphUtil.getNSequence(alg.readID, alg.readLength);
				  	

					myAss.currentReadCount ++;
					myAss.currentBaseCount += read.length();
	
					List<BDPath> paths=myAss.simGraph.uniqueBridgesFinding(read, hits);
					if(paths!=null) {
						usefulness=true;
					    paths.stream().forEach(p->myAss.simGraph.reduceUniquePath(p));
					}
			  	}
			  	responseObserver.onNext(ResponseAssembly.newBuilder().setUsefulness(usefulness).setReadId(request.getReadId()).build());
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
