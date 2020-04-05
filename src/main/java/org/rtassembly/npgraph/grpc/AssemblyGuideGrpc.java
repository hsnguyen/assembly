package org.rtassembly.npgraph.grpc;

import static io.grpc.MethodDescriptor.generateFullMethodName;
import static io.grpc.stub.ClientCalls.asyncBidiStreamingCall;
import static io.grpc.stub.ClientCalls.asyncClientStreamingCall;
import static io.grpc.stub.ClientCalls.asyncServerStreamingCall;
import static io.grpc.stub.ClientCalls.asyncUnaryCall;
import static io.grpc.stub.ClientCalls.blockingServerStreamingCall;
import static io.grpc.stub.ClientCalls.blockingUnaryCall;
import static io.grpc.stub.ClientCalls.futureUnaryCall;
import static io.grpc.stub.ServerCalls.asyncBidiStreamingCall;
import static io.grpc.stub.ServerCalls.asyncClientStreamingCall;
import static io.grpc.stub.ServerCalls.asyncServerStreamingCall;
import static io.grpc.stub.ServerCalls.asyncUnaryCall;
import static io.grpc.stub.ServerCalls.asyncUnimplementedStreamingCall;
import static io.grpc.stub.ServerCalls.asyncUnimplementedUnaryCall;

/**
 * <pre>
 * Interface exported by the server.
 * </pre>
 */
@javax.annotation.Generated(
    value = "by gRPC proto compiler (version 1.28.0)",
    comments = "Source: npgraph_service.proto")
public final class AssemblyGuideGrpc {

  private AssemblyGuideGrpc() {}

  public static final String SERVICE_NAME = "assembly.AssemblyGuide";

  // Static method descriptors that strictly reflect the proto.
  private static volatile io.grpc.MethodDescriptor<org.rtassembly.npgraph.grpc.RequestAssembly,
      org.rtassembly.npgraph.grpc.ResponseAssembly> getGetAssemblyContributionMethod;

  @io.grpc.stub.annotations.RpcMethod(
      fullMethodName = SERVICE_NAME + '/' + "GetAssemblyContribution",
      requestType = org.rtassembly.npgraph.grpc.RequestAssembly.class,
      responseType = org.rtassembly.npgraph.grpc.ResponseAssembly.class,
      methodType = io.grpc.MethodDescriptor.MethodType.UNARY)
  public static io.grpc.MethodDescriptor<org.rtassembly.npgraph.grpc.RequestAssembly,
      org.rtassembly.npgraph.grpc.ResponseAssembly> getGetAssemblyContributionMethod() {
    io.grpc.MethodDescriptor<org.rtassembly.npgraph.grpc.RequestAssembly, org.rtassembly.npgraph.grpc.ResponseAssembly> getGetAssemblyContributionMethod;
    if ((getGetAssemblyContributionMethod = AssemblyGuideGrpc.getGetAssemblyContributionMethod) == null) {
      synchronized (AssemblyGuideGrpc.class) {
        if ((getGetAssemblyContributionMethod = AssemblyGuideGrpc.getGetAssemblyContributionMethod) == null) {
          AssemblyGuideGrpc.getGetAssemblyContributionMethod = getGetAssemblyContributionMethod =
              io.grpc.MethodDescriptor.<org.rtassembly.npgraph.grpc.RequestAssembly, org.rtassembly.npgraph.grpc.ResponseAssembly>newBuilder()
              .setType(io.grpc.MethodDescriptor.MethodType.UNARY)
              .setFullMethodName(generateFullMethodName(SERVICE_NAME, "GetAssemblyContribution"))
              .setSampledToLocalTracing(true)
              .setRequestMarshaller(io.grpc.protobuf.ProtoUtils.marshaller(
                  org.rtassembly.npgraph.grpc.RequestAssembly.getDefaultInstance()))
              .setResponseMarshaller(io.grpc.protobuf.ProtoUtils.marshaller(
                  org.rtassembly.npgraph.grpc.ResponseAssembly.getDefaultInstance()))
              .setSchemaDescriptor(new AssemblyGuideMethodDescriptorSupplier("GetAssemblyContribution"))
              .build();
        }
      }
    }
    return getGetAssemblyContributionMethod;
  }

  /**
   * Creates a new async stub that supports all call types for the service
   */
  public static AssemblyGuideStub newStub(io.grpc.Channel channel) {
    io.grpc.stub.AbstractStub.StubFactory<AssemblyGuideStub> factory =
      new io.grpc.stub.AbstractStub.StubFactory<AssemblyGuideStub>() {
        @java.lang.Override
        public AssemblyGuideStub newStub(io.grpc.Channel channel, io.grpc.CallOptions callOptions) {
          return new AssemblyGuideStub(channel, callOptions);
        }
      };
    return AssemblyGuideStub.newStub(factory, channel);
  }

  /**
   * Creates a new blocking-style stub that supports unary and streaming output calls on the service
   */
  public static AssemblyGuideBlockingStub newBlockingStub(
      io.grpc.Channel channel) {
    io.grpc.stub.AbstractStub.StubFactory<AssemblyGuideBlockingStub> factory =
      new io.grpc.stub.AbstractStub.StubFactory<AssemblyGuideBlockingStub>() {
        @java.lang.Override
        public AssemblyGuideBlockingStub newStub(io.grpc.Channel channel, io.grpc.CallOptions callOptions) {
          return new AssemblyGuideBlockingStub(channel, callOptions);
        }
      };
    return AssemblyGuideBlockingStub.newStub(factory, channel);
  }

  /**
   * Creates a new ListenableFuture-style stub that supports unary calls on the service
   */
  public static AssemblyGuideFutureStub newFutureStub(
      io.grpc.Channel channel) {
    io.grpc.stub.AbstractStub.StubFactory<AssemblyGuideFutureStub> factory =
      new io.grpc.stub.AbstractStub.StubFactory<AssemblyGuideFutureStub>() {
        @java.lang.Override
        public AssemblyGuideFutureStub newStub(io.grpc.Channel channel, io.grpc.CallOptions callOptions) {
          return new AssemblyGuideFutureStub(channel, callOptions);
        }
      };
    return AssemblyGuideFutureStub.newStub(factory, channel);
  }

  /**
   * <pre>
   * Interface exported by the server.
   * </pre>
   */
  public static abstract class AssemblyGuideImplBase implements io.grpc.BindableService {

    /**
     */
    public void getAssemblyContribution(org.rtassembly.npgraph.grpc.RequestAssembly request,
        io.grpc.stub.StreamObserver<org.rtassembly.npgraph.grpc.ResponseAssembly> responseObserver) {
      asyncUnimplementedUnaryCall(getGetAssemblyContributionMethod(), responseObserver);
    }

    @java.lang.Override public final io.grpc.ServerServiceDefinition bindService() {
      return io.grpc.ServerServiceDefinition.builder(getServiceDescriptor())
          .addMethod(
            getGetAssemblyContributionMethod(),
            asyncUnaryCall(
              new MethodHandlers<
                org.rtassembly.npgraph.grpc.RequestAssembly,
                org.rtassembly.npgraph.grpc.ResponseAssembly>(
                  this, METHODID_GET_ASSEMBLY_CONTRIBUTION)))
          .build();
    }
  }

  /**
   * <pre>
   * Interface exported by the server.
   * </pre>
   */
  public static final class AssemblyGuideStub extends io.grpc.stub.AbstractAsyncStub<AssemblyGuideStub> {
    private AssemblyGuideStub(
        io.grpc.Channel channel, io.grpc.CallOptions callOptions) {
      super(channel, callOptions);
    }

    @java.lang.Override
    protected AssemblyGuideStub build(
        io.grpc.Channel channel, io.grpc.CallOptions callOptions) {
      return new AssemblyGuideStub(channel, callOptions);
    }

    /**
     */
    public void getAssemblyContribution(org.rtassembly.npgraph.grpc.RequestAssembly request,
        io.grpc.stub.StreamObserver<org.rtassembly.npgraph.grpc.ResponseAssembly> responseObserver) {
      asyncUnaryCall(
          getChannel().newCall(getGetAssemblyContributionMethod(), getCallOptions()), request, responseObserver);
    }
  }

  /**
   * <pre>
   * Interface exported by the server.
   * </pre>
   */
  public static final class AssemblyGuideBlockingStub extends io.grpc.stub.AbstractBlockingStub<AssemblyGuideBlockingStub> {
    private AssemblyGuideBlockingStub(
        io.grpc.Channel channel, io.grpc.CallOptions callOptions) {
      super(channel, callOptions);
    }

    @java.lang.Override
    protected AssemblyGuideBlockingStub build(
        io.grpc.Channel channel, io.grpc.CallOptions callOptions) {
      return new AssemblyGuideBlockingStub(channel, callOptions);
    }

    /**
     */
    public org.rtassembly.npgraph.grpc.ResponseAssembly getAssemblyContribution(org.rtassembly.npgraph.grpc.RequestAssembly request) {
      return blockingUnaryCall(
          getChannel(), getGetAssemblyContributionMethod(), getCallOptions(), request);
    }
  }

  /**
   * <pre>
   * Interface exported by the server.
   * </pre>
   */
  public static final class AssemblyGuideFutureStub extends io.grpc.stub.AbstractFutureStub<AssemblyGuideFutureStub> {
    private AssemblyGuideFutureStub(
        io.grpc.Channel channel, io.grpc.CallOptions callOptions) {
      super(channel, callOptions);
    }

    @java.lang.Override
    protected AssemblyGuideFutureStub build(
        io.grpc.Channel channel, io.grpc.CallOptions callOptions) {
      return new AssemblyGuideFutureStub(channel, callOptions);
    }

    /**
     */
    public com.google.common.util.concurrent.ListenableFuture<org.rtassembly.npgraph.grpc.ResponseAssembly> getAssemblyContribution(
        org.rtassembly.npgraph.grpc.RequestAssembly request) {
      return futureUnaryCall(
          getChannel().newCall(getGetAssemblyContributionMethod(), getCallOptions()), request);
    }
  }

  private static final int METHODID_GET_ASSEMBLY_CONTRIBUTION = 0;

  private static final class MethodHandlers<Req, Resp> implements
      io.grpc.stub.ServerCalls.UnaryMethod<Req, Resp>,
      io.grpc.stub.ServerCalls.ServerStreamingMethod<Req, Resp>,
      io.grpc.stub.ServerCalls.ClientStreamingMethod<Req, Resp>,
      io.grpc.stub.ServerCalls.BidiStreamingMethod<Req, Resp> {
    private final AssemblyGuideImplBase serviceImpl;
    private final int methodId;

    MethodHandlers(AssemblyGuideImplBase serviceImpl, int methodId) {
      this.serviceImpl = serviceImpl;
      this.methodId = methodId;
    }

    @java.lang.Override
    @java.lang.SuppressWarnings("unchecked")
    public void invoke(Req request, io.grpc.stub.StreamObserver<Resp> responseObserver) {
      switch (methodId) {
        case METHODID_GET_ASSEMBLY_CONTRIBUTION:
          serviceImpl.getAssemblyContribution((org.rtassembly.npgraph.grpc.RequestAssembly) request,
              (io.grpc.stub.StreamObserver<org.rtassembly.npgraph.grpc.ResponseAssembly>) responseObserver);
          break;
        default:
          throw new AssertionError();
      }
    }

    @java.lang.Override
    @java.lang.SuppressWarnings("unchecked")
    public io.grpc.stub.StreamObserver<Req> invoke(
        io.grpc.stub.StreamObserver<Resp> responseObserver) {
      switch (methodId) {
        default:
          throw new AssertionError();
      }
    }
  }

  private static abstract class AssemblyGuideBaseDescriptorSupplier
      implements io.grpc.protobuf.ProtoFileDescriptorSupplier, io.grpc.protobuf.ProtoServiceDescriptorSupplier {
    AssemblyGuideBaseDescriptorSupplier() {}

    @java.lang.Override
    public com.google.protobuf.Descriptors.FileDescriptor getFileDescriptor() {
      return org.rtassembly.npgraph.grpc.AssemblyGuideProto.getDescriptor();
    }

    @java.lang.Override
    public com.google.protobuf.Descriptors.ServiceDescriptor getServiceDescriptor() {
      return getFileDescriptor().findServiceByName("AssemblyGuide");
    }
  }

  private static final class AssemblyGuideFileDescriptorSupplier
      extends AssemblyGuideBaseDescriptorSupplier {
    AssemblyGuideFileDescriptorSupplier() {}
  }

  private static final class AssemblyGuideMethodDescriptorSupplier
      extends AssemblyGuideBaseDescriptorSupplier
      implements io.grpc.protobuf.ProtoMethodDescriptorSupplier {
    private final String methodName;

    AssemblyGuideMethodDescriptorSupplier(String methodName) {
      this.methodName = methodName;
    }

    @java.lang.Override
    public com.google.protobuf.Descriptors.MethodDescriptor getMethodDescriptor() {
      return getServiceDescriptor().findMethodByName(methodName);
    }
  }

  private static volatile io.grpc.ServiceDescriptor serviceDescriptor;

  public static io.grpc.ServiceDescriptor getServiceDescriptor() {
    io.grpc.ServiceDescriptor result = serviceDescriptor;
    if (result == null) {
      synchronized (AssemblyGuideGrpc.class) {
        result = serviceDescriptor;
        if (result == null) {
          serviceDescriptor = result = io.grpc.ServiceDescriptor.newBuilder(SERVICE_NAME)
              .setSchemaDescriptor(new AssemblyGuideFileDescriptorSupplier())
              .addMethod(getGetAssemblyContributionMethod())
              .build();
        }
      }
    }
    return result;
  }
}
