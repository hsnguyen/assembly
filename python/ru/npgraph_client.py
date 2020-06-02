import grpc
import npgraph_service_pb2
import npgraph_service_pb2_grpc
import mappy as mp
import logging
import time

def run():

    a = mp.Aligner("/home/sonhoanghguyen/Projects/readuntil/simulation/npgraph/test/assembly_graph.fasta")  # load or build index
    if not a: raise Exception("ERROR: failed to load/build index")


    with grpc.insecure_channel('localhost:2105') as channel:
        stub = npgraph_service_pb2_grpc.AssemblyGuideStub(channel)
        print("Connected with server at localhost:2105")

        for name, seq, qual in mp.fastx_read("/home/sonhoanghguyen/Projects/readuntil/simulation/npgraph/test/E_coli_K-12_MG1655_good_long.fastq.gz"):
            #1. make request
            request = npgraph_service_pb2.RequestAssembly()
            request.read_id = name
            for hit in a.map(seq): # traverse alignments
                #print("{}\t{}\t{}\t{}".format(hit.ctg, hit.r_st, hit.r_en, hit.cigar_str))
                request.hits_list.append(npgraph_service_pb2.AlignmentMsg(query_name=name,query_length=len(seq),query_start=hit.q_st,query_end=hit.q_en,strand=hit.strand>0,target_name=hit.ctg,target_length=hit.ctg_len,target_start=hit.r_st,target_end=hit.r_en,quality=hit.mapq,score=hit.mlen))

            #2. get and print response
            if len(request.hits_list)>0:
                try:
                    start_time = time.time()
                    response = stub.GetAssemblyContribution(request)
                    print("{}: {} in {:.5f} seconds".format(response.read_id, response.usefulness, time.time()-start_time))
                except grpc.RpcError as e:
                    print("{}: errorcode={}".format(request.read_id, str(e.code())))
                    continue
            else:
                print("{}: unmapped!".format(request.read_id))
                continue

if __name__ == "__main__":
    logging.basicConfig()
    run()
