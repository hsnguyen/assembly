package org.rtassembly.npgraph;

import java.io.File;
import java.lang.invoke.MethodHandles;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.log4j.Logger;

import japsa.bio.np.ErrorCorrection;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;

//Module for consensus bridging of contigs when they're not connected in the assembly graph
public class ConsensusCaller {
    private static final Logger logger = Logger.getLogger(MethodHandles.lookup().lookupClass().getSimpleName());
	
	volatile HashMap<String, List<Sequence>> bridgingReads;
	volatile HashMap<String, Sequence> consensusReads;
	private Set<String> connectedPairs; //save connected node pairs in the assembly graph so that we don't need to save the in-betwwen long-read
	String msa;
	ConsensusCaller(){
		consensusReads = new HashMap<String, Sequence>();
		bridgingReads = new HashMap<String, List<Sequence>>();
		connectedPairs = new HashSet<>();
		msa="kalign"; //default because of popularity, recommended one is spoa
	}
	ConsensusCaller(String msa){
		this();
		this.msa=msa;
	}
	public void setConsensusMSA(String aligner){
		aligner=aligner.toLowerCase().trim();
		if(aligner.startsWith("spoa"))
			msa="spoa";
		else if(aligner.startsWith("poa"))
			msa="poa";
		else if(aligner.startsWith("kalign3"))
			msa="kalign3";
		else if(aligner.startsWith("kalign"))
			msa="kalign";
		else
			msa="none";

	}
	public String getConsensusMSA(){return msa;}
	public void addBridgingRead(String id, Sequence seq){
		//no need to add more reads if consensus sequence is already determined
		if(consensusReads.containsKey(id))
			return;
		else if(!bridgingReads.containsKey(id))
			bridgingReads.put(id, new ArrayList<Sequence>());
		
		getBridgingReadList(id).add(seq);
		if(getBridgingReadsNumber(id) >= BDGraph.GOOD_SUPPORT) {
			//TODO: check sequences length
			setConsensusSequence(id);
		}
	}
	public List<Sequence> getBridgingReadList(String id){
		return bridgingReads.get(id);
	}
	public int getBridgingReadsNumber(String id){
		return getBridgingReadList(id)==null?-1:getBridgingReadList(id).size();
	}
	public Sequence getConsensus(String id, boolean force){
		if(force && getBridgingReadsNumber(id)>=BDGraph.MIN_SUPPORT)
			setConsensusSequence(id); //do it even not enough bridging reads 
		
		return consensusReads.get(id);
	}
	
	//Scan the list and check for the length consistency & remove the odd-one out
	//If none is removed -> good, return 1
	//If several are removed -> return 0
	//If more than 20% are removed -> bad list, return -1 and not to considered for consensus calling 
	private int scanElementReads(String id) {
		if(!bridgingReads.containsKey(id))
			return -1;
		int size=bridgingReads.get(id).size();
		double average=bridgingReads.get(id).stream().mapToInt(Sequence::length).average().getAsDouble();
		bridgingReads.get(id).removeIf(s->GraphUtil.approxCompare(average, s.length())!=0);
		if(bridgingReads.size()==0)
			return -1;
		else if(bridgingReads.size()==size)
			return 1;
		else if(bridgingReads.size() > (1-BDGraph.R_TOL*size))
			return 0;
		else {
			bridgingReads.remove(id);
			return -1;
		}
		
	}
	
	private void setConsensusSequence(String id){
		Sequence consensus=null;
		try {
			ErrorCorrection.msa=msa;
			consensus=ErrorCorrection.consensusSequence(getBridgingReadList(id), AlignedRead.tmpFolder+File.separator+id);
		} catch (Exception e) {
			logger.debug("Invalid consensus calling:\n {} Pick first read for the consensus of bridge " + id, e);
			consensus=getBridgingReadList(id).get(0);
		}
		consensusReads.put(id, consensus);
		bridgingReads.remove(id); //don't need anymore, release the entry to free (not much) memory!
		
	}
	
	//TODO: when to save (important) and what to save!!!
	public void saveBridgingReadsFromAlignments(AlignedRead read){
		for(AlignedRead r:read.splitAtPotentialAnchors()){
			String id=r.getEndingsID();
			if(consensusReads.containsKey(id) || connectedPairs.contains(id))
				continue;
			Alignment first=r.getFirstAlignment(), last=r.getLastAlignment();
			HashMap<String,Integer> shortestMap = BDGraph.getShortestTreeFromNode(last.node, !last.strand, Math.abs(last.readStart-first.readEnd));
			String key=first.node.getId()+ (first.strand?"o":"i");
			if(shortestMap.containsKey(key))
				connectedPairs.add(id);
			else{
				saveCorrectedSequenceInBetween(r);
			}
		}
			
	}
	//Save the long read sequence between 2 unique nodes' alignments: start and end 
	//with the aligned parts are replaced by corresponding reference parts (of Illumina data)
	//return length of non-Illumina bases from the success written sequence
	private void saveCorrectedSequenceInBetween(AlignedRead read){
		Alignment 	start = read.getFirstAlignment(), 
					end = read.getLastAlignment();
		BDNode 	fromContig = start.node,
				toContig = end.node;
		if(		fromContig.getId().compareTo(toContig.getId()) > 0 
			|| (fromContig==toContig && start.strand==false && end.strand==false)) //to agree with BDEdge.createID()
		{
			read.reverse();
			start = read.getFirstAlignment();
			end = read.getLastAlignment();
			fromContig = start.node;
			toContig = end.node;
			
		}
		//1. Appending k-flanking sequence from the start node
		Sequence flank0=((Sequence)fromContig.getAttribute("seq"))
						.subSequence(	(int)fromContig.getNumber("len")-BDGraph.getKmerSize(), 
										(int)fromContig.getNumber("len"));
		
		if(!start.strand){
			flank0=((Sequence)fromContig.getAttribute("seq"))
					.subSequence(0, BDGraph.getKmerSize());
			flank0=Alphabet.DNA.complement(flank0);
		}
		
		SequenceBuilder seqBuilder = new SequenceBuilder(Alphabet.DNA(), 1024*1024,  read.readSequence.getName());
		seqBuilder.append(flank0);
		
		//2. Filling the sequence in-between
		int posReadEnd = start.getReadPositionAtReferencePosition(start.strand?(int)fromContig.getNumber("len"):1);
		int posReadFinal = end.getReadPositionAtReferencePosition(end.strand?1:(int)toContig.getNumber("len"));
		
		for (Alignment record:read.alignments){
			BDNode contig = record.node;
			if (posReadEnd >= posReadFinal)
				break; 


			if (record.readAlignmentEnd() <= posReadEnd)
				continue;				
		
			if (record.readAlignmentStart() > posReadEnd){
				//Really need to fill in using read information
				int newPosReadEnd = Math.min(posReadFinal, record.readAlignmentStart());
				if (newPosReadEnd > posReadEnd){
					seqBuilder.append(read.readSequence.subSequence(posReadEnd-1, newPosReadEnd-1)); //subsequence is 0-index
					posReadEnd = newPosReadEnd;
					
				}
				if (posReadEnd >= posReadFinal)
					continue;//Done

				//Now get information on the contig from start
				if (contig == toContig)
					break;
				if (record.strand){
					int refLeft = record.refStart;
					int refRight = record.refEnd;

					if (posReadFinal <= record.readAlignmentEnd()){
						refRight = record.getReferencePositionAtReadPosition(posReadFinal); 
						posReadEnd = posReadFinal;
					}else{
						posReadEnd = record.readAlignmentEnd();
					}
					if(refLeft > refRight)
						continue;
			
					seqBuilder.append(((Sequence)contig.getAttribute("seq")).subSequence(refLeft-1, refRight-1));

				}else{//neg strand
					int refRight = record.refStart;
					int refLeft = record.refEnd;

					if (posReadFinal <= record.readAlignmentEnd()){
						refLeft = record.getReferencePositionAtReadPosition(posReadFinal); 
						posReadEnd = posReadFinal;
					}else{
						posReadEnd = record.readAlignmentEnd();
					}
					if(refLeft < refRight)
						continue;
					
				
					seqBuilder.append(Alphabet.DNA.complement(((Sequence)contig.getAttribute("seq")).subSequence(refRight-1, refLeft-1)));
				}
			}else{//Now get information on the contig from start
				if (contig == toContig)
					break;
				if (record.strand){
					int refLeft = record.getReferencePositionAtReadPosition(posReadEnd);						
					int refRight = record.refEnd;

					if (posReadFinal <= record.readAlignmentEnd()){
						refRight = record.getReferencePositionAtReadPosition(posReadFinal); 
						posReadEnd = posReadFinal;
					}else{
						posReadEnd = record.readAlignmentEnd();
					}
					if(refLeft > refRight)
						continue;
					
					
					seqBuilder.append(((Sequence)contig.getAttribute("seq")).subSequence(refLeft-1, refRight-1));
				}else{//neg strand						
					int refLeft = record.getReferencePositionAtReadPosition(posReadEnd);		
					int refRight = record.refStart;

					if (posReadFinal <= record.readAlignmentEnd()){
						refRight = record.getReferencePositionAtReadPosition(posReadFinal); 
						posReadEnd = posReadFinal;
					}else{
						posReadEnd = record.readAlignmentEnd();
					}
					if(refLeft < refRight)
						continue;
					
					seqBuilder.append(Alphabet.DNA.complement(((Sequence)contig.getAttribute("seq")).subSequence(refRight - 1, refLeft-1)));
				}
			}
		}
		
		//3. Make sure the k-flanking sequence of the ending node is included at last.
		
		Sequence flank1=((Sequence)toContig.getAttribute("seq"))
						.subSequence(0, BDGraph.getKmerSize());
		if(!end.strand){
				flank1=((Sequence)toContig.getAttribute("seq"))
						.subSequence(	(int)toContig.getNumber("len")-BDGraph.getKmerSize(), 
										(int)toContig.getNumber("len"));
				flank1=Alphabet.DNA.complement(flank1);
		}
		if(posReadEnd>posReadFinal){
			//scan for overlap
			int overlap=GraphUtil.overlap(flank0, flank1);
			if(overlap>0)
				seqBuilder.append(flank1.subSequence(overlap, flank1.length()));
		}else{
			//just append
			seqBuilder.append(flank1);
		}
		
		
		String key=read.getEndingsID();
		logger.debug("Save read "+read.readSequence.getName()+" to bridge "+key+": "
						+fromContig.getId()+(start.strand?"+":"-")+" -> "+toContig.getId() + (end.strand?"+":"-"));		Sequence seq=seqBuilder.toSequence();
		addBridgingRead(key, seq);
	}
}
