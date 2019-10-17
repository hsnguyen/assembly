package org.rtassembly.npgraph;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import japsa.seq.Alphabet;
import japsa.seq.FastaReader;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
//Module for consensus bridging of contigs when they're not connected in the assembly graph
public class ConsensusCaller {
	private static final Logger LOG = LoggerFactory.getLogger(ConsensusCaller.class);
	
	HashMap<String, List<Sequence>> bridgingReads;

	ConsensusCaller(){
		bridgingReads = new HashMap<String, List<Sequence>>();
	}
	/* 
	 * Adapt from japsa without need to output fasta input file
	 */
	public Sequence consensusSequence(String prefix) throws IOException, InterruptedException{
		Sequence consensus = null;	
		String msa = BDGraph.MSA;
		String 	faoFile = AlignedRead.tmpFolder+File.separator+prefix+"_"+msa+".fasta";
		File consFile = new File(faoFile);
		if(!consFile.exists()){
			String faiFile=AlignedRead.tmpFolder+File.separator+prefix+".fasta";
			//1.0 Check faiFile?
			if(!Files.isReadable(Paths.get(faiFile))){
				return null;
			}
			FastaReader faiReader =  new FastaReader(faiFile);
			int count=0, minGap=Integer.MAX_VALUE;
			Sequence seq=null;
			while((seq=faiReader.nextSequence(Alphabet.DNA()))!=null){
				//if there is no valid MSA detected, output the sequence with the least non-Illumina bases
				if(!msa.startsWith("kalign")){
					String[] toks=seq.getDesc().split("=");
					try{
						int curGap=Integer.parseInt(toks[1]);
						if( curGap < minGap){
							minGap=curGap;
							consensus=seq;
						}		
					}catch(NumberFormatException e){
						if(HybridAssembler.VERBOSE)
							LOG.info("Pattern %s=%d not found in the header of input FASTA file {}", faiFile);
						if(consensus==null || consensus.length() > seq.length())
							consensus=seq;
					}
				}
				count++;
			}
			
			faiReader.close();
			if(count<BDGraph.MIN_SUPPORT) //at least 3 of the good read count is required
				return null;
			
			//2.0 Run multiple alignment. 
			// For now, only kalign were chosen due to its speedy operation
			String cmd  = "";

			if (msa.startsWith("kalign"))
				cmd = "kalign -gpo 60 -gpe 10 -tgpe 0 -bonus 0 -q -i " + faiFile	+ " -o " + faoFile;
//			else if (msa.startsWith("poa"))
//				cmd = "poa -read_fasta " + faiFile + " -clustal " + faoFile + " -hb blosum80.mat";
//			else if (msa.startsWith("muscle"))
//				cmd = "muscle -in " + faiFile + " -out " + faoFile + " -maxiters 5 -quiet";				
//			else if (msa.startsWith("clustal")) 
//				cmd = "clustalo --force -i " + faiFile + " -o " + faoFile;
//			else if (msa.startsWith("msaprobs"))
//				cmd = "msaprobs -o " + faoFile + " " + faiFile;
//			else if (msa.startsWith("mafft"))
//				cmd = "mafft_wrapper.sh  " + faiFile + " " + faoFile;
			else{
				if(HybridAssembler.VERBOSE)
					LOG.info("Unknown msa function {}, pick the best read as consensus!", msa);
				consensus.setName(prefix+"_consensus");
				Files.deleteIfExists(Paths.get(faiFile));
				Files.deleteIfExists(Paths.get(faoFile));
				return consensus;
			}
			
			if(HybridAssembler.VERBOSE)
				LOG.info("Running " + cmd);
			Process process = Runtime.getRuntime().exec(cmd);
			process.waitFor();
			if(HybridAssembler.VERBOSE)
				LOG.info("Done " + cmd);
			Files.deleteIfExists(Paths.get(faiFile));

		}
				


		//3.0 Read in multiple alignment
		ArrayList<Sequence> seqList = new ArrayList<Sequence>();
		{
			SequenceReader msaReader = FastaReader.getReader(faoFile);
			Sequence nSeq = null;
			while ((nSeq = msaReader.nextSequence(Alphabet.DNA())) != null) {
				seqList.add(nSeq);						
			}
			msaReader.close();
		}

		//4.0 get consensus			
		{
			int [] coef = new int[seqList.size()];			
			for (int y = 0; y < seqList.size(); y++){
				coef[y] = 1;				
			}
			//TODO: combine error profiles?? 
			int [] counts = new int[6];

			SequenceBuilder sb = new SequenceBuilder(Alphabet.DNA(), seqList.get(0).length());
			for (int x = 0; x < seqList.get(0).length();x++){		
				Arrays.fill(counts, 0);
				for (int y = 0; y < seqList.size(); y++){					
					byte base = seqList.get(y).getBase(x);
					if (base >= 6) 
						counts[4] += coef[y];//N
					else
						counts[base] += coef[y];
				}//for y
				int maxIdx = 0;
				for (int y = 1; y < counts.length; y++){
					if (counts[y] > counts[maxIdx])
						maxIdx = y;					
				}//for y
				if (maxIdx < Alphabet.DNA.GAP){//not a gap
					sb.append((byte)maxIdx);
				}//if
			}//for x
			sb.setName(prefix+"_consensus");
			if(HybridAssembler.VERBOSE)
				LOG.info(sb.getName() + "  " + sb.length());
			consensus = sb.toSequence();
		}
		
		Files.deleteIfExists(Paths.get(faoFile));
		return consensus;
	}
	
	
	//no check: use with care
	public void saveReadToDisk(AlignedRead read){
		if(BDGraph.getReadsNumOfBrg(read.getEndingsID())>=BDGraph.MAX_LISTING)
			return;
		
		if(read.getEFlag()>=3)
			saveCorrectedSequenceInBetween(read);
		else
			read.splitAtPotentialAnchors().forEach(r->saveCorrectedSequenceInBetween(r));
			
	}
	//Save the long read sequence between 2 unique nodes' alignments: start and end 
	//with the aligned parts are replaced by corresponding reference parts (of Illumina data)
	//return length of non-Illumina bases from the success written sequence
	public int saveCorrectedSequenceInBetween(AlignedRead read){
		Alignment 	start = read.getFirstAlignment(), 
					end = read.getLastAlignment();
		BDNode 	fromContig = start.node,
				toContig = end.node;
		int noise=0;
		if(		fromContig.getId().compareTo(toContig.getId()) > 0 
			|| (fromContig==toContig && start.strand==false && end.strand==false)) //to agree with BDEdge.createID()
		{
			read.reverse();
			start = read.getFirstAlignment();
			end = read.getLastAlignment();
			fromContig = start.node;
			toContig = end.node;
			
		}
//		sortAlignment();
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
//		int posReadEnd   = start.readAlignmentEnd();
//		int posReadFinal = end.readAlignmentStart();// I need as far as posReadFinal
		//extrapolating...
		int posReadEnd = start.getReadPositionAtReferencePosition(start.strand?(int)fromContig.getNumber("len"):1);
		int posReadFinal = end.getReadPositionAtReferencePosition(end.strand?1:(int)toContig.getNumber("len"));
		
		for (Alignment record:read.alignments){
			BDNode contig = record.node;
//			System.out.println("Current alignment record: " + record.toString());
//			System.out.println("posReadEnd=" + posReadEnd + " posReadFinal=" + posReadFinal);
//			System.out.println("Read " + readSequence.getName() + " length=" + readSequence.length());
			if (posReadEnd >= posReadFinal)
				break; 


			if (record.readAlignmentEnd() <= posReadEnd)
				continue;				
		
			if (record.readAlignmentStart() > posReadEnd){
				//Really need to fill in using read information
				int newPosReadEnd = Math.min(posReadFinal, record.readAlignmentStart());
				if (newPosReadEnd > posReadEnd){
					seqBuilder.append(read.readSequence.subSequence(posReadEnd-1, newPosReadEnd-1)); //subsequence is 0-index
					noise+=newPosReadEnd-posReadEnd;
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
			}//FIXME: scan for overlap between suggested references instead of blindly go based on alignments
			else{//Now get information on the contig from start
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
		//TODO: double-check this case
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
		if(HybridAssembler.VERBOSE)
			LOG.info("Save read %s to bridge %s: %s -> %s", read.readSequence.getName(), key,
																	fromContig.getId() + (start.strand?"+":"-"),
																	toContig.getId() + (end.strand?"+":"-")
																	);

		String fileName=AlignedRead.tmpFolder+File.separator+key+".fasta";
		try {
			SequenceOutputStream out = new SequenceOutputStream(new FileOutputStream(fileName,true));
			Sequence seq=seqBuilder.toSequence();
			seq.setDesc("noise="+noise);
			seq.writeFasta(out);
			out.close();
			BDGraph.addBrg2ReadsNum(key);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return -1;
		}
		
		return noise;
	}
}
