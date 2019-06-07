/*****************************************************************************
 * Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/**************************     REVISION HISTORY    **************************
 * 31/12/2014 - Minh Duc Cao: Created ReadFilling.java
 * 20/08/2018 - Son Nguyen: Adapted to AlignedRead.java                                       
 *  
 ****************************************************************************/

package org.rtassembly.npgraph;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceBuilder;
import japsa.seq.SequenceOutputStream;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;


public class AlignedRead{
	public static String tmpFolder=System.getProperty("usr.dir")+File.separator+"npGraph_tmp"; //folder to save spanning reads of the bridge

	/**
	 * The read sequence
	 */
	Sequence readSequence;		
	private boolean sorted = false;
	
	//This is only used in uniqueBridgesFinding()
	private int eFlag=0; // 0: both ends are from non-unique nodes; 1: start node is unique; 2: end node is unique; 3: both ends are from unique nodes
	
	ArrayList<Alignment> alignments;	

	public AlignedRead(Sequence read, ArrayList<Alignment> alignmentList){
		readSequence = read;
		alignments = alignmentList;
	}
	public AlignedRead(Sequence read, Alignment alg){
		this(read, new ArrayList<Alignment>());
		append(alg);
	}
	
	public Sequence getReadSequence(){
		return readSequence;
	}

	public ArrayList<Alignment> getAlignmentRecords(){
		return alignments;
	}
	public void append(Alignment alg){
		alignments.add(alg);
	}
	public void appendAll(ArrayList<Alignment> algs){
		alignments.addAll(algs);
	}
	public Alignment getFirstAlignment(){
		if(alignments==null || alignments.isEmpty())
			return null;
		return alignments.get(0);
	}
	public Alignment getLastAlignment(){
		if(alignments==null || alignments.size() <= 1)
			return null;
		return alignments.get(alignments.size()-1);
	}
	
	public void setEFlag(int flag){
		eFlag=flag;
	}
	public int getEFlag(){
		return eFlag;
	}

	
	public void reverse(){
		//return an (conceptually the same) read filling with the a reverse read
		Sequence revRead = Alphabet.DNA.complement(readSequence);
		revRead.setName("REV"+readSequence.getName());
		ArrayList<Alignment> revAlignments = new ArrayList<Alignment>(); 
		
		for (Alignment alignment:alignments)
			revAlignments.add(0,alignment.reverseRead());
		
		readSequence=revRead;
		alignments=revAlignments;

	}
	
	public void sortAlignment(){
		if (!sorted){				
			Collections.sort(alignments);
			sorted = true;
		}
	}
	public String getEndingsID() {
		if(alignments==null || alignments.isEmpty())
			return "-,-";
		else{
			BDEdgePrototype tmp = new BDEdgePrototype(getFirstAlignment(),getLastAlignment());
			return tmp.toString();
			
		}
	}
	public String getCompsString() {
		String retval = "";
		if(alignments!=null && !alignments.isEmpty()){
			for(Alignment alg:alignments)
				retval+=alg.node.getId()+ (alg.strand?"+":"-");
			
		}
		return retval;
	}
	//get ScaffoldVector a->b from two alignments of *this*
	public ScaffoldVector getVector(Alignment a, Alignment b) {
		if(a==b)
			return new ScaffoldVector();
		
		int 	alignedReadLen = Math.abs(a.readEnd - a.readStart) + Math.abs(b.readEnd - b.readStart),
				alignedRefLen = Math.abs(a.refEnd - a.refStart) + Math.abs(b.refEnd - b.refStart);
		double rate = 1.0 * alignedRefLen/alignedReadLen;		

		int alignP = (int) ((b.readStart - a.readStart) * rate);
		int alignD = (a.strand == b.strand)?1:-1;

		//(rough) relative position from ref_b (contig of b) to ref_a (contig of a) in the assembled genome
		int gP = (alignP + (a.strand ? a.refStart:-a.refStart) - (b.strand?b.refStart:-b.refStart));
		if (!a.strand)
			gP = -gP;	

		return new ScaffoldVector(gP, alignD);	
	}
	
	//Split an AlignedRead at a specific alignment record.
	//E.g. <...a,b,c,d...> split at (c) becoming: <...a,b,c> <c,d,...>
	public ArrayList<AlignedRead> split(Alignment cutPoint){
		ArrayList<AlignedRead> retval = new ArrayList<AlignedRead>();
		//recursive function?
		if(alignments!=null){
			int idx=alignments.indexOf(cutPoint);
			if(idx>=0){
				AlignedRead left = new AlignedRead( readSequence.subSequence(0, cutPoint.readAlignmentEnd()+1), 
													new ArrayList<>(alignments.subList(0, idx+1))),
							right = new AlignedRead(readSequence.subSequence(cutPoint.readAlignmentStart(), readSequence.length()), 
													new ArrayList<>(alignments.subList(idx, alignments.size())));
				retval.add(left);
				retval.add(right);
			}				
		}			
		return retval;
	}
	
	//Save the long read sequence between 2 unique nodes' alignments: start and end 
	//with the aligned parts are replaced by corresponding reference parts (of Illumina data)

	public void saveCorrectedSequenceInBetween(){
		Alignment 	start = getFirstAlignment(), 
					end = getLastAlignment();
		BDNode 	fromContig = start.node,
				toContig = end.node;
		
		if(SimpleBinner.getBinIfUniqueNow(fromContig)==null || SimpleBinner.getBinIfUniqueNow(toContig)==null)
			return;
	
		
		SequenceBuilder seqBuilder = new SequenceBuilder(Alphabet.DNA5(), 1024*1024,  readSequence.getName());
		if(start.readAlignmentStart() > end.readAlignmentEnd())
			reverse();
		
		sortAlignment();
		//1. Appending (k-1)-flanking sequence from the start node
		Sequence flank0=((Sequence)fromContig.getAttribute("seq"))
						.subSequence(	(int)fromContig.getNumber("len")-BDGraph.getKmerSize()+1, 
										(int)fromContig.getNumber("len"));
		
		if(start.strand)
			flank0=Alphabet.DNA.complement(flank0);

		seqBuilder.append(flank0);
		
		//2. Filling the sequence in-between
//		int posReadEnd   = start.readAlignmentEnd();
//		int posReadFinal = end.readAlignmentStart();// I need as far as posReadFinal
		//extrapolating...
		int posReadEnd = start.getReadPositionAtReferencePosition(start.strand?(int)fromContig.getNumber("len"):1);
		int posReadFinal = end.getReadPositionAtReferencePosition(end.strand?1:(int)toContig.getNumber("len"));
		
		for (Alignment record:alignments){
			BDNode contig = record.node;
			System.out.println("Current alignment record: " + record.toString());
			System.out.println("posReadEnd=" + posReadEnd + " posReadFinal=" + posReadFinal);
			System.out.println("Read " + readSequence.getName() + " length=" + readSequence.length());
			if (posReadEnd >= posReadFinal)
				break; 


			if (record.readAlignmentEnd() <= posReadEnd)
				continue;				
		
			if (record.readAlignmentStart() > posReadEnd){
				//Really need to fill in using read information
				int newPosReadEnd = Math.min(posReadFinal, record.readAlignmentStart());
				if (newPosReadEnd > posReadEnd){
					seqBuilder.append(readSequence.subSequence(posReadEnd-1, newPosReadEnd-1)); //subsequence is 0-index
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

				}else{//neg strain
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
			}//if record.readAlignmentStart() > posReadEnd
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
		
		//3. Make sure the (k-1)-flanking sequence of the ending node is included at last.
		
		Sequence flank1=((Sequence)toContig.getAttribute("seq"))
						.subSequence(0, BDGraph.getKmerSize()-1);
		if(end.strand)
				flank1=Alphabet.DNA.complement(flank1);
		
		if(posReadEnd>posReadFinal){
			//scan for overlap
			int overlap=GraphUtil.overlap(flank0, flank1);
			if(overlap>0)
				seqBuilder.append(flank1.subSequence(overlap, flank1.length()));
		}else{
			//just append
			seqBuilder.append(flank1);
		}
		
		//seqBuilder.toSequence();
		String fileName=tmpFolder+File.separator+BDEdge.createID(fromContig, toContig, start.strand, !end.strand)+".fasta";
		try {
			SequenceOutputStream out = new SequenceOutputStream(new FileOutputStream(fileName,true));
			seqBuilder.toSequence().writeFasta(out);
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			System.exit(1);
		}
	}
	
}