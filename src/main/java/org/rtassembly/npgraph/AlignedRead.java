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

import java.util.ArrayList;
import java.util.Collections;


public class AlignedRead{
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
	
	//Return the long read sequence between 2 unique nodes' alignments: start and end 
	//with the aligned parts are replaced by corresponding reference parts (of Illumina data)
	//TODO: check if including 2 flanking sequences increases the poa performance

	public Sequence getCorrectedSequence(Alignment start, Alignment end){
		assert alignments.contains(start)&&alignments.contains(end):"Ending alignments not belong to the read!";
		BDNode 	fromContig = start.node,
				toContig = end.node;
		SequenceBuilder seqBuilder = new SequenceBuilder(Alphabet.DNA5(), 1024*1024,  getEndingsID());
		if(start.readAlignmentStart() > end.readAlignmentEnd())
			reverse();
		//1. Appending (k-1)-flanking sequence from the start node
		
		
		//2. Filling the sequence in-between
		int posReadEnd   = start.readAlignmentEnd();
		int posReadFinal = end.readAlignmentStart();// I need as far as posReadFinal
		// locate the last position being extended...
		if(posReadEnd >= posReadFinal ){
			//TODO: Scan for overlap, append sufficient DNA and return;
		}
		
		for (Alignment record:alignments){
			BDNode contig = record.node;
//			if (contig.getIndex() == fromContig.getIndex())
//				continue;

			if (posReadEnd >= posReadFinal -1)
				break; 


			if (record.readAlignmentEnd() <= posReadEnd)
				continue;				
		
			if (record.readAlignmentStart() > posReadEnd){
				//Really need to fill in using read information
				int newPosReadEnd = Math.min(posReadFinal - 1, record.readAlignmentStart() -1);
				if (newPosReadEnd > posReadEnd){
					seqBuilder.append(readSequence.subSequence(posReadEnd, newPosReadEnd));
					posReadEnd = newPosReadEnd;
					
				}
				if (posReadEnd + 1 >= posReadFinal)
					continue;//Done

				//Now get information on the contig from start
				if (contig.getIndex() == toContig.getIndex())
					continue;//tandem
				if (record.strand){
					int refLeft = record.refStart;
					int refRight = record.refEnd;

					if (posReadFinal <= record.readAlignmentEnd()){
						refRight = record.positionOnRef(posReadFinal) -1; 
						posReadEnd = posReadFinal -1;
					}else{
						posReadEnd = record.readAlignmentEnd();
					}
					if(refLeft > refRight)
						continue;
			
					seqBuilder.append(((Sequence)contig.getAttribute("seq")).subSequence(refLeft - 1, refRight));

				}else{//neg strain
					int refRight = record.refStart;
					int refLeft = record.refEnd;

					if (posReadFinal <= record.readAlignmentEnd()){
						refLeft = record.positionOnRef(posReadFinal) + 1; 
						posReadEnd = posReadFinal -1;
					}else{
						posReadEnd = record.readAlignmentEnd();
					}
					if(refLeft < refRight)
						continue;
					
				
					seqBuilder.append(Alphabet.DNA.complement(((Sequence)contig.getAttribute("seq")).subSequence(refRight - 1, refLeft)));
				}
			}//if record.readAlignmentStart() > posReadEnd
			else{//Now get information on the contig from start
				if (contig.getIndex() == toContig.getIndex())
					continue;//tandem
				if (record.strand){
					int refLeft = record.positionOnRef(posReadEnd) + 1;						
					int refRight = record.refEnd;

					if (posReadFinal <= record.readAlignmentEnd()){
						refRight = record.positionOnRef(posReadFinal) -1; 
						posReadEnd = posReadFinal -1;
					}else{
						posReadEnd = record.readAlignmentEnd();
					}
					if(refLeft > refRight)
						continue;
					
					
					seqBuilder.append(((Sequence)contig.getAttribute("seq")).subSequence(refLeft - 1, refRight));
				}else{//neg strand						
					int refLeft = record.positionOnRef(posReadEnd) + 1;		
					int refRight = record.refStart;

					if (posReadFinal <= record.readAlignmentEnd()){
						refRight = record.positionOnRef(posReadFinal) + 1; 
						posReadEnd = posReadFinal -1;
					}else{
						posReadEnd = record.readAlignmentEnd();
					}
					if(refLeft < refRight)
						continue;
					
					seqBuilder.append(Alphabet.DNA.complement(((Sequence)contig.getAttribute("seq")).subSequence(refRight - 1, refLeft)));
				}
			}
		}
		
		//3. Appending the k-flanking sequence of the ending node
		
		
		
		return seqBuilder.toSequence();
		
	}
}