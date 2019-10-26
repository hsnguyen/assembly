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
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


public class AlignedRead{
	private static final Logger LOG = LoggerFactory.getLogger(AlignedRead.class);

	public static int PSEUDO_ID=1;
	public static String tmpFolder=System.getProperty("usr.dir")+File.separator+"npGraph_tmp"; //folder to save spanning reads of the bridge
	//The query long read sequence
	Sequence readSequence;
	private boolean sorted = false;
	
	//This is only used in uniqueBridgesFinding()
	// 0: both ends are from non-unique nodes; 
	// 1: start node is unique; 
	// 2: end node is unique; 
	// 3: both ends are from unique nodes
	private int eFlag=0; 
	
	ArrayList<Alignment> alignments;	

	public AlignedRead(Sequence read, ArrayList<Alignment> alignmentList){
		readSequence = read;
		alignments = alignmentList;
	}
	public AlignedRead(Sequence read, Alignment alg){
		this(read, new ArrayList<Alignment>());
		append(alg);
	}
	
	public Sequence getLongReadSequence(){ //the actual read used as query	
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
		readSequence = Alphabet.DNA.complement(readSequence);
		ArrayList<Alignment> revAlignments = new ArrayList<Alignment>(); 
		
		for (Alignment alignment:alignments)
			revAlignments.add(0,alignment.reverseRead());
		
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
		else
			return BDEdge.createID(getFirstAlignment().node, getLastAlignment().node, getFirstAlignment().strand, !getLastAlignment().strand);
		
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
	
	public ArrayList<AlignedRead> splitAtPotentialAnchors(){
		if(alignments==null||alignments.isEmpty())
			return null;
		ArrayList<AlignedRead> retval = new ArrayList<AlignedRead>();
		ArrayList<Alignment> curList=new ArrayList<>();
		Alignment start=null;
		for(int i=0; i<alignments.size();i++){
			Alignment curAlg=alignments.get(i);
			curList.add(curAlg);
			if(SimpleBinner.isPotentialAnchorNode(curAlg.node)){
				if(start!=null){
					retval.add(new AlignedRead(readSequence,curList));
				}
				start=curAlg;
				curList=new ArrayList<>();
				curList.add(curAlg);

			}
		}
		return retval;
	}
	
}