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
 * 31/12/2014 - Minh Duc Cao: Created AlignmentRecord.java
 * 20/08/2018 - Son Nguyen: Adapted to Alignment.java                                       
 *  
 ****************************************************************************/
package org.rtassembly.npgraph;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import japsa.seq.PAFRecord;
import japsa.seq.Sequence;

/*
 * Class holding alignment of a read to reference from a SAMRecord
 * Only keep aligned coordinates, not the query read itself
 */
public class Alignment implements Comparable<Alignment> {
//	public final static int OVERHANG_THRES=500; 
	public final static int GOOD_QUAL=60; 

	public static int MIN_QUAL=10; 

	public int quality;

	public String readID;
	BDNode node; //keep an instance of contig reference here in case it's removed from BDGraph

	public int refStart, refEnd;  //1-based position on ref of the start and end of the alignment
	
	//Position on read of the start and end of the alignment (using the direction of read) 
	//readStart map to refStart, readEnd map to refEnd. 
	//readStart < readEnd if strand = true, else readStart > readEnd
	public int readStart = 0, readEnd = 0;	
	
	//read length
	public int readLength = 0;		

	//TODO: take into account NM (edit distance) and AS (alignment score) tag from SAM file
	public int score = 0; //just alignment length for now
	
	public boolean strand = true;//positive
	public boolean prime = true;//primary alignment
	public boolean goodMargin = false;
	public boolean useful = false;
	
	Cigar cigar=null;
	
	public Alignment(BDNode node){ this.node=node;}
	public Alignment(PAFRecord paf, BDNode node) {
		this(node);
		readID=paf.qname;
		quality=paf.qual;
		score=paf.score;
		readLength=paf.qlen;
		refStart=paf.tstart; refEnd=paf.tend;
		strand=paf.strand;
		cigar=paf.getCigar();
		if (!strand){
			readStart=paf.qend; 
			readEnd=paf.qstart;
		}else{
			readStart=paf.qstart;
			readEnd=paf.qend;
		}
		
		//these temporary variable to determine usefulness
		int readLeft = paf.qstart -1;
		int readRight = readLength - paf.qend;

		int refLeft = refStart - 1;
		int refRight = ((Sequence) node.getAttribute("seq")).length() - refEnd;
		int overhangTolerance = (int) Math.min(BDGraph.A_TOL, BDGraph.R_TOL*node.getNumber("len"));
		if (
				(readLeft < overhangTolerance || refLeft < overhangTolerance)
				 && (readRight  < overhangTolerance || refRight < overhangTolerance)
				 && Math.min(refLeft,refRight) < overhangTolerance
			)
			goodMargin=true;
		
		if	(	goodMargin
				&& prime //TODO: there are useful secondary alignment!!!
//				&& alignLength > BDGraph.getKmerSize() //FIXME: 
				&& quality >= MIN_QUAL
			)
			useful = true;
	}
	
	public Alignment(SAMRecord sam, BDNode node) {
		this(node);
//		readID = Integer.parseInt(sam.getReadName().split("_")[0]);
		readID = sam.getReadName();
		quality = sam.getMappingQuality();
		prime=!sam.getNotPrimaryAlignmentFlag();
		this.node = node;

		refStart = sam.getAlignmentStart();
		refEnd = sam.getAlignmentEnd();
		
		cigar = sam.getCigar();			
		boolean enterAlignment = false;						
		//////////////////////////////////////////////////////////////////////////////////
		
		for (final CigarElement e : cigar.getCigarElements()) {
			final int  length = e.getLength();
			switch (e.getOperator()) {
			case H :
			case S :					
			case P : //pad is a kind of clipped
				if (enterAlignment)
					readEnd = readLength;
				readLength += length;
				break; // soft clip read bases
			case I :	                	
			case M :					
			case EQ :
			case X :
				if (!enterAlignment){
					readStart = readLength + 1;
					enterAlignment = true;
				}
				readLength += length;
				break;
			case D :
			case N :
				if (!enterAlignment){
					readStart = readLength + 1;
					enterAlignment = true;
				}
				break;				
			default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
			}//case
		}//for
		if (readEnd == 0)
			readEnd = readLength;
		//these temporary variable to determine usefulness
		int readLeft = readStart -1;
		int readRight = readLength - readEnd;

		int refLeft = refStart - 1;
		int refRight = ((Sequence) node.getAttribute("seq")).length() - refEnd;
		
		try{
			score = sam.getIntegerAttribute("AS");
		}catch(RuntimeException e){
			System.out.println("AS tag not found for score, use alignment length instead!");
			score =  refEnd + 1 - refStart;
		}
		
		if (sam.getReadNegativeStrandFlag()){			
			strand = false;
			//need to convert the alignment position on read the correct direction 
			readStart = 1 + readLength - readStart;
			readEnd = 1 + readLength - readEnd;
		}
		
		int overhangTolerance = (int) Math.min(BDGraph.A_TOL, BDGraph.R_TOL*node.getNumber("len"));
		if (
				(readLeft < overhangTolerance || refLeft < overhangTolerance)
				 && (readRight  < overhangTolerance || refRight < overhangTolerance)
				 && Math.min(refLeft,refRight) < overhangTolerance
			)
			goodMargin=true;
		
		if	(	goodMargin
				&& prime //TODO: there are useful secondary alignment!!!
//				&& alignLength > BDGraph.getKmerSize() //FIXME: 
				&& quality >= MIN_QUAL
			)
			useful = true;

	}
	
	
	public int readAlignmentStart(){
		return Math.min(readStart,readEnd);
	
	}
	
	public int readAlignmentEnd(){
		return Math.max(readStart,readEnd);
	}

	public Alignment reverseRead(){
		Alignment revAlign = new Alignment(node);
		revAlign.readID = readID;
		revAlign.refStart = refStart;
		revAlign.refEnd = refEnd;
		
		revAlign.readLength = readLength;
		revAlign.readStart = readLength - readStart + 1;//1-index
		revAlign.readEnd = readLength - readEnd + 1;//1-index
		revAlign.strand = !strand;
		revAlign.useful = useful;			
		revAlign.score = score;
		revAlign.cigar = cigar;//https://www.biostars.org/p/289583/
		revAlign.quality = quality;

		return revAlign;
	}
	public String toString() {
		return node.getAttribute("name")  
				+ ": " + refStart 
				+ " -> " + refEnd
				+ " / " + ((Sequence) node.getAttribute("seq")).length()
				+ " map to "
				+ readID
				+ ": " + readStart 
				+ " -> " + readEnd
				+ " / " + readLength				
				+ ", strand: " + (strand?"+":"-")
				+ ", prime: " + (prime?"yes":"no")
				+ ", margin: " + (goodMargin?"good":"bad")
				+ ", qual=: " + quality;
	}
	
	/**
	 * Re-implement what already in htsjdk (2.10.1) SAMRecord because we remove SAMRecord from Alignment
	 * Return the position on the reference that corresponds to a given position on the read.
	 *  
	 * @param posInRead
	 * @param record
	 * @return
	 */
	public int getReferencePositionAtReadPosition(int readLookingPosition){
		int location=-1;
			
//		if (readLookingPosition < readAlignmentStart() || readLookingPosition > readAlignmentEnd())
//			return 1;

		int posOnRead = readStart, endOfRead=readEnd;
		int posOnRef = refStart;
		if (!strand) {
			readLookingPosition = readLength + 1 - readLookingPosition; // use direction of ref (forward)
			posOnRead = readLength + 1 - posOnRead;
			endOfRead = readLength + 1 - endOfRead;
		}

		if( readLookingPosition < posOnRead) {
			location = refStart + readLookingPosition - posOnRead;
		}else if(readLookingPosition > endOfRead) {
			location = refEnd + readLookingPosition - posOnRead;
		}else {
			if(cigar==null){ //approximate
				location = (int) (posOnRef + (readLookingPosition - posOnRead)*Math.abs((refEnd-refStart)*1.0/(readEnd-readStart)));	
			}else{	
				boolean found=false;
				for (final CigarElement e : cigar.getCigarElements()) {
					if(found)
						break;
					final int  length = e.getLength();
					switch (e.getOperator()) {
					case H :
					case S :					
					case P :
						break; // ignore pads and clips
					case I :				
						//insert
						if (posOnRead + length < readLookingPosition){
							posOnRead += length;				
						}else{
							location = posOnRef;
							found=true;
						}
						break;
					case M ://match or mismatch				
					case EQ://match
					case X ://mismatch
						if (posOnRead + length < readLookingPosition){
							posOnRead += length;
							posOnRef += length;
						}else{
							location = posOnRef + readLookingPosition - posOnRead;
							found=true;
						}
						break;
					case D :
						posOnRef += length;
						break;
					case N :	
						posOnRef += length;
						break;								
					default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
					}//casse
				}//for		
			}
		
		}
		location=location>refStart?location:refStart;
		location=location<refEnd?location:refEnd;
		
		return location;
//		return location<=refEnd&&location>=refStart?location:-1;
	}
	/**
	 * Return the position on the read that corresponds to a given position
	 * on reference with extrapolate if neccessary.
	 *  
	 * @param posInRef
	 * @param record
	 * @return
	 */
	public int getReadPositionAtReferencePosition(int refLookingPosition){
		// read htsjdk.samtools.* API
		int location=-1;

		if(refLookingPosition < refStart) {
			location = strand?readStart+refLookingPosition-refStart:readStart-refLookingPosition+refStart;
		}else if(refLookingPosition > refEnd) {
			location = strand?readEnd+refLookingPosition-refEnd:readEnd-refLookingPosition+refEnd;
		}
		else{
			// current coordinate on sense/anti-sense read that has the same direction as the ref contig
			int posOnRead = strand?readStart:readLength-readStart+1;
			 // current position on ref 
			int posOnRef = refStart;
			
			if(cigar==null){ //approximate
				location = (int) (posOnRead + (refLookingPosition - posOnRef)*Math.abs(1.0*(readEnd-readStart)/(refEnd-refStart)));
			}else{
				boolean found=false;
				for (final CigarElement e : cigar.getCigarElements()) {
					if(found)
						break;
					final int  length = e.getLength();
					switch (e.getOperator()) {
					case H :
					case S :					
					case P :
						break; // ignore pads and clips
					case I :			
						posOnRead += length;
						break;	
					case M ://match or mismatch				
					case EQ://match
					case X ://mismatch
						if (posOnRef + length < refLookingPosition){
							posOnRef += length;
							posOnRead += length;
						}else{
							location = refLookingPosition + posOnRead - posOnRef;
							found=true;
						}
						break;
					case D :
					case N :	
						//delete
						if (posOnRef + length < refLookingPosition){
							posOnRef += length;				
						}else{
							location = posOnRead;
							found=true;
						}
						break;	
					default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
					}//casse
				}//for		
			}
			//convert back to coordinate based on read direction
			location = strand?location:readLength-location+1; //1-index
		}

		location=location>1?location:1;
		location=location<readLength?location:readLength;
		return location;
//		return location>=readAlignmentStart()&&location<=readAlignmentEnd()?location:-1;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Comparable#compareTo(java.lang.Object)
	 */
	@Override
	public int compareTo(Alignment o) {			
		return readAlignmentStart() - o.readAlignmentStart();
	}
}
