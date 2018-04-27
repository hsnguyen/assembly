package org.rtassembly.concatemer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import org.rtassembly.scaffold.AlignmentRecord;
import org.rtassembly.scaffold.Contig;
import org.rtassembly.scaffold.ContigBridge;
import org.rtassembly.scaffold.ReadFilling;
import org.rtassembly.scaffold.ScaffoldVector;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.FastqSequence;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.HTSUtilities;

public class SequenceChopper {
	Sequence ref;
	SequenceChopper(String refFile) throws IOException{
		SequenceReader reader = SequenceReader.getReader(refFile);
		ref=reader.nextSequence(Alphabet.DNA()); //ref only have 1 sequence
		reader.close();
	}
	// Chop the concatemers (long nanopore reads) based on the reference
	public void chopLongReadsBasedOnRef(String samFile, String outFile, double qual) throws IOException {
		SequenceOutputStream fout = SequenceOutputStream.makeOutputStream(outFile);
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader=SamReaderFactory.makeDefault().open(new File(samFile));
		SAMRecordIterator iter = reader.iterator();

		String readID = "";
		AlignmentRecord myRec = null;
		ArrayList<AlignmentRecord> samList = null;// alignment record of the same read;		
		
		Contig refContig=new Contig(0, ref);//adapt to the AlignmentRecord API
		FastqSequence fullSeq = null;
		while (iter.hasNext()) {
			SAMRecord rec = iter.next();
			if(rec.getReadUnmappedFlag())
				continue;

			myRec = new AlignmentRecord(rec, refContig);
			
			if (!readID.equals(myRec.readID)){
				//chop and print seq ...
				if(samList!=null && samList.size() >= 1)
					chopAndPrint(fullSeq, samList, fout);
				//re-initialize
				samList = new ArrayList<AlignmentRecord>();
				readID = myRec.readID;	
				fullSeq=new FastqSequence(Alphabet.DNA5(), rec);
			}
//			if (rec.getMappingQuality() < qual || !myRec.useful)	
			if (rec.getMappingQuality() < qual || Math.abs(myRec.readEnd-myRec.readStart) < .8*ref.length()) //80% error rate of nanopore data		
				continue;		
			
			samList.add(myRec);

		}// while
		
		//FastqSequence seq.print(fout);
		fout.close();
		iter.close();
		reader.close();
	}
	
	private void chopAndPrint(FastqSequence read, ArrayList<AlignmentRecord> records, SequenceOutputStream out) {
		if(read==null || records == null || records.size()==0)
			return;
		System.out.println("===========================================================");

		int[] 	start = new int[records.size()], 
				stop = new int[records.size()];
		int count=0;
		for(AlignmentRecord rec:records) {
			start[count]=mapToRead(0, rec);
			stop[count]=mapToRead(ref.length()-1, rec);
			count++;
		}
		System.out.println("Read " + read.getName() + " length=" + read.length() + ":");
		for(int i=0;i<count;i++) {
			System.out.printf("\t%d -> %d\n", start[i], stop[i]);
		}
		
	}
	static int mapToRead(int posOnRef, AlignmentRecord record){
		// read htsjdk.samtools.* API
		int location = -1;
		
		if ((posOnRef - record.refStart)*(posOnRef - record.refEnd) >= 0){
			if (Math.abs(posOnRef-record.refStart) > Math.abs(posOnRef-record.refEnd))
				location = record.strand?record.readEnd+posOnRef-record.refEnd:record.readEnd-posOnRef+record.refEnd;			
			else
				location = record.strand?record.readStart+posOnRef-record.refStart:record.readStart-posOnRef+record.refStart;
		}
		else{
			// current coordinate on read, followed the reference contig's direction
			int posOnRead = record.strand?record.readStart:record.readLength-record.readStart+1;
			 // current position on ref 
			int pos = record.refStart;
			
			for (final CigarElement e : record.getCigars()) {
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
					if (pos + length < posOnRef){
						pos += length;
						posOnRead += length;
					}else{
						location = posOnRef + posOnRead - pos;
					}
					break;
				case D :
				case N :	
					//delete
					if (pos + length < posOnRef){
						pos += length;				
					}else{
						location = posOnRead;
					}
					break;	
				default : throw new IllegalStateException("Case statement didn't deal with cigar op: " + e.getOperator());
				}//casse
			}//for		
			//convert back to coordinate based on read direction
			location = record.strand?location:record.readLength-location+1;
		}
		
//		System.out.println( "Contig (ref): " + record.getContig().getName() + " Read: " + record.readID + " Strand: " + record.strand);   
//		System.out.println( "\tOn contig: " + record.refStart	+ " -> " + record.refEnd
//				+ " Len: " + record.getContig().length() + " Cut point: " + posOnRef);
//		System.out.println( "\tOn read: " + record.readStart	+ " -> " + record.readEnd
//				+ " Len: " + record.readLength + " Alleged cut point: " + location);
		
		location=location>0?location:0;
		location=location<record.readLength?location:record.readLength-1;
		
//		System.out.println( "\tOn read: " + record.readStart	+ " -> " + record.readEnd
//				+ " Len: " + record.readLength + " Final cut point: " + location);
		return location;
	}
//	private ScaffoldVector getVector(ReadFilling readSequence, AlignmentRecord a, AlignmentRecord b) {
//
//		// Rate of aligned lengths: ref/read (illumina contig/nanopore read)
//		int 	alignedReadLen = Math.abs(a.readEnd - a.readStart) + Math.abs(b.readEnd - b.readStart),
//				alignedRefLen = Math.abs(a.refEnd - a.refStart) + Math.abs(b.refEnd - b.refStart);
//		double rate = 1.0 * alignedRefLen/alignedReadLen;		
//		int alignP = (int) ((b.readStart - a.readStart) * rate);
//		int alignD = (a.strand == b.strand)?1:-1;
//		//(rough) relative position from ref_b (contig of b) to ref_a (contig of a) in the assembled genome
//		int gP = (alignP + (a.strand ? a.refStart:-a.refStart) - (b.strand?b.refStart:-b.refStart));
//		if (!a.strand)
//			gP = -gP;
//
//		return new ScaffoldVector(gP, alignD);
//	}
	
	public static void main(String[] args) throws IOException {
		SequenceChopper tony_tony_chopper = new SequenceChopper("/home/s_hoangnguyen/Projects/plant_virus/ref/CaMV.fasta");
		tony_tony_chopper.chopLongReadsBasedOnRef("/home/s_hoangnguyen/Projects/plant_virus/alignment/barcode08_pass_mm2.sam", "/home/s_hoangnguyen/Projects/plant_virus/test.out", 10);

	}

}
