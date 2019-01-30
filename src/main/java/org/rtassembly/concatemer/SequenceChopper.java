package org.rtassembly.concatemer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;

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
import org.rtassembly.scaffold.AlignmentRecord;
import org.rtassembly.scaffold.Contig;

public class SequenceChopper {
	Sequence ref;
	SequenceChopper(String refFile) throws IOException{
		SequenceReader reader = SequenceReader.getReader(refFile);
		ref=reader.nextSequence(Alphabet.DNA()); //ref only have 1 sequence
		reader.close();
	}
	// Chop the concatemers (long nanopore reads) based on the reference
	public void chopLongReadsBasedOnRef(String samFile, String outFile, double qual) throws IOException {
		HashMap<Integer,ArrayList<String>> count2Reads = new HashMap<>(); // concatemer information
		ArrayList<Integer> allLength = new ArrayList<Integer>();
		
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
				if(samList!=null && samList.size() >= 1) {
					int isomerCount=chopAndPrint(fullSeq, samList, fout);
					if(count2Reads.get(isomerCount)==null) {
						ArrayList<String> tmp = new ArrayList<String>();
						tmp.add(readID);
						count2Reads.put(isomerCount, tmp);
					}else
						count2Reads.get(isomerCount).add(readID);
				}
				//re-initialize
				samList = new ArrayList<AlignmentRecord>();
				readID = myRec.readID;	
				fullSeq=new FastqSequence(Alphabet.DNA5(), rec);
//				fullSeq.print(fout);
				allLength.add(fullSeq.length());
			}
//			if (rec.getMappingQuality() < qual || !myRec.useful)	
			if (rec.getMappingQuality() < qual || Math.abs(myRec.readEnd-myRec.readStart) < .8*ref.length()) //80% error rate of nanopore data		
				continue;		
			
			samList.add(myRec);

		}// while
		
		ArrayList<Integer> counts = new ArrayList<Integer>(count2Reads.keySet());
		Collections.sort(counts);
		System.out.println("====================================================");
		System.out.println("====================================================");
		int cov=0;
		for(int c:counts) {
			System.out.printf("Number of %d-concatemer is %d\n",c,count2Reads.get(c).size());
			for(String id:count2Reads.get(c))
				System.out.println(id);
			cov+=c*count2Reads.get(c).size();
			System.out.println();
		}
		System.out.printf("Total=%dX\n",cov);
		fout.close();
		iter.close();
		reader.close();
		
		System.out.println("====================================================");
		System.out.println("## Length of all reads mapped:");
		Collections.sort(allLength);
		for(int l:allLength)
			System.out.printf("%d ", l);
	}
	
	private int chopAndPrint(FastqSequence read, ArrayList<AlignmentRecord> records, SequenceOutputStream out) throws IOException {
		if(read==null || records == null || records.size()==0)
			return 0;
		System.out.println("===========================================================");

		Sequence motif = new Sequence(Alphabet.DNA5(), "TGGTATCAGAGC", "template");
		CrossCorrelationScanner scanner = new CrossCorrelationScanner(motif);
				
		int[] 	start = new int[records.size()], 
				stop = new int[records.size()];
		
		int[] 	start2 = new int[records.size()], 
				stop2 = new int[records.size()];
		
		int count=0;
		for(AlignmentRecord rec:records) {
			start[count]=mapToRead(0, rec);
			stop[count]=mapToRead(ref.length()-1, rec);
			
			start2[count] = scanner.scanAround(read, start[count], rec.strand);
			stop2[count] = scanner.scanAround(read, stop[count], rec.strand);
			
			count++;
		}
		System.out.println("Read " + read.getName() + " length=" + read.length() + ":");
		for(int i=0;i<count;i++) {
			System.out.printf("\t%d -> %d After scan: %d -> %d\n", start[i], stop[i], start2[i], stop2[i]);
			FastqSequence isomer=read.subSequence(Math.min(start2[i], stop2[i]), Math.max(start2[i], stop2[i]));
			isomer.print(out);
		}
		return count;
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
//		SequenceChopper tony_tony_chopper = new SequenceChopper("/home/s_hoangnguyen/Projects/plant_virus/ref/BSMYV.fasta");
		SequenceChopper tony_tony_chopper = new SequenceChopper("/home/s_hoangnguyen/Projects/plant_virus/ref/CaMV.fasta");

		tony_tony_chopper.chopLongReadsBasedOnRef("/home/s_hoangnguyen/Projects/plant_virus/alignment/barcode08_pass_mm2.sam", "/home/s_hoangnguyen/Projects/plant_virus/barcode08_isomers.fastq", 1);

	}

}
