package org.rtassembly.concatemer;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Set;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import japsa.seq.Alphabet;
import japsa.seq.FastqSequence;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;

import org.rtassembly.npscarf.AlignmentRecord;
import org.rtassembly.npscarf.Contig;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

public class ReferenceBasedChopper {
	private static final Logger LOG = LoggerFactory.getLogger(ReferenceBasedChopper.class);

	HashMap<String, Contig> refMap = new HashMap<String, Contig>();
	ReferenceBasedChopper(String refFile) throws IOException{
		SequenceReader reader = SequenceReader.getReader(refFile);
		Sequence ref;
		int id=1;
		while((ref=reader.nextSequence(Alphabet.DNA()))!=null){
			refMap.put(ref.getName(), new Contig(id++, ref));
		}
		reader.close();
	}
	// Chop the concatemers (long nanopore reads) based on the reference
	public void chopLongReadsBasedOnRef(String samFile, String outFolder, double qual) throws IOException {
		HashMap<Integer,ArrayList<String>> count2Reads = new HashMap<>(); // concatemer information
		ArrayList<Integer> allLength = new ArrayList<Integer>();
		
		HashMap<String, SequenceOutputStream> outMap = new HashMap<>();
		Set<String> refs = refMap.keySet();
		for(String r:refs){
			String fileName=outFolder+File.separator+r+".fastq";
			SequenceOutputStream fout = SequenceOutputStream.makeOutputStream(fileName);
			outMap.put(r, fout);
		}
		
		
		SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
		SamReader reader;
		if ("-".equals(samFile))
			reader = SamReaderFactory.makeDefault().open(SamInputResource.of(System.in));
		else
			reader = SamReaderFactory.makeDefault().open(new File(samFile));
		
		SAMRecordIterator iter = reader.iterator();
		
		
		String readID = "";
		AlignmentRecord algRec = null;
		ArrayList<AlignmentRecord> alignmentList = null;// alignment record of the same read;		
		
		Contig refContig;//adapt to the AlignmentRecord API
		FastqSequence fullSeq = null;
		while (iter.hasNext()) {
			SAMRecord samRecord = iter.next();
			if(samRecord.getReadUnmappedFlag())
				continue;
			refContig=refMap.get(samRecord.getReferenceName());
			if(refContig==null){
				LOG.error("Cannot find {} from reference list!", samRecord.getReferenceName());
				continue;
			}
				
			algRec = new AlignmentRecord(samRecord, refContig);
			
			if (!readID.equals(algRec.readID)){
				//chop and print seq ...
				if(alignmentList!=null && alignmentList.size() >= 1) {
					int isomerCount=chopAndPrint(fullSeq, refContig, alignmentList, outMap.get(samRecord.getReferenceName()));
					if(count2Reads.get(isomerCount)==null) {
						ArrayList<String> tmp = new ArrayList<String>();
						tmp.add(readID);
						count2Reads.put(isomerCount, tmp);
					}else
						count2Reads.get(isomerCount).add(readID);
				}
				//re-initialize
				alignmentList = new ArrayList<AlignmentRecord>();
				readID = algRec.readID;	
				fullSeq=new FastqSequence(Alphabet.DNA(), samRecord);
//				fullSeq.print(fout);
				allLength.add(fullSeq.length());
			}
//			if (rec.getMappingQuality() < qual || !myRec.useful)	
			if (samRecord.getMappingQuality() < qual || Math.abs(algRec.readEnd-algRec.readStart) < .8*refContig.length()) //80% error rate of nanopore data		
				continue;		
			
			alignmentList.add(algRec);

		}// while
		
		ArrayList<Integer> counts = new ArrayList<Integer>(count2Reads.keySet());
		Collections.sort(counts);
		int cov=0;
		for(int c:counts) {
			LOG.info("\nNumber of {}-concatemer is {}",c,count2Reads.get(c).size());
			for(String id:count2Reads.get(c))
				LOG.info(id);
			cov+=c*count2Reads.get(c).size();
		}
		LOG.info("Total={}X\n",cov);
		for(SequenceOutputStream os:outMap.values()){
			os.close();
		}
		iter.close();
		reader.close();
		
		//Print mapped reads length to draw histogram
//		System.out.println("====================================================");
//		System.out.println("## Length of all reads mapped:");
//		Collections.sort(allLength);
//		for(int l:allLength)
//			System.out.printf("%d ", l);
	}
	
	private int chopAndPrint(FastqSequence read, Contig ref, ArrayList<AlignmentRecord> records, SequenceOutputStream out) throws IOException {
		if(read==null || records == null || records.size()==0)
			return 0;
		System.out.println("===========================================================");

		Sequence motif = new Sequence(Alphabet.DNA(), "TGGTATCAGAGC", "template");
		NKDFChopper scanner = new NKDFChopper(motif);
				
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

	
	public static void main(String[] args) throws IOException {
//		SequenceChopper tony_tony_chopper = new SequenceChopper("/home/s_hoangnguyen/Projects/plant_virus/ref/BSMYV.fasta");
		ReferenceBasedChopper tony_tony_chopper = new ReferenceBasedChopper("/home/sonhoanghguyen/Projects/concatemers/ref.fasta");
		tony_tony_chopper.chopLongReadsBasedOnRef("/home/sonhoanghguyen/Projects/concatemers/all_pass_mm2.bam", "/home/sonhoanghguyen/Projects/concatemers/monomers/", 1);

	}

}
