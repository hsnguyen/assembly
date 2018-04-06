package org.rtassembly.concatemer;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;

public class CrossCorrelationScanner{
	Sequence template;
	public CrossCorrelationScanner(Sequence template) {
		// TODO Auto-generated constructor stub
		this.template=template;

	}
	
	public ArrayList<Integer> scan(Sequence query) {
		ArrayList<Integer> retval = new ArrayList<>();
		int 	l1=template.length(),
				l2=query.length();
		/*
		 * 				0			l1-1
		 * 				|			|
		 * template	: 	-------------
		 * 							|
		 * 							i
		 * 							|
		 * query	:			--------------------------------|-------------|	
		 * 						|								|	
		 *						0								l2-1		l1+l2-2
		 *
		 *	Overlap: template (l1-i,l1) <=> query (q_start,i)
		 */
		
		for(int i=0;i<l1+l2-1;i++) {
			int 	q_start=(i>=l1?i-l1+1:0), //coordinate on query of the overlap starting point
					q_end=(i<=l2-1?i:l2-1),	 //coordinate on query of the overlap ending point	
					t_start=(i>=l1?0:l1-i-1), //coordinate on template of the overlap starting point
					t_end=(i<=l2-1?l1-1:l1+l2-i-2); //coordinate on template of the overlap ending point
			int score=0;
			while(q_start<=q_end && t_start<=t_end)
				score+=(query.charAt(q_start++)==template.charAt(t_start++)?1:0);
			retval.add(score);
		}
		
		return retval;
	}
	
	public static void main(String args[]) throws IOException {
		/***********************************************************************/
		String 	inputFileName = "/home/s_hoangnguyen/Projects/plant_virus/barcode08_pass.fastq.gz",
				outputFileName = "/home/s_hoangnguyen/Projects/plant_virus/barcode08.signal";
		
		
		Sequence motif = new Sequence(Alphabet.DNA5(), "TGGTATCAGAGC", "template");
		CrossCorrelationScanner scanner = new CrossCorrelationScanner(motif),
								scanner_rev = new CrossCorrelationScanner(Alphabet.DNA().complement(motif));
		SequenceReader reader = SequenceReader.getReader(inputFileName);
		Sequence seq;
		PrintWriter writer = new PrintWriter(new FileWriter(outputFileName)); 
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
//			System.out.println("");
			ArrayList<Integer> 	result = scanner.scan(seq),
								result_ret = scanner_rev.scan(seq);
			writer.println(">"+seq.getName());
			for(int sig:result)				
				writer.print(sig+" ");
			writer.println();
			for(int sig:result_ret)				
				writer.print(sig+" ");
			writer.println();
		}
		writer.close();
	}
}
