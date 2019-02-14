/*
 * Reference-free concatemers detection on base-called sequence
 * Using Normalized Kronecker delta function(NKDF)  as metric for
 * auto-correlation scanner.
 */
package org.rtassembly.concatemer;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import japsa.seq.Alphabet;
import japsa.seq.Alphabet.DNA;
import japsa.seq.Sequence;
import japsa.seq.SequenceReader;

public class NKDFChopper{
	public static final int SCAN_WINDOW=100;
	Sequence template;
	
	public NKDFChopper(Sequence template) {
		// TODO Auto-generated constructor stub
		this.template=template;

	}
	
	public ArrayList<Double> scan(Sequence query) {
		ArrayList<Double> retval = new ArrayList<>();
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
			double score=0.0;
			int overlap=q_end-q_start;
			if(overlap==0){
				retval.add(0.0);
			}else{
				while(q_start<=q_end && t_start<=t_end)
					score+=(query.charAt(q_start++)==template.charAt(t_start++)?1:0);
				retval.add(score/overlap); //normalized it
			}
		}
		
		return retval;
	}
	public int scanAround(Sequence query, int coordinate, boolean strand) {
		int retval=-1;
		int l=query.length();
		int lower=(coordinate-SCAN_WINDOW)>0?(coordinate-SCAN_WINDOW):0,
			upper=(coordinate+SCAN_WINDOW)<l?(coordinate+SCAN_WINDOW):(l-1);
		int newcoord = coordinate-lower+1;
		Sequence target=query.subSequence(lower, upper);
		if(strand) {
			target=DNA.complement(target);
			newcoord=target.length()-newcoord;
		}
		ArrayList<Double> spectrum=scan(target);
		double max_score=0;
		int min_dist=l;
		for(int i=0;i<spectrum.size();i++) {
			int distance=Math.abs(newcoord-i);
			
			if(spectrum.get(i)>max_score || (spectrum.get(i) == max_score && distance < min_dist)) {
				max_score=spectrum.get(i);
				retval=lower+i;
				min_dist = distance;
			}
		}
//		System.out.printf("Max and closest hit at %d, distance=%d, score=%d\n",retval, min_dist,max_score);
		return retval;
	}
	
//  //calculate auto-correlation using FFT on DNA model sequence {1,-1,i,-i}}. 
//	//copy to NSDFChopper to test>Not used!!!
//  public double[] nuclAutoCorrelationFFT(double[] nuclSeq) throws IOException {
//      int n = nuclSeq.length/2; //nuclSeq length is even!
//      double[] 	x2 = nuclSeq,	//input signal (complex)
//      			ac2 = new double[2*n], 	//autocorrelation in freq domain
//      			retval = new double[n];	//result (time domain)
//      
//      System.out.println("Length n=" + x2.length);
//      DoubleFFT_1D fft = new DoubleFFT_1D(n);
//      fft.complexForward(x2);
//      
//      //convolute in time domain = multiple in freq domain
//		PrintWriter writer = new PrintWriter(new FileWriter(DATA+"concat7_dft.signal.xls"));
//		writer.print(name+" ");
//		
//      for (int i = 0; i < 2*n-1; i += 2) {
//          ac2[i] = sqr(x2[i]) + sqr(x2[i+1]);
//          ac2[i+1] = 0;
//          writer.printf("%.5f\n",ac2[i]); //normalized it
//      }
//      
//		writer.close();
//		
//      DoubleFFT_1D ifft = new DoubleFFT_1D(n); 
//      ifft.complexInverse(ac2, true);
//      for(int i=0;i<n;i++)
//      	retval[i]=ac2[2*i];
//      
//      return retval;
//
//  }
  
	public static void main(String args[]) throws IOException {
		/***********************************************************************/
		String 	inputFileName = "/home/sonhoanghguyen/Projects/concatemers/data/jie/B02_A4_allpass.fastq",
				outputFileName = "/home/sonhoanghguyen/Projects/concatemers/data/jie/B02_A4_allpass.signal";
		
		
		SequenceReader reader = SequenceReader.getReader(inputFileName);
		Sequence seq;
		PrintWriter writer = new PrintWriter(new FileWriter(outputFileName)); 
		int count=1;
		while ((seq = reader.nextSequence(Alphabet.DNA())) != null){
			String name=seq.getName().split(" ")[0];
			System.out.println("Read number " + count++ + ": " + name + " length= " + seq.length() + "...");
			NKDFChopper scanner = new NKDFChopper(seq);
			ArrayList<Double> 	result = scanner.scan(seq); //scan to itself
			writer.print(name+" ");
			for(double sig:result)				
				writer.printf("%.5f ",sig);
			writer.println();
		}
		
		writer.close();
	}
}
