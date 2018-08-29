package org.rtassembly.experiment;

import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jtransforms.fft.DoubleFFT_1D;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;

public class RawConcatemer {
	String name;
	short[] signal;
	public RawConcatemer(String f5File){
		try {			
			IHDF5Reader reader = HDF5Factory.openForReading(f5File);
			//find data here: "Raw/Reads/Read_17553/Signal"
			String pathToReadData="Raw/Reads/"+reader.getGroupMembers("Raw/Reads").get(0);
			name=reader.string().getAttr(pathToReadData,"read_id");
			signal=reader.int16().readArray(pathToReadData+"/Signal");
			reader.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public ArrayList<Double> scan() {
		ArrayList<Double> retval = new ArrayList<>();
		int 	l=signal.length;
		System.out.println("length="+l);
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
		
		for(int i=0;i<2*l-1;i++) {
			System.out.print("n="+i+"/"+(2*l)+"...");
			int 	q_start=(i>=l?i-l+1:0), //coordinate on query of the overlap starting point
					q_end=(i<=l-1?i:l-1),	 //coordinate on query of the overlap ending point	
					t_start=(i>=l?0:l-i-1), //coordinate on template of the overlap starting point
					t_end=(i<=l-1?l-1:2*l-i-2); //coordinate on template of the overlap ending point
			double score=0.0;
			int overlap=q_end-q_start;
			if(overlap==0){
				retval.add(0.0);
			}else{
				while(q_start<=q_end && t_start<=t_end)
					score+=1.0/(1.0+Math.abs(signal[q_start++]-signal[t_start++]));
				retval.add(score/overlap); //normalized it
			}
			System.out.println("f(n)="+(score/overlap));
		}
		
		return retval;
	}
	
	public void printCrossCorrelation(String outFile) throws IOException {
		PrintWriter writer = new PrintWriter(new FileWriter(outFile)); 
		writer.print(name+" ");

		int 	l=signal.length;
//		System.out.println("length="+l);
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
		
		for(int i=0;i<l-1;i++) {
			System.out.print("n="+i+"/"+(2*l)+"...");
			int 	q_start=(i>=l?i-l+1:0), //coordinate on query of the overlap starting point
					q_end=(i<=l-1?i:l-1),	 //coordinate on query of the overlap ending point	
					t_start=(i>=l?0:l-i-1), //coordinate on template of the overlap starting point
					t_end=(i<=l-1?l-1:2*l-i-2); //coordinate on template of the overlap ending point
			double score=0.0, value=0.0;
			int overlap=q_end-q_start;
			
			if(overlap>0){
				while(q_start<=q_end && t_start<=t_end)
					score+=1.0/(1.0+Math.abs(signal[q_start++]-signal[t_start++]));
				
				value=score/overlap;
			}
			writer.printf("%.5f ",value); //normalized it

//			System.out.println("f(n)="+(score/overlap));
		}
		writer.println();
		writer.close();
	}
	
	
    void print(String msg, double [] x) {
        System.out.println(msg);
        for (double d : x) System.out.println(d);
    }

    /**
     * This is a "wrapped" signal processing-style autocorrelation. 
     * For "true" autocorrelation, the data must be zero padded.  
     */
    public void bruteForceAutoCorrelation(double [] x, double [] ac) {
        Arrays.fill(ac, 0);
        int n = x.length;
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                ac[j] += x[i] * x[(n + i - j) % n];
            }
        }
    }

    private double sqr(double x) {
        return x * x;
    }
    public void smooth(double[] in, double[] out, int window){
    	double beg=in[0], end=in[window];
    	for(int i=0;i<window;i++){
    		out[0]+=in[i];
    	}
    	out[0]/=window;
    	for(int i=1;i<in.length;i++){
    		beg=in[i-1];
    		end=in[(i+window-1)%in.length];
    		out[i]=out[i-1]+ (end-beg)/window;
    	}
    }
    public void fftAutoCorrelation(double [] x, double [] ac) {
        int n = x.length;
        // Assumes n is even.
        DoubleFFT_1D fft = new DoubleFFT_1D(n);
        fft.realForward(x);
        ac[0] = sqr(x[0]);
        // ac[0] = 0;  // For statistical convention, zero out the mean 
        ac[1] = sqr(x[1]);
        for (int i = 2; i < n-1; i += 2) {
            ac[i] = sqr(x[i]) + sqr(x[i+1]);
            ac[i+1] = 0;
        }
        DoubleFFT_1D ifft = new DoubleFFT_1D(n); 
        ifft.realInverse(ac, true);
        // For statistical convention, normalize by dividing through with variance
        for (int i = 1; i < n; i++)
            ac[i] /= ac[0];
        ac[0] = 1;
    }

    void test() {
        double [] data = new double [signal.length];
        SummaryStatistics stats=new SummaryStatistics();
        long curTime=System.currentTimeMillis();

        for (int j=0;j<signal.length;j++) {
            data[j] = (double)signal[j];
            stats.addValue(signal[j]);
        }
        double 	mean=stats.getMean(),
        		std=stats.getStandardDeviation();
        for (int j=0;j<data.length;j++) {
            data[j] = (data[j]-mean)/std;
        }
        System.out.printf("Done init arrays in %d secs\n", (System.currentTimeMillis()-curTime)/1000);
        curTime=System.currentTimeMillis();

        double [] ac1 = new double [data.length];
        double [] ac2 = new double [data.length];
        bruteForceAutoCorrelation(data, ac1);
        System.out.printf("Done brute force autocorrelation in  %d secs\n", (System.currentTimeMillis()-curTime)/1000);
        curTime=System.currentTimeMillis();
        
        fftAutoCorrelation(data, ac2);
        double [] ac2_smt = new double [data.length];
        smooth(ac2, ac2_smt, 1000);
        
        System.out.printf("Done FFT autocorrelation in  %d secs\n", (System.currentTimeMillis()-curTime)/1000);
        curTime=System.currentTimeMillis();
//        Print to file
        try {
    		PrintWriter writer = new PrintWriter(new FileWriter("/home/sonhoanghguyen/Projects/concatemers/data/raw/concat7_fft.signal"));
    		writer.print(name+" ");
    		for(double value:ac2)
    			writer.printf("%.5f ",value); //normalized it
    		writer.close();
    		
    		PrintWriter writer2 = new PrintWriter(new FileWriter("/home/sonhoanghguyen/Projects/concatemers/data/raw/concat7_fft_1kstep.txt")); 
    		for(double value:ac2_smt)
    			writer2.printf("%.5f ",value); //normalized it
    		writer2.close();
        } catch (IOException e) {
        		throw new RuntimeException(e);
        } 
        
        System.out.printf("Done writing to files in %d secs\n", (System.currentTimeMillis()-curTime)/1000);
        curTime=System.currentTimeMillis();
    }
    
//    /* Find autocorrelation peaks */
//    public List<Integer> findPeaks() {
//        List<Integer> peaks = new ArrayList<>();
//        int max = 1;
//        maxACF = 0;
//        if (correlations.length > 1) {
//            boolean positive = (correlations[1] > correlations[0]);
//            for (int i = 2; i < correlations.length; i++) {
//                if (!positive && correlations[i] > correlations[i - 1]) {
//                    max = i;
//                    positive = !positive;
//                } else if (positive && correlations[i] > correlations[max]) {
//                    max = i;
//                } else if (positive && correlations[i] < correlations[i - 1]) {
//                    if (max > 1 && correlations[max] > ACF_THRESH) {
//                        peaks.add(max);
//                        if (correlations[max] > maxACF) { maxACF = correlations[max]; }
//                    }
//                    positive = !positive;
//                }
//            }
//        }
//        return peaks;
//    }
	public static void main(String[] args){
		RawConcatemer concat=new RawConcatemer("/home/sonhoanghguyen/Projects/concatemers/data/raw/imb17_013486_20171130__MN17279_sequencing_run_20171130_Ha_BSV_CaMV1_RBarcode_35740_read_17553_ch_455_strand.fast5");
//		try {
//			concat.printCrossCorrelation("/home/sonhoanghguyen/Projects/concatemers/data/raw/concat7.signal");
//		} catch (IOException e) {
//			e.printStackTrace();
//		} 
		concat.test();

	}

}
