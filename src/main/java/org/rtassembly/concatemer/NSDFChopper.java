/*
 * Reference-free concatemers detection on base-called sequence
 * Using  Normalized  Square  Difference Function (NSDF)  as metric for
 * auto-correlation scanner.
 ** Based on http://www.cs.otago.ac.nz/tartini/papers/A_Smarter_Way_to_Find_Pitch.pdf
 */
package org.rtassembly.concatemer;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.stream.DoubleStream;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jtransforms.fft.DoubleFFT_1D;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import japsa.seq.Sequence;


public class NSDFChopper {
	String name;
	short[] signal;
	
	public NSDFChopper() {};
	public NSDFChopper(String f5File){
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
	
//	public ArrayList<Double> scan() {
//		ArrayList<Double> retval = new ArrayList<>();
//		int 	l=signal.length;
//		System.out.println("length="+l);
//		/*
//		 * 				0			l1-1
//		 * 				|			|
//		 * template	: 	-------------
//		 * 							|
//		 * 							i
//		 * 							|
//		 * query	:			--------------------------------|-------------|	
//		 * 						|								|	
//		 *						0								l2-1		l1+l2-2
//		 *
//		 *	Overlap: template (l1-i,l1) <=> query (q_start,i)
//		 */
//		
//		for(int i=0;i<2*l-1;i++) {
//			System.out.print("n="+i+"/"+(2*l)+"...");
//			int 	q_start=(i>=l?i-l+1:0), //coordinate on query of the overlap starting point
//					q_end=(i<=l-1?i:l-1),	 //coordinate on query of the overlap ending point	
//					t_start=(i>=l?0:l-i-1), //coordinate on template of the overlap starting point
//					t_end=(i<=l-1?l-1:2*l-i-2); //coordinate on template of the overlap ending point
//			double score=0.0;
//			int overlap=q_end-q_start;
//			if(overlap==0){
//				retval.add(0.0);
//			}else{
//				while(q_start<=q_end && t_start<=t_end)
//					score+=1.0/(1.0+Math.abs(signal[q_start++]-signal[t_start++]));
//				retval.add(score/overlap); //normalized it
//			}
//			System.out.println("f(n)="+(score/overlap));
//		}
//		
//		return retval;
//	}
//	
//	public void printCrossCorrelation(String outFile) throws IOException {
//		PrintWriter writer = new PrintWriter(new FileWriter(outFile)); 
//		writer.print(name+" ");
//
//		int 	l=signal.length;
////		System.out.println("length="+l);
//		/*
//		 * 				0			l1-1
//		 * 				|			|
//		 * template	: 	-------------
//		 * 							|
//		 * 							i
//		 * 							|
//		 * query	:			--------------------------------|-------------|	
//		 * 						|								|	
//		 *						0								l2-1		l1+l2-2
//		 *
//		 *	Overlap: template (l1-i,l1) <=> query (q_start,i)
//		 */
//		
//		for(int i=0;i<l-1;i++) {
//			System.out.print("n="+i+"/"+(2*l)+"...");
//			int 	q_start=(i>=l?i-l+1:0), //coordinate on query of the overlap starting point
//					q_end=(i<=l-1?i:l-1),	 //coordinate on query of the overlap ending point	
//					t_start=(i>=l?0:l-i-1), //coordinate on template of the overlap starting point
//					t_end=(i<=l-1?l-1:2*l-i-2); //coordinate on template of the overlap ending point
//			double score=0.0, value=0.0;
//			int overlap=q_end-q_start;
//			
//			if(overlap>0){
//				while(q_start<=q_end && t_start<=t_end)
//					score+=1.0/(1.0+Math.abs(signal[q_start++]-signal[t_start++]));
//				
//				value=score/overlap;
//			}
//			writer.printf("%.5f ",value); //normalized it
//
////			System.out.println("f(n)="+(score/overlap));
//		}
//		writer.println();
//		writer.close();
//	}
//	
//	
//    void print(String msg, double [] x) {
//        System.out.println(msg);
//        System.out.print("[ ");
//        for (double d : x) System.out.printf("%.2f ",d);
//        System.out.println("]");
//    }
//
//    /**
//     * This is a "wrapped" signal processing-style autocorrelation. 
//     * For "true" autocorrelation, the data must be zero padded.  
//     */
//    public double[] bruteForceAutoCorrelation(double [] x) {
//    	double[] retval = new double[x.length];
//        Arrays.fill(retval, 0);
//        int n = x.length;
//        for (int j = 0; j < n; j++) {
//            for (int i = 0; i < n; i++) {
////                retval[j] += x[i] * x[(n + i - j) % n];
//            	if(i>=j)
//            		retval[j] += x[i] * x[i-j];
//            }
//        }
//        return retval;
//    }

    private double sqr(double x) {
        return x * x;
    }
    //Running average
    //TODO: need another low pass filter
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
    //calculate auto-correlation using FFT on DNA sequence {1,-1,i,-i}}
    public double[] nuclAutoCorrelationFFT(Sequence nuclSeq) {
        int n = nuclSeq.length();
        double[] 	x2 = nuclSeq.seq2sig(),	//input signal (complex)
        			ac2 = new double[2*n], 	//autocorrelation in freq domain
        			retval = new double[n];	//result (time domain)
        
        
        DoubleFFT_1D fft = new DoubleFFT_1D(2*n);
        fft.complexForward(x2);
        
        //convolute in time domain = multiple in freq domain
        for (int i = 0; i < 2*n-1; i += 2) {
            ac2[i] = sqr(x2[i]) + sqr(x2[i+1]);
            ac2[i+1] = 0;
        }
        DoubleFFT_1D ifft = new DoubleFFT_1D(2*n); 
        ifft.realInverse(ac2, true);
        for(int i=0;i<n;i++)
        	retval[i]=ac2[i];
        
        return retval;

    }
    
    //calculate r'(t) = \sum{x_j*x_{j+t}} for real signal of raw data
    public double[] signalAutoCorrelationFFT(double [] rawSeq) {
        int n = rawSeq.length;
        double[] 	x2 = Arrays.copyOf(rawSeq, 2*n),
        			ac2 = new double[2*n],
        			retval = new double[n];
        DoubleFFT_1D fft = new DoubleFFT_1D(2*n);
        fft.realForward(x2);
        ac2[0] = sqr(x2[0]);
        ac2[1] = sqr(x2[1]);
        for (int i = 2; i < 2*n-1; i += 2) {
            ac2[i] = sqr(x2[i]) + sqr(x2[i+1]);
            ac2[i+1] = 0;
        }
        DoubleFFT_1D ifft = new DoubleFFT_1D(2*n); 
        ifft.realInverse(ac2, true);
        for(int i=0;i<n;i++)
        	retval[i]=ac2[i];
        
        return retval;

    }
    
    //calculate m'(t) = \sum{x_j^2+x_{j+t}^2} for real signal of raw data
    public double[] lagSquareSum(double[] x) {
    	double[] 	xsqr = new double[x.length],
    				retval = new double[x.length];
    	
    	Arrays.parallelSetAll(xsqr, i->sqr(x[i]));
    	double ssqr = DoubleStream.of(xsqr).parallel().sum();
    	
    	retval[0]=2*ssqr;
    	for(int i=1;i<x.length;i++) {
    		retval[i]=retval[i-1]-xsqr[i-1]-xsqr[x.length-i];
    	}
    	
    	return retval;
    }

    void printNSDFSignal() {
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
//        double [] ac1 = bruteForceAutoCorrelation(data);
//        System.out.printf("Done brute force autocorrelation in  %d secs\n", (System.currentTimeMillis()-curTime)/1000);
//        curTime=System.currentTimeMillis();
        
        double [] 	r = signalAutoCorrelationFFT(data),
        			m = lagSquareSum(data),
        			n = new double[data.length];
        Arrays.parallelSetAll(n, i->2*r[i]/m[i]) ;        
        
        double [] n1k = new double [data.length];
        smooth(n, n1k, 20000);
        
        System.out.printf("Done FFT autocorrelation in  %d secs\n", (System.currentTimeMillis()-curTime)/1000);
        curTime=System.currentTimeMillis();
//        Print to file
        try {
    		PrintWriter writer = new PrintWriter(new FileWriter(DATA+"concat7_mpm.signal.xls"));
    		writer.print(name+" ");
    		for(double value:n)
    			writer.printf("%.5f\n",value); //normalized it
    		writer.close();
    		
    		PrintWriter writer2 = new PrintWriter(new FileWriter(DATA+"concat7_mpm_20kstep.xls")); 
    		for(double value:n1k)
    			writer2.printf("%.5f\n",value); //normalized it
    		writer2.close();
    		
            System.out.printf("Done writing to files in %d secs\n", (System.currentTimeMillis()-curTime)/1000);
            curTime=System.currentTimeMillis();
            
        } catch (IOException e) {
        		throw new RuntimeException(e);
        } 
        

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
    
//    void test() {
//        double [] 	data = {2,3,-1};
//        double [] ac1 = bruteForceAutoCorrelation(data);
//        print("Brute force ACF: ",ac1);
//        
//        double [] ac2 = fftAutoCorrelation(data);
//        print("FFT ACF: ",ac2);
////        double [] ac2_smt = new double [data.length];
////        smooth(ac2, ac2_smt, 1000);
//        print("test", lagSquareSum(data));
//        
//    }
    static String DATA="/home/sonhoanghguyen/Projects/concatemers/data/raw/";
	public static void main(String[] args){
		NSDFChopper concat=new NSDFChopper(DATA+"imb17_013486_20171130__MN17279_sequencing_run_20171130_Ha_BSV_CaMV1_RBarcode_35740_read_17553_ch_455_strand.fast5");
//		try {
//			concat.printCrossCorrelation(DATA+"concat7.signal");
//		} catch (IOException e) {
//			e.printStackTrace();
//		} 
		concat.printNSDFSignal();

	}

}
