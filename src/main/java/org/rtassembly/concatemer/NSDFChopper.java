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
import japsa.seq.Alphabet;
import japsa.seq.FastqReader;
import japsa.seq.Sequence;


public class NSDFChopper {
	String name;
	short[] signal;
	
	public NSDFChopper(Sequence seq) {
//		signal=seq.seq2sig();
		signal=new short[seq.length()];
		for(int i=0;i<seq.length();i++)
			signal[i]=seq.getBase(i);
		if(signal==null) //unsupported DNA sequence (A,C,G,T only)
			System.exit(1);
		else
			name=seq.getName();
	}
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

    //calculate r'(t) = \sum{x_j*x_{j+t}} for real signal of raw data
    public double[] signalAutoCorrelationFFT(double [] rawSeq) {
        int n = rawSeq.length;
        double cutFreq=100;
        double[] 	x2 = Arrays.copyOf(rawSeq, 2*n),
        			ac2 = new double[2*n],
        			retval = new double[n];
        DoubleFFT_1D fft = new DoubleFFT_1D(2*n);
        fft.realForward(x2);
        ac2[0] = sqr(x2[0]);
        ac2[1] = sqr(x2[1]);
        for (int i = 2; i < 2*n-1; i += 2) {
        	if(i<cutFreq)
        		ac2[i] = sqr(x2[i]) + sqr(x2[i+1]);
        	else
        		ac2[i] = 0;//simple
            ac2[i+1] = 0;
        }
        DoubleFFT_1D ifft = new DoubleFFT_1D(2*n); 
        ifft.realInverse(ac2, true);
        for(int i=0;i<n;i++)
        	retval[i]=ac2[i];
        
        return retval;

    }
    
    /************************************************************
     *  Try to use window function instead of zeroing freq domain 
     *  Not work yet...
     ************************************************************/
	private double sinc(double x) {
		if (x == 0)
			return 1;
		return Math.sin(Math.PI*x)/(Math.PI*x);
	}
    //calculate r'(t) = \sum{x_j*x_{j+t}} for real signal of raw data.
    public double[] signalAutoCorrelationFFT2(double [] signal) {
        int n = signal.length;
        double cutFreq=100.0/(double)n;
        double[] 	filter = new double[n],
        			filterFFT = new double[n*2],
        			signalFFT = new double[n*2],
        			signalFiltered = new double[n];
        
        DoubleFFT_1D fft = new DoubleFFT_1D(n);
    
		for (int i = 0; i < n; i++) {
			// see http://en.wikipedia.org/wiki/Sinc_filter
			double sincFilter = 2*cutFreq*sinc(2*cutFreq*(i-(n-1)/2.0)/(double)(n-1));
			// applying a Blackman window
			filter[i] = (0.42-0.5*Math.cos((2*Math.PI*i)/(double)(n-1))+0.08*Math.cos((4*Math.PI*i)/(double)(n-1))) * sincFilter;
		}

		// copying time domain filter data to the FFT buffer
		for (int i = 0; i < n; i++){
			filterFFT[2*i] = filter[i];
			filterFFT[2*i+1] = 0;
		}


		fft.complexForward(filterFFT);
		
		for (int i = 0; i < signal.length; i++){
			signalFFT[2 * i] = signal[i];
			signalFFT[2*i+1] = 0;
		}
		// calculating the fft of the data
		fft.complexForward(signalFFT);
		// pointwise multiplication of the filter and audio data in the frequency domain
		for (int i = 0; i < filterFFT.length; i += 2) {
			double temp = signalFFT[i] * filterFFT[i] - signalFFT[i+1]*filterFFT[i+1];
			signalFFT[i+1] = signalFFT[i] * filterFFT[i+1] + signalFFT[i+1] * filterFFT[i]; // imaginary part
			signalFFT[i] = temp; // real part
			System.out.println(signalFFT[i] + "," + signalFFT[i+1] + ",");
		}

		fft.complexInverse(signalFFT, true);
        
		for(int i=0;i<n;i++)
			signalFiltered[i]=signalFFT[2*i];
        return signalFiltered;

    }
    /******************************************************************/
    
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

    void printRawSignal(int window) throws IOException {
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
        
        double [] 	r = signalAutoCorrelationFFT(data), //TODO:signalAutoCorrelationFFT2
    				m = lagSquareSum(data),
    				n = new double[data.length];
        Arrays.parallelSetAll(n, i->2*r[i]/m[i]) ;        
        double [] n1k = new double [data.length];
        smooth(n, n1k, window);
        
        System.out.printf("Done FFT autocorrelation in  %d secs\n", (System.currentTimeMillis()-curTime)/1000);
        curTime=System.currentTimeMillis();
//        Print to file
		PrintWriter writer = new PrintWriter(new FileWriter(DATA+"concat7_mpm.signal.xls"));
		writer.print(name+" ");
		for(double value:n)
			writer.printf("%.5f\n",value); //normalized it
		writer.close();
		
		PrintWriter writer2 = new PrintWriter(new FileWriter(DATA+"concat7_mpm_"+window+".xls")); 
		for(double value:n1k)
			writer2.printf("%.5f\n",value); //normalized it
		writer2.close();
		
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
    

//    static String DATA="/home/sonhoanghguyen/Projects/concatemers/data/raw/";
    static String DATA="/home/sonhoanghguyen/Projects/concatemers/data/test/";
	public static void main(String[] args){
		try {
			FastqReader reader = new FastqReader(DATA+"concat7.fastq");
			Sequence seq=reader.nextSequence(Alphabet.DNA4());
			NSDFChopper concat = new NSDFChopper(seq);
//			NSDFChopper concat=new NSDFChopper(DATA+"imb17_013486_20171130__MN17279_sequencing_run_20171130_Ha_BSV_CaMV1_RBarcode_35740_read_17553_ch_455_strand.fast5");
			concat.printRawSignal(100);
			
			reader.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}


	}

}
