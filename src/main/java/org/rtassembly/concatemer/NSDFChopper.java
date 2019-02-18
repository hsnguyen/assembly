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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.stream.DoubleStream;

import org.apache.commons.math3.stat.descriptive.SummaryStatistics;
import org.jtransforms.fft.DoubleFFT_1D;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import ch.systemsx.cisd.hdf5.HDF5Factory;
import ch.systemsx.cisd.hdf5.IHDF5Reader;
import japsa.seq.Alphabet;
import japsa.seq.FastqReader;
import japsa.seq.Sequence;


public class NSDFChopper {
    private static final Logger LOG = LoggerFactory.getLogger(NSDFChopper.class);

	String name;
	short[] signal;
	DoubleFFT_1D fft;
    double cutFreq=100.0;
    static int MIN_MONOMER=2000;
    static int WINDOW=2; //running average window
    static double THRES=.5; //ignore peaks smaller than half of the max
    static double ERROR=.2; //error rate of nanopore data

    double [] r, m, n; //follow notation in McLeod Pitch Method
    List<Integer> chopper;
    public NSDFChopper(){}
	public NSDFChopper(Sequence seq) {
		this();
		setData(seq);
	}
	
	public NSDFChopper(String f5File){
		this();
		setData(f5File);
	}
	
	public void setData(Sequence seq) {
//		signal=seq.seq2sig();
		signal=new short[seq.length()];
		for(int i=0;i<seq.length();i++){
			switch(seq.getBase(i)){
				case 0: //A <=> -2
					signal[i]=-2;
					break;
				case 1: //C <=> -1
					signal[i]=-1;
					break;
				case 2: //G <=> 1
					signal[i]=1;
					break;
				case 3: //T <=> 2
					signal[i]=2;
					break;
				default:
					LOG.error("Invalid DNA character (only ACGT)!");
					System.exit(1);
			}
		
		}
		
		if(signal==null) //unsupported DNA sequence (A,C,G,T only)
			System.exit(1);
		else
			name=seq.getName();
		
        r = new double[signal.length];
    	m = new double[signal.length];
    	n = new double[signal.length];
    	chopper=new ArrayList<Integer>();
	}
	
	public void setData(String f5File){
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
		
        r = new double[signal.length];
    	m = new double[signal.length];
    	n = new double[signal.length];
    	chopper=new ArrayList<Integer>();
	}
	
    private double sqr(double x) {
        return x * x;
    }
    
    //calculate r'(t) = \sum{x_j*x_{j+t}} for real signal of raw data
    public void signalAutoCorrelationFFT(double[] rawSeq, double[] retval) {
        int n = rawSeq.length;
        double[] 	signalFFT = Arrays.copyOf(rawSeq, 2*n),
        			acFFT = new double[2*n];
        fft = new DoubleFFT_1D(2*n);
        fft.realForward(signalFFT);
        
        
        
        acFFT[0] = sqr(signalFFT[0]);
        acFFT[1] = sqr(signalFFT[1]);
        for (int i = 2; i < 2*n-1; i += 2) {
        	acFFT[i] = sqr(signalFFT[i]) + sqr(signalFFT[i+1]);
            acFFT[i+1] = 0;
        }

        fft.realInverse(acFFT, true);
        for(int i=0;i<n;i++)
        	retval[i]=acFFT[i];
        
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
    public void lowPassFilter(double[] signal) {
        int n = signal.length;
        double[] 	filter = new double[n],
        			filterFFT = new double[n*2],
        			signalFFT = new double[n*2];
        
        fft = new DoubleFFT_1D(n);
        double ft=cutFreq/(double)n; //normalised transition frequency
		for (int i = 0; i < n; i++) {
			// see http://en.wikipedia.org/wiki/Sinc_filter
			double sincFilter = 2*ft*sinc(2*ft*(i-(n-1)/2.0));
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
			signalFFT[2*i] = signal[i];
			signalFFT[2*i+1] = 0;
		}
		// calculating the fft of the data
		fft.complexForward(signalFFT);
		// pointwise multiplication of the filter and audio data in the frequency domain
		for (int i = 0; i < filterFFT.length; i += 2) {
			double temp = signalFFT[i] * filterFFT[i] - signalFFT[i+1]*filterFFT[i+1];
			signalFFT[i+1] = signalFFT[i] * filterFFT[i+1] + signalFFT[i+1] * filterFFT[i]; // imaginary part
			signalFFT[i] = temp; // real part
//			System.out.println(i/2 + " : " + signalFFT[i] + "," + signalFFT[i+1] + ",");
		}

		fft.complexInverse(signalFFT, true);
        
		for(int i=0;i<n;i++)
			signal[(i<(n-1)/2)?(i+(n+1)/2):(i-(n-1)/2)]=signalFFT[2*i]; //phase response of recursive filter: http://www.dspguide.com/ch19/4.htm

    }
    /******************************************************************/
    
    //calculate m'(t) = \sum{x_j^2+x_{j+t}^2} for real signal of raw data
    public void lagSquareSum(double[] x, double[] retval) {
    	double[] 	xsqr = new double[x.length];
    	
    	Arrays.parallelSetAll(xsqr, i->sqr(x[i]));
    	double ssqr = DoubleStream.of(xsqr).parallel().sum();
    	
    	retval[0]=2*ssqr;
    	for(int i=1;i<x.length;i++) {
    		retval[i]=retval[i-1]-xsqr[i-1]-xsqr[x.length-i];
    	}
    	
    }

    void concatemersDetection() throws IOException {
    	LOG.info("====================================================================");
    	/******************************************
    	 * 1. Calculate NSDF by FFT
    	 *****************************************/
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
               
        signalAutoCorrelationFFT(data, r);
        lagSquareSum(data,m);
        Arrays.parallelSetAll(n, i->2*r[i]/m[i]) ;        
//        Print to file
//		printArrayToFile(n,DATA+"mpm.signal.xls");
		lowPassFilter(n);
//		printArrayToFile(n, DATA+"mpm_lpf.xls"); 
	
        LOG.info("Done NSDF calculation in {} secs", (System.currentTimeMillis()-curTime)/1000);
        curTime=System.currentTimeMillis();
            
    	/******************************************
    	 * 2. Find peaks to chop
    	 *****************************************/
        findPeaks();
        LOG.info("Done peak picking in {} secs", (System.currentTimeMillis()-curTime)/1000);
        if(chopper==null||chopper.isEmpty())
        	LOG.info("Processing read {} length={}: not a concatemer!", name, signal.length);
        else
        	LOG.info("Processing read {} length={}: {}-concatemers!",name, signal.length,chopper.size());
        
    }
    
    @SuppressWarnings("unused")
	private void printArrayToFile(double[] input,  String filename){
		try {
			PrintWriter writer = new PrintWriter(new FileWriter(filename));
			for(double value:input)
				writer.printf("%.6f\n",value);
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

    }

    /* Find autocorrelation peaks */
    private void findPeaks() {
        List<Integer> 	peaks = new ArrayList<>(),
        				candidates = new ArrayList<>();
//        HashMap<Integer,Double> candidates = new HashMap<>();
        int localMaxIdx = -1;
        double 	localMax=0.0, globalMax = 0.0;
                
        if (n.length > 1) {
            boolean positive = false; //positively sloped zero crossing
            
            for (int i = MIN_MONOMER; i < n.length*(1-1/cutFreq); i++) {
            	if(n[i-1] <= 0 && n[i] > 0){
            		if(positive && localMaxIdx>0){
//            			candidates.put(maxIdx, localMax);
            			candidates.add(localMaxIdx);
            			if(globalMax<localMax){
            				globalMax=localMax;
            			}
            		}
            		
            		positive=true;
            		localMaxIdx=-1;
            		localMax=0.0;
            	}
            	
            	if(n[i] > localMax){
            		localMaxIdx=i;
            		localMax=n[i];
            	}
            		
            }
        }
        List<Integer> prominentPeaks = new ArrayList<>();
        double sum=0.0;
//		System.out.println("Candidates peaks:");

        for(int i:candidates){
        	if(n[i] > THRES*globalMax){
        		prominentPeaks.add(i);
            	sum+=n[i];
//        		System.out.println(i);
        	}
        }
        
        //loop over prominent peaks to find chopping coordinates
        int period=1;
    	for(int i=0;i<prominentPeaks.size();i++){
    		period=prominentPeaks.get(i);//try this as the monomer length
    		int sIdx=i;
    		//searchNextPeak(period, i, prominenPeaks)
    		peaks=new ArrayList<>();
    		peaks.add(period);
//    		System.out.printf("...adding %d to 0\n", period);
			int k=1; //peaks.get(k) = coord of k-th monomer
			int prevPeakCoord=peaks.get(k-1);
			while(++sIdx<prominentPeaks.size()){
				int coord=prominentPeaks.get(sIdx);
				int distance=coord-prevPeakCoord;
				if(Math.abs(distance-period) < period*ERROR){ //legal coord for next peak
					if(peaks.size()<k+1){
//						System.out.printf("...adding %d to %d\n", coord, k);
						peaks.add(coord);
					}else if(n[peaks.get(k)] < n[coord]){
//						System.out.printf("...replacing %d to %d\n", coord, k);
						peaks.set(k, coord);
					}

				}else if(distance > period){
					if(peaks.size()<k+1){
//						System.out.printf("Interrupted at %d: distance=%d > period=%d\n",coord, distance,period);
						break; //not found next consecutive peak
					}
					else{
						period=peaks.get(k)-prevPeakCoord;
						prevPeakCoord=peaks.get(k++);
						sIdx--;
					}
				}
				
			}

    		//check if sum of spikes from peaks is greater than half of all candidates
			double psum=peaks.stream().map(t->n[t]).reduce(0.0, Double::sum);
//			System.out.println("psum=" + psum + " sum="+sum);
			if(psum > sum*(1-THRES) && peaks.size()==n.length/peaks.get(0)){ //only the found peaks are significant ones
				break;
			}
			else 
				continue;
    	}    
//        System.out.printf("Detected %d-concatemers!\n", peaks.size());
//        for(int i:peaks)
//        	System.out.println(i);
    	chopper=peaks;
    }
    

	public static void main(String[] args){
		NSDFChopper tony = new NSDFChopper();
		HashMap<Integer, ArrayList<String>> histogram = new HashMap<>();
		try {
			FastqReader reader = new FastqReader(args[0]);
			Sequence seq=null;
			while((seq=reader.nextSequence(Alphabet.DNA4()))!=null){
				tony.setData(seq);
				tony.concatemersDetection();
				if(tony.chopper!=null && !tony.chopper.isEmpty()){
					int key=tony.chopper.size();
					ArrayList<String> clist=histogram.get(key);
					if(clist==null){
						clist=new ArrayList<>();
						histogram.put(key, clist);
					}
					clist.add(tony.name);
						
				}
			}

			reader.close();
			
			PrintWriter writer = new PrintWriter(new FileWriter(args[1]));
			for(int key:histogram.keySet()){
				writer.printf(">%d\n",key);
				for(String read:histogram.get(key))
					writer.printf("%s\n",read);
			}
			writer.close();
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

//        /**********************************************************
//         * Print GC to investigate why it takes so long to finish up
//         */
//        long totalGarbageCollections = 0;
//        long garbageCollectionTime = 0;
//
//        for(GarbageCollectorMXBean gc :
//            ManagementFactory.getGarbageCollectorMXBeans()) {
//
//            long count = gc.getCollectionCount();
//
//            if(count >= 0) {
//                totalGarbageCollections += count;
//            }
//
//            long time = gc.getCollectionTime();
//
//            if(time >= 0) {
//                garbageCollectionTime += time;
//            }
//        }
//
//        System.out.println("Total Garbage Collections: "
//            + totalGarbageCollections);
//        System.out.println("Total Garbage Collection Time (ms): "
//            + garbageCollectionTime);
        
	}

}
