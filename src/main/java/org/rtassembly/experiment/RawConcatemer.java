package org.rtassembly.experiment;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

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
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		RawConcatemer concat=new RawConcatemer("/home/sonhoanghguyen/Projects/concatemers/data/raw/imb17_013486_20171130__MN17279_sequencing_run_20171130_Ha_BSV_CaMV1_RBarcode_35740_read_17553_ch_455_strand.fast5");
		ArrayList<Double> retval = concat.scan();
		PrintWriter writer = new PrintWriter(new FileWriter("/home/sonhoanghguyen/Projects/concatemers/data/raw/concat7.signal")); 
		writer.print(concat.name+" ");
		for(double sig:retval)				
			writer.printf("%.5f ",sig);
		writer.println();
		
		writer.close();
	}

}
