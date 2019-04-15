package org.rtassembly;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import org.rtassembly.concatemer.NSDFChopper;

import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.seq.SequenceReader;
import japsa.util.CommandLine;


public class ConcatChopperCmd extends CommandLine{
	public ConcatChopperCmd(){
		super();

		addString("seq", null, "Name of the base-called FASTQ/FASTA input file.");
		addString("raw", null, "Name of the raw signal FAST5 input file.");
		addString("output", "chopper", "Name of the output folder.");
				
		addInt("maxK", 100, "Maximum number of monomer copy in a read. Use for cutoff frequency in low-pass filter.");
		addInt("minL", 2000, "Minimum length of monomer or genome size of interests.");
		
		addStdHelp();
	}
	
	public static void main(String[] args) throws IOException{
		CommandLine cmdLine = new ConcatChopperCmd();		
		args = cmdLine.stdParseLine(args);

		/***********************************************************************/
		String 	seqInput = cmdLine.getStringVal("seq"),
				rawInput = cmdLine.getStringVal("raw"),
				outputDir = cmdLine.getStringVal("output");

		File outDir=new File(outputDir);
		if(!outDir.isDirectory() && !outDir.mkdirs()){
			System.err.println("Error creating output directory "+ outputDir);
			System.exit(1);
		}
		
		NSDFChopper.CUTOFF_FREQ = cmdLine.getIntVal("maxK");
		NSDFChopper.MIN_MONOMER=cmdLine.getIntVal("minL");

		NSDFChopper tony = new NSDFChopper();
		HashMap<Integer, ArrayList<String>> histogram = new HashMap<>();
		if(seqInput!=null){
			try {
				SequenceReader reader = SequenceReader.getReader(seqInput);
				Sequence seq=null;
				SequenceOutputStream outFile = SequenceOutputStream.makeOutputStream(outputDir+File.separator+"all_monomers.fastq");
				while((seq=reader.nextSequence(Alphabet.DNA4()))!=null){
					tony.setData(seq);
					tony.concatemersDetection();
					ArrayList<Integer> chopper=new ArrayList<>(tony.chopper());
					if(chopper!=null && !chopper.isEmpty()){
						int key=chopper.size();
						ArrayList<String> clist=histogram.get(key);
						if(clist==null){
							clist=new ArrayList<>();
							histogram.put(key, clist);
						}
						clist.add(tony.getName());
						
						
						//chop and print
						int prev=0, next;
						for(int i=0;i<chopper.size();i++){
							next=chopper.get(i);
							Sequence isomer=seq.subSequence(prev,next);
							isomer.setName(seq.getName()+" chopping="+(i+1)+":"+prev+"-"+next);
							isomer.print(outFile);
							
							prev=next;
						}
					}
				}
	
				reader.close();
				outFile.close();
				
				PrintWriter writer = new PrintWriter(new FileWriter(outputDir+File.separator+"count.txt"));
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
		}else if(rawInput!=null){
			System.err.println("To be done...");
		}else{
			System.err.println("Please provide input!\n" + cmdLine.usageString());
		}
	}
	
}
