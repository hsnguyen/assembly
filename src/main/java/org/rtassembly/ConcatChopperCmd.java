package org.rtassembly;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.rtassembly.concatemer.NSDFChopper;

import japsa.bio.np.ErrorCorrection;
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
		addString("msa", "kalign", "Name of the MSA for consensus calling.");
				
		addInt("maxK", 100, "Maximum number of monomer copy in a read. Use for cutoff frequency in low-pass filter.");
		addInt("minL", 2000, "Minimum length of monomer or genome size of interests.");
		addInt("consensus", 10, "Minimum concatemers for consensus calling a read.");

		addStdHelp();
	}
	
	public static void main(String[] args) throws IOException{
		CommandLine cmdLine = new ConcatChopperCmd();		
		args = cmdLine.stdParseLine(args);

		/***********************************************************************/
		String 	seqInput = cmdLine.getStringVal("seq"),
				rawInput = cmdLine.getStringVal("raw"),
				outputDir = cmdLine.getStringVal("output"),
				msa = cmdLine.getStringVal("msa");
		int cThres = cmdLine.getIntVal("consensus");

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
						//chop and print consensus sequence if concatemer is long enough (10-concatemer)
						if(key>=cThres) {
							String id=tony.getName().split("\\s+")[0];
							int prev=0, next;
							List<Sequence> seqList=new ArrayList<>();
							for(int i=0;i<chopper.size();i++){
								next=chopper.get(i);
								seqList.add(seq.subSequence(prev,next));
								prev=next;
							}
							try {
								ErrorCorrection.msa=msa;
								Sequence consensus=ErrorCorrection.consensusSequence(seqList, id);
								consensus.setName(id+"_c"+key+"_consensus");
								consensus.writeFasta(outputDir+File.separator+id+"_c"+key+"_consensus.fasta");
	
							} catch (InterruptedException e) {
								// TODO Auto-generated catch block
								e.printStackTrace();
								continue;
							}
						}
						
					}
				}
	
				reader.close();
				
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
