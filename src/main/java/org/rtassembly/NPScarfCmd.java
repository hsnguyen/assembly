/*****************************************************************************
 * Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/**************************     REVISION HISTORY    **************************
 * 30/06/2014 - Minh Duc Cao: Created                                        
 *  
 ****************************************************************************/

package org.rtassembly;

import org.rtassembly.scaffold.ContigBridge;
import org.rtassembly.scaffold.RealtimeScaffolding;
import org.rtassembly.scaffold.ScaffoldGraph;
import org.rtassembly.scaffold.ScaffoldGraphDFS;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import japsa.seq.SequenceReader;
import japsa.util.CommandLine;
import japsa.util.Logging;

import java.io.BufferedReader;
import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author sonnguyen, minhduc
 * 
 */
public class NPScarfCmd extends CommandLine{
	private static final Logger LOG = LoggerFactory.getLogger(NPScarfCmd.class);

	public NPScarfCmd(){
		super();
		addString("seqFile", null, "Name of the assembly file (sorted by length)",true);

		addString("input", "-", "Name of the input file, - for stdin", true);
		addString("format", "sam", "Format of the input: fastq/fasta or sam/bam", true);
		addBoolean("index", true, "Whether to index the contigs sequence by the aligner or not.");
		addString("bwaExe", "bwa", "Path to bwa");
		addInt("bwaThread", 4, "Theads used by bwa");
		addBoolean("long", false, "Whether report all sequences, including short/repeat contigs (default) or only long/unique/completed sequences.");
		addString("assembler", "spades", "Name of the assembler used for Illumina assembly: SPAdes (default) or ABySS.");
		addString("spadesDir", null, "Name of the output folder by SPAdes: assembly graph and paths will be used for better gap-filling.");
		addString("prefix", "out", "Prefix for the output files");	
		addString("genes", null , "Realtime annotation: name of annotated genes in GFF 3.0 format");
		addString("resistGene", null , "Realtime annotation: name of antibiotic resistance gene fasta file");
		addString("insertSeq", null , "Realtime annotation: name of IS fasta file");
		addString("oriRep", null, "Realtime annotation: name of fasta file containing possible origin of replication");
		//addInt("marginThres", 1000, "Margin threshold: to limit distance to the contig's ends of the alignment used in bridging."); 
		addInt("minContig", 300, "Minimum contigs length that are used in scaffolding."); 
		addInt("maxRepeat", 7500, "Maximum length of repeat in considering species."); 

		addDouble("cov", 0, "Expected average coverage of Illumina, <=0 to estimate");
		addInt("qual", 1, "Minimum quality");
		addInt("support", 1, "Minimum supporting long read needed for a link between markers");

		addBoolean("realtime", false, "Process in real-time mode. Default is batch mode (false)");
		addInt("read", 50,  "Minimum number of reads between analyses");		
		addInt("time", 10,   "Minimum number of seconds between analyses");
		addBoolean("verbose", false, "Turn on debugging mode");

		addStdHelp();		

	} 	
	//static boolean hardClip = false;

	public static void main(String[] args) throws 
	IOException, InterruptedException {
		CommandLine cmdLine = new NPScarfCmd();		
		args = cmdLine.stdParseLine(args);

		/***********************************************************************/
		String prefix = cmdLine.getStringVal("prefix");
		//String bamFile = cmdLine.getStringVal("bamFile");

		String input = cmdLine.getStringVal("input");
		String bwaExe = cmdLine.getStringVal("bwaExe");
		int bwaThread = cmdLine.getIntVal("bwaThread");
		String format = cmdLine.getStringVal("format").toLowerCase();

		String sequenceFile = cmdLine.getStringVal("seqFile"),
				spadesFolder = cmdLine.getStringVal("spadesDir"),

				genesFile = cmdLine.getStringVal("genes"),
				resistFile = cmdLine.getStringVal("resistGene"),
				isFile = cmdLine.getStringVal("insertSeq"),
				oriFile = cmdLine.getStringVal("oriRep");


		File 	graphFile = new File(spadesFolder+"/assembly_graph.fastg"),
				pathFile = new File(spadesFolder+"/contigs.paths");
		String assembler = cmdLine.getStringVal("assembler").toLowerCase();

		if(assembler.equals("abyss")){
			ScaffoldGraph.assembler=0b01;
			//graphFile = ;
			//pathFile = ;
		}

		if (format.startsWith("fastq") ||
				format.startsWith("fasta") ||
				format.startsWith("fq") ||
				format.startsWith("fa")){
			try{
				ProcessBuilder pb = new ProcessBuilder(bwaExe).redirectErrorStream(true);
				Process process =  pb.start();
				//BWA process doesn't produce gzip-compressed output
				BufferedReader bf = SequenceReader.openInputStream(process.getInputStream());


				String line;
				String version = "";
				Pattern versionPattern = Pattern.compile("^Version:\\s(\\d+\\.\\d+\\.\\d+).*");
				Matcher matcher=versionPattern.matcher("");
				
				while ((line = bf.readLine())!=null){				
					matcher.reset(line);
					if (matcher.find()){
					    version = matcher.group(1);
					    break;//while
					}
					
									
				}	
				bf.close();
				
				if (version.length() == 0){
					LOG.error(bwaExe + " is not the right path to bwa. bwa is required");
					System.exit(1);
				}else{
					LOG.info("bwa version: " + version);
					if (version.compareTo("0.7.11") < 0){
						LOG.error(" Require bwa of 0.7.11 or above");
						System.exit(1);
					}
				}
				
				//run indexing 
				if(cmdLine.getBooleanVal("index")){
					LOG.info("bwa index running...");
					ProcessBuilder pb2 = new ProcessBuilder(bwaExe,"index",sequenceFile);
					Process indexProcess =  pb2.start();
					indexProcess.waitFor();
					LOG.info("bwa index finished!");
				}
			}catch (IOException e){
				System.err.println(e.getMessage());
				System.exit(1);
			}

		}else if (format.startsWith("sam") || format.startsWith("bam")){
			// no problem
		}else{
			LOG.error("Unrecognized format: " + format);
			System.exit(1);
		}

		if(spadesFolder !=null){
			if(ScaffoldGraph.assembler==0b00  && graphFile.exists() && pathFile.exists()){
				LOG.info("===> Use assembly graph and path from SPAdes!");				
			}else if (ScaffoldGraph.assembler==0b01){
				File f = new File(spadesFolder);
				File[] matchingFiles = f.listFiles(new FilenameFilter() {
				    public boolean accept(File dir, String name) {
				        return name.endsWith("contigs.dot");
				    }
				});
				if(matchingFiles.length != 1){
					LOG.info("Failed to looking for an unique *-contigs.dot file in " + spadesFolder + " . Proceeding without assembly graph...");
					spadesFolder=null;
				} else{
					graphFile=matchingFiles[0];
					LOG.info("===> Use assembly graph from short-read assembler ABySS: " + graphFile);
				}
				
			}
			
		}
		else{
			LOG.info("Not found any legal assembly output folder, assembly graph thus not included!");
			spadesFolder=null;
		}
		
		

		int 	//marginThres = cmdLine.getIntVal("marginThres"),
		minContig = cmdLine.getIntVal("minContig"),
		minSupport = cmdLine.getIntVal("support"),
		maxRepeat = cmdLine.getIntVal("maxRepeat");
		//if(marginThres < 0)
		//	LOG.exit("Marginal threshold must not be negative", 1);
		if(minContig <= 0) {
			LOG.error("Minimum contig length has to be positive");
			System.exit(1);

		}if(minSupport <= 0) {
			LOG.error("Minimum supporting reads has to be positive");
			System.exit(1);
		}
		if(maxRepeat <= 0) {
			LOG.error("Maximal possible repeat length has to be positive", 1);
		}


		ScaffoldGraph.minContigLength = minContig;
		ScaffoldGraph.minSupportReads = minSupport;	
		ScaffoldGraph.maxRepeatLength = maxRepeat;		
		//ScaffoldGraph.marginThres = marginThres;
		ScaffoldGraph.verbose = cmdLine.getBooleanVal("verbose");
		ScaffoldGraph.reportAll = !cmdLine.getBooleanVal("long");

		double cov = cmdLine.getDoubleVal("cov");
		int qual = cmdLine.getIntVal("qual");
		if(qual < 0) {
			LOG.error("Phred score of quality has to be positive");
			System.exit(1);
		}

		int number = cmdLine.getIntVal("read"),
				time = cmdLine.getIntVal("time");

		if(number <= 0) {
			LOG.error("Number of reads has to be positive");
			System.exit(1);
		}
		if(time < 0) {
			LOG.error("Sleeping time must not be negative");
			System.exit(1);
		}
		/**********************************************************************/

		ScaffoldGraph graph;
		boolean rt = cmdLine.getBooleanVal("realtime");
		ContigBridge.relaxFilling();
		if(rt){
			RealtimeScaffolding rtScaffolding = new RealtimeScaffolding(sequenceFile, genesFile, resistFile, isFile, oriFile, "-");

			graph = rtScaffolding.graph;
			if(prefix != null)
				graph.prefix = prefix;
			if(spadesFolder!=null)
				synchronized(graph){
					graph.readMore(spadesFolder+"/assembly_graph.fastg",spadesFolder+"/contigs.paths");
				}
			if (cov <=0)
				cov = ScaffoldGraph.estimatedCov;

			rtScaffolding.scaffolding(input, number, time, cov/1.6, qual, format, bwaExe, bwaThread, sequenceFile);

		}
		else{
			graph = new ScaffoldGraphDFS(sequenceFile, genesFile, resistFile, isFile, oriFile);
			if(spadesFolder!=null)
				graph.readMore(spadesFolder+"/assembly_graph.fastg",spadesFolder+"/contigs.paths");

			if (cov <=0)
				cov = ScaffoldGraph.estimatedCov;

			graph.makeConnections(input, cov / 1.6, qual, format, bwaExe, bwaThread, sequenceFile);

			graph.connectBridges();
			if(prefix != null)
				graph.prefix = prefix;

			ContigBridge.forceFilling();
			graph.printSequences(true,true);
		}

	}
}
