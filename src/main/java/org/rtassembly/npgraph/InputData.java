package org.rtassembly.npgraph;

import java.lang.invoke.MethodHandles;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import javafx.beans.property.BooleanProperty;
import javafx.beans.property.SimpleBooleanProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.beans.property.StringProperty;

public class InputData {
    private static final Logger logger = LogManager.getLogger(MethodHandles.lookup().lookupClass());

	public ShortReadInputs srIn;
	public LongReadInputs lrIn;
	
	public InputData() {
		srIn=new ShortReadInputs();
		lrIn=new LongReadInputs();
	}
	public boolean validateShortReadInput() {return srIn.validate();}
	public boolean validateLongReadInput() {return lrIn.validate();}
	
	public final void setUseSPAdesPath(boolean owr) {srIn.useSPAdesPath.set(owr);}
	public final boolean getUseSPAdesPath() {return srIn.useSPAdesPath.get();}
	public BooleanProperty useSPAdesPathProperty() {return srIn.useSPAdesPath;}
		
	public final void setAligner(String tool) {	lrIn.aligner.set(tool);}
	public final String getAligner() {return lrIn.aligner.get();}
	public StringProperty alignerProperty(){return lrIn.aligner;}
	
	public final void setAlignerOpts(String setting) {lrIn.alignerOpt.set(setting);}
	public final String getAlignerOpts() {return lrIn.alignerOpt.get();}
	public StringProperty alignerOptProperty(){return lrIn.alignerOpt;}
	
	public final void setMSA(String tool) {	lrIn.msa.set(tool);}
	public final String getMSA() {return lrIn.msa.get();}
	public StringProperty msaProperty(){return lrIn.msa;}
	
	public final void setBinReadsInput(String brInput) {srIn.binReadsInput.set(brInput);}
	public final String getBinReadsInput() {return srIn.binReadsInput.get();}
	public StringProperty binReadsInputProperty() {return srIn.binReadsInput;}
	
	public final void setShortReadsInput(String srInput) {srIn.shortReadsInput.set(srInput);}
	public final String getShortReadsInput() {return srIn.shortReadsInput.get();}
	public StringProperty shortReadsInputProperty() {return srIn.shortReadsInput;}
	
	public final void setLongReadsInput(String lrInput) {lrIn.longReadsInput.set(lrInput);
	}
	public final String getLongReadsInput() {return lrIn.longReadsInput.get();}
	public StringProperty longReadsInputProperty() {return lrIn.longReadsInput;}
	
	public final void setShortReadsInputFormat(String srInputFormat) {srIn.shortReadsInputFormat.set(srInputFormat);}
	public final String getShortReadsInputFormat() {return srIn.shortReadsInputFormat.get();}
	public StringProperty shortReadsInputFormatProperty() {return srIn.shortReadsInputFormat;}
	
	public final void setLongReadsInputFormat(String lrInputFormat) {
		lrInputFormat=lrInputFormat.toLowerCase();
		if(	lrInputFormat.contains("fasta") || lrInputFormat.contains("fa") || lrInputFormat.contains("fna")
			|| lrInputFormat.contains("fastq") || lrInputFormat.contains("fq"))
			lrIn.longReadsInputFormat.set("fasta/fastq");
		else if(lrInputFormat.contains("sam") || lrInputFormat.contains("bam"))
			lrIn.longReadsInputFormat.set("sam/bam");
		else if(lrInputFormat.contains("paf"))
			lrIn.longReadsInputFormat.set("paf");
	}
	public final String getLongReadsInputFormat() {return lrIn.longReadsInputFormat.get();}
	public StringProperty longReadsInputFormatProperty() {return lrIn.longReadsInputFormat;}

	
	/*
	 * Class represent input data regarding Illumina short-read	
	 */
	class ShortReadInputs{
	    private BooleanProperty useSPAdesPath;
	    private StringProperty 	shortReadsInput, 
								binReadsInput, 
								shortReadsInputFormat;
	    ShortReadInputs(){
	    	useSPAdesPath=new SimpleBooleanProperty(false);
	    	shortReadsInput=new SimpleStringProperty("");
	    	binReadsInput=new SimpleStringProperty("");
	    	shortReadsInputFormat=new SimpleStringProperty("");
	    	
	        shortReadsInput.addListener((observable, oldValue, newValue) -> 
				{
					String fn = ((String)observable.getValue()).toLowerCase();
					if(	fn.endsWith(".fastg")) 
						setShortReadsInputFormat("fastg");
					else if(fn.endsWith(".gfa"))
						setShortReadsInputFormat("gfa");
				}	 

	        );
	        
	        shortReadsInputFormat.addListener((observable, oldValue, newValue) -> 
				{
					if(!getShortReadsInput().toLowerCase().endsWith(newValue))
						setShortReadsInput("");
				}	 

	        );
	    }
	    
	    boolean validate() {
			if(!getShortReadsInputFormat().equals("fastg") && !getShortReadsInputFormat().equals("gfa")){
				logger.error("Please specify a correct format of graph file!");
				return false;
			}
				
			if(!GraphUtil.checkFile(getShortReadsInput()))
				return false;
			
			return true;
	    }
	}
	/*
	 * Class represent input data involving ONT long-read	
	 */
	class LongReadInputs{
	    private StringProperty 	aligner,
	    						alignerOpt,
	    						msa,
								longReadsInput,
								longReadsInputFormat;
	    LongReadInputs(){
	    	aligner=new SimpleStringProperty("");
	    	alignerOpt=new SimpleStringProperty("");
	    	msa=new SimpleStringProperty("");
	    	longReadsInput=new SimpleStringProperty("");
	    	longReadsInputFormat=new SimpleStringProperty("");
	    	
	        longReadsInput.addListener( (observable, oldValue, newValue) -> 
	    		{
					String fn = ((String)observable.getValue()).toLowerCase();
					if(	fn.endsWith(".fasta") || fn.endsWith(".fa") || fn.endsWith("fna")
						|| fn.endsWith(".fastq") || fn.endsWith(".fq")
						|| fn.endsWith(".fasta.gz") || fn.endsWith(".fa.gz") || fn.endsWith("fna.gz")
						|| fn.endsWith(".fastq.gz") || fn.endsWith(".fq.gz") 
						) 
						setLongReadsInputFormat("fasta/fastq");
					else if(fn.endsWith(".sam") || fn.endsWith(".bam")) 
						setLongReadsInputFormat("sam/bam");
					else if(fn.endsWith(".paf")) 
						setLongReadsInputFormat("paf");	
	    		}	 
	        );
	        
	        longReadsInputFormat.addListener((observable, oldValue, newValue) -> 
				{
					String oldFile=getLongReadsInput().toLowerCase();
					if(oldFile.equals("-"))
						return;
					if(	newValue.equals("fasta/fastq") 
								&& !oldFile.endsWith(".fasta") && !oldFile.endsWith(".fa") && !oldFile.endsWith("fna")
								&& !oldFile.endsWith(".fastq") && !oldFile.endsWith(".fq")
								&& !oldFile.endsWith(".fasta.gz") && !oldFile.endsWith(".fa.gz") && !oldFile.endsWith("fna.gz")
								&& !oldFile.endsWith(".fastq.gz") && !oldFile.endsWith(".fq.gz") 
								) 
						setLongReadsInput("");
							
					if(newValue.equals("sam/bam") && !oldFile.endsWith(".sam") && !oldFile.endsWith(".bam"))
						setLongReadsInput("");
					
					if(newValue.equals("paf") && !oldFile.endsWith(".paf"))
						setLongReadsInput("");
				}	 

	        );
	        
	        aligner.addListener( (observable, oldValue, newValue) ->
	        	{
					String aligner=(String)observable.getValue();
					if(aligner.toLowerCase().endsWith("minimap2"))
						setAlignerOpts("-t4 -k15 -w5");
					else if (aligner.toLowerCase().endsWith("bwa"))
						setAlignerOpts("-t4 -k11 -W20 -r10 -A1 -B1 -O1 -E1 -L0 -a -Y");			
				}	 

	        );
	    }
	    
	    boolean validate() {	    	
			//accept the case when no long read data is provided. Just output simplified assembly graph then.
			if(getLongReadsInput().isEmpty()) {
				logger.warn("No long read data is provided. Only output the simplified assembly graph and stop!");
				return true;
			}
			
			//long-read file format must be known if provided
			if (getLongReadsInput().isEmpty()) {
				logger.warn("Please specify the long-read input format!");
				return false;
			}else if(getLongReadsInputFormat().equals("fasta/fastq")){
				if(getAligner().isEmpty() || getAlignerOpts().isEmpty()) {
					logger.warn("Please specify aligner and/or alignment arguments!");
					return false;
				}
			}else if(!getLongReadsInputFormat().equals("sam/bam")
					&& !getLongReadsInputFormat().equals("paf") ){
				logger.error("Please specify a correct alignment format (BAM/SAM/PAF)!");
				return false;
			}
			
			//check file valid if not input from stdin
			if(!getLongReadsInput().equals("-") && !GraphUtil.checkFile(getLongReadsInput()))
				return false;
						
			
	    	return true;
	    }
	}
}