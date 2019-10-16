#!/bin/bash
###Must have spades (megahit), metabat(maxbin), bwa, minimap2, samtools in the $PATH
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <Illumina_PE1.fastq.gz> <Illumina_PE2.fastq.gz>"
	exit 1	
fi

SDIR=`dirname $0`
FDIR=`dirname $1`/`date +%s` #an unique folder in the same level as input data to output intermediate results
mkdir -p ${FDIR}

# Path to npGraph jar file
# npgraph="${SDIR}/../target/assembly-0.1.0-SNAPSHOT.jar"

PE1=$1
PE2=$2

ITE=3 #number of stages
## After each iteration, $i/: 
# assembly_graph_with_scaffolds.gfa is the original graph
# npgraph_assembly.gfa is the simplified graph
# npgraph_assembly.bin is binning info for the simplified graph 
for i in `seq ${ITE}`
do
	
	## Short read assembly using metaSPAdes
	spades.py --meta -k 55,77,99,107,127 -1 ${PE1} -2 ${PE2} -o ${FDIR}/${i}
	## Run npgraph without long reads to have the significant contigs, reassembly reads that are unaligned to this significant contigs
	awk '/^S/{print ">"$2; print $3;}' ${FDIR}/${i}/assembly_graph_with_scaffolds.gfa | fold > ${FDIR}/${i}/assembly_graph_with_scaffolds.fasta
	
	bwa index ${FDIR}/${i}/assembly_graph_with_scaffolds.fasta
	bwa mem -t8 ${FDIR}/${i}/assembly_graph_with_scaffolds.fasta ${PE1} ${PE2} | samtools view -hF 2304 - | samtools sort -@8 -o ${s}/metaSPAdes/bin.bam - \
	&& jgi_summarize_bam_contig_depths --outputDepth ${FDIR}/${i}/orig_cov.txt ${FDIR}/${i}/bin.bam \
	&& rm ${FDIR}/${i}/bin.bam
	
	metabat  --saveCls --noBinOut --inFile ${FDIR}/${i}/assembly_graph_with_scaffolds.fasta --abdFile ${FDIR}/${i}/orig_cov.txt --outFile ${FDIR}/${i}/assembly_graph_with_scaffolds.bin 

	java -cp ${npgraph} org.rtassembly.NPGraphCmd -si ${FDIR}/${i}/assembly_graph_with_scaffolds.gfa -sb ${FDIR}/${i}/assembly_graph_with_scaffolds.bin -output ${FDIR}/${i}/ > ${FDIR}/${i}/npgraph.log 2>&1
	
	bwa index ${FDIR}/${i}/npgraph_assembly.fasta
	bwa mem -t8 ${FDIR}/${i}/npgraph_assembly.fasta ${PE1} ${PE2} | samtools view -bh - > ${FDIR}/${i}/alignments.bam
	
	## extract unaligned read pairs
	samtools view -u  -f 4 -F 264 ${FDIR}/${i}/alignments.bam  > ${FDIR}/${i}/tmps1.bam
	samtools view -u -f 8 -F 260 ${FDIR}/${i}/alignments.bam  > ${FDIR}/${i}/tmps2.bam
	samtools view -u -f 12 -F 256 ${FDIR}/${i}/alignments.bam > ${FDIR}/${i}/tmps3.bam
	samtools merge -u - ${FDIR}/${i}/tmps[123].bam | samtools sort -@8 -n -o ${FDIR}/${i}/unmapped.bam -
	samtools fastq -n -1 ${FDIR}/${i}/unmapped_reads_1.fastq -2 ${FDIR}/${i}/unmapped_reads_2.fastq ${FDIR}/${i}/unmapped.bam
	
	## binning the simplified graph
	samtools view -hF 2304 ${FDIR}/${i}/alignments.bam | samtools sort -@8 -o ${FDIR}/${i}/alignments_sorted.bam \
	&& jgi_summarize_bam_contig_depths --outputDepth ${FDIR}/${i}/sim_cov.txt ${FDIR}/${i}/alignments_sorted.bam && rm ${FDIR}/${i}/alignments_sorted.bam
	metabat  --saveCls --noBinOut --inFile ${FDIR}/${i}/ngraph_assembly.fasta --abdFile  ${FDIR}/${i}/sim_cov.txt --outFile ${FDIR}/${i}/ngraph_assembly.bin	

	PE1=${FDIR}/${i}/unmapped_reads_1.fastq
	PE2=${FDIR}/${i}/unmapped_reads_2.fastq

done

## Now combined all the graphs with their corresponding binning files

