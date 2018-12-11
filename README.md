# *npGraph* - Resolve assembly graph in real-time using nanopore data
This is another real-time scaffolder beside [npScarf](https://github.com/mdcao/npScarf). Instead of using contig sequences as pre-assemblies, this tool is able to work on assembly graph (from [SPAdes](http://cab.spbu.ru/software/spades/)). 
The batch algorithm has been implemented in hybrid assembler module of [Unicycler](https://github.com/rrwick/Unicycler) and others.

<p align="center">
  <img src="http://drive.google.com/uc?export=view&id=10sZZk9QzjOiI_P3qxX47eHkKkaIUN_qB" alt="npGraph"/>
</p>

## Introduction
*npScarf* is the real-time hybrid assembler that use the stream of long reads to bridge the Illumina contigs together, expecting to give more complete genome sequences while the sequencing process is still ongoing. The pipeline has been applied sucessfully for microbial genomics and even bigger data sets. However, due to its greedy approach over the noisy data, it is difficult to eliminate all mis-assemblies without further pre-processing and parameter tuning. To help prevent this issue, the assembly graph - bulding block graph structure for the contigs - should be used as the source for bridging algorithm. 
This approach can give better accuracy, but as the trade-off, are more computational expensive and challenging to adapt in the real-time mode.

A (rather simple at the moment) Graphical User Interface is implemented for better interaction. [GraphStream](http://graphstream-project.org/), a dynamic graph library for Java, has been employed for such task. At the moment, only a few among sea of amazing features of this library had been used, leaving many room for visualiazation improvements.
## Quick installation guide
The tool is included in a Java project that can be built with maven2 by following command:
```
 mvn clean package
```
to generate a jar file containing **npGraph** and several other modules.

The code has been developed with *Java 1.8.0_144* that enables lambda expression and JavaFX, so equal or later version is expected to build & run the tool properly.

## Documentation
For detail options of the commandline interface:
```
java -cp target/assembly-0.0.1-SNAPSHOT.jar org.rtassembly.NPGraphCmd -h

Usage: 
Options:
  --si=s          Name of the short-read assembly file.
                  (default='')
  --sf=s          Format of the assembly input file. Accepted format are FASTG, GFA
                  (default='')
  --li=s          Name of the long-read data input file, - for stdin.
                  (default='')
  --lf=s          Format of the long-read data input file. This may be FASTQ/FASTA (MinION reads) or SAM/BAM (aligned with the assembly graph already)
                  (default='')
  --output=s      Name of the output folder.
                  (default='/tmp/')
  --sb=s          Name of the metaBAT file for binning information (experimental).
                  (default='')
  --aligner=s     Aligner tool that will be used, either minimap2 or bwa
                  (default='')
  --algPath=s     Absolute path to the binary aligner file
                  (default='')
  --algOpt=s      Settings used by aligner to align long reads to the contigs
                  (default='')
  --overwrite     Whether to overwrite or reuse the intermediate file
                  (default='true')
  --spaths        Whether to use SPAdes contigs.paths for bridging.
                  (default='false')
  --qual=i        Minimum quality of alignment to considered
                  (default='1')
  --dfs=i         Number of DFS steps to search
                  (default='15')
  --gui           Whether using GUI or not.
                  (default='false')
  --help          Display this usage and exit
                  (default='false')


```

or invoke the GUI of the module and control it interactively by
```
java -cp target/assembly-0.0.1-SNAPSHOT.jar org.rtassembly.NPGraphCmd -gui
```
A proper combination of command line and GUI can provide an useful streaming pipeline that copes well with MinION output data. This practice allows the assembly to take place abreast of nanopore sequencing run.

### Input
The mandatory inputs are short-read pre-assemblies file (*-si*) and long-read data (*-li*).
The pre-assemblies must be given as the assembly graph outputted from SPAdes in either FASTG or GFA file (normally assembly_graph.fastg or assembly_graph.gfa).

The long-read data will be used for bridging and can be given as DNA sequences (FASTA/FASTQ format, possible .gz) or alignment records (SAM/BAM). If the sequences are given, then it's mandatory to have either BWA-MEM or minimap2 installed in your system to do the alignment between long reads and the pre-assemblies. Alternative option is to use your favourite aligner and provide SAM/BAM to *npGraph*. *npGraph* will try to guess the format of the inputs based on the extensions, but sometimes you'll have to specify it yourself (e.g. when "-" is provided to read from *stdin*).

It is important to emphasis the quality of the assembly graph to the final results. [Unicycler](https://github.com/rrwick/Unicycler) pre-process the graph data by running SPAdes with multiple *kmer* options to chose the best one. Unfortunatly, *npGraph* doesn't employ such technique thus if the graph is not good, you should do the task for yourself before running the tool. Normally, 60X Illumina MiSeq data would give decent SPAdes assembly graph. The better the assembly graph is, the more complete and accurate assembly you'll get.
It doesn't do neither any polishing or other exhaustive post-processing for the final assembly assuming the quality is equivalent to the short-read data which is good enough.

### Output
The tool generate the assembly in FASTA sequences. I'll output GFA file if needed.
Also, the stats of the assembly process is reported as the program is running.

### GUI
The GUI includes the dashboard for control the settings of the program and a separate window for the assembly graph.
From the second window, the colored vertices are unique contigs while the white ones are either unknown or repeats. The number of colors (other than white) indicates number of populations (e.g. chromosome vs plasmids, or different bins in metagenomics).

More features would be added later.
### Note
* GUI consumes memory, considering increase heap size (e.g. -Xmx16G). If the number of vertices in the assembly graph greater than 1000, you shouldn't show the graph in real-time if running on the normal desktop.
* If aligner is used along side, there is more resource occupied. Considering separate alignment and npGraph+GUI on different machines communicating through network socket e.g. by Japsa utility [jsa.util.streamServer](https://japsa.readthedocs.io/en/latest/tools/jsa.util.streamServer.html) and [jsa.util.streamClient](https://japsa.readthedocs.io/en/latest/tools/jsa.util.streamClient.html)
* Most suitable for bacterial data (assembly graph not too complicated). Consider using *npScarf* for bigger data set.
## Reference
Publication on the way... of procrastination.

## License
Similar to [Japsa](https://github.com/mdcao/japsa) project, tools included in this repo is available under BSD-like license.
