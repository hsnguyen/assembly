# *npGraph* - Resolve assembly graph in real-time using long read data
This is another real-time scaffolder beside [npScarf](https://github.com/mdcao/npScarf). Instead of using contig sequences as pre-assemblies, this tool is able to work on assembly graph (from [SPAdes](http://cab.spbu.ru/software/spades/)). 
The batch algorithm has been implemented in hybrid assembler module of [Unicycler](https://github.com/rrwick/Unicycler) and others.

### Introduction


A (rather simple at the moment) Graphical User Interface is implemented for better interaction. [GraphStream](http://graphstream-project.org/) has been employed for such task.
### Quick installation guide
The tool is included in a Java project that can be built with maven2 by following command:
```
 mvn clean package
```
to generate a jar file containing **npGraph** and several other modules.

The code has been developed with platform running *Java 1.8.0_144* that enables lambda expression and JavaFX, so equal or later version is expected to build & run the tool properly.

### Documentation
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

A proper combination of command line and GUI can provide an useful streaming fashion that copes well with MinION output data.
This practice allows the assembly to take place abreast of nanopore sequencing run.

### Reference
Publication on the way (or not!)

### License
Similar to [Japsa]{https://github.com/mdcao/japsa} project, tools included in this repo is available under BSD-like license.
