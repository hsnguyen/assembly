# *npConcatemer* - Rapid method to detect RCA concatemeric signal from viral genomes using Nanopore reads.

## Introduction
Small genomes (around ten kilobase pairs), such as viral or bacterial plasmid, can be deep sequenced by using MinION and Rolling Circle Amplification (RCA). In which, concatemeric molecules generated from RCA can be sequenced directly by the thumbnail sequencer without using restriction enzymes to physically divide them into monomers.

*npConcatemer* is a tool to assist computational analysis on the MinION data of concatemeric long reads in an attempt to produce viral genome assembly in an efficient way.

More details to come...

## Documentation
For detail options of the commandline interface:
```
java -cp target/assembly-x.x.x-SNAPSHOT.jar org.rtassembly.ConcatChopperCmd -h

Usage: 
Options:
  --seq=s         Name of the base-called FASTQ/FASTA input file.
                  (default='null')
  --raw=s         Name of the raw signal FAST5 input file.
                  (default='null')
  --output=s      Name of the output folder.
                  (default='chopper')
  --maxK=i        Maximum number of monomer copy in a read. Use for cutoff frequency in low-pass filter.
                  (default='100')
  --minL=i        Minimum length of monomer or genome size of interests.
                  (default='2000')
  --help          Display this usage and exit
                  (default='false')

```

### Input

### Output

## Reference
Publication on the way.

