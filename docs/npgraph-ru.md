# *npGraphServer* - gRPC server to run selective sequencing aiming at complete assembly

## Introduction
ReadUntil protocol for targeted sequencing requires cricial real-time operations on the streaming output of Nanopore devices.
Here, we implemented *npGraphServer* that provides gRPC interface to invoke from ReadUntil (now ReadFish) python code.

## Documentation
For detail options of the commandline interface:
```
java -cp target/assembly-x.x.x.jar org.rtassembly.NPGraphServerCmd -h

Usage: 
Options:
  --port=i        Port number for which the server will listen to.
                  (default='2105')
  --si=s          Name of the short-read assembly file.
                  (default='')
  --sf=s          Format of the assembly input file. Accepted format are FASTG, GFA
                  (default='')
  --output=s      Output folder for temporary files and the final assembly npgraph_assembly.fasta
                  (default='/tmp/')
  --sb=s          Name of the metaBAT file for binning information (experimental).
                  (default='')
  --overwrite     Whether to overwrite or reuse the intermediate file
                  (default='true')
  --sp            Whether to use SPAdes contigs.paths for bridging.
                  (default='false')
  --qual=i        Minimum quality of alignment to considered
                  (default='10')
  --mincov=i      Minimum number of reads spanning a confident bridge
                  (default='3')
  --maxcov=i      Cut-off number of reads spanning a confident bridge
                  (default='20')
  --depth=i       Limit depth for searching path between 2 neighbors
                  (default='300')
  --anchor=i      Lowerbound of an anchor contig's length.
                  (default='1000')
  --unique=i      Lowerbound of a unique contig's length.
                  (default='10000')
  --expect=i      The readuntil protocol will expect a read to span this length for an unblock/proceed decision.
                  (default='10000')
  --time=i        Time interval to considered for real-time assembly.
                  (default='10')
  --read=i        Read interval to considered for real-time assembly.
                  (default='50')
  --gui           Whether using GUI or not.
                  (default='false')
  --keep          Whether to keep extremely-low-coveraged contigs.
                  (default='false')
  --help          Display this usage and exit
                  (default='false')


```
The provided functionality is similar to *npGraph* except no consensus is supported since MSA is slow for this time-critical server. Also PAF is the only alignment format used for rRPC communications since its light-weighted property.

### How to
To use the pipeline, first you need to familiarize yourself with [*ReadFish*](https://github.com/LooseLab/readfish). A slight modification of this repo to operate with *npGraphServer* is [*here*](https://github.com/hsnguyen/ru). So after install original ReadUntil code, you'll have to update (from the virtualenv) by
```
pip install -U git+https://github.com/hsnguyen/ru@master
```


After finishing all installation, you can start to test running selective sequencing with *npGraphServer*.Beside the Guppy basecall server, you'll need to start a *npGraphServer* as well for the run.
```
java -cp target/assembly-x.x.x.jar org.rtassembly.NPGraphServerCmd -si assembly_graph.[fastg|gfa] -output output/ -read 1 -expect 20000 -gui
```
Follow the procedure similar to *ReadFish* but using the [*npgraph_selection.toml*](https://github.com/hsnguyen/ru/blob/master/examples/npgraph_selection.toml) instead.

## Reference



