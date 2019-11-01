# A JAVA project for streaming assembly methods using Nanopore data
This Java project aims to develop tools for streaming assembly. In particular, assembly pipelines were migrated here from their code base in [Japsa](https://github.com/mdcao/npScarf) for better developing and maintenance. 

## Quick installation guide
### Requirements
* Linux operating system; not tested for Mac, Window. 
* Java 11+.
* For *npGraph* ,if the pipeline requires an aligner (e.g. raw sequences in FASTA/FASTQ are provided instead of SAM/BAM), [minimap2](https://github.com/lh3/minimap2) (recommended) or [bwa](https://github.com/lh3/bwa) (later than 0.7.11) must be included.

### Install
After cloning the project, the tool is included and can be built with maven2 by following command:
```
 mvn clean package
```
to generate a JAR file containing application modules (target/assembly-x.x.x-SNAPSHOT.jar).

Or you can download directly the JAR file from a release version without having to compile the source.
### Docker
User can build an image from the Dockerfile that also includes [bwa](https://github.com/lh3/bwa), [minimap2](https://github.com/lh3/minimap2) for aligners and [kalign3](https://github.com/TimoLassmann/kalign), [spoa](https://github.com/rvaser/spoa) for MSA-based consensus calling.
```
docker build -t japsa .
```
and then run the container (*npGraph* with GUI, *minimap2* and *kalign v3* by default)
```
docker run --rm -it -e DISPLAY -v $HOME/.Xauthority:/home/developer/.Xauthority -v <local_data_folder>:/data --net=host japsa
```
or overide the default parameters
```
docker run --rm -it -v <local_data_folder>:/data japsa org.rtassembly.NPGraphCmd ...
```
## Modules
For now two modules are on development:
* [*npGraph*](docs/npgraph.md): streaming hybrid assembly using Nanopore data
* [*npConcatemer*](docs/npconcatemer.md): using signal processing to detect concatemeric reads for viral genomes.

## License
Similar to [Japsa](https://github.com/mdcao/japsa) project, tools included in this repo is available under BSD-like license.
