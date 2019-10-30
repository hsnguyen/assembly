# A Java project for streaming assembly methods using Nanopore data
This Java project aims to develop tools for streaming assembly, as an extension of [Japsa](https://github.com/mdcao/npScarf). In particular, assembly pipeline had been migrated here from its code base in [Japsa](https://github.com/mdcao/npScarf) for more convenient developing and maintaining the modules. In the future, the project might include additional analysis pipelines for Nanopore data as well.

## Quick installation guide
### Requirements
* Linux operating system; not tested for Mac, Window. 
* Java 11+.
* For [*npGraph*](docs/npgraph.md) ,if the pipeline requires an aligner (e.g. raw sequences in FASTA/FASTQ are provided instead of SAM/BAM), [minimap2](https://github.com/lh3/minimap2) (recommended) or [bwa](https://github.com/lh3/bwa) (later than 0.7.11) must be included.

### Install by maven
After cloning the project, the tool can be built with maven2 (3.6.0) by following command:
```
git clone https://github.com/hsnguyen/assembly.git

mvn clean package
```
to generate a JAR file containing application modules (target/assembly-x.x.x-SNAPSHOT.jar).

Or you can download directly the JAR file from a release version without having to compile the source.
### Docker
User can build an image directly from the Dockerfile that also includes [bwa](https://github.com/lh3/bwa), [minimap2](https://github.com/lh3/minimap2) for aligners and [kalign3](https://github.com/TimoLassmann/kalign), [spoa](https://github.com/rvaser/spoa) for MSA-based consensus calling.
```
docker build -t npgraph .
```
The image is also made available on DockerHub as well
```
docker pull nguyenhoangson/npgraph
```

After having the docker image, one can run [*npGraph*](docs/npgraph.md) with GUI, *minimap2* and *kalign v3* by default
```
docker run --rm -it -e DISPLAY -v $HOME/.Xauthority:/home/developer/.Xauthority -v <local_data_folder>:/data --net=host npgraph
```
You might need to run ```xhost +``` to disable X server access control on your local machine before running [*npGraph*](docs/npgraph.md) with GUI by above command.

If you want to disable GUI mode, overide the default executable behaviour, or invoke another module than [*npGraph*](docs/npgraph.md) then provide appropriate parameters, e.g.
```
docker run --rm -it -v <local_data_folder>:/data npgraph org.rtassembly.NPGraphCmd ...
```


## Modules
For now two modules are on development:
* [*npGraph*](docs/npgraph.md): streaming hybrid assembly using Nanopore data
* [*npConcatemer*](docs/npsignal.md): using signal processing to detect concatemeric reads for viral genomes (experimental)

## License
Similar to [Japsa](https://github.com/mdcao/japsa) project, tools included in this repo is available under BSD-like license.
