# A JAVA project for streaming assembly methods using Nanopore data
This Java project aims to develop tools for streaming assembly. In particular, assembly pipelines were migrated here from their code base in [Japsa](https://github.com/mdcao/npScarf) for better developing and maintenance. 

## Quick installation guide
### Requirements
* Linux operating system; not tested for Mac, Window. 
* The project has been developed with *Oracle Java 1.8.0_144* that enables lambda expression and JavaFX included, so similar or later JVM version is expected to build & run the tool properly. 

If JavaFX is not found from your compiler (e.g. OpenJDK or Oracle Java later version that might not include it in the future), you need to download jfxrt.jar and specify it from POM. So it's better to have Java 8 with JavaFX bundled.
* For *npGraph* ,if the pipeline requires an aligner (e.g. raw sequences in FASTA/FASTQ are provided instead of SAM/BAM), [minimap2](https://github.com/lh3/minimap2) (recommended) or [bwa](https://github.com/lh3/bwa) (later than 0.7.11) must be included.

### Install
After cloning the project, the tool is included and can be built with maven2 by following command:
```
 mvn clean package
```
to generate a JAR file containing application modules (target/assembly-x.x.x-SNAPSHOT.jar).

Or you can download directly the JAR file from a release version without having to compile the source.

## Modules
For now two modules are on development:
* [*npGraph*](docs/npgraph.md): streaming hybrid assembly using Nanopore data
* [*npSignal*](docs/npsignal.md): using signal processing to detect concatemeric reads for viral genomes.

More to come in the future. Please support.
## License
Similar to [Japsa](https://github.com/mdcao/japsa) project, tools included in this repo is available under BSD-like license.
