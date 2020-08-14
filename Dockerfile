# Dockerfile
#
# BUILD STAGE: BWA, MINIMAP2, KALIGN, NPGRAPH 
#
FROM maven:3.6.0-jdk-11-slim AS build
# set non-interactive mode
ENV DEBIAN_FRONTEND noninteractive

# Install dependencies
RUN	apt-get update && \
	apt-get install -y --no-install-recommends apt-utils && \
  	apt-get install --yes git \
	python3 \
	python3-pkg-resources \
	build-essential \
	gcc-multilib \
	dh-autoreconf \
	cmake \
	zlib1g-dev 
#Minimap2
WORKDIR /build
RUN git clone https://github.com/lh3/minimap2.git 
WORKDIR /build/minimap2 
RUN git checkout v2.17
RUN make
RUN mkdir -p /build/bin/ && cp -p minimap2 /build/bin/
#BWA
WORKDIR /build
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /build/bwa
RUN git checkout v0.7.17
RUN make
RUN cp -p bwa /build/bin/
#KALIGN3
WORKDIR /build
RUN git clone https://github.com/TimoLassmann/kalign.git
WORKDIR /build/kalign
RUN git checkout 3.1
RUN ./autogen.sh --prefix=/build/
RUN make && make install
#ABPOA
WORKDIR /build
RUN git clone https://github.com/yangao07/abPOA.git
WORKDIR /build/abPOA
RUN git checkout v1.0.3
RUN make
RUN cp bin/abpoa /build/bin/
#SPOA
WORKDIR /build
RUN git clone --recursive https://github.com/rvaser/spoa
WORKDIR /build/spoa
RUN git checkout 3.0.1
WORKDIR /build/spoa/build
RUN cmake -DCMAKE_BUILD_TYPE=Release -Dspoa_build_executable=ON ..
RUN make
RUN cp bin/spoa /build/bin/
#NPGRAPH
WORKDIR /build/app
ADD pom.xml .
RUN ["/usr/local/bin/mvn-entrypoint.sh", "mvn", "verify", "clean", "--fail-never"]
ADD . .
RUN ["mvn", "package"]
RUN cp -p target/assembly-0.2.1-beta.jar /build/bin/
#READFISH
#Waiting for code update for MinKNOW core 4.04

#
# Package stage
#
FROM openjdk:11-jre-slim
RUN apt-get update && apt-get install libgtk-3-0 libglu1-mesa -y && apt-get update

COPY --from=build /build/bin /usr/local/bin
EXPOSE 8080
ENTRYPOINT ["java","-cp","/usr/local/bin/assembly-0.2.1-beta.jar"]
CMD ["org.rtassembly.NPGraphCmd","--msa=kalign3","--aligner=minimap2","--gui"]
