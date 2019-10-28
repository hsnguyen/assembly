# Dockerfile

FROM  phusion/baseimage:0.9.17

# Use ubuntu 16.04 base image
FROM ubuntu:16.04
MAINTAINER  Son Hoang Nguyen <s.hoangnguyen@imb.uq.edu.au>
# set non-interactive mode
ENV DEBIAN_FRONTEND noninteractive

# Install dependencies
RUN apt update && \
    apt install --yes git \
    python3 \
    python3-pkg-resources \
    build-essential \
	gcc-multilib \
	apt-utils \
    zlib1g-dev 
#Minimap2
WORKDIR /tmp
RUN git clone https://github.com/lh3/minimap2.git 
WORKDIR /tmp/minimap2 
RUN git checkout v2.17
RUN make
RUN cp -p minimap2 /usr/local/bin

#BWA
WORKDIR /tmp
RUN git clone https://github.com/lh3/bwa.git
WORKDIR /tmp/bwa
RUN git checkout v0.7.17
RUN make
RUN cp -p bwa /usr/local/bin




# Cleanup
RUN rm -rf /tmp/minimap2 /tmp/bwa
RUN apt-get clean
RUN apt remove --purge --yes git build-essential gcc-multilib apt-utils zlib1g-dev && \
    apt autoremove --purge --yes
