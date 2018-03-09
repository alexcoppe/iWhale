FROM ubuntu:17.10

ENV bwa_version 0.7.17

ADD  http://downloads.sourceforge.net/project/bio-bwa/bwa-${bwa_version}.tar.bz2 /tmp/
COPY checkFiles.py /usr/bin/checkFiles.py
COPY SConstruct /tmp/

RUN apt-get update \
        && apt-get install -y apt-utils \
        && apt-get install -y python scons bzip2 make gcc zlib1g-dev \
        && cd /tmp/ && tar xjvf bwa-${bwa_version}.tar.bz2 \
        && cd /tmp/bwa-${bwa_version} \
        && make \
        && mv /tmp/bwa-${bwa_version}/bwa /usr/bin \
        && rm /tmp/bwa-${bwa_version}.tar.bz2 && rm -rf /tmp/bwa-${bwa_version}


WORKDIR /working

ENTRYPOINT ["python","/usr/bin/checkFiles.py","-s","/tmp"]
