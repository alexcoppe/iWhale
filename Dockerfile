FROM ubuntu:17.10

ENV bwa_version 0.7.17
ENV picard_version 2.17.11

ADD  http://downloads.sourceforge.net/project/bio-bwa/bwa-${bwa_version}.tar.bz2 /tmp/
ADD https://github.com/broadinstitute/picard/releases/download/${picard_version}/picard.jar /tmp/
ADD http://javadl.oracle.com/webapps/download/AutoDL?BundleId=230532_2f38c3b165be4555a1fa6e98c45e0808 /tmp/java.tar.gz

COPY checkFiles.py /usr/bin/checkFiles.py
COPY SConstruct /tmp/
COPY config.py /tmp/

RUN apt-get update \
        && apt-get install -y apt-utils \
        && apt-get install -y python scons bzip2 make gcc zlib1g-dev \
        && cd /tmp/ && tar xjvf bwa-${bwa_version}.tar.bz2 \
        && cd /tmp/bwa-${bwa_version} \
        && make \
        && mv /tmp/bwa-${bwa_version}/bwa /usr/bin \
        && cd /tmp/ && tar xvzf java.tar.gz \
        && rm /tmp/bwa-${bwa_version}.tar.bz2 && rm -rf /tmp/bwa-${bwa_version} && rm /tmp/java.tar.gz


WORKDIR /working

ENTRYPOINT ["python","/usr/bin/checkFiles.py","-s","/tmp"]
