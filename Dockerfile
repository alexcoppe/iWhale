FROM ubuntu:17.10

ENV bwa_version 0.7.17
ENV picard_version 2.17.11
ENV gatk4_version 4.0.2.1
ENV gatk3_version 3.8-1
ENV PATH "$PATH:/tmp/jre1.8.0_161/bin/"

ADD  http://downloads.sourceforge.net/project/bio-bwa/bwa-${bwa_version}.tar.bz2 /tmp/
ADD https://github.com/broadinstitute/picard/releases/download/${picard_version}/picard.jar /tmp/
ADD http://javadl.oracle.com/webapps/download/AutoDL?BundleId=230532_2f38c3b165be4555a1fa6e98c45e0808 /tmp/java.tar.gz
ADD https://github.com/broadinstitute/gatk/releases/download/${gatk4_version}/gatk-${gatk4_version}.zip /tmp/
ADD "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=${gatk3_version}-0-gf15c1c3ef" /tmp/gatk3.bz2
ADD https://software.broadinstitute.org/gatk/download/auth?package=M1 /tmp/mutect.zip

COPY checkFiles.py /usr/bin/checkFiles.py
COPY SConstruct /tmp/
COPY config.py /tmp/

RUN apt-get update \
        && apt-get install -y apt-utils \
        && apt-get install -y python scons bzip2 make gcc zlib1g-dev unzip bedtools \
        && cd /tmp/ && tar xjvf bwa-${bwa_version}.tar.bz2 \
        && cd /tmp/bwa-${bwa_version} \
        && make \
        && mv /tmp/bwa-${bwa_version}/bwa /usr/bin \
        && cd /tmp/ && tar xvzf java.tar.gz \
        && unzip gatk-${gatk4_version}.zip \
        && mv gatk-${gatk4_version} gatk4 \
        && tar xvjf gatk3.bz2 \
        && mv GenomeAnalysisTK-${gatk3_version}-0-gf15c1c3ef gatk3 \
        && unzip /tmp/mutect.zip \
        && rm /tmp/bwa-${bwa_version}.tar.bz2 && rm -rf /tmp/bwa-${bwa_version} && rm /tmp/java.tar.gz && rm /tmp/gatk-${gatk4_version}.zip && rm /tmp/gatk3.bz2 && rm /tmp/mutect.zip


WORKDIR /working

ENTRYPOINT ["python","/usr/bin/checkFiles.py","-s","/tmp"]
