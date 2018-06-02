FROM ubuntu:16.04

ENV bwa_version 0.7.17
ENV picard_version 2.17.11
ENV gatk4_version 4.0.2.1
ENV gatk3_version 3.8-1
ENV PATH "$PATH:/tmp/jre1.8.0_161/bin/"
ENV strelka2_version 2.9.2
ENV varscan_version 2.4.2
ENV snpeff_version 4_3t

ADD  http://downloads.sourceforge.net/project/bio-bwa/bwa-${bwa_version}.tar.bz2 /tmp/
ADD https://github.com/broadinstitute/picard/releases/download/${picard_version}/picard.jar /tmp/
ADD http://javadl.oracle.com/webapps/download/AutoDL?BundleId=230532_2f38c3b165be4555a1fa6e98c45e0808 /tmp/java.tar.gz
ADD https://github.com/broadinstitute/gatk/releases/download/${gatk4_version}/gatk-${gatk4_version}.zip /tmp/
ADD "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=${gatk3_version}-0-gf15c1c3ef" /tmp/gatk3.bz2
ADD https://software.broadinstitute.org/gatk/download/auth?package=M1 /tmp/mutect.zip
ADD http://compgen.bio.unipd.it/downloads/java-7-oracle.tar.gz /tmp/java7.tar.gz
ADD https://github.com/Illumina/strelka/releases/download/v${strelka2_version}/strelka-${strelka2_version}.release_src.tar.bz2 /tmp/
ADD https://github.com/dkoboldt/varscan/releases/download/${varscan_version}/VarScan.v${varscan_version}.jar /tmp/
ADD https://raw.githubusercontent.com/alexcoppe/varscan_accessories/master/vs_format_converter.py /tmp/
ADD https://downloads.sourceforge.net/project/snpeff/snpEff_v${snpeff_version}_core.zip /tmp/


COPY checkFiles.py /usr/bin/checkFiles.py
COPY SConstruct /tmp/
COPY configuration.py /tmp/
COPY Scons_variant_calling /tmp/

RUN apt-get update \
    && apt-get install -y apt-utils \
    && apt-get install -y python scons bzip2 make gcc zlib1g-dev unzip bedtools g++ samtools  \
    && cd /tmp/ && tar xjvf bwa-${bwa_version}.tar.bz2 \
    && cd /tmp/bwa-${bwa_version} \
    && make \
    && mv /tmp/bwa-${bwa_version}/bwa /usr/bin \
    && cd /tmp/ && tar xvzf java.tar.gz \
    && tar xvzf java7.tar.gz \
    && unzip gatk-${gatk4_version}.zip \
    && mv gatk-${gatk4_version} gatk4 \
    && tar xvjf gatk3.bz2 \
    && mv GenomeAnalysisTK-${gatk3_version}-0-gf15c1c3ef gatk3 \
    && unzip /tmp/mutect.zip \
    && tar -xjf strelka-${strelka2_version}.release_src.tar.bz2 \
    && mkdir build && cd build \
    && ../strelka-${strelka2_version}.release_src/configure --jobs=4 --prefix=/tmp/strelka \
    && make -j4 install \
    && cd /tmp/ && mv VarScan.v${varscan_version}.jar VarScan.jar \
    && cd /tmp/ \
    && unzip snpEff_v${snpeff_version}_core.zip \
    && rm /tmp/bwa-${bwa_version}.tar.bz2 && rm -rf /tmp/bwa-${bwa_version} && rm /tmp/java.tar.gz && rm /tmp/gatk-${gatk4_version}.zip && rm /tmp/gatk3.bz2 && rm /tmp/mutect.zip  && rm snpEff_v${snpeff_version}_core.zip  \
    && rm /tmp/java7.tar.gz && rm -rf strelka-${strelka2_version}.release_src.tar.bz2 strelka-${strelka2_version}.release_src

WORKDIR /working

ENTRYPOINT ["python","/usr/bin/checkFiles.py","-s","/tmp"]
