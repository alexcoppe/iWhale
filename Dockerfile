FROM ubuntu:17.10

RUN apt-get update \
        && apt-get install -y apt-utils \
        && apt-get install -y python scons

WORKDIR /working

ENTRYPOINT ["scons"]
