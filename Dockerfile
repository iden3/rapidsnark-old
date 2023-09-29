FROM ubuntu:jammy

RUN apt-get update && apt-get upgrade -y

# Update the package list and install necessary dependencies
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt install -y nodejs npm cmake build-essential pkg-config libssl-dev libgmp-dev libsodium-dev nasm awscli git tar

RUN npm install -g yarn npx

# Clone rapidsnark
COPY . ./rapidsnark
WORKDIR /rapidsnark

RUN git submodule init
RUN git submodule update
RUN npm install
RUN npx task createFieldSources
RUN npx task buildPistache
RUN npx task buildProver
RUN chmod +x /rapidsnark/build/prover