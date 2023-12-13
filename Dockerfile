# Utilisez l'image de base Ubuntu
FROM ubuntu:latest

RUN apt-get update -y

# Installez les bibliothèques CBLAS et CLAPACK
RUN apt-get install -y libblas-dev liblapack-dev

WORKDIR /app

# Copie
COPY . .


CMD ["/bin/bash"]

