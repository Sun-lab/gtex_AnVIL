docker build -t bioconductor_docker_trecase:devel .

docker image tag bioconductor_docker_trecase:devel sunway1999/bioconductor_trecase:0.1

docker image push sunway1999/bioconductor_trecase:0.1

 