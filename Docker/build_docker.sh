
docker build -t bioconductor_trecase:0.1 .

docker image tag bioconductor_trecase:0.1 sunway1999/bioconductor_trecase:0.1

docker image push sunway1999/bioconductor_trecase:0.1

