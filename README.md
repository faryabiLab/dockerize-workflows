![image](https://github.com/user-attachments/assets/546e29be-0802-4272-843f-786062e37367)
# Dockerized Genomic Processing Pipelines
Genomic data processing pipelines written in Workflow Definition Language (WDL), making use of Docker containers to ensure reproducability and limit time spent setting up a compute environment. The Docker containers used have been curated for usage in the Faryabi Lab processing pipelines, and can be found on [Dockerhub](https://hub.docker.com/u/faryabilab).

## Setup
Clone this repository to your local machine:
```
git clone https://github.com/faryabiLab/dockerize-workflows
```
#### Dependencies
* [Cromwell WDL engine jarfile](https://github.com/broadinstitute/cromwell)
* [JDK](https://www.oracle.com/java/technologies/downloads/)
* [Docker](https://www.docker.com/products/personal/)

#### Usage documentation available [here](docs/USAGE.md).
#### Tips for troubleshooting common errors can be found [here](TROUBLESHOOTING.md).
