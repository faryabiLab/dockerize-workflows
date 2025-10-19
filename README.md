![Docker](https://img.shields.io/badge/docker-image-blue?logo=docker)
![Size](https://img.shields.io/github/repo-size/faryabilab/dockerize-workflows)
![Cromwell](https://img.shields.io/badge/Cromwell-Workflow%20Ready-4A90E2?logo=scala)
![WDL](https://img.shields.io/badge/WDL-Validated-ff69b4?logo=googlecloud)
# Dockerized Genomic Processing Pipelines
Genomic data processing pipelines written in Workflow Definition Language (WDL), making use of Docker containers to ensure reproducability and limit time spent setting up a compute environment. The Docker containers used have been curated for usage in the Faryabi Lab processing pipelines, and can be found on [Dockerhub](https://hub.docker.com/u/faryabilab).

## Setup
Clone this repository to your local machine:
```
git clone https://github.com/faryabiLab/dockerize-workflows
```
#### Dependencies
* [Cromwell WDL engine jarfile](https://github.com/broadinstitute/cromwell)
* [JDK](https://www.oracle.com/java/technologies/downloads/) 8+
* [Docker](https://www.docker.com/products/personal/)

#### Usage documentation available [here](docs/USAGE.md).
#### Tips for troubleshooting common errors can be found [here](docs/TROUBLESHOOTING.md).
