#!/bin/bash
# Create and run a new MySQL Docker container

NAME="PipelineDatabase"
ROOTPWD="@noah1234"
USER="pipeline"
PWD="Run@pipelines9061"
IMG="mysql/mysql-server:latest"

if docker ps --format '{{.Names}}' | grep -q "^${NAME}$"; then
    echo "Container ${NAME} is already running — stopping it..."
    docker stop "${NAME}"
fi

if docker ps -a --format '{{.Names}}' | grep -q "^${NAME}$"; then
    echo "Container ${NAME} exists — removing it..."
    docker rm "${NAME}"
fi

docker run \
-p 52000:3306 \
--name "${NAME}" \
-e MYSQL_ROOT_PASSWORD="${ROOTPWD}" \
-e MYSQL_DATABASE="${NAME}" \
-e MYSQL_USER="${USER}" \
-e MYSQL_PASSWORD="${PWD}" \
-d "${IMG}"

echo -e "\n\nStarting container, please wait ~30 seconds.\nCheck container status with 'docker ps'"
