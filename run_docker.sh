#!/bin/bash
set -eo pipefail

docker build -t python_medaka -f Dockerfile .
docker run -i -v $(pwd):/home/docker/model -t python_medaka
