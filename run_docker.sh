#!/bin/bash
set -eo pipefail

docker build -t python_tabak -f Dockerfile .
docker run -i -v $(pwd):/home/docker/tabak -t python_tabak
