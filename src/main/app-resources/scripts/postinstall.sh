#!/bin/bash

/opt/anaconda/bin/conda install -y --file /application/dependencies/R/packages.list
/opt/anaconda/bin/conda create --name cairo-env r-cairo
