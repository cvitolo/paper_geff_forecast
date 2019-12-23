#!/bin/bash
# This script is used to manipulate GEFF data.

module load cdo
cdo -select,date=1980-01-01,2016-12-31 /scratch/rd/nen/perClaudia/erai/fwi.nc /hugetmp/publications/FireForecasting/data/fwi.nc
