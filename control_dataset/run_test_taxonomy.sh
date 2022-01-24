#!/bin/bash

# needs haybaler to have been run first
cd haybaler_output/
cp ../*.sh ../*.py .
bash run_haybaler_tax.sh
