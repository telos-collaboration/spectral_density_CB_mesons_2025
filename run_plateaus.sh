#!/bin/bash

cd plateaus

snakemake --core 6 --use-conda assets/plots/test.pdf
