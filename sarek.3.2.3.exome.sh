#!/bin/bash

nextflow run nf-core/sarek --input xaa.csv --outdir test1 --genome GATK.GRCh38 -profile docker --wes TRUE --aligner bwa-mem --save_mapped True --tools haplotypecaller,merge --max_cpus 64 --max_memory 40.GB --max_time 480.h --igenomes_base /disk1/oneomics_analysis/references -resume -r 3.2.3 --intervals xgen-exome-hyb-panel-v2-targets-hg38.bed
