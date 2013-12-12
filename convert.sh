#!/bin/bash
        for i in `seq 1 10`; do
            echo item: $i
            python convertformat.py -i nexus fasta ancestorT${i}.nex
        done
