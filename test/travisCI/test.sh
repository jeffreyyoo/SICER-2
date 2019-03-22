#!/bin/bash

sicer SICER -t ./test/treatment_1.bed -c ./test/control_1.bed -s hg38 --wig_output

sicer SICER -t ./test/treatment_2.bed -c ./test/control_2.bed -s hg38 --wig_output

sicer_df SICER -t ./test/treatment_1.bed ./test/treatment_2.bed -c ./test/control_1.bed ./test/control_2.bed -s hg38 --wig_output

if python3 ./test/travisCI/compare.py; then
	echo "Test success"
	exit 0
else
	echo "Test failed"
	exit 1
fi
