#!/bin/bash

sicer SICER -t treatment.bed -c control.bed -s hg38 --opt_output

if python3 compare.py; then
	echo "Test success"
	exit 0
else
	echo "Test failed"
	exit 1
fi 