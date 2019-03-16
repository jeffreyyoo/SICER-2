#!/bin/bash

sicer SICER -t ./test/treatment.bed -c ./test/control.bed -s hg38 --opt_output

if python3 ./test/travisCI/compare.py; then
	echo "Test success"
	exit 0
else
	echo "Test failed"
	exit 1
fi 