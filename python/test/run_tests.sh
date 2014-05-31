#!/bin/bash 
# 
# run_tests.sh
# 
# Runs all unit tests. 

tests=( test_*.py )
for test in ${tests[@]}; do
	python "$test"
	if [ $? -ne 0 ]; then
		break
	fi
done

