#!/bin/bash

for file in *.fasta;do
	echo "$file"
	grep "^>" "$file" >> headers.txt
done

for file in *.gb;do
	echo "$file"
	grep "LOCUS" "$file" >> headers.txt
done
