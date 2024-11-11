#!/bin/bash

cd ../../data

while read -r accession; do
  echo "Downloading accession: $accession"
  prefetch $accession
  fasterq-dump $accession | gzip
done < accessions.txt
