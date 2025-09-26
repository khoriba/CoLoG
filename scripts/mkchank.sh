#!/usr/bin/env bash

files=($(ls -1 *_R1*.gz | perl -pe 's/_.+//' | perl -pe 's/\.fasta//'));

for FILE in "${files[@]}"
  do
      mkdir "${FILE}" && mv "${FILE}"_*.gz "${FILE}"
  done
