#!/bin/bash
LIST=($(ls -1 -d * | grep -v "list" | grep -v ".txt"));
RESULT_FILE="result-summary-table.txt";

touch ${RESULT_FILE};

for LINE in "${LIST[@]}"
do 
  if [ -s "${LINE}/result-summary.txt" ]; then
    cat "${LINE}"/result-summary.txt | perl -pe 's/ : /\t/g' | perl -pe 's/$/\thoge/' | HV-change.pl | head -n 1 >> ${RESULT_FILE}
    break
  fi
done

for LINE in "${LIST[@]}"
do 
  if [ -s "${LINE}/result-summary.txt" ]; then
    cat "${LINE}/result-summary.txt" | perl -pe 's/ : /\t/g'| perl -pe 's/$/\thoge/' | HV-change.pl | awk 'NR==2 {print $0}' >> ${RESULT_FILE};
  fi
done
