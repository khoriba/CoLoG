#!/bin/bash
LIST=($(awk -F'\\t' '{print $1}' list.txt));
OUTF=$(basename `pwd`);
for DIRNAME in "${LIST[@]}"
 do
    bestseq=""
    if [ -f "${DIRNAME}/${DIRNAME}.draft.fasta"  ];then
	bestseq="${DIRNAME}/${DIRNAME}.draft.fasta"
    fi
    if [ -f "${DIRNAME}/${DIRNAME}.a5miseq-complete.fasta"  ];then
        bestseq="${DIRNAME}/${DIRNAME}.a5miseq-complete.fasta"
    fi
    if [ -f "${DIRNAME}/${DIRNAME}.skesa-complete.fasta"  ];then
        bestseq="${DIRNAME}/${DIRNAME}.skesa-complete.fasta"
    fi
    if [ -f "${DIRNAME}/${DIRNAME}.vcf-consensus-complete.fasta"  ];then
        bestseq="${DIRNAME}/${DIRNAME}.vcf-consensus-complete.fasta"
    fi
    #echo "$bestseq"
    if [ $bestseq != "" ];then
     	cat $bestseq >> temp.mfa
    fi
 done
seqkit replace -p '^(\S+)' -r '{kv}' -k list.txt -o RENAMED_${OUTF}-rebind.fasta temp.mfa
rm temp.mfa
