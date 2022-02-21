#!/usr/local/bin/bash

# @author T. Paysan-Lafosse
# @brief this script concatains fasta files and generates Eukaryotic microbiome clusters
#usage: bsub -q production-rh74 -M 600000 -R "rusage[mem=600000]" -oo clustering.log -J cluster_euk -Pbigmem -n 16 ./cluster.sh config.ini

if [[ $# -lt 1 ]]; then
    echo "Illegal number of parameters. Usage: ./cluster.sh config_file"
    exit 2
fi

CONFIG_FILE=$1
if [[ -s $CONFIG_FILE ]];then
    . $CONFIG_FILE
else
    echo "Config file ${CONFIG_FILE} not found"
    exit
fi

SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

set -e

mkdir -p $cluster_dir
fasta_file="${cluster_dir}/euk.fasta"
if [[ ! -s $fasta_file  ]];then
    echo "Combining fasta files"
    `cat $output_fasta/* > $fasta_file`
else
    echo "Fasta files already combined"
fi

#Clustering
mmseqs="/nfs/production/interpro/metagenomics/peptide_db/mmseqs"
cd $cluster_dir

DB="${cluster_dir}/euk_seqs"
cluster_file="${DB}.cluster_seq.fa"

if [[ ! -s $cluster_file ]];then
    echo "Starting clustering"
    $mmseqs createdb $fasta_file $DB.mmseqs
    $mmseqs linclust $DB.mmseqs $DB.cluster "${SUBDIR}/tmp/" --min-seq-id 0.5 -c 0.75 --cov-mode 0 --threads 16
    $mmseqs result2repseq $DB.mmseqs $DB.cluster $DB.cluster_rep
    $mmseqs result2flat $DB.mmseqs $DB.mmseqs $DB.cluster_rep $DB.cluster_rep.fa
    $mmseqs createtsv $DB.mmseqs $DB.mmseqs $DB.cluster $DB.cluster.tsv
    $mmseqs createseqfiledb $DB.mmseqs $DB.cluster $DB.cluster_seq
    $mmseqs result2flat $DB.mmseqs $DB.mmseqs $DB.cluster_seq $DB.cluster_seq.fa
    echo "Clustering completed successfully"
fi

if [[ -s $cluster_file ]];then
    cd $mgnify_dir
    # run get_stats.py
    python get_stats.py -i $cluster_file
fi

stat_file="${cluster_file}_percent_euk"
if [[ -s $stat_file ]];then
    # next step: run generate_alignments.py
    python generate_alignments.py -i $stat_file -f Pfam-E -n 100 -b 1
fi