from Bio import SeqIO
import re
from collections import defaultdict


def bwa_runner(orfs_file, genes_file, bwa):
    import os
    os.system(bwa + ' index ' + orfs_file)
    os.system(bwa + ' mem {0} {1} -t {2} > align.sam'.format(orfs_file, genes_file, '15'))


def read_clusters(orfs_file):
    clusters = defaultdict(set)
    for orf in SeqIO.parse(orfs_file, 'fasta'):
        m = re.findall(r'(cluster\d+)_(ORF_\d+)', orf.id)[0]
        if m:
            clusters[m[0]].add(m[1])
    return clusters


def alignments_stats(clusters):
    import pysam
    i = 0
    clusters_count = defaultdict(int)
    aligned_genes = set()
    for rec in pysam.Samfile('align.sam', 'rb'):
        total_match = sum([i[1] for i in rec.cigar if i[0] == 0])
        total_len = sum([i[1] for i in rec.cigar])
        if rec.cigarstring is not None:
            m = re.findall(r'(cluster\d+)_ORF_\d+', rec.reference_name)[0]
            if total_match/total_len > 0.98:
                i += 1
                clusters_count[m] = clusters_count.get(m, 0) + 1
                aligned_genes.add(rec.qname)
    print(str(sum(map(bool, clusters_count.values()))) + '/' + str(len(clusters)))
    print(len(aligned_genes), i)


def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('orfs', type=str, help='clustered ORFs in fasta')
    parser.add_argument('genes', type=str, help='genes in fasta')
    return parser.parse_args()


def get_config():
    import configparser
    config = configparser.ConfigParser()
    config.read('conf.cfg')
    return config['TOOLS']['bwa']


if __name__ == '__main__':
    bwa_exe = get_config()
    args = get_args()
    bwa_runner(args.orfs, args.genes, bwa_exe)
    alignments_stats(read_clusters(args.orfs))