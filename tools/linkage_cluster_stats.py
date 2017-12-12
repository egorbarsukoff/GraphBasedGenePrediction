from main.ORFsClustering import Clusters
import re
from collections import defaultdict


def bwa_runner(orfs_file, genes_file, bwa):
    import os
    # os.system(bwa + ' index ' + orfs_file)
    os.system(bwa + ' mem {0} {1} -t {2} > align.sam'.format(orfs_file, genes_file, '15'))


def alignments_stats(oc):
    import pysam
    i = 0
    clusters_count = defaultdict(int)
    for rec in pysam.Samfile('align.sam', 'rb'):
        total_match = sum([i[1] for i in rec.cigar if i[0] == 0])
        total_len = sum([i[1] for i in rec.cigar])
        if rec.cigarstring is not None:
            if total_match/total_len > 0.98:
                i += 1
                clusters_count[oc[int(re.findall(r'ORF_(\d+)', rec.reference_name)[0])]] += 1
    return sum(map(bool, clusters_count.values()))


def clusters_plot(a):
    import numpy as np
    import matplotlib as mpl
    mpl.use('Agg')
    import matplotlib.pyplot as plt
    number_of_clusters = []
    good_clusters = []
    for i in np.arange(0, 1, 0.05):
        a.clustering(i)
        number_of_clusters.append(len(a.clusters))
        good_clusters.append(alignments_stats(a.orfs_clusters))
        a.reset()
    plt.plot(np.arange(0, 1, 0.05), number_of_clusters, label='Clusters')
    plt.plot(np.arange(0, 1, 0.05), good_clusters, label='Good clusters')
    plt.xlabel('Non-intersection ratio')
    plt.ylabel('N')
    plt.savefig('plot.png')


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
    # bwa_runner('/home/ebarsukov/python/genepred/data/meta/1k7d/ORFs.fasta', '/home/ebarsukov/Bacterials/meta/scattered_genes/sc_genes_seq.fasta', get_config())
    a = Clusters.from_file('/home/ebarsukov/python/genepred/data/meta/1k7d/linkage_matrix.npy')
    # clusters_plot(a)
    a.clustering(-0.02)
    print(alignments_stats(a.orfs_clusters))
