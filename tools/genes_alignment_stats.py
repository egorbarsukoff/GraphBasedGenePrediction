import pysam

def alignment(orfs_file, genes_file, bwa):
    import os
    os.system(bwa + ' index ' + orfs_file)
    os.system(bwa + ' mem {0} {1} -t {2} > align.sam'.format(orfs_file, genes_file, '15'))


def alignment_stats():
    stats = {}
    i = 0
    for rec in pysam.Samfile('align.sam', 'rb'):
        total_match = sum([i[1] for i in rec.cigar if i[0] == 0])
        total_len = sum([i[1] for i in rec.cigar])
        if rec.cigarstring is not None:
            if total_match/total_len > 0.98:
                i += 1
    print(i)

def get_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('orfs', type=str, help='ORFs in fasta')
    parser.add_argument('genes', type=str, help='genes in fasta')
    return parser.parse_args()


def get_config():
    import configparser
    config = configparser.ConfigParser()
    config.read('conf.cfg')
    return config['TOOLS']['bwa']

if __name__ == '__main__':
    args = get_args()
    bwa = get_config()
    alignment(args.orfs, args.genes, bwa)
    alignment_stats()