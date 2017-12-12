import argparse
import configparser
import os

from main.remove_orfs import remove_orfs
from main.find_partial_genes import find_partial_genes
from main.functions import log
from main.ORFFinderInGraph import ORFsInGraph
from main.GFAgraph import GFAgraph
from main.Clustering import clustering


def get_config():
    configs = configparser.ConfigParser()
    configs.read('conf.cfg')
    return configs


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('graph', type=str, help='Assembly graph in GFA format')
    parser.add_argument('contigs', type=str, help='Contigs sequnces in fasta format')
    parser.add_argument('outdir', type=str, help='Output folder')
    return parser.parse_args()


if __name__ == '__main__':

    log('__init_log__')
    args = get_args()
    config = get_config()
    try:
        os.mkdir(args.outdir)
    except FileExistsError:
        pass
    try:
        os.mkdir(args.outdir+'temp/')
    except FileExistsError:
        pass

    log('Start ORFs search')
    a = ORFsInGraph(GFAgraph().load_graph(args.graph), int(config['SEARCH SETTINGS']['minlen']),
                    int(config['SEARCH SETTINGS']['maxcyc']), int(config['SEARCH SETTINGS']['maxdepth']))
    a.find_orfs()
    a.write_orfs(args.outdir+'temp/')
    log('Find partial genes in contigs')
    find_partial_genes(config['TOOLS']['prodigal'], args.contigs, args.outdir+'temp/')
    log('Find suitable ORFs')
    remove_orfs(config['TOOLS']['bwa'], args.outdir + 'temp/' + 'ORFs.fasta', args.outdir + 'temp/' +
                'partial_genes.fasta', args.outdir + 'temp/', 1)
    log('Clustering...')
    clustering(args.outdir+'temp/cleaned_ORFs.fasta', args.outdir + 'temp/' + 'ORFs.path', args.graph,
               config['SEARCH SETTINGS']['clustering_setting'], args.outdir + 'temp/')
    os.rename(args.outdir + 'temp/cleaned_ORFs.fasta', args.outdir + 'cleaned_ORFs.fasta')
    os.rename(args.outdir + 'temp/ORFs.path', args.outdir + 'ORFs.path')
    os.rename(args.outdir + 'temp/clustering_results.txt', args.outdir + 'clustering_results.txt')
    log('Done. Results saved in {0}'.format(args.outdir))