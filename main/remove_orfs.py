import os
import re

import pysam
from Bio import SeqIO


def bwa_runner(bwa_exe, orfs_seq, partial_genes, output, threads):
    os.system('{0} index {1}'.format(bwa_exe, partial_genes))
    os.system('{0} mem {2} {1} -t {3} > {4}'.format(bwa_exe, orfs_seq, partial_genes, threads,
                                                    output+'ORFs_on_genes_align.sam'))


def genes_alignments(outdir):
    """
    Make list of ORFs that contain partial genes using the alignment all ORFs to partial genes
    :param outdir: directory with ORFs_on_genes_align.sam
    :return: list with ORFs name
    """
    good_orfs = set()
    c = 0
    align_file = pysam.Samfile(outdir+'ORFs_on_genes_align.sam', 'rb')
    for rec in align_file:
        c += 1
        if rec.cigarstring is not None:
            gene_in_contig_len = int(re.findall(r'_len_(\d+)', rec.reference_name)[0])
            total_match = sum([i[1] for i in rec.cigar if i[0] == 0])
            if total_match / gene_in_contig_len > 0.95:
                good_orfs.add(rec.qname)
    return list(good_orfs)


def write_genes(good_orfs, outdir):
    records = []
    orfs = SeqIO.to_dict(SeqIO.parse(outdir+'ORFs.fasta', 'fasta'))
    for orf in good_orfs:
        records.append(orfs[orf])
    SeqIO.write(records, outdir+'cleaned_ORFs.fasta', 'fasta')


def remove_orfs(bwa_exe, orfs, genes, outdir, threads):
    bwa_runner(bwa_exe, orfs, genes, outdir, threads)
    write_genes(genes_alignments(outdir), outdir)


