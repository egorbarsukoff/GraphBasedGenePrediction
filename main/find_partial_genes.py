import subprocess

from BCBio import GFF
from Bio import Seq
from Bio import SeqIO
from Bio import SeqRecord


def remove_500bp(contigs, output):
    """
    Remove short contigs
    :param contigs: contigs .fasta file
    :param output: output folder
    :return:
    """
    long_contigs = []
    for i in SeqIO.parse(contigs, "fasta"):
        if len(i.seq) > 500:
            long_contigs.append(i)
    SeqIO.write(long_contigs, output+'long_contigs.fasta', 'fasta')

def prodigal_runner(prodigal_exe, outdir):
    subprocess.call([prodigal_exe, '-i', outdir+'long_contigs.fasta', '-o', outdir+'prodigal_output.gff',
                     '-p', 'meta', '-f', 'gff'])


def partial_genes_extract(contigs_path, outdir):
    """
    Extract partial genes annotated in .gff file from contigs
    :param contigs_path: path to contig file in .fasta
    :param outdir: directory with prodigal_output.gff
    :return: list of SeqRecord
    """
    contigs = SeqIO.to_dict(SeqIO.parse(contigs_path, "fasta"))
    partial_genes = []
    with open(outdir+'prodigal_output.gff') as f:
        for rec in GFF.parse(f, base_dict=contigs):
            for feature in rec.features:
                extract_seq = str(feature.extract(rec.seq))
                if feature.qualifiers['partial'][0] in ['11', '01', '10'] and len(extract_seq) > 110:
                    subrecord = SeqRecord.SeqRecord(Seq.Seq(extract_seq),
                                                    id='gene_'+str(len(partial_genes))+'_len_'+str(len(extract_seq)) +
                                                    '_conf_'+str(feature.qualifiers['conf'][0]),
                                                    description='patrial gene')
                    partial_genes.append(subrecord)
    return partial_genes


def write_seq(outdir, sequences):
    with open(outdir+'partial_genes.fasta', "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")


def find_partial_genes(prodigal_exe, contigs, outdir):
    remove_500bp(contigs, outdir)
    prodigal_runner(prodigal_exe, outdir)
    write_seq(outdir, partial_genes_extract(contigs, outdir))



