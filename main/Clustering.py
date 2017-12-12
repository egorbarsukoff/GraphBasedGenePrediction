import re

import numpy as np
import scipy.cluster
import scipy.sparse
import scipy.spatial.distance as ssd

from .GFAgraph import GFAgraph
from .functions import log


class Clustering:
    def __init__(self, graph, orfs):
        self.graph = graph
        self.orfs = orfs
        self.orfs_dict = {}
        self.distances = None
        self.linkage_matrix = None
        self.clusters = dict()


    @classmethod
    def from_files(cls, seq_path, paths_path, graph_path):
        """
        Create object from graph, ORFs files
        :param seq_path: path to .fasta file with ORFs sequences
        :param paths_path: path to .path file with ORF's paths
        :param graph_path: path to graph
        :return: LinkageCluster
        """
        from Bio import SeqIO
        graph = GFAgraph().load_graph(graph_path)
        orfs = []
        with open(seq_path) as seq_file, open(paths_path) as path_file:
            paths_dict = {}
            orfid = ''
            for ind, line in enumerate(path_file.read().split('\n')):
                if ind % 2 == 0:
                    try:
                        orfid = re.findall(r'(ORF_\d+)', line)[0]
                    except:
                        continue
                else:
                    paths_dict[orfid] = [int(re.findall(r'^(\d+),', line)[0])] + re.findall('(\d+)([+-]),', line) + \
                                        [int(re.findall(r',(\d+)$', line)[0])]
            for rec in SeqIO.parse(seq_file, 'fasta'):
                orfs.append((str(rec.seq), cls.add_edges_frame(paths_dict[rec.id], graph), rec.id))
        return cls(graph, orfs)

    def save_linkage_matrix(self, output):
        np.savetxt(output, self.linkage_matrix)

    @staticmethod
    def add_edges_frame(edges, graph):
        """
        Calculate frames for every edge
        :param edges: ORF's path (list)
        :return: list of edge's frames
        """
        frames = []
        for ind, edge in enumerate(edges[1:-1]):
            if ind == 0:
                if edge[1] == '-':
                    frames.append((len(graph.segments[edge[0]]) - edges[0] - 1) % 3)
                else:
                    frames.append(edges[0] % 3)
            else:
                frames.append((frames[-1] + len(graph.segments[edges[ind][0]]) - graph.overlap_size) % 3)
        return [edges[0]] + [edges[i+1] + tuple([frames[i]]) for i in range(len(frames))] + [edges[-1]]

    def distance(self, orf1, orf2):
        """
        Distance between 0(absilutely close ORFs) and 1(almost not overlapping) or 100(absolutely not overlapping)
        :param orf1: ORF's seq
        :param orf2:
        :return: distance >0
        """
        if orf1 == orf2:
            return 0
        dist = 0
        edges_intersect = set(orf1[1][1:-1]).intersection(set(orf2[1][1:-1]))
        if not edges_intersect:
            return 100
        for orf in [orf1, orf2]:
            for ind, edge in enumerate(orf[1]):
                if edge not in edges_intersect:
                    if ind not in [0, 1, len(orf[1]) - 2, len(orf[1]) - 1]:
                        dist += len(self.graph.segments[edge[0]])
                    elif ind == 1:
                        if '+' == orf[1][ind][1]:
                            dist += len(self.graph.segments[edge[0]]) - orf[1][0]
                        else:
                            dist += orf[1][0] - 1
                    elif ind == len(orf[1]) - 2:
                        if '+' == orf[1][ind][1]:
                            dist += orf[1][-1]
                        else:
                            dist += len(self.graph.segments[edge[0]]) - orf[1][-1]
        if dist > 100000:
            print('w')
        return dist/(len(orf1[0]) + len(orf2[0]))

    def create_distance_matrix(self):
        """
        Create distance matrix for self.orfs
        :return:
        """
        self.distances = np.ndarray(shape=(len(self.orfs), len(self.orfs)))
        log('Calculate {0} distances...'.format(int(len(self.orfs) * (len(self.orfs) + 1) / 2)))
        count = 0
        for ind1, orf1 in enumerate(self.orfs):
            for ind2, orf2 in enumerate(self.orfs):
                if count % 1000000 == 0:
                    log('Processed {0} distances'.format(count))
                if ind1 == ind2:
                    break
                if ind1 > ind2:
                    self.distances[ind1, ind2] = self.distance(orf1, orf2)
                    self.distances[ind2, ind1] = self.distances[ind1, ind2]
                count += 1

    def create_linkage_matrix(self):
        self.linkage_matrix = scipy.cluster.hierarchy.linkage(ssd.squareform(self.distances), method='complete')

    def make_clusters(self, distance):
        """
        Perform clustering with given threshold. self.cluster[i] will be contain number of cluster for self.orfs[i]
        :param distance: distance threshold [0, 1]
        :return:
        """
        for i, o in enumerate(scipy.cluster.hierarchy.fcluster(self.linkage_matrix, distance, 'distance')):
            if self.clusters.get(o) is None:
                self.clusters[o] = set()
            self.clusters[o].add(i)

    def output_clustering(self, output):
        """
        Save clustering results in .txt.
        :param output:
        :return:
        """
        with open(output, 'w') as f:
            for i, cluster in enumerate(self.clusters.values()):
                f.write('Cluster {}\n'.format(i))
                for orf in cluster:
                    f.write('{},'.format(self.orfs[orf][2]))
                f.write('\n')



def clustering(orfs, paths, graph, distance, output):
    obj = Clustering.from_files(orfs, paths, graph)
    obj.create_distance_matrix()
    log('Clustering...')
    obj.create_linkage_matrix()
    obj.save_linkage_matrix(output+'linkage.txt')
    obj.linkage_matrix = np.loadtxt(output+'linkage.txt')
    obj.make_clusters(distance)
    obj.output_clustering(output+'clustering_results.txt')