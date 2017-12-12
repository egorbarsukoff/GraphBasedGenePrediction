import argparse
import re

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .functions import log, invert


class ORFsInGraph:
    def __init__(self, graph, minlen=500, cycle_multiplicity=3, depth=14):
        """
        :param graph: GFAgrpah object
        :param minlen: minimum len of search
        :param cycle_multiplicity: max passed of one edge
        :param depth: max depth
        """
        self.graph = graph
        self.minlen = minlen
        self.cycle_multiplicity = cycle_multiplicity
        '''
        ORFs represents as dict with seq as key and path as value
        Path is a list with int in first and last pos (pos of start/and) and (edge, orientation) in another places
        '''
        self.orfs = {}
        self.max_depth = depth

    # Not works in _find_stop_codons
    starts = ['ATG', 'GTG', 'TTG']
    stops = ['TAA', 'TGA', 'TAG']

    def write_orfs(self, outdir):
        """
        Write ORFs from self.orfs and create files ORFs.[fasta, path] with sequences and paths
        :param outdir: saves path (with "/" in the end)
        """
        with open(outdir + 'ORFs.fasta', 'w') as f, open(outdir + 'ORFs.path', 'w') as path:
            counter = 0
            for o in self.orfs:
                SeqIO.write(SeqRecord(Seq(o), id='ORF_{0}'.format(counter), description=self.orfs[o][1]), f, 'fasta')
                path.write('ORF_{0} '.format(counter) + 'max_edge_len: {0}\n'.format(self.max_edge_len(o)) +
                           str(self.orfs[o][0][0]) + ',' + ''.join([i[0] + i[1] + ',' for i in self.orfs[o][0][1:-1]]) +
                           str(self.orfs[o][0][-1]) + '\n')
                counter += 1

    def get_orfs(self):
        """
        Return processed ORFs
        :return orfs: list of tuples: (seq, path)
        """
        return list(self.orfs.items())

    def max_edge_len(self, orf):
        """
        Give len of max part of sequences that contains in a single edge
        :param orf: seq of ORF
        :return int: len of max part of sequences that contains in a single edge
        """
        max_len = 0
        for edge in self.orfs[orf][0][2:-2]:
            if len(self.graph.segments[edge[0]]) > max_len:
                max_len = len(self.graph.segments[edge[0]])
        if self.orfs[orf][0][1][1] == '+':
            if max_len < len(self.graph.segments[self.orfs[orf][0][1][0]]) - self.orfs[orf][0][0]:
                max_len = len(self.graph.segments[self.orfs[orf][0][1][0]]) - self.orfs[orf][0][0]
        elif self.orfs[orf][0][1][1] == '-':
            if max_len < self.orfs[orf][0][0] + 1:
                max_len = self.orfs[orf][0][0] + 1
        if self.orfs[orf][0][-2][1] == '+':
            if max_len < self.orfs[orf][0][-1] + 1:
                max_len = self.orfs[orf][0][-1] + 1
        elif self.orfs[orf][0][-2][1] == '-':
            if max_len < len(self.graph.segments[self.orfs[orf][0][-2][0]]) - self.orfs[orf][0][-1] + 1:
                max_len = len(self.graph.segments[self.orfs[orf][0][-2][0]]) - self.orfs[orf][0][-1] + 1
        return max_len

    def _find_stop_codons(self):
        """
        Find all stop-codons if a graph
        :return (string, int, string): list of tuples (edge, locate, orientation)
        """
        codons = []

        def add_codon(start, orient):
            codons.append((segment, start, orient))

        for segment in self.graph.segments:
            for it in re.finditer('TAA|TGA|TAG', self.graph.segments[segment]):
                add_codon(it.start(), '+')
            for it in re.finditer('CTA|TCA|TTA', self.graph.segments[segment]):
                add_codon(it.end(), '-')
        return codons

    def _continuation(self, segment):
        """
        Make reverse compliment (if orientation is '-') and cut the overlapping part of the given edge
        :param (name of segment, orientation)
        :return : requered string
        """
        name, orient = segment
        if orient == '+':
            return self.graph.segments[name][:-self.graph.overlap_size]
        if orient == '-':
            return self.graph.rc_segments[name][:-self.graph.overlap_size]

    def _read_str(self, context):
        """
        Return the countination sequence from given context
        :param context: DFSContext object
        :return string: contination string
        """
        if context.depth != 1:
            s = self._continuation(context.pos) + context.shift
        else:
            if context.pos[2] == '+':
                s = self.graph.segments[context.pos[0]][:context.pos[1] + 3]
            else:
                s = self.graph.rc_segments[context.pos[0]][:len(self.graph.segments[context.pos[0]]) - (context.pos[1]
                                                                                                        - 3)]
        return s

    @staticmethod
    def _reverse_orf_read(s):
        """
        Reads given string back until encounters stop-codon. Also write last encounter of start-codon
        :param s: string
        :return t, t_part_orf, fl, is_stop:
            t - pos of end reading
            t_part_orf - pos of current ORF's start
            fl - True if start codon was read
            is_stop - True if the reading was interrupted by stop codon
        """
        t = len(s)
        t_part_orf = t
        fl = False
        while True:  # do..while
            if s[t - 3:t] in ORFsInGraph.starts:
                t_part_orf = t - 3
                fl = True
            t -= 3
            if (s[t - 3:t] in ORFsInGraph.stops) or t < 3:
                break
        return t, t_part_orf, fl, s[t - 3:t] in ORFsInGraph.stops

    class DFSContext:
        def __init__(self, orf, cand_orf, path, cand_path, pos, shift, used, depth):
            """
            :param orf: (string) biggest ORF in current path
            :param cand_orf: (string) current path sequence
            :param path: (list) path of orf
            :param cand_path: (list) path of cand_orf
            :param pos: current pos (name of edge, side of edge)
            :param shift: (string) 0-2 symbols from previous edge that not read in _revese_orf_reading
            :param used: dict that stored numbers of edges pass
            :param depth: current search depth
            """
            self.orf = orf
            self.cand_orf = cand_orf
            self.path = path
            self.cand_path = cand_path
            self.pos = pos
            self.shift = shift
            self.used = used
            self.depth = depth
            self.is_stop = False

        @classmethod
        def new_context(cls, pos):
            """
            Create start context from finded stop-codon
            :param pos: (edge name, start point, orintation)
            :return: DFSContext object
            """
            return cls('', '', [], [pos[1]], pos, '', {}, 1)

        def update_context(self, obj):
            """
            Do context transformations, update all context field
            :param obj: ORFsInGraph object(for functions access)
            :return:
            """
            # s - current work string
            s = obj._read_str(self)
            # t - pos of end reading
            # t_part_orf - pos of current ORF's start
            # fl - True if start codon was read
            # is_stop - True if the reading was interrupted by stop codon
            t, t_part_orf, fl, self.is_stop = ORFsInGraph._reverse_orf_read(s)
            if self.depth == 1:
                # (edge, pos, orientation) from first point to (edge, orientation)
                self.pos = self.pos[0], self.pos[2]
            self.cand_path = [self.pos] + self.cand_path  # Update cand_path by current pos
            if fl:  # if new ORF was found
                self.orf = s[t_part_orf:] + self.cand_orf  # Save the biggest current ORF
                self.path = [t_part_orf if self.pos[1] == '+' else len(s) - t_part_orf - len(self.shift)] + \
                             self.cand_path  # Save path of the biggest current ORF and add start pos
            self.cand_orf = s[t:] + self.cand_orf  # Update cand_orf by current work string
            self.shift = s[:t]  # New shift
            self.used[self.pos[0:2]] = self.used.get(self.pos, 0) + 1  # Update used edges
            self.pos = self.pos[0], invert(self.pos[1])  # Go to opposide side of edge

        def next_context(self, c):
            """
            Create next context for DFS
            :param c: next (edge, orientation) in path
            :return: DFSContext object
            """
            return ORFsInGraph.DFSContext(self.orf, self.cand_orf, self.path, self.cand_path.copy(),
                                          c, self.shift, self.used.copy(), self.depth + 1)

    def _add_orf(self, context, reason):
        """
        add finded ORF in self.orfs
        :param context: DFSContext where ORF search stopped
        :param reason: reason for stopped
        :return:
        """
        if len(context.orf) > self.minlen and context.orf not in self.orfs and len(context.path) > 3:
            self.orfs[context.orf] = (context.path + [], reason)
            self.orfs[context.orf] += tuple([self.max_edge_len(context.orf)])

    def _check_and_add(self, context):
        """
        Check finded ORF for requerment and call _add_orf with necessary reason
        :param context: DFSContext
        :return:
        """
        if context.depth > self.max_depth:
            self._add_orf(context, 'max_depth')
        elif context.is_stop:
            self._add_orf(context, 'stop')
        elif self.graph.dead_end(context.pos):
            self._add_orf(context, 'dead_end')
        else:
            return False
        return True

    def find_orfs(self):
        """
        Start DFS search from all stop-codons
        :return:
        """

        def dfs(context):
            """
            DFS search for ORF
            :param context: DFSContext object
            :return:
            """
            context.update_context(self)
            if self._check_and_add(context):
                return
            flag = True  # Stay True if DFS is not continue (max cycle cutoff here)
            for c in self.graph.connections[(context.pos[0], context.pos[1])]:
                if context.used.get(c, 0) <= self.cycle_multiplicity:
                    dfs(context.next_context(c))
                    flag = False
            if flag:
                self._add_orf(context, 'max_cyc')

        log('Searching stop codons...')
        stops = self._find_stop_codons()
        log('Done. Total found {0} stop codons'.format(len(stops)))
        log('Starts DFS from {0} points'.format(len(stops)))
        count = 0
        for stop in stops:
            dfs(ORFsInGraph.DFSContext.new_context(stop))
            count += 1
            if count % 100000 == 0:
                log('Processed ' + str(count) + ' of ' + str(len(stops)) + ' stop codons')
        log('Done. Total found {0} ORFs'.format(len(self.orfs)))


def get_args():
    parser = argparse.ArgumentParser(description='find scattered ORFs in graph')
    parser.add_argument('input', type=str, help='graph in GFA1 format')
    parser.add_argument('-o', '--outdir', type=str, default='.', help='folder for output')
    parser.add_argument('-d', '--depth', type=int, default=14, help='max depth of search. default: 14')
    parser.add_argument('-ml', '--minlen', type=int, default=300, help='min lenght of stored ORF. default: 300')
    parser.add_argument('-c', '--maxcyc', type=int, default=3, help='max number of passes of a edge. default: 3')
    return parser.parse_args()


if __name__ == '__main__':
    log('__init_log__')
    args = get_args()
    log('Search ORFs in {0} with: \nMax depth of search: {1}\nMin lenght of ORF: {2}\nMax number of passes of a \
edge: {3}'.format(args.outdir, args.depth, args.minlen, args.maxcyc))
    a = ORFsInGraph(GFAgraph().load_graph(args.input), args.minlen, args.maxcyc, args.depth)
    a.find_orfs()
    a.write_orfs(args.outdir)
    log('Results saved in {0}/orfs.[fasta, paths]'.format(args.outdir))
