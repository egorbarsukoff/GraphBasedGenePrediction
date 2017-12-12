from .functions import log, rev_comp, invert


class GFAgraph:
    def __init__(self):
        self.segments = {}
        self.connections = {}
        # '+' - right side of segment
        # '-' - left side of segment
        # connection ('1', '+') - ('2', '-') report that right side of '1' connected with left side of '2'
        self.rc_segments = {}
        self.overlap_size = None  # Overlap size must be equal for all links and SIGAR must be like (\d+)M

    def load_graph(self, filein):
        """
        Load graph functions
        :param filein: file with GFA graph
        :return: GFAgraph object
        """
        log('Loading graph...')
        self._gfa_parse(filein)
        log('Compute reverse compliments...')
        self.create_rc_seg()
        log('Graph loaded. Total {0} nodes and {1} links'.format(len(self.segments),
                                                                 int(len(self.connections) / 2)))
        return self

    def create_rc_seg(self):
        for s in self.segments.items():
            self.rc_segments[s[0]] = rev_comp(s[1])

    def _gfa_parse(self, file):
        """
        Parse links and segments if .gfa and fill fields
        :param file: file with GFA graph
        :return:
        """
        with open(file) as f:
            for line in f.readlines():
                line = line.rstrip('\n').split('\t')
                if line[0] == 'S':
                    self.segments[line[1]] = line[2]
                elif line[0] == 'L':
                    if self.connections.get((line[3], invert(line[4]))) is None:
                        self.connections[(line[3], invert(line[4]))] = set()
                    if self.connections.get((line[1], line[2])) is None:
                        self.connections[(line[1], line[2])] = set()
                    self.connections[(line[3], invert(line[4]))].add((line[1], line[2]))
                    self.connections[(line[1], line[2])].add((line[3], invert(line[4])))
                    if self.overlap_size is None:
                        self.overlap_size = int(line[5][:-1])

    def dead_end(self, pos):
        """
        True if pos is not dead end, otherwise False
        :param pos: current pos
        :return: bool
        """
        return self.connections.get(pos) is None
