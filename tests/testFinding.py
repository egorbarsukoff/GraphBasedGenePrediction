import unittest

from GFAgraph import GFAgraph
from ORFFinderInGraph import ORFsInGraph


class TestSimpleFunctions(unittest.TestCase):
    def setUp(self):
        graph = GFAgraph()
        graph.segments = {'1': 'A1A2A3A4A5A6A7A8A9T1T', '2': 'T1T2T3T4T5T6T7T8T9'}
        graph.rc_segments = {'1': 'A1A9T8T7T6T5T4T3T2T1T', '2': '9A8A7A6A5A4A3A2A1A'}
        graph.connections = {('1', '+'): {('2', '-')}, ('1', '-'): {('2', '-')},
                             ('2', '-'): {('1', '-'), ('1', '+')}}
        graph.overlap_size = 3
        self.obj1 = ORFsInGraph(graph)

    def test_continuation_direct(self):
        self.assertEqual('A1A2A3A4A5A6A7A8A9', self.obj1._continuation(('1', '+')))

    def test_continuation_reverse(self):
        self.assertEqual('A1A9T8T7T6T5T4T3T2', self.obj1._continuation(('1', '-')))

    def test_read_str_direct_not_first(self):
        self.assertEqual(self.obj1._read_str(ORFsInGraph.DFSContext('', '', [], [], ('1', '+'), 'AA', {}, 2)),
                         'A1A2A3A4A5A6A7A8A9AA')

    def test_read_str_reverse_not_first(self):
        self.assertEqual(self.obj1._read_str(ORFsInGraph.DFSContext('', '', [], [], ('1', '-'), 'AA', {}, 2)),
                         'A1A9T8T7T6T5T4T3T2AA')

    def test_read_str_direct_first(self):
        self.assertEqual(self.obj1._read_str(ORFsInGraph.DFSContext('', '', [], [], ('1', 5, '+'), '', {}, 1)),
                         'A1A2A3A4')

    def test_read_str_revese_first(self):
        self.assertEqual(self.obj1._read_str(ORFsInGraph.DFSContext('', '', [], [], ('1', 5, '-'), '', {}, 1)),
                         'A1A9T8T7T6T5T4T3T2T')

    def test_reverse_orf_read_stop_in_the_beggining(self):
        self.assertEqual(
            self.obj1._reverse_orf_read(
                'TTGACCCTCTACCCTTGTGGGAGAGGGCCTTGCGCAGCAAGGGGTGAGGGGGTCTTTAAACAGGTAAGGGCGGAGGGG'), (0, 0, True, False))

    def test_reverse_orf_read_stop_simple(self):
        self.assertEqual(
            self.obj1._reverse_orf_read(
                'TAGACCCTCTACCCTTGTGGGAGAGGGCCTTGCGCAGCAAGGGGTGAGGGGGTCTTTAAACAGGTAAGGGCGGAGGGG'), (3, 78, False, True))


class TestStopCodonSearch(unittest.TestCase):
    def setUp(self):
        graph = GFAgraph()
        graph.segments = {'1': 'TAGTGATAATAG', '2': 'CTATCATTACTA'}
        graph.rc_segments = {'1': 'CTATTATCACTA', '2': 'TAGTAATGATAG'}
        graph.connections = {('1', '+'): {('2', '+')}, ('2', '+'): {('1', '+')}}
        graph.overlap_size = 3
        self.obj2 = ORFsInGraph(graph)

    def test_stop_codon_search_simple(self):
        ans = [('2', 3, '-'), ('2', 6, '-'), ('2', 9, '-'),('1', 0, '+'), ('1', 3, '+'), ('1', 6, '+'), ('1', 9, '+')]
        output = self.obj2._find_stop_codons()
        for a in ans:
            self.assertIn(a, output)

    def test_stop_codon_search_overlaps(self):
        op_ans = [('2', 12, '-'), ('1', 9, '+')]
        output = self.obj2._find_stop_codons()
        self.assertTrue(op_ans[0] in output or op_ans[1] in output)


class TestMaxEdgeLen(unittest.TestCase):
    def setUp(self):
        graph = GFAgraph()
        graph.segments = {'1': '123456789', '2': '123456', '3': '1234567', '4': '12345678910'}
        self.obj = ORFsInGraph(graph)

    def test_max_in_middle(self):
        self.obj.orfs = {'6789123456123456712345': ([5, ('1', '+'), ('2', '+'),('3', '+'),('4', '+'), 5], '')}
        self.assertEqual(self.obj.max_edge_len('6789123456123456712345'), 7)

    def text_max_in_beggining_direst(self):
        self.obj.orfs = {'23456789123456123456712345': ([1, ('1', '+'), ('2', '+'), ('3', '+'), ('4', '+'), 5], '')}
        self.assertEqual(self.obj.max_edge_len('23456789123456123456712345'), 8)

    def test_max_in_beggining_revese(self):
        self.obj.orfs = {'987654321123456123456712345': ([8, ('1', '-'), ('2', '+'), ('3', '+'), ('4', '+'), 5], '')}
        self.assertEqual(self.obj.max_edge_len('987654321123456123456712345'), 9)

    def test_max_in_end_direct(self):
        self.obj.orfs = {'6789123456123456712345678910': ([5, ('1', '+'), ('2', '+'),('3', '+'),('4', '+'), 10], '')}
        self.assertEqual(self.obj.max_edge_len('6789123456123456712345678910'), 11)

    def test_max_in_end_revese(self):
        self.obj.orfs = {'678912345612345670198765432': ([5, ('1', '+'), ('2', '+'), ('3', '+'), ('4', '-'), 2], '')}
        self.assertEqual(self.obj.max_edge_len('678912345612345670198765432'), 10)


if __name__ == '__main__':
    unittest.main()