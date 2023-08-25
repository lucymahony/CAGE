# Unit tests for triad_CAGE_peak_comparisons.py

import unittest
import pandas as pd


from traid_CAGE_peak_comparisons import  make_gene_dict, check_for_cage_peak, score_transcript,\
    check_for_dominant_cage_peak, describe_pattern_true_and_false


class TestTriadCAGE(unittest.TestCase):
    def test_make_gene_dict(self):
        test_triad_table = pd.DataFrame({'A': ['TraesCS7A02G360600'],
                                         'B': ['TraesCS7B02G267100'],
                                         'D': ['TraesCS7D02G362400']})
        expected_gene_dict = {
            'TraesCS7A02G360600': ['7A', '+', 534003132],
            'TraesCS7B02G267100': ['7B', '-', 490352447],
            'TraesCS7D02G362400': ['7D', '-', 465601777]}

        self.assertDictEqual(make_gene_dict(test_triad_table), expected_gene_dict)

    def test_check_for_cage_peak(self):
        test_cage_peaks = pd.DataFrame({'group': [1, 1, 1, 1, 1],
                                        'group_name': ['IS', 'IS', 'IS', 'IS', 'IS'],
                                        'seqnames': ['chr1A', 'chr1A', 'chr1A', 'chr1A', 'chr1A'],
                                        'start': [10, 20, 30, 40, 50],
                                        'end': [15, 25, 35, 45, 55],
                                        'width': [1, 1, 1, 1, 1],
                                        'strand': ['+', '+', '+', '+', '+'],
                                        'cluster': [1, 1, 1, 1, 1],
                                        'nr_ctss': [1, 1, 1, 1, 1],
                                        'dominant_ctss': [1, 1, 1, 1, 1],
                                        'tpm': [1, 1, 1, 1, 1],
                                        'tpm.dominant_ctss': [1, 1, 1, 1, 1],
                                        'min_density': [1, 1, 1, 1, 1],
                                        'max_density': [1, 1, 1, 1, 1]})
        # True tests
        self.assertTrue(check_for_cage_peak('chr1A', '+', 1, 13, test_cage_peaks))
        self.assertTrue(check_for_cage_peak('chr1A', '+', 15, 19, test_cage_peaks))
        self.assertTrue(check_for_cage_peak('chr1A', '+', 15, 25, test_cage_peaks))
        self.assertTrue(check_for_cage_peak('chr1A', '+', 16, 20, test_cage_peaks))
        # False tests
        self.assertFalse(check_for_cage_peak('chr1A', '+', 1, 9, test_cage_peaks))
        self.assertFalse(check_for_cage_peak('chr1A', '+', 16, 19, test_cage_peaks))
        self.assertFalse(check_for_cage_peak('chr1A', '-', 1, 13, test_cage_peaks))

    def test_check_for_dominant_cage_peak(self):
        test_cage_peaks = pd.DataFrame({'group': [1, 1, 1, 1, 1],
                                       'group_name': ['IS', 'IS', 'IS', 'IS', 'IS'],
                                       'seqnames': ['chr1A', 'chr1A', 'chr1A', 'chr1A', 'chr1A'],
                                       'start': [10, 20, 30, 40, 50],
                                       'end': [15, 25, 35, 45, 55],
                                       'width': [1, 1, 1, 1, 1],
                                       'strand': ['+', '+', '+', '+', '+'],
                                       'cluster': [1, 1, 1, 1, 1],
                                       'nr_ctss': [1, 1, 1, 1, 1],
                                       'dominant_ctss': [1, 1, 1, 1, 1],
                                       'tpm': [1, 1, 1, 1, 1],
                                       'tpm.dominant_ctss': [1, 1, 1, 1, 1],
                                       'min_density': [1, 1, 1, 1, 1],
                                       'max_density': [1, 1, 1, 1, 1]})
        self.assertTrue(check_for_dominant_cage_peak('chr1A', '+', 0, 2, test_cage_peaks))
        self.assertTrue(check_for_dominant_cage_peak('chr1A', '+', 1, 2, test_cage_peaks))
        self.assertTrue(check_for_dominant_cage_peak('chr1A', '+', 0, 1, test_cage_peaks))
        self.assertFalse(check_for_cage_peak('chr1A', '+', 2, 9, test_cage_peaks))
    def test_score_transcript(self):
        test_gene_dict = {
            'TraesCS7A02G360600': ['chr7A', '+', 534879418]}
        test_transcript = 'TraesCS7A02G360600'
        test_cage_peaks = pd.DataFrame({'group': [3, 4, 4, 4],
                                        'group_name': ['RO', 'SP', 'SP', 'SP'],
                                        'seqnames': ['chr7A', 'chr7A', 'chr7A', 'chr7A'],
                                        'start': [534396840, 534879418, 534396840, 534879994],
                                        'end': [534396855, 534879586, 534396855, 534880039],
                                        'width': [16, 169, 16, 46],
                                        'strand': ['+', '+', '+', '+'],
                                        'cluster': [35263, 37195, 38151, 38152],
                                        'nr_ctss': [2, 2, 2, 4],
                                        'dominant_ctss': [534879418, 534879419, 534396841, 534880039],
                                        'tpm': [1, 1, 1, 1],
                                        'tpm.dominant_ctss': [1, 1, 1, 1],
                                        'min_density': [1, 1, 1, 1],
                                        'max_density': [1, 1, 1, 1]})
        expected = score_transcript(test_transcript, test_gene_dict, 1500, 0, test_cage_peaks)
        self.assertTrue(expected)

    def test_describe_pattern_true_and_false(self):
        data = pd.DataFrame({
            'A': [True, False, True, False, True, True],
            'B': [True, True, False, True, False, True],
            'D': [True, True, False, True, False, True]
        })
        actual = describe_pattern_true_and_false(data)
