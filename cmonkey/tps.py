# vi: sw=4 ts=4 et:
"""tps.py - cMonkey tps specific module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import os
import organism
import scoring
import microarray
import stringdb
import network as nw
import motif
import meme
import datamatrix as dm
import membership as memb
import util
import cmonkey_run

MULTIPROCESSING=True
#MULTIPROCESSING=False

CHECKPOINT_INTERVAL = 100
CHECKPOINT_FILE = None
CACHE_DIR = 'tpscache'

THESAURUS_FILE = 'tps/tps.synonyms.gz'
STRING_LINKS = 'tps/string_links.tps.tab'

NUM_CLUSTERS = 400
ROW_WEIGHT = 6.0
NUM_ITERATIONS = 2000
#NETWORK_SCORE_INTERVAL = 7
#MOTIF_SCORE_INTERVAL = 10
MAX_CLUSTER_ROWS = 50

SEQUENCE_TYPES = ['upstream','downstream']
SEARCH_DISTANCES = {'upstream': (0, 450),'downstream': (0,600)}
SCAN_DISTANCES = {'upstream': (0, 450),'downstream': (0,600)}
UPSTREAM_FILE = 'tps/tps.upstream.-400.50.csv'
DOWNSTREAM_FILE = 'tps/tps.downstream.-100.500.csv'
SEQ_FILENAMES = {'upstream': UPSTREAM_FILE, 'downstream': DOWNSTREAM_FILE }

MAX_MOTIF_WIDTH = 18
ADD_SET_ENRICHMENT = False
ADD_MEME = True
ADD_WEEDER = False
WEEDER_SEQ_TYPE = 'upstream'

"""these are the default meme iterations ("meme.iters") in the R version"""
MOTIF_ITERS = range( 400, 1200, 100 ) + \
              range( 1250, 1500, 50 ) + \
              range( 1525, 1800, 25 ) + \
              range( 1810, max( NUM_ITERATIONS, 1820 ), 10 )

#MOTIF_START_ITERATION = 600
#MOTIF_UPDATE_INTERVAL = 10
#MOTIF_COMPUTE_INTERVAL = 100

mode = 'normal'
#mode = 'debug'
#mode = 'short'
if mode == 'debug':
    NUM_ITERATIONS = 200
    MOTIF_ITERS = [2,100,200]
    NETWORK_SCORE_INTERVAL = 5
    NUM_CLUSTERS = 100

if mode == 'short':
    NUM_ITERATIONS = 500
    MOTIF_ITERS = [100] + range(200,500,50)

def motif_iterations(iteration):
    return iteration in MOTIF_ITERS

#def network_iterations(iteration):
#    return iteration > 0 and iteration % NETWORK_SCORE_INTERVAL == 0

class TpsCMonkeyRun(cmonkey_run.CMonkeyRun):

    def __init__(self, organism_code, ratio_matrix, num_clusters):
        cmonkey_run.CMonkeyRun.__init__(self, organism_code, ratio_matrix, num_clusters=num_clusters)
        self.__organism = None
        self['row_scaling'] = ROW_WEIGHT
        self['string_file'] = STRING_LINKS
        self['cache_dir'] = CACHE_DIR
        self['sequence_types'] = SEQUENCE_TYPES
        self['search_distances'] = SEARCH_DISTANCES
        self['scan_distances'] = SCAN_DISTANCES
        self['memb.max_cluster_rows_allowed'] = MAX_CLUSTER_ROWS
        self['num_iterations'] = NUM_ITERATIONS
        self['multiprocessing'] = MULTIPROCESSING
        self.CHECKPOINT_INTERVAL = CHECKPOINT_INTERVAL

    def organism(self):
        if self.__organism == None:
            self.__organism = self.make_tps()
        return self.__organism

    def make_tps(self):
        """returns a tps organism object"""
        nw_factories = [stringdb.get_network_factory2(self['string_file'], 1.0)]
        return organism.GenericOrganism(
            'tps', THESAURUS_FILE, nw_factories,
            seq_filenames=SEQ_FILENAMES,
            search_distances=self['search_distances'],
            scan_distances=self['scan_distances'])

    def meme_suite(self, seqtype):
        """upstream meme suite"""
        background_file = meme.global_background_file(
            self.organism(), self.ratio_matrix.row_names, seqtype, use_revcomp=True)
        return meme.MemeSuite430(max_width=MAX_MOTIF_WIDTH,
                                 background_file=background_file)

    def make_meme_scoring(self, seqtype, meme_suite, sequence_filters):
        motif_scaling_fun = scoring.get_default_motif_scaling(self['num_iterations'])
        return motif.MemeScoringFunction(
            self.organism(), self.membership(),
            self.ratio_matrix, meme_suite,
            seqtype=seqtype,
            sequence_filters=sequence_filters,
            scaling_func=motif_scaling_fun,
            num_motif_func=motif.default_nmotif_fun,
            #update_in_iteration=scoring.schedule(MOTIF_START_ITERATION, MOTIF_UPDATE_INTERVAL),
            #motif_in_iteration=scoring.schedule(MOTIF_START_ITERATION, MOTIF_COMPUTE_INTERVAL),
            update_in_iteration=motif_iterations,
            motif_in_iteration=motif_iterations,
            config_params=self.config_params)

    def make_weeder_scoring(self, seqtype, meme_suite, sequence_filters):
        motif_scaling_fun = scoring.get_default_motif_scaling(self['num_iterations'])
        return motif.WeederScoringFunction(
            self.organism(), self.membership(), self.ratio_matrix,
            meme_suite,
            seqtype=seqtype,
            sequence_filters=sequence_filters,
            scaling_func=motif_scaling_fun,
            num_motif_func=motif.default_nmotif_fun,
            update_in_iteration=motif_iterations,
            motif_in_iteration=motif_iterations,
            config_params=self.config_params)

    def make_network_scoring(self, scaling_fun):
        return nw.ScoringFunction(self.organism(),
                                  self.membership(),
                                  self.ratio_matrix,
                                  scaling_func=scaling_fun,
                                  run_in_iteration=scoring.schedule(1, 7),
                                  config_params=self.config_params)

    def make_row_scoring(self):
        """returns the row scoring function"""
        sequence_filters = []
        network_scaling_fun = scoring.get_default_network_scaling(self['num_iterations'])
        meme_scoring = None
        weeder_scoring = None

        row_scoring = microarray.RowScoringFunction(
            self.membership(), self.ratio_matrix,
            lambda iteration: self['row_scaling'],
            config_params=self.config_params)
        scoring_funcs = [row_scoring,
                         self.make_network_scoring(network_scaling_fun)]

        if ADD_MEME:
            meme_scoring = self.make_meme_scoring('upstream',
                                                  self.meme_suite('upstream'),
                                                  sequence_filters)

        if ADD_WEEDER:
            weeder_scoring = self.make_weeder_scoring(WEEDER_SEQ_TYPE,
                                                      self.meme_suite(WEEDER_SEQ_TYPE),
                                                      sequence_filters)
        if ADD_MEME and ADD_WEEDER:
            scoring_funcs.append(scoring.ScoringFunctionCombiner(
                    self.membership(),
                    [meme_scoring, weeder_scoring],
                    scaling_func=lambda iteration: 0.5,
                    config_params=self.config_params))
        else:
            if ADD_MEME:
                scoring_funcs.append(meme_scoring)
            if ADD_WEEDER:
                scoring_funcs.append(weeder_scoring)

        return scoring.ScoringFunctionCombiner(self.membership(), scoring_funcs,
                                               config_params=self.config_params)

if __name__ == '__main__':
    print('cMonkey (Python port) (c) 2011, Institute for Systems Biology')
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    if len(sys.argv) < 2:
        print('Usage: ./tps.sh <ratio-file> [checkpoint-file]')
    else:
        if len(sys.argv) > 2:
            CHECKPOINT_FILE = sys.argv[2]
        matrix_factory = dm.DataMatrixFactory([dm.nochange_filter, dm.center_scale_filter])
        ratios = sys.argv[1]
        if mode == 'debug': ratios = 'tps_ratios.1000x4.tsv'
        infile = util.DelimitedFile.read(ratios, has_header=True, quote='\"')
        matrix = matrix_factory.create_from(infile)
        cmonkey_run = TpsCMonkeyRun('tps', matrix, NUM_CLUSTERS)
        if CHECKPOINT_FILE and os.path.exists(CHECKPOINT_FILE): cmonkey_run.run_from_checkpoint(CHECKPOINT_FILE)
        else: cmonkey_run.run()
