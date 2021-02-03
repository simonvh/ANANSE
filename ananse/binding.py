#!/usr/bin/env python

# Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
#
# This module is free software. You can redistribute it and/or modify it under
# the terms of the MIT License, see the file COPYING included with this
# distribution.

"""Predict TF binding network"""

# Python imports
import os
import pickle
from tqdm import tqdm
import warnings
from tempfile import NamedTemporaryFile
from loguru import logger

import numpy as np
import pandas as pd
from sklearn.preprocessing import minmax_scale
import dask.dataframe as dd
from scipy.stats import rankdata

from pybedtools import BedTool
from genomepy import utils
from gimmemotifs.scanner import Scanner
from gimmemotifs.motif import read_motifs
from gimmemotifs.utils import as_fasta, pfmfile_location

from ananse import __file__
from ananse.utils import mytmpdir, set_width

warnings.filterwarnings("ignore")


def clear_tfs(motifs2factors, tffile, include_notfs=False, rm_curated=True):
    """
    filter unreal TFs from motif database

    Arguments:
        motifs2factors {[type]} -- [motifs2factors]
        tffile {[type]} -- [real tfs]

    Returns:
        [type] -- [motifs2factors]
    """
    ft = pd.read_csv(motifs2factors, sep="\t")

    ft['Factor'] = ft.Factor.str.upper()
    if rm_curated:
        # "Curated" is manually curated or direct evidence for binding. For instance a ChIP-seq predicted motif is an N in this column
        ft = ft.loc[ft.Curated == "Y"]
    if not include_notfs:
        tfs = pd.read_csv(tffile, header = None)[0].tolist()
        ft = ft.loc[ft.Factor.isin (tfs)]
    # replace T to TBXT
    ft = ft.replace("T", "TBXT")
    ft = ft.replace("t", "tbxt")

    ft.rename(columns={"Factor": "factor"}, inplace=True)
    return ft


class Binding(object):
    def __init__(self, ncore=1, genome="hg38", gene_bed=None, pfmfile=None, include_notfs=False, rm_curated=True, etype="hg38H3K27ac", tffile=None):
        self.ncore = ncore
        self.genome = genome

        # dream_model.txt is the logistic regression model.
        package_dir = os.path.dirname(__file__)
        # self.etype = etype

        # TODO: instead of etype, expose model selection an input variable, with these 2 as mentioned included.
        self.model = os.path.join(package_dir, "db", "dream_model_p300.pickle")
        # if self.genome == "hg38" and self.etype == "hg38H3K27ac":
        #     self.model = os.path.join(package_dir, "db", "dream_model_h3k27ac.pickle")
        # elif self.etype == "p300" or self.etype == "ATAC":
        #     self.model = os.path.join(package_dir, "db", "dream_model_p300.pickle")
        # else:
        #     raise TypeError("""The input enhancer data type should hg38H3K27ac, p300 or ATAC.
        #     It is not possible set -e to hg38H3K27ac if the genome is not hg38.
        #     Please provide a enhancer type with -e argument. By default is hg38H3K27ac.""")

        # filter tfs?
        self.include_notfs = include_notfs
        # remove curated?
        self.rm_curated = rm_curated

        # load real tfs
        self.tffile = tffile
        if self.tffile is None:
            self.tffile = os.path.join(package_dir, "db", "tfs.txt")
        # self.tffile = "db/tfs.txt"

        # Motif information file
        self.pfmfile = pfmfile_location(pfmfile)
        self.motifs2factors = self.pfmfile.replace(".pfm", ".motif2factors.txt")

        self.filtermotifs2factors = clear_tfs(self.motifs2factors, self.tffile, self.include_notfs, self.rm_curated)

    @staticmethod
    def get_peak_scores(enhancer_regions_bed, outfile):
        """
        accepts a BED4 file with a scoring metric in the 4th column (originally RPKM).
        """
        # When we built model, the peak intensity was ranked and scaled.
        peaks = pd.read_table(enhancer_regions_bed, names=["chrom", "start", "end", "peak_score"])
        peaks["peak"] = (
            peaks["chrom"].astype(str)
            + ":"
            + peaks["start"].astype(str)
            + "-"
            + peaks["end"].astype(str)
        )
        # TODO: IMPORTANT: log10 transform may not be OK for every scoring factor. Expose as variable?
        peaks["log10_peak_score"] = np.log10(peaks["peak_score"] + 1)
        peaks["Scaled_peak_score"] = minmax_scale(peaks["log10_peak_score"])
        peaks["Ranked_peak_score"] = minmax_scale(rankdata(peaks["log10_peak_score"]))

        # TODO: original column names ------v
        # add = peaks["peakRPKM"][peaks["peakRPKM"] > 0].min()
        # peaks["log10_peakRPKM"] = np.log10(peaks["peakRPKM"] + add)
        # peaks["peakRPKMScale"] = minmax_scale(peaks["log10_peakRPKM"])
        # peaks["peakRPKMRank"] = minmax_scale(rankdata(peaks["log10_peakRPKM"]))

        cols = ["peak", "peak_score", "log10_peak_score", "Scaled_peak_score", "Ranked_peak_score"]
        peaks[cols].to_csv(outfile, sep="\t", index=False)

    def get_motif_scores(self, enhancer_regions_bed, outfile):
        """
        Scan for TF binding motifs in potential enhancer regions.
        """
        pfmscorefile = outfile  # NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
        seqs = [s.split(" ")[0] for s in as_fasta(enhancer_regions_bed, genome=self.genome).ids]

        s = Scanner(ncpus=self.ncore)
        s.set_motifs(self.pfmfile)
        s.set_genome(self.genome)
        s.set_threshold(threshold=0.0)

        # generate GC background index
        _ = s.best_score([], zscore=True, gc=True)

        with open(self.pfmfile) as f:
            motifs = read_motifs(f)

        # Run 10k peaks per scan.
        chunksize = 10000

        open(outfile, 'w').close()
        with tqdm(total=len(seqs)) as pbar:
            for chunk in range(0, len(seqs), chunksize):
                pfm_score = []
                chunk_seqs = seqs[chunk:chunk+chunksize]
                # We are using GC-normalization for motif scanning as many enhancer binding regions are GC-enriched.
                chunk_scores = s.best_score(chunk_seqs, zscore=True, gc=True)
                for seq, scores in zip(chunk_seqs, chunk_scores):
                    for motif, score in zip(motifs, scores):
                        pfm_score.append([motif.id, seq, score])
                    pbar.update(1)
                pfm_score = pd.DataFrame(pfm_score, columns=["motif", "enhancer", "zscore"])
                pfm_score = pfm_score.set_index("motif")
                pfm_score["zscoreRank"] = minmax_scale(rankdata(pfm_score["zscore"]))

                # When we built model, rank and minmax normalization was used.
                cols = ["enhancer", "zscore", "zscoreRank"]
                write_header = True if chunk == 0 else False
                # TODO: default mode was 'w'. does that not overwrite every loop iteration? Set to 'a' (append)
                pfm_score[cols].to_csv(pfmscorefile, sep="\t", header=write_header, mode='a')

    def get_binding_score(self, pfm, peak):
        """Infer TF binding score from motif z-score and peak intensity.

        Arguments:
            pfm {[type]} -- [motif scan result]
            peak {[type]} -- [peak intensity]

        Returns:
            [type] -- [the predicted tf binding table]
        """

        # Load model
        with open(self.model, "rb") as f:
            clf = pickle.load(f)

        ft = self.filtermotifs2factors

        r = pfm.merge(peak, left_on="enhancer", right_on="peak")[
            ["motif", "enhancer", "zscore", "log10_peakRPKM"]
        ]
        r = r.merge(ft, left_on="motif", right_on="Motif")
        r = r.groupby(["factor", "enhancer"])[["zscore", "log10_peakRPKM"]].mean()
        r = r.dropna().reset_index()

        table = r.compute(num_workers=self.ncore)
        table["binding"] = clf.predict_proba(table[["zscore", "log10_peakRPKM"]])[:, 1]

        return table

    def run_binding(self, peak_bed, outfile):
        intermediate_dir = os.path.join(os.path.dirname(outfile), "intermediate")
        utils.mkdir_p(intermediate_dir)

        logger.info("Peak initialization")
        filter_bed = os.path.join(intermediate_dir, "filter.bed")
        if not os.path.exists(filter_bed) or os.stat(filter_bed).st_size == 0:
            set_width(self.genome, peak_bed, filter_bed)

        logger.info("Motif scan")
        pfm_weight = os.path.join(intermediate_dir, "pfm_weights.tsv")
        if not os.path.exists(pfm_weight) or os.stat(pfm_weight).st_size == 0:
            self.get_motif_scores(filter_bed, pfm_weight)

        logger.info("Predicting TF binding sites")
        peak_weight = os.path.join(intermediate_dir, "peak_weights.tsv")
        if not os.path.exists(peak_weight) or os.stat(peak_weight).st_size == 0:
            self.get_peak_scores(filter_bed, peak_weight)

        pfm = dd.read_csv(pfm_weight, sep="\t")
        peak = dd.read_csv(peak_weight, sep="\t", blocksize=200e6)
        table = self.get_binding_score(pfm, peak)

        logger.info("Save results")
        table.to_csv(outfile, sep="\t", index=False)
