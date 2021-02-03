import os
import logging

from genomepy import Genome, utils
from ananse import exceptions
from ananse.utils import set_width

# use logger settings defined in __init__
logger = logging.getLogger(__name__)


class Enhancer(object):
    """
    This implementation:
    Create a BED3+1 file, the 4th column containing a scoring method.
    The output width of the regions is normalized to 200 bp, centered around a summit if possible.
    Regions that don't fit after normalization are shifted/shrunk.

    Quans implementation:
    The scoring column was created by counting the number of BAM reads intersecting with the BED peaks regions

    narrowPeak example:                                   v-- 7: signal           v-- 9: summit
    Chr1    481     681     name1        37      .       4.05358 5.95365 3.71434 101
    Chr1    85875   86082   name2        27      .       3.52564 4.84918 2.72004 166

    broadPeak example:                                                                                          v-- 13: signal (?)
    chr1	11012359	11013547	name1	11	.	11012359	11013547	0	4	1,682,362,1	0,8,819,1187	0	0	0
    chr1	11059288	11061047	name2	70	.	11059288	11061047	0	4	1,251,432,1	0,69,1083,1758	0	0	0
    """
    def __init__(self, genome, peak, output, signal_column=6, summit_column=-1):
        self.genome = genome
        self.peak = peak
        self.output = output
        self.outdir = os.path.dirname(self.output)
        self.signal_column = signal_column
        self.summit_column = summit_column

    @staticmethod
    def bedn2bed4(bed_in, bed_out, signal_column=6):
        """
        Converts the input BED to a BED3+1
        columns: chromosome, start, end and score_column

        broadPeak/narrowPeak column 7: "signalValue - Measurement of overall (usually, average) enrichment for the region."
        """
        with open(bed_in) as old, open(bed_out, "w") as new:
            for line in old:
                line = line.split()
                newline = "\t".join(line[0:3] + [line[signal_column]]) + "\n"
                new.write(newline)

    def run_enhancer(self):
        utils.mkdir_p(os.path.join(self.outdir, "intermediate"))
        intermediate_bed = os.path.join(self.outdir, "intermediate", "wide_bed4.bed")

        set_width(genome=self.genome, bed_in=self.peak, bed_out=intermediate_bed, summit_col=self.summit_column)
        self.bedn2bed4(intermediate_bed, self.output, self.signal_column)

        if os.path.exists(self.output):
            utils.rm_rf(os.path.join(self.outdir, "intermediate"))

# def bam_cov(bamfile, peakfile, coverage_bed):
#     """
#     Identify and score enhancer sites by quantifying overlap in enhancer data
#     """
#     bed = BedTool(peakfile)
#     bed.multi_bam_coverage(bams=[bamfile], output=coverage_bed)


# def set_peak_size(genome, peak, seqlen=200):
#     """
#     normalize peaks width to SEQLEN bp (centered around the summit)
#     shifts the peak to best fit on the chromosome
#     """
#     half_seqlen = int(seqlen/2)
#     chrom_sizes = Genome(genome).sizes
#
#     npeaks = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#     with open(peak) as peakfile, open(npeaks.name, "w") as npeakfile:
#         for line in peakfile:
#             a = line.split()
#             chrm = a[0]
#             start = int(a[1])
#             summit = int(a[9])
#
#             if chrm not in chrom_sizes.keys():
#                 logger.error(f"Peak file {os.path.basename(peak)} "
#                              f"contains a chromosome notation not present in the genome: {chrm}")
#
#             chrm_len = chrom_sizes[chrm]
#             if chrm_len <= seqlen:
#                 npeakfile.write(f"{chrm}\t{0}\t{chrm_len}\n")
#             else:
#                 # adjust the summit for the chromosome boundaries
#                 nsummit = start+summit
#                 nsummit = max(nsummit, 0 + half_seqlen)
#                 nsummit = min(nsummit, chrm_len - half_seqlen)
#                 npeakfile.write(f"{chrm}\t{nsummit-half_seqlen}\t{nsummit+half_seqlen}\n")
#
#     return npeaks.name


# def shrink_bed(genome, bed_in, bed_out, seqlen=200, skip_missing_chrm=True):
#     """
#     normalize BED region width to SEQLEN bp
#     shifts the peak to best fit on the chromosome
#     """
#     half_seqlen = seqlen // 2
#     chrom_sizes = Genome(genome).sizes
#
#     with open(bed_in) as old, open(bed_out, "w") as new:
#         for line in old:
#             line = line.split()
#             chrm = str(line[0])
#             start = int(line[1])
#             end = int(line[2])
#             rest = line[3:]
#             summit = (start + end) // 2
#
#             if chrm not in chrom_sizes.keys():
#                 if skip_missing_chrm:
#                     continue
#                 else:
#                     exceptions.error(f"File {os.path.basename(bed_in)} "
#                                      f"contains a chromosome notation not present in the genome: {chrm}")
#
#             chrm_len = chrom_sizes[chrm]
#             if chrm_len <= seqlen:
#                 nstart = str(0)
#                 nend = str(chrm_len)
#             else:
#                 # adjust the summit for the chromosome boundaries
#                 nsummit = start+summit
#                 nsummit = max(nsummit, 0 + half_seqlen)
#                 nsummit = min(nsummit, chrm_len - half_seqlen)
#
#                 nstart = str(nsummit-half_seqlen)
#                 nend = str(nsummit+half_seqlen)
#
#             new.write("\t".join([chrm, nstart, nend] + rest) + "\n")
#
#     return bed_out


# 1: set peak width depending on data type
# 2: score the peaks using the bam coverage
# 3: "quantile normalize"
# 4: set output peak widths to 200 bp.


# class Enhancer(object):
#     def __init__(self, genome, bam, peak, output):
#         # package_dir = os.path.dirname(__file__)
#         self.genome = genome
#         self.bam = bam
#         self.peak = peak
#         self.output = output
#         self.outdir = os.path.dirname(self.output)
    #     self.peak_2k = os.path.join(package_dir, "db", "enhancer_peak_2000.bed")
    #     self.peak_rank = os.path.join(package_dir, "db", "peak_rank_hg38_h3k27ac.txt")
    #
    # def quantile_normalize(self, bed_input, quantile_bed):
    #     # TODO: what is going on here! Qnorm a single file with a random other file?
    #
    #     # TODO: this just loads the file to a list
    #     rank = []
    #     with open(self.peak_rank) as p:
    #         for i in p:
    #             rank.append(float(i[:-1]))
    #
    #     # TODO: this replaces the multiBamCov score with the scores from the list (as long as the list is longer)
    #     bed = pd.read_csv(bed_input, header=None, sep="\t")
    #     t = np.searchsorted(np.sort(bed[3]), bed[3])
    #     # bed[1] = [int(i)+900 for i in bed[1].tolist()]
    #     # bed[2] = [int(i)-900 for i in bed[2].tolist()]
    #     bed[3] = [rank[i] for i in t]                   # TODO: replaces the score
    #
    #     #quantile_bed = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
    #     bed.to_csv(quantile_bed, sep="\t", header=False, index=False)
    #     #return quantile_bed.name
    #
    # def run_intersect(self, quantile_bed, intersected_bed):
    #     """
    #     only keep the shrunk peaks that still intersect the original peaks
    #     # TODO: I think that's the purpose?
    #     """
    #     # TODO: switch to BedTool.intersect()
    #     cmd = f"bedtools intersect -a {quantile_bed} -b {self.peak} -wa > {intersected_bed}"
    #     process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
    #     process.wait()
    #     return self.output

    # def run_enhancer(self):
    #     utils.mkdir_p(os.path.join(self.outdir, "intermediate"))
    #     intermediate_bed = os.path.join(self.outdir, "intermediate", "wide_bed4.bed")
    #
    #     bedn2bed4(self.peak, intermediate_bed)
    #     shrink_bed(self.genome, intermediate_bed, self.output)
    #
    #     utils.rm_rf(os.path.join(self.outdir, "intermediate"))


        # cov_bed = os.path.join(self.outdir, "intermediate", "multi_bam_cov.bed")
        # if not os.path.exists(cov_bed) or os.stat(cov_bed).st_size == 0:
        #     logger.info("multi_bam_cov.bed is being generated")
        #     bam_cov(self.bam, self.peak_2k, cov_bed)
        #     logger.info("multi_bam_cov.bed completed")
        #
        # quantile_bed = os.path.join(self.outdir, "intermediate", "quantile_normalized.bed")
        # if not os.path.exists(quantile_bed) or os.stat(quantile_bed).st_size == 0:
        #     logger.info("quantile_normalized.bed is being generated")
        #     self.quantile_normalize(cov_bed, quantile_bed)
        #     logger.info("quantile_normalized.bed completed")
        #
        # intersected_bed = os.path.join(self.outdir, "intermediate", "intersected.bed")
        # if not os.path.exists(intersected_bed) or os.stat(intersected_bed).st_size == 0:
        #     logger.info("intersected.bed is being generated")
        #     self.run_intersect(quantile_bed, intersected_bed)
        #     logger.info("intersected.bed completed")
        #
        # if not os.path.exists(self.output) or os.stat(self.output).st_size == 0:
        #     logger.info("enhancer data is being generated")
        #     shrink_bed(self.genome, intersected_bed, self.output, 200)
        #     logger.info("enhancer data normalized!")


# class P300Enhancer(object):
#     def __init__(self, genome, bam, peak, output):
#         package_dir = os.path.dirname(__file__)
#         self.genome = genome
#         self.bam = bam
#         self.peak = peak
#         self.output = output
#         self.peak_rank = os.path.join(package_dir, "db", "peak_rank.txt")
#
#     def quantileNormalize(self, bed_input):
#         enhancer_number = pd.read_csv(self.peak, sep="\t", header=None).shape[0]
#         rank = pd.read_csv(self.peak_rank, header=None).sample(n=enhancer_number, random_state=1).sort_values(0)[0].tolist()
#
#         bed = pd.read_csv(bed_input, header=None, sep="\t")
#         t = np.searchsorted(np.sort(bed[3]), bed[3])
#         bed[3] = [rank[i] for i in t]
#         # bed[1] = [int(i) for i in bed[1].tolist()]
#         # bed[2] = [int(i) for i in bed[2].tolist()]
#
#         quantile_bed = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         bed.to_csv(quantile_bed, sep="\t", header=False, index=False)
#         return quantile_bed.name
#
#     def run_enhancer(self):
#         # peak200 = set_peak_size(self.genome, self.peak, 200)
#         peak200 = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         shrink_bed(self.genome, self.peak, peak200, 200)
#         bed_cov = bam_cov(self.bam, peak200)
#         quantile_bed = self.quantileNormalize(bed_cov)
#
#         cmd = f"mv {quantile_bed} {self.output}"
#         process = sp.Popen(cmd, shell=True, stdout=sp.PIPE)
#         process.wait()
#
#
# class AtacEnhancer(object):
#     def __init__(self, bam, peak, output, genome="hg38"):
#         package_dir = os.path.dirname(__file__)
#         self.genome = genome
#         self.bam = bam
#         self.peak = peak
#         self.output = output
#         self.peak_rank = os.path.join(package_dir, "db", "peak_rank.txt")
#
#     def quantileNormalize(self, peak, bed_input, bed_output):
#         enahcer_number = pd.read_csv(peak, sep="\t", header=None).shape[0]
#         rank = pd.read_csv(self.peak_rank, header=None).sample(n = enahcer_number, random_state = 1).sort_values(0)[0].tolist()
#
#         bed = pd.read_csv(bed_input, header=None, sep="\t")
#         t = np.searchsorted(np.sort(bed[3]), bed[3])
#         bed[3] = [rank[i] for i in t]
#         # bed[1] = [int(i)+900 for i in bed[1].tolist()]
#         # bed[2] = [int(i)-900 for i in bed[2].tolist()]
#
#         bed.to_csv(bed_output, sep="\t", header=False, index=False)
#         return bed_output
#
#     def run_enhancer(self):
#         peak2000 = set_peak_size(self.genome, self.peak, 2000)
#         bed_cov = bam_cov(self.bam, peak2000)
#         quantile_bed = self.quantileNormalize(self.peak, bed_cov, self.output)
#         shrink_bed(self.genome, quantile_bed, self.output, 200)





"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""






# #!/usr/bin/env python
# 
# # Copyright (c) 2009-2019 Quan Xu <qxuchn@gmail.com>
# #
# # This module is free software. You can redistribute it and/or modify it under
# # the terms of the MIT License, see the file COPYING included with this
# # distribution.
# 
# 
# # TODO: samples with more than 108431 peaks
# # TODO: bam.bai file
# 
# import os
# 
# import numpy as np
# import pandas as pd
# from tempfile import NamedTemporaryFile
# import subprocess
# from pybedtools import BedTool
# from genomepy import Genome
# 
# from ananse import mytmpdir
# import ananse
# 
# class Enhancer(object):
#     def __init__(self, bam_input, epeak, bed_output, genome="hg38"):
#         
#         package_dir = os.path.dirname(ananse.__file__)
# 
#         self.genome = genome
#         self.bam_input = bam_input
#         self.epeak = epeak
#         self.bed_output = bed_output
#         self.peak_2k =  os.path.join(package_dir, "db", "enhancer_peak_2000.bed")
#         self.peak_rank =  os.path.join(package_dir, "db", "peak_rank_hg38_h3k27ac.txt")
# 
#     def runCov(self, bam_input):
#         covfile = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         covcmd = f"multiBamCov -bams {bam_input} -bed {self.peak_2k} > {covfile.name}"
#         process = subprocess.Popen(covcmd, shell=True, stdout=subprocess.PIPE)
#         process.wait()
#         return covfile.name
# 
#     def quantileNormalize(self, bed_input):
#         rank=[]
#         with open(self.peak_rank) as p:
#             for i in p:
#                 rank.append(float(i[:-1]))
# 
#         bed = pd.read_csv(bed_input, header=None, sep="\t")
#         t = np.searchsorted(np.sort(bed[3]), bed[3])
#         bed[3] = [rank[i] for i in t]
#         bed[1] = [int(i)+900 for i in bed[1].tolist()]
#         bed[2] = [int(i)-900 for i in bed[2].tolist()]
#         
#         quantile_bed = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         bed.to_csv(quantile_bed, sep="\t", header=False, index=False)
#         return quantile_bed.name
# 
#     def runIntersect(self, epeak, bed_input, bed_output):
#         intercmd = f"bedtools intersect -a {bed_input} -b {epeak} -wa > {bed_output}"
#         process = subprocess.Popen(intercmd, shell=True, stdout=subprocess.PIPE)
#         process.wait()
# 
#     def run_enhancer(self, bed_input, epeak, bed_output):
#         bed_input = self.runCov(self.bam_input)
#         quantile_bed = self.quantileNormalize(bed_input)
#         self.runIntersect(self.epeak, quantile_bed, bed_output)
#      
# class P300Enhancer(object):
#     def __init__(self, bam_input, epeak, bed_output, genome="hg38"):
#         
#         package_dir = os.path.dirname(ananse.__file__)
# 
#         self.genome = genome
#         self.bam_input = bam_input
#         self.epeak = epeak
#         self.bed_output = bed_output
#         self.peak_rank =  os.path.join(package_dir, "db", "peak_rank.txt")
# 
#     def set_peak_size(self, peak_bed, seqlen=200):
# 
#         """
#         set all input peaks to 200bp
#         Arguments:
#             peak_bed {[bed]} -- [input peak bed file]
#
#         Keyword Arguments:
#             seqlen {int} -- [peak length] (default: {200})
#
#         Returns:
#             [type] -- [200bp peak file]
#         """
#         gsizedic = Genome(self.genome).sizes
#
#         peaks = BedTool(peak_bed)
#         fl2 = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#
#         for peak in peaks:
#             if peak.length < seqlen or peak.length > seqlen:
#                 # get the summit and the flanking low and high sequences
#                 summit = (peak.start + peak.end) // 2
#                 start, end = summit - seqlen // 2, summit + seqlen // 2
#             else:
#                 start, end = peak.start, peak.end
#
#             # remove seq which langer than chromosome length or smaller than 0
#             if start > 0 and end < int(gsizedic[peak.chrom]):
#                 fl2.write(f"{peak.chrom}\t{start}\t{end}\n")
#
#         return fl2.name
#
#     def mk_peak(self, epeak):
#         epeak200 = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         with open(epeak) as peakfile, open(epeak200.name,"w") as npeakfile:
#             for line in peakfile:
#                 a=line.split()
#                 chrm=a[0]
#                 start=int(a[1])
#                 summit=int(a[9])
#                 nsummit=start+summit
#                 if nsummit<100:
#                     nsummit=100
#                 npeakfile.write(f"{chrm}\t{nsummit-100}\t{nsummit+100}\n")
#         return epeak200.name
#
#     def quantileNormalize(self, epeak, bed_input, bed_output):
#         enahcer_number = pd.read_csv(epeak, sep="\t", header=None).shape[0]
#         rank = pd.read_csv(self.peak_rank, header=None).sample(n = enahcer_number, random_state = 1).sort_values(0)[0].tolist()
#
#         bed = pd.read_csv(bed_input, header=None, sep="\t")
#         t = np.searchsorted(np.sort(bed[3]), bed[3])
#         bed[3] = [rank[i] for i in t]
#         bed[1] = [int(i) for i in bed[1].tolist()]
#         bed[2] = [int(i) for i in bed[2].tolist()]
#
#         bed.to_csv(bed_output, sep="\t", header=False, index=False)
#
#     def runCov(self, bam_input, clear_epeak200):
#         covfile = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         covcmd = f"multiBamCov -bams {bam_input} -bed {clear_epeak200} > {covfile.name}"
#         process = subprocess.Popen(covcmd, shell=True, stdout=subprocess.PIPE)
#         process.wait()
#         return covfile.name
#
#     def run_enhancer(self, bed_input, epeak, bed_output):
#         epeak200 = self.mk_peak(epeak)
#         clear_epeak200 = self.set_peak_size(epeak200)
#         # os.system(f"cp {clear_epeak200} ./")
#         bed_cov = self.runCov(self.bam_input, clear_epeak200)
#         quantile_bed = self.quantileNormalize(epeak, bed_cov, bed_output)
#
#
# class AtacEnhancer(object):
#     def __init__(self, bam_input, epeak, bed_output, genome="hg38"):
#
#         package_dir = os.path.dirname(ananse.__file__)
#
#         self.genome = genome
#         self.bam_input = bam_input
#         self.epeak = epeak
#         self.bed_output = bed_output
#         self.peak_rank =  os.path.join(package_dir, "db", "peak_rank.txt")
#
#     def set_peak_size(self, peak_bed, seqlen=200):
#         """set all input peaks to 200bp
#
#         Arguments:
#             peak_bed {[bed]} -- [input peak bed file]
#
#         Keyword Arguments:
#             seqlen {int} -- [peak length] (default: {200})
#
#         Returns:
#             [type] -- [200bp peak file]
#         """
#         gsizedic = Genome(self.genome).sizes
#
#         peaks = BedTool(peak_bed)
#         fl2 = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#
#         for peak in peaks:
#
#             if peak.length < seqlen or peak.length > seqlen:
#                 # get the summit and the flanking low and high sequences
#                 summit = (peak.start + peak.end) // 2
#                 start, end = summit - seqlen // 2, summit + seqlen // 2
#             else:
#                 start, end = peak.start, peak.end
#             # remove seq which langer than chromosome length or smaller than 0
#             if start > 0 and end < int(gsizedic[peak.chrom]):
#                 fl2.write(f"{peak.chrom}\t{start}\t{end}\n")
#
#         # return npeaks
#         return fl2.name
#
#     def mk_peak(self, epeak):
#         epeak200 = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         with open(epeak) as peakfile, open(epeak200.name,"w") as npeakfile:
#             for line in peakfile:
#                 a=line.split()
#                 chrm=a[0]
#                 start=int(a[1])
#                 summit=int(a[9])
#                 nsummit=start+summit
#                 if nsummit<100:
#                     nsummit=100
#                 npeakfile.write(f"{chrm}\t{nsummit-100}\t{nsummit+100}\n")
#         return epeak200.name
#
#     def quantileNormalize(self, epeak, bed_input, bed_output):
#         enahcer_number = pd.read_csv(epeak, sep="\t", header=None).shape[0]
#         rank = pd.read_csv(self.peak_rank, header=None).sample(n = enahcer_number, random_state = 1).sort_values(0)[0].tolist()
#
#         bed = pd.read_csv(bed_input, header=None, sep="\t")
#         t = np.searchsorted(np.sort(bed[3]), bed[3])
#         bed[3] = [rank[i] for i in t]
#         bed[1] = [int(i)+900 for i in bed[1].tolist()]
#         bed[2] = [int(i)-900 for i in bed[2].tolist()]
#
#         bed.to_csv(bed_output, sep="\t", header=False, index=False)
#
#     def runCov(self, bam_input, clear_epeak2k):
#         covfile = NamedTemporaryFile(mode="w", dir=mytmpdir(), delete=False)
#         covcmd = f"multiBamCov -bams {bam_input} -bed {clear_epeak2k} > {covfile.name}"
#         process = subprocess.Popen(covcmd, shell=True, stdout=subprocess.PIPE)
#         process.wait()
#         return covfile.name
#
#     def run_enhancer(self, bed_input, epeak, bed_output):
#         epeak200 = self.mk_peak(epeak)
#         clear_epeak2k = self.set_peak_size(epeak200, 2000)
#         bed_cov = self.runCov(self.bam_input, clear_epeak2k)
#         quantile_bed = self.quantileNormalize(epeak, bed_cov, bed_output)
#
#
