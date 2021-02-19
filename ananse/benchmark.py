import urllib
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    roc_auc_score,
)
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re
import sys
import os

from loguru import logger

logger.remove()
logger.add(
    sys.stderr, format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | {level} | {message}"
)

def download_trrust_reference(outfile):
    """Download TRRUST GRN reference.

    Parameters
    ----------
    outfile : str
        Name of output file.
    """
    TRRUST_URL = "https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv"

    edges = []
    with urllib.request.urlopen(
        TRRUST_URL,
    ) as f:
        for line in f.readlines():
            tf, target, regtype, pmid = line.decode().strip().split("\t")
            # Just skip repression for now
            if regtype in ["Activation", "Unknown"]:
                edges.append([tf, target, 1])
    edges = pd.DataFrame(edges, columns=["tf", "target", "interaction"])
    edges.to_csv(outfile, sep="\t", index=False)


def fix_columns(df):
    """Make sure network has a tf and a target column."""

    df.columns = df.columns.str.lower()
    df = df.rename(
        columns={
            "source": "tf",
            "source_target": "tf_target",
            "target_gene": "target",
        }
    )

    if "tf_target" in df.columns:
        df[["tf", "target"]] = df["tf_target"].str.split("_", expand=True).iloc[:, :2]
        df = df.drop(columns=["tf_target"])

    if not "tf" in df.columns:
        raise ValueError("Expect a column named 'source' or 'tf'")
    if not "target" in df.columns:
        raise ValueError("Expect a column named 'target' or 'target_gene'")
    return df

def get_tfs():
    valid_factors = pd.read_excel(
        "https://www.biorxiv.org/content/biorxiv/early/2020/12/07/2020.10.28.359232/DC1/embed/media-1.xlsx",
        engine="openpyxl",
        sheet_name=1,
    )
    valid_factors = valid_factors.loc[
        valid_factors["Pseudogene"].isnull(), "HGNC approved gene symbol"
    ].values
    valid_factors = [f for f in valid_factors if f != "EP300"]
    return valid_factors


def prepare_reference_network(network, filter_tfs=True):
    """Generate reference network.

    This network contains all possible edges, based on the TFs
    and the target genes in the input. TFs are optionally filtered
    to contain only validated TFs.

    Returns
    -------
        DataFrame with column `"interaction"` having 1 for a validated
        edge and 0 otherwise.
    """
    if isinstance(network, pd.DataFrame):
        df = network.reset_index()
    elif isinstance(network, str):
        if network.endswith("feather"):
            df = pd.read_feather(network)
        else:
            df = pd.read_table(network)
    else:
        raise ValueError("Unknown network type, need DataFrame or filename.")

    df = fix_columns(df)

    interaction_column = None
    for col in df.columns:
        if col in ["tf", "target"]:
            continue
        vals = df[col].unique()

        if len(vals) in [1, 2] and 1 in vals:
            interaction_column = col
            break

    tfs = set(df["tf"].unique())
    if filter_tfs:
        valid_tfs = set(get_tfs())
        tfs = list(tfs.intersection(valid_tfs))

    targets = df["target"].unique()

    logger.info(
        f"{os.path.split(network)[-1]} reference - {len(tfs)} TFs, {len(targets)} targets, {df.shape[0]} edges."
    )

    total = []
    for tf in tfs:
        for target in targets:
            total.append([tf, target])
    total = pd.DataFrame(total, columns=["tf", "target"]).set_index(["tf", "target"])

    if interaction_column is not None:
        logger.info(f"Using '{interaction_column}' as interaction column.")
        df = df.set_index(["tf", "target"])[[interaction_column]].rename(
            columns={interaction_column: "interaction"}
        )
    else:
        logger.info("No column with 1 found, assuming all lines are positive edges.")
        df = df.set_index(["tf", "target"])
        df["interaction"] = 1

    return total.join(df[["interaction"]]).fillna(0)

def read_network(fname, name=None):
    """Read a network file.

    Will make sure the resulting DataFrame has a "tf", "factor" MultiIndex.
    Will rename and reprocess columns if necessary.

    Parameters
    ----------
    fname : str
        Filename.
    name : str, optional
        Name for interaction column.
    
    Returns
    -------
    pandas.DataFrame
        DataFrame containing edges.
    """
    network = fname
    if fname.endswith("feather"):
        df = pd.read_feather(network)
    else:
        df = pd.read_table(network)
    df = fix_columns(df)
    df = df.set_index(["tf", "target"])

    # Assuming last column is the edge weight
    df = df.iloc[:, [-1]]
    if name is not None:
        df.columns = [name]

    return df

def _read_dorothea_reference():
    dorothea = pd.read_table(
        "/ceph/rimlsfnwi/data/moldevbio/share_moldevbio/heeringen/ananse/benchmark/dorothea/entire_database.txt"
    )
    cols = [
        "is_evidence_chip_seq",
        "is_evidence_curated",
        "is_evidence_inferred",
        "is_evidence_tfbs",
    ]
    dorothea = dorothea.set_index(["tf", "target"])[cols]
    for col in cols:
        dorothea[col] = dorothea[col].astype(int)
    dorothea["dorothea"] = np.any(dorothea[cols] == 1, 1).astype(int)

    dorothea = dorothea.reset_index()

    tfs = set(dorothea["tf"].unique())
    valid_tfs = set(get_tfs())
    tfs = list(tfs.intersection(valid_tfs))
    targets = dorothea["target"].unique()
    logger.info(
        f"Dorothea reference - {len(tfs)} TFs, {len(targets)} targets, {dorothea.shape[0]} edges."
    )

    total = []
    for tf in tfs:
        for target in targets:
            total.append([tf, target])
    total = pd.DataFrame(total, columns=["tf", "target"]).set_index(["tf", "target"])

    dorothea = dorothea.set_index(["tf", "target"])
    dorothea = total.join(dorothea)
    dorothea = dorothea.fillna(0)
    return dorothea


def _read_enrichr_perturbation_reference():
    """Uses the TF perturbations from Enrichr[1,2] to create reference edges.

    Targets are defined by up- or down-regulated gened from the following sets:
    Up: INDUCTION, ACTIVATION, OE.
    Down: KD, KO, INACTIVATION, DEPLETION, SIRNA, SHRNA, KNOCKOUT, DELETION INHIBITION.

    The TF and targets in the DataFrame consists of the Cartesian product of
    all TFs and target genes that occur in the set.

    Returns
    -------
        DataFrame with tf-target edges.

    References
    ----------
    .. [1] Chen EY, Tan CM, Kou Y, Duan Q, Wang Z, Meirelles GV, Clark NR, Ma'ayan A.
       "Enrichr: interactive and collaborative HTML5 gene list enrichment analysis
       tool." BMC Bioinformatics. 2013;128(14)

    .. [2] Kuleshov MV, Jones MR, Rouillard AD, Fernandez NF, Duan Q, Wang Z,
       Koplev S, Jenkins SL, Jagodnik KM, Lachmann A, McDermott MG, Monteiro CD,
       Gundersen GW, Ma'ayan A. "Enrichr: a comprehensive gene set enrichment
       analysis web server 2016 update." Nucleic Acids Research. 2016; gkw377.
    """
    infile = "/ceph/rimlsfnwi/data/moldevbio/share_moldevbio/heeringen/ananse/benchmark/enrichr/TF_Perturbations_Followed_by_Expression.txt"
    p = re.compile("(\w+)\s+(\w+)\s+(.+)\s+(\w+)")
    all_info = []
    edges = []
    with open(infile) as f:
        for line in f:
            vals = line.strip().split("\t")
            m = re.search(p, vals[0])
            # print(m.groups())
            all_info.append(m.groups(0))
            if (
                m.group(2) in ["INDUCTION", "ACTIVATION", "OE"] and m.group(4) == "UP"
            ) or (
                m.group(2)
                in [
                    "KD",
                    "KO",
                    "INACTIVATION",
                    "DEPLETION",
                    "SIRNA",
                    "SHRNA",
                    "KNOCKOUT",
                    "DELETION",
                    "INHIBITION",
                ]
                and m.group(4) == "DOWN"
            ):
                tf = m.group(1)
                for target in vals[2:]:
                    edges.append([tf, target])
    all_info = pd.DataFrame(all_info, columns=["tf", "exp", "info", "up_down"])

    perturb_df = pd.DataFrame(edges, columns=["tf", "target"])
    tfs = set(perturb_df["tf"].unique())
    targets = perturb_df["target"].unique()

    logger.info(
        f"TF perturbation reference - {len(tfs)} TFs, {len(targets)} targets, {perturb_df.shape[0]} edges."
    )

    perturb_df["experiments"] = 1
    perturb_df = perturb_df.groupby(["tf", "target"]).count()
    perturb_df["interaction"] = 1
    perturb_df.columns = ["perturb_experiments", "perturb_interaction"]

    valid_tfs = set(get_tfs())
    tfs = list(tfs.intersection(valid_tfs))

    total = []
    for tf in tfs:
        for target in targets:
            total.append([tf, target])
    total = pd.DataFrame(total, columns=["tf", "target"]).set_index(["tf", "target"])

    perturb_df = total.join(perturb_df).fillna(0)
    return perturb_df

def read_reference(fname):
    if fname.lower() == "dorothea":
        return _read_dorothea_reference()
    if fname.lower() == "perturbation":
        return _read_enrichr_perturbation_reference()
    return prepare_reference_network(fname)

def _validate_files(fnames):
    file_error = False
    for fname in fnames:
        if not os.path.exists(fname):
            print(f"file {fname} does not exist")
            file_error = True
    if file_error:
        raise ValueError("One or more files not found!")


def read_networks(network_dict):
    """Read predicted networks.

    Input is a dictionary with name as key and filename as value.
    """
    # Validate files first
    _validate_files(network_dict.values())

    df = pd.DataFrame({"tf": [], "target": []}).set_index(["tf", "target"])

    for name, fname in network_dict.items():
        logger.info(f"Reading {name}")
        tmp = read_network(fname, name=name)
        logger.info(f"Merging {name}")
        df = df.join(tmp, how="outer")

    return df
