from collections import defaultdict
# Suppress ete3's invalid escape sequence warnings
import warnings
import re
warnings.filterwarnings(
    "ignore",
    category=SyntaxWarning,
    module="ete3"
)
from ete3 import Tree
from operator import itemgetter
from typing import Dict, List, Union, Any, TypedDict, Optional, Set, Tuple
import math
import multiprocessing
import os
import glob
import pickle
import resource
import shutil
import statistics
import subprocess
import sys
import traceback
import zlib
from Bio import SeqIO
from hmmlearn.base import BaseHMM
from scipy.stats import pearsonr
os.environ['KMP_WARNINGS'] = 'off'
import numpy
import pandas as pd
import tqdm


from zol import data_dictionary, util, zol

#: Number of target genomes processed together in identify_gc_instances. DIAMOND hit summaries are written as one
#: (zstd-compressed) file per block of this many genomes (sorted by name), so each block reads only its own hits.
SAMPLE_BLOCK_SIZE = 1000
#: Upper bound on files kept open at once when sharding DIAMOND hits or writing per-block hit summaries.
MAX_OPEN_FILES = 128
#: DIAMOND blastp output columns used by fai's homology search.
FAI_DIAMOND_COLUMNS = ["qseqid", "sseqid", "pident", "evalue", "bitscore", "qlen", "slen"]
DIAMOND_SUMMARY_HEADER = "Locus Tag\tHomolog Group\tSample\tBest Bitscore\tE-value\tIdentity\tQuery Length\tSubject Length\n"


def subset_genbank_for_query_locus(
    full_gw_genbank,
    locus_genbank,
    locus_proteins,
    reference_contig,
    reference_start,
    reference_end,
    log_object,
) -> None:
    """
    Description:
    This function extracts a locus specific GenBank from a reference genome based on user provided coordinates.
    ********************************************************************************************************************
    Parameters:
    - full_gw_genbank: Reference genome GenBank file that should have CDS features.
    - locus_genbank: The locus-specific GenBank file to create.
    - locus_proteins: The locus-specific protein FASTA file to create.
    - reference_contig: The identifier of the reference scaffold / contig with the locus.
    - reference_start: The starting position of the locus in the reference genome.
    - reference_end: The ending position of the locus in the reference genome.
    - log_object: A logging object.
    ********************************************************************************************************************
    """
    try:
        util.create_genbank(
            full_gw_genbank,
            locus_genbank,
            reference_contig,
            reference_start,
            reference_end,
        )

        with open(locus_proteins, "w") as locus_proteins_handle:
            with open(locus_genbank) as olg:
                for rec in SeqIO.parse(olg, "genbank"):
                    for feature in rec.features:
                        if feature.type != "CDS":
                            continue
                        lt = feature.qualifiers.get("locus_tag")[0]
                        prot_seq = feature.qualifiers.get("translation")[0]
                        locus_proteins_handle.write(
                            f">{lt}\n" + str(prot_seq) + "\n"
                        )

    except Exception as e:
        sys.stderr.write(
            "Issues subsetting locus-specific GenBank from genome-wide GenBank.\n"
        )
        log_object.error(
            "Issues subsetting locus-specific GenBank from genome-wide GenBank."
        )
        util.report_cause(fallback_traceback=True)
        sys.exit(1)


def parse_homolog_group_matrix(orthogroup_matrix_file, log_object) -> Dict[str, Any]:
    """
    Description:
    This function parses an ortholog group matrix and returns a dictionary mapping protein names / locus_tags to ortholog
    groups.
    ********************************************************************************************************************
    Parameters:
    - orthogroup_matrix_file: The ortholog group vs sample matrix file, where cells correspond to locus tag identifiers.
    - log_object: A logging object.
    ********************************************************************************************************************
    Return:
    - protein_to_hg: A dictionary mapping proteins names / locus_tags to ortholog groups.
    ********************************************************************************************************************
    """
    protein_to_hg: Dict[str, Any] = {}
    try:
        with open(orthogroup_matrix_file) as oog:
            for i, line in enumerate(oog):
                line = line.strip()
                ls = line.split("\t")
                if i == 0:
                    continue
                hg = ls[0]
                for lts in ls[1:]:
                    for lt in lts.split(", "):
                        protein_to_hg[lt] = hg
    except Exception as e:
        sys.stderr.write("Issues mapping proteins to ortholog groups.\n")
        log_object.error("Issues mapping proteins to ortholog groups.")
        util.report_cause(fallback_traceback=True)
        sys.exit(1)
    return protein_to_hg


def collapse_proteins_using_diamond_cluster(protein_fasta, nr_protein_fasta, log_object, diamond_params="--approx-id 95.0 --mutual-cover 90.0", threads=1, memory=4) -> Dict[str, Any]:
    """
    Description:
    This function collapses query proteins using diamond cluster. To ensure that criteria are met 
    diamond recluster is also applied. 
    ********************************************************************************************************************
    Parameters:
    - protein_fasta: A file with the original set of protein queries in FASTA form.
    - nr_protein_fasta: The path to the non-redundant (nr) protein queries after diamond cluster representative selection 
                        in FASTA format.
    - log_object: A logging object.
    - diamond_params: The parameters to use for diamond cluster. Default is "--approx-id 95.0 --mutual-cover 90.0".
    - threads: The number of threads to use for diamond cluster. Default is 1.
    - memory: The maximum memory to use in GB.
    ********************************************************************************************************************
    Return:
    - protein_to_hg: A dictionary mapping proteins names / locus_tags to representative proteins.
    ********************************************************************************************************************
    """
    protein_to_hg: Dict[str, Any] = {}
    try:
        diamond_cluster_file = nr_protein_fasta + "_clusters.tsv"
        diamond_refined_cluster_file = nr_protein_fasta + "_refined_clusters.tsv"

        util.diamond_cluster_and_recluster(protein_fasta, diamond_cluster_file, diamond_refined_cluster_file, 
                                           log_object, threads=threads, diamond_params=diamond_params, 
                                           cluster_method="cluster", skip_recluster=False, mem=memory)

        prefix = ".".join(protein_fasta.split("/")[-1].split(".")[:-1])
        
        # Parse diamond clustering output
        reps = set([])
        with open(diamond_refined_cluster_file) as occf:
            for line in occf:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    continue
                rep_id = parts[0]
                mem_id = parts[1]
                
                protein_to_hg[prefix + "|" + mem_id] = rep_id
                reps.add(rep_id)

        # Create representative fasta file to use as query
        # and add self-mappings for cluster representatives
        with open(nr_protein_fasta, "w") as outf_handle:
            with open(protein_fasta) as fasta_handle:
                for record in SeqIO.parse(fasta_handle, "fasta"):
                    protein_id = record.id
                    if protein_id not in reps:
                        protein_to_hg[prefix + "|" + protein_id] = protein_id                
                        outf_handle.write(f">{protein_id}\n{str(record.seq)}\n")
                    else:
                        outf_handle.write(f">{protein_id}\n{str(record.seq)}\n")

    except Exception as e:
        msg = f"Issues running diamond cluster or recluster."
        log_object.error(msg)
        sys.stderr.write(msg + "\n")
        util.report_cause(fallback_traceback=True)
        sys.exit(1)
    return protein_to_hg


def parse_coords_from_genbank(genbanks, log_object) -> Dict[str, Dict[str, Any]]:
    """
    Description:
    This function gathers the coordinates for CDS features from multiple GenBank files.
    ********************************************************************************************************************
    Parameters:
    - genbanks: A set of GenBank files.
    - log_object: A logging object.
    ********************************************************************************************************************
    Return:
    - comp_gene_info: A 3 tier dicitonary with the first key corresponding to the gene cluster prefix / name, the second
                      key corresponding to the locus tag identifier, and the third is information pertaining to CDS
                      locations (scaffold, start position, end position).
    ********************************************************************************************************************
    """
    # get coords for individual genes
    comp_gene_info = defaultdict(dict)
    try:
        for gbk in genbanks:
            prefix = ".".join(gbk.split("/")[-1].split(".")[:-1])
            gene_locations = util.parse_gbk(
                gbk, prefix, log_object, use_either_lt_or_pi=True
            )
            # Type assertion to help the type checker
            assert gene_locations is not None
            for g in gene_locations:  # type: ignore
                comp_gene_info[prefix][g] = gene_locations[g]  # type: ignore
    except Exception as e:
        sys.stderr.write(
            "Issues with parsing coordinates of genes in GenBanks!\n"
        )
        log_object.error(
            "Issues with parsing coordinates of genes in GenBanks!"
        )
        util.report_cause(fallback_traceback=True)
        sys.exit(1)
    return dict(comp_gene_info)


def gen_consensus_sequences(
    genbanks, outdir, log_object, threads=1, aligner='famsa') -> List[str]:
    """
    Description:
    This function performs zol-based ortholog group determination across CDS protein sequences from query GenBank files
    provided and then determines consensus protein sequences for each ortholog group identified.
    ********************************************************************************************************************
    Parameters:
    - genbanks: A set of query GenBank files.
    - outdir: The output directory / workspace where to write intermediate files / analyses.
    - log_object: A logging object.
    - threads: The number of threads to use.
    - aligner: The aligner to use for MSA. Options are 'famsa' (default), 'muscle', or 'muscle-super5'.
    ********************************************************************************************************************
    Return:
    - orthogroup_matrix_file: An ortholog group vs sample matrix file, where cells correspond to locus tag identifiers.
    - consensus_prot_seqs_faa: Path to a FASTA file containing consensus protein sequences for ortholog groups
                               determined.
    ********************************************************************************************************************
    """
    # determine orthologs
    outdir = os.path.abspath(outdir) + '/'
    prot_dir = outdir + "CDS_Protein/"
    nucl_dir = outdir + "CDS_Nucleotide/"
    og_dir = outdir + "Determine_Orthogroups/"
    ortho_matrix_file = og_dir + "Orthogroups.tsv"
    util.setup_ready_directory([prot_dir, nucl_dir])
    try:
        for gbk in genbanks:
            try:
                prefix = ".".join(gbk.split("/")[-1].split(".")[:-1])
                proteins, nucleotides, upstream_regions = (
                    util.parse_genbank_for_cds_proteins_and_dna(gbk, log_object, use_either_lt_or_pi=True)  # type: ignore
                )
                protein_outf = prot_dir + f"{prefix}.faa"
                with open(protein_outf, "w") as protein_handle:
                    for lt in proteins:
                        protein_handle.write(f">{prefix}|{lt}\n{proteins[lt]}\n")

            except Exception as e:
                msg = f"Issues with parsing the GenBank {gbk}\n"
                sys.stderr.write(msg + "\n")
                util.report_cause(fallback_traceback=True)
                log_object.error(msg)
                sys.exit(1)

        fo_cmd = [
            "findOrthologs.py",
            "-p",
            prot_dir,
            "-o",
            og_dir,
            "-c",
            str(threads)]
        util.run_cmd_via_subprocess(fo_cmd, log_object=log_object,
                                    check_files = [ortho_matrix_file])

        if not (util.check_core_homolog_groups_exist(ortho_matrix_file)):
            log_object.error(f"Issue validating core ortholog group(s) exist: {ortho_matrix_file}")
            sys.stderr.write(traceback.format_exc())
            sys.exit(1)

    except Exception as e:
        sys.stderr.write("Issues with determining ortholog / ortholog groups!\n")
        util.report_cause(fallback_traceback=True)
        log_object.error("Issues with determining ortholog / ortholog groups!")
        sys.exit(1)

    proc_dir = outdir + "Homolog_Group_Processing/"
    hg_prot_dir = proc_dir + "OG_Protein_Sequences/"
    hg_nucl_dir = proc_dir + "OG_Nucleotide_Sequences/"
    prot_algn_dir = proc_dir + "OG_Protein_Alignments/"
    phmm_dir = proc_dir + "OG_Profile_HMMs/"
    cons_dir = proc_dir + "OG_Consensus_Sequences/"
    consensus_prot_seqs_faa = outdir + "OG_Consensus_Seqs.faa"
    util.setup_ready_directory(
        [proc_dir, prot_algn_dir, phmm_dir, cons_dir, hg_prot_dir]
    )
    zol.partition_sequences_by_homolog_groups(
        ortho_matrix_file,
        prot_dir,
        nucl_dir,
        hg_prot_dir,
        hg_nucl_dir,
        log_object,
    )
    zol.create_protein_alignments(
        hg_prot_dir,
        prot_algn_dir,
        log_object,
        aligner=aligner,
        threads=threads,
    )
    zol.create_profile_hmms_and_consensus_seqs(
        prot_algn_dir, phmm_dir, cons_dir, log_object, threads=1
    )
    with open(consensus_prot_seqs_faa, "w") as consensus_prot_seqs_handle:
        for f in os.listdir(cons_dir):
            with open(cons_dir + f) as ocf:
                for rec in SeqIO.parse(ocf, "fasta"):
                    consensus_prot_seqs_handle.write(f">{f.split('.cons.faa')[0]}\n{rec.seq}\n")

    return [ortho_matrix_file, consensus_prot_seqs_faa]


def load_target_genome_info(
    target_annotation_information,
    target_genomes_pkl_dir,
    diamond_results,
    valid_tg_samples,
    log_object,
    min_genes_per_scaffold=1,
    min_distinct_genes_per_scaffold=1,
) -> Dict[str, Any]:
    """
    Description:
    Load information pertaining to CDS locations from target genomes from sample pickle files created during prepTG.
    ********************************************************************************************************************
    Parameters:
    - target_annotation_information: A set / dictionary of target genome sample identifiers.
    - target_genomes_pkl_dir: The directory where pickle files on target genome CDS location information are stored.
    - diamond_results: The results of the DIAMOND blastp search.
    - valid_tg_samples: The valid target genome samples.
    - log_object: A logging object.
    - lowmem_mode: Whether to only regard CDS found on scaffolds with hits from the gene cluster queries.
    - min_genes_per_scaffold: The minimum number of ORFs hit by query proteins for the scaffold to be considered and
                              CDS information for it stored.
    - min_distinct_genes_per_scaffold: The minimum number of distinct genes per scaffold to be considered for a 
                                       gene-cluster segment.
    ********************************************************************************************************************
    Return:
    - target_genome_info: A dictionary of dictionaries:
            - gene_locations: A multi - tiered dictionary with location information for CDS features.
            - scaffold_genes: A double dictionary with keys corresponding to samples and secondary keys to scaffolds and
                              values to CDS features found on them.
            - boundary_genes: A dictionary with keys corresponding to samples and values to sets of CDS features nearby
                              scaffold edges.
            - gene_id_to_order: A multi-tiered dictionary mapping CDS identifiers to their order along scaffolds.
            - gene_order_to_id: A multi-tiered dictionary mapping gene order indices along scaffolds to CDS identifiers.
    ********************************************************************************************************************
    """
    gene_locations: Dict[str, Any] = {}
    scaffold_genes: Dict[str, Any] = {}
    boundary_genes: Dict[str, Any] = {}
    gene_id_to_order: Dict[str, Any] = {}
    gene_order_to_id: Dict[str, Any] = {}

    try:
        for sample in target_annotation_information:
            if not sample in valid_tg_samples or not sample in diamond_results or \
                len(diamond_results[sample]['lts']) < min_genes_per_scaffold or \
                len(diamond_results[sample]['hgs']) < min_distinct_genes_per_scaffold:
                continue

            sample_info_pkl_file = target_genomes_pkl_dir + sample + ".pkl"
            genbank_info = None
            try:
                if not (os.path.isfile(sample_info_pkl_file)):
                    raise AssertionError()
                with open(sample_info_pkl_file, "rb") as handle:
                    genbank_info = pickle.load(handle)
            except Exception as e:
                sys.stderr.write(
                    f"Could not find pickle file with CDS coordinate information for sample {sample}\n"
                )
                log_object.error(
                    f"Could not find pickle file with CDS coordinate information for sample {sample}"
                )
                util.report_cause(fallback_traceback=True)
                sys.exit(1)

            gene_locs = genbank_info[0]
            scaff_genes = genbank_info[1]
            bound_genes = genbank_info[2]
            gito = genbank_info[3]
            goti = genbank_info[4]

            valid_scaffolds = set([])
            gene_locs_filt: Dict[str, Any] = {}
            scaff_genes_filt = defaultdict(set)
            bound_genes_filt = set([])
            gito_filt = defaultdict(dict)
            goti_filt = defaultdict(dict)
            for scaffold in scaff_genes:
                if (
                    len(
                        scaff_genes[scaffold].intersection(
                            diamond_results[sample]
                        )
                    )
                    >= min_genes_per_scaffold
                ):
                    valid_scaffolds.add(scaffold)
                    scaff_genes_filt[scaffold] = scaff_genes[scaffold]
                    gito_filt = gito[scaffold]
                    goti_filt = goti[scaffold]
                    for lt in scaff_genes[scaffold]:
                        gene_locs_filt[lt] = gene_locs[lt]
                    bound_genes_filt = bound_genes_filt.union(
                        scaff_genes[scaffold].intersection(bound_genes)
                    )
            gene_locations[sample] = gene_locs_filt
            scaffold_genes[sample] = scaff_genes_filt
            boundary_genes[sample] = bound_genes_filt
            gene_id_to_order[sample] = gito_filt
            gene_order_to_id[sample] = goti_filt

    except Exception as e:
        log_object.error(
            "Issues with parsing CDS location information from genes of target genomes."
        )
        sys.stderr.write(
            "Issues with parsing CDS location information from genes of target genomes.\n"
        )
        util.report_cause(fallback_traceback=True)
        sys.exit(1)

    target_genome_gene_info = {
        "gene_locations": gene_locations,
        "scaffold_genes": scaffold_genes,
        "boundary_genes": boundary_genes,
        "gene_id_to_order": gene_id_to_order,
        "gene_order_to_id": gene_order_to_id,
    }
    return target_genome_gene_info


def run_diamond_blastp(
    target_concat_genome_db,
    query_fasta,
    diamond_work_dir,
    log_object,
    diamond_sensitivity="very-sensitive",
    diamond_block_size=None,
    evalue_cutoff=1e-3,
    threads=1,
    k_hits="0",
    compute_query_coverage=False,
) -> None:
    """
    Description:
    This functions runs DIAMOND blastp analysis in fai for assessing homology of target genome proteins to query
    proteins.
    ********************************************************************************************************************
    Parameters:
    - target_concat_genome_db: The path to the concatenated FASTA file comprising proteomes from all target genomes.
    - query_fasta: The path to the query / reference FASTA file.
    - diamond_work_dir: The workspace directory where DIAMOND blastp alignment and intermediate files should be written.
    - log_object: A logging object.
    - diamond_sensitivity: The sensitivity mode to run DIAMOND blastp with.
    - evalue_cutoff: The maximum E - value cutoff to regard an alignment to a target genome protein as homologous to a
                     query protein.
    - threads: The number of threads to use.
    - k_hits: The number of hits to return.
    - compute_query_coverage: Whether to compute the query coverage - used for the simple BLASTp approach of abon,
                              atpoc, and apos.
    ********************************************************************************************************************
    Return:
    - diamond_results_file: Path to the DIAMOND blastp resulting alignment file.
    ********************************************************************************************************************
    """

    diamond_results_file = diamond_work_dir + "DIAMOND_Results.txt"
    try:
        diamond_blastp_cmd = [
            "diamond",
            "blastp",
            "--ignore-warnings",
            "--threads",
            str(threads),
            "--query",
            query_fasta,
            "--db",
            target_concat_genome_db,
            "--outfmt",
            "6",
            "qseqid",
            "sseqid",
            "pident",
            "length",
            "mismatch",
            "gapopen",
            "qstart",
            "qend",
            "sstart",
            "send",
            "evalue",
            "bitscore",
            "qlen",
            "slen", 
            "-k" + k_hits,
            "--out",
            diamond_results_file,
            "--evalue",
            str(evalue_cutoff),
        ]
        if compute_query_coverage:
            diamond_blastp_cmd = ["diamond", "blastp", "--ignore-warnings", "--threads",
                                  str(threads), "--query", query_fasta, "--db", target_concat_genome_db,
                                  "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore",
                                  "qlen", "slen", "qcovhsp", "-k" + k_hits, "--out", diamond_results_file,
                                  "--evalue", str(evalue_cutoff)]

        if diamond_sensitivity != "default":
            diamond_blastp_cmd.append("--" + diamond_sensitivity)

        if diamond_block_size is not None:
            diamond_blastp_cmd.append("--block-size")
            diamond_blastp_cmd.append(str(diamond_block_size))

        util.run_cmd_via_subprocess(diamond_blastp_cmd, log_object=log_object, 
                                    check_files = [diamond_results_file])

    except Exception as e:
        msg = "Issues with running DIAMOND blastp."
        log_object.error(msg)
        sys.stderr.write(msg + "\n")
        util.report_cause(fallback_traceback=True)
        sys.exit(1)

    return diamond_results_file

# TypedDict for kq_top_hits structure
class KqTopHitInfo(TypedDict):
    hg_list: List[str]
    best_bitscore: float

def parse_tg_faa_file(tg_faa_file, log_object):
    """
    Description:
    This function parses the target genome FASTA file and returns a dictionary of cluster information.
    ********************************************************************************************************************
    Parameters:
    - tg_faa_file: The path to the target genome FASTA file.
    - log_object: A logging object.
    ********************************************************************************************************************
    Return:
    - rep_prot_to_nonreps: A dictionary of members in cluster with key corresponding to the representative protein.
    ********************************************************************************************************************
    """
    try:
        rep_prot_to_nonreps = defaultdict(set)
        with open(tg_faa_file) as occf:
            for rec in SeqIO.parse(occf, "fasta"):
                cluster_rep = rec.id
                rep_prot_to_nonreps[cluster_rep] = set([cluster_rep])
        return rep_prot_to_nonreps 
    except Exception as e:
        msg = f"Issues parsing target genome FASTA file:\n{tg_faa_file}.\n"
        sys.stderr.write(msg + "\n")
        log_object.info(msg)
        util.report_cause()
        sys.exit(1)

def sample_blocks(samples, block_size=SAMPLE_BLOCK_SIZE) -> List[List[str]]:
    """ Split target genome sample names (sorted) into the blocks processed together by identify_gc_instances. """
    return [list(block) for block in util.chunks(sorted(samples), block_size)]


def block_summary_file(summary_dir, block_index) -> str:
    """ Path of the (zstd-compressed) DIAMOND hit summary for one block of target genomes. """
    return f"{summary_dir}Block_{block_index}.tsv.zst"


def max_open_files() -> int:
    """
    Number of files fai may keep open at once when sharding DIAMOND hits or writing block summaries: at most
    MAX_OPEN_FILES, and comfortably below the process's open-file limit (e.g. 256 by default on macOS).
    """
    soft_limit = 256
    try:
        soft_limit, _ = resource.getrlimit(resource.RLIMIT_NOFILE)
        if soft_limit == resource.RLIM_INFINITY or soft_limit <= 0:
            soft_limit = 1024
    except Exception:
        pass
    return max(1, min(MAX_OPEN_FILES, (soft_limit - 64) // 2))


def run_diamond_blastp_into_shards(
    target_concat_genome_db,
    query_fasta,
    diamond_work_dir,
    num_shards,
    log_object,
    diamond_sensitivity="very-sensitive",
    diamond_block_size=None,
    evalue_cutoff=1e-3,
    threads=1,
    diamond_index_chunks=None,
) -> List[str]:
    """
    Description:
    Runs DIAMOND blastp of the query proteins against the target genomes database and streams its output directly
    into num_shards zstd-compressed shard files, so the full (potentially very large) alignment table is never written
    to disk uncompressed or held in memory. Hits are assigned to shards by a hash of the subject (representative)
    protein, so every hit to a given representative - and therefore to every member of its cluster - ends up in the
    same shard, which lets each shard be reduced independently with bounded memory.
    ********************************************************************************************************************
    Parameters:
    - target_concat_genome_db: The DIAMOND database of target genome proteins (from prepTG).
    - query_fasta: The path to the query / reference FASTA file.
    - diamond_work_dir: The workspace directory where shard files and the DIAMOND log are written.
    - num_shards: The number of shard files to split hits into.
    - log_object: A logging object.
    - diamond_sensitivity: The sensitivity mode to run DIAMOND blastp with.
    - diamond_block_size: DIAMOND block size (optional).
    - evalue_cutoff: The maximum E-value cutoff to regard an alignment as homologous.
    - threads: The number of threads to use.
    - diamond_index_chunks: Number of chunks for DIAMOND's seed index processing (-c; optional).
    ********************************************************************************************************************
    Return:
    - shard_files: Paths to the zstd-compressed shard files (columns as in FAI_DIAMOND_COLUMNS).
    ********************************************************************************************************************
    """
    diamond_cmd = ["diamond", "blastp", "--ignore-warnings", "--threads", str(threads), "--query", query_fasta,
                   "--db", target_concat_genome_db, "--outfmt", "6"] + FAI_DIAMOND_COLUMNS + \
                  ["-k0", "--evalue", str(evalue_cutoff)]
    if diamond_sensitivity != "default":
        diamond_cmd.append("--" + diamond_sensitivity)
    if diamond_block_size is not None:
        diamond_cmd += ["--block-size", str(diamond_block_size)]
    if diamond_index_chunks is not None:
        diamond_cmd += ["--index-chunks", str(diamond_index_chunks)]

    shard_files = [f"{diamond_work_dir}DIAMOND_Hits_Shard_{i}.tsv.zst" for i in range(num_shards)]
    diamond_log_file = diamond_work_dir + "DIAMOND_blastp.log"
    log_object.info(f"Running DIAMOND blastp, streaming hits into {num_shards} shard(s):\n{' '.join(diamond_cmd)}")
    try:
        shard_handles = [util.open_zstd_text_writer(shard_file) for shard_file in shard_files]
        try:
            with open(diamond_log_file, "w") as diamond_log_handle:
                # stderr goes to a file (not a pipe) so DIAMOND can never block on it while we read stdout
                proc = subprocess.Popen(diamond_cmd, stdout=subprocess.PIPE, stderr=diamond_log_handle,
                                        text=True, bufsize=1024 * 1024)
                assert proc.stdout is not None
                for line in proc.stdout:
                    fields = line.split("\t", 2)
                    if len(fields) < 3:
                        continue
                    shard = zlib.crc32(fields[1].encode()) % num_shards
                    shard_handles[shard].write(line)
                return_code = proc.wait()
        finally:
            for shard_handle in shard_handles:
                shard_handle.close()

        if return_code != 0:
            with open(diamond_log_file) as diamond_log_handle:
                diamond_stderr = diamond_log_handle.read()[-5000:]
            raise RuntimeError(f"DIAMOND blastp exited with code {return_code}:\n{diamond_stderr}")
    except Exception as e:
        msg = f"Issues with running DIAMOND blastp: {e}"
        log_object.error(msg)
        sys.stderr.write(msg + "\n")
        util.report_cause(fallback_traceback=True)
        sys.exit(1)

    return shard_files


def process_diamond_blastp_shards(
    shard_files, target_annot_information, rep_prot_to_nonreps, summary_dir, log_object,
    block_size=SAMPLE_BLOCK_SIZE, delete_shards=True
) -> Dict[str, Dict[str, Set[str]]]:
    """
    Description:
    Processes sharded DIAMOND blastp results (from run_diamond_blastp_into_shards): finds the best hit(s) - by
    bitscore, keeping ties - for each target (representative) protein, expands representatives to all members of
    their cluster (when prepTG clustered the database with DIAMOND linclust), and writes the hits as one
    zstd-compressed summary per block of target genomes (see sample_blocks / block_summary_file).

    Memory is bounded by one shard at a time, and at most max_open_files() block summaries are open at once; if there
    are more blocks than that, summaries are written in several passes over the (compressed) shards.
    ********************************************************************************************************************
    Parameters:
    - shard_files: Paths to the zstd-compressed DIAMOND hit shards.
    - target_annot_information: Dictionary with target genome sample names as keys.
    - rep_prot_to_nonreps: Dictionary mapping cluster representative proteins to tuples of all cluster members (or
                           None if the database was not clustered).
    - summary_dir: Directory where per-block hit summaries are written.
    - log_object: A logging object.
    - block_size: The number of target genomes per block (must match identify_gc_instances).
    - delete_shards: Whether to delete the shard files once processed.
    ********************************************************************************************************************
    Return:
    - diamond_results: A dictionary with sample identifiers as keys and values as a dictionary with the following keys:
            - lts: A set of locus tags.
            - hgs: A set of homolog groups.
    ********************************************************************************************************************
    """
    diamond_results: Dict[str, Dict[str, Set[str]]] = defaultdict(lambda: {'lts': set(), 'hgs': set()})
    try:
        util.setup_ready_directory([summary_dir], delete_if_exist=True)
        blocks = sample_blocks(target_annot_information.keys(), block_size)
        sample_to_block = {sample: i for i, block in enumerate(blocks) for sample in block}

        for block_group in util.chunks(list(range(len(blocks))), max_open_files()):
            block_handles = {i: util.open_zstd_text_writer(block_summary_file(summary_dir, i)) for i in block_group}
            try:
                for block_handle in block_handles.values():
                    block_handle.write(DIAMOND_SUMMARY_HEADER)
                for shard_file in shard_files:
                    # best hit(s) per subject (representative) protein: [best bitscore, [hit fields, ...]]
                    best_hits: Dict[str, List[Any]] = {}
                    with util.open_text_file(shard_file) as shard_handle:
                        for line in shard_handle:
                            hg, rep_lt, identity, evalue, bitscore, qlen, slen = line.rstrip("\n").split("\t")
                            bitscore_value = float(bitscore)
                            hit = (hg, bitscore, evalue, identity, qlen, slen)
                            current_best = best_hits.get(rep_lt)
                            if current_best is None or bitscore_value > current_best[0]:
                                best_hits[rep_lt] = [bitscore_value, [hit]]
                            elif bitscore_value == current_best[0]:
                                current_best[1].append(hit)

                    for rep_lt, (_, hits) in best_hits.items():
                        members = (rep_lt,)
                        if rep_prot_to_nonreps is not None and rep_lt in rep_prot_to_nonreps:
                            members = rep_prot_to_nonreps[rep_lt]
                        for member in members:
                            sample, lt = member.split("|", 1)
                            block_index = sample_to_block.get(sample)
                            if block_index is None or block_index not in block_handles:
                                continue
                            block_handle = block_handles[block_index]
                            diamond_results[sample]['lts'].add(lt)
                            for hg, bitscore, evalue, identity, qlen, slen in hits:
                                diamond_results[sample]['hgs'].add(hg)
                                block_handle.write(f"{lt}\t{hg}\t{sample}\t{bitscore}\t{evalue}\t{identity}\t{qlen}\t{slen}\n")
                    del best_hits
            finally:
                for block_handle in block_handles.values():
                    block_handle.close()

        samples_without_hits = sorted(set(target_annot_information.keys()).difference(diamond_results.keys()))
        if len(samples_without_hits) > 0:
            msg = (f"Warning: Did not detect homology whatsoever at requested E-value for {len(samples_without_hits)} "
                   f"target genome(s) - they are listed in the log file.")
            sys.stderr.write(msg + "\n")
            log_object.warning(msg + "\n" + "\n".join(samples_without_hits))

        if delete_shards:
            for shard_file in shard_files:
                os.remove(shard_file)
    except Exception as e:
        log_object.error("Issues with processing of DIAMOND blastp results.")
        sys.stderr.write("Issues with processing of DIAMOND blastp results.\n")
        util.report_cause(fallback_traceback=True)
        sys.exit(1)

    return dict(diamond_results)


def run_and_process_diamond_blastp(
    target_concat_genome_db,
    query_fasta,
    target_annot_information,
    rep_prot_to_nonreps,
    diamond_work_dir,
    summary_dir,
    log_object,
    diamond_sensitivity="very-sensitive",
    diamond_block_size=None,
    evalue_cutoff=1e-3,
    threads=1,
    block_size=SAMPLE_BLOCK_SIZE,
    diamond_index_chunks=None,
) -> Dict[str, Dict[str, Set[str]]]:
    """
    Runs fai's DIAMOND blastp search and processes the results into per-block hit summaries (see
    run_diamond_blastp_into_shards and process_diamond_blastp_shards). Roughly one shard is used per block of target
    genomes, capped by max_open_files(). Returns the diamond_results dictionary of process_diamond_blastp_shards.
    """
    num_blocks = max(1, math.ceil(len(target_annot_information) / block_size))
    num_shards = min(max_open_files(), num_blocks)
    shard_files = run_diamond_blastp_into_shards(
        target_concat_genome_db, query_fasta, diamond_work_dir, num_shards, log_object,
        diamond_sensitivity=diamond_sensitivity, diamond_block_size=diamond_block_size,
        evalue_cutoff=evalue_cutoff, threads=threads, diamond_index_chunks=diamond_index_chunks
    )
    return process_diamond_blastp_shards(
        shard_files, target_annot_information, rep_prot_to_nonreps, summary_dir, log_object, block_size=block_size
    )


def map_key_proteins_to_homolog_groups(
    query_fasta, key_protein_queries_fasta, work_dir, log_object, threads=1, skip_key_protein_validation=False
) -> Set[str]:
    """
    Description:
    Function to align key protein queries FASTA to general query ortholog FASTA file and determine which general
    query ortholog groups / non-redundant protein sequences correspond to the key ones provided by the user. The best
    hit for key proteins to general query sequences will be taken and assigned by bitscore (assuming an E-value
    theshold of 1e-10 is met).
    ********************************************************************************************************************
    Parameters:
    - query_fasta: A FASTA file of query ortholog group (consensus) sequences.
    - key_protein_queries_fasta: A FASTA file of key sequences.
    - work_dir: The workspace directory where to perform mapping alignment and write intermediate files.
    - log_object: A logging object.
    - threads: The number of threads to use.
    - skip_key_protein_validation: Whether to skip validation of key protein queries against query proteins.
    ********************************************************************************************************************
    Return:
    - key_hgs: A set of "key" designated query ortholog groups.
    ********************************************************************************************************************
    """
    key_hgs = set([])
    try:
        work_dir = os.path.abspath(work_dir) + '/'
        search_res_dir = work_dir + "Alignment_Results/"
        util.setup_ready_directory([search_res_dir])

        query_dmnd_db = query_fasta.split(".faa")[0] + ".dmnd"
        align_result_file = (
            work_dir + "Key_Proteins_to_General_Queries_Alignment.txt"
        )
        diamond_cmd = ["diamond", "makedb", "--ignore-warnings", "--in", query_fasta,
                       "-d", query_dmnd_db, "--threads", str(threads), ";", "diamond",
                       "blastp", "--ignore-warnings", "--threads", str(threads),
                       "--very-sensitive", "--query", key_protein_queries_fasta,
                       "--db", query_dmnd_db, "--outfmt", 
                       "6", "qseqid", "sseqid", "pident", "bitscore", "qcovhsp", "scovhsp", "--out",
                       align_result_file, "--evalue", "1e-5"]

        util.run_cmd_via_subprocess(diamond_cmd, log_object=log_object, check_files = [align_result_file])

        kq_top_hits: Dict[str, KqTopHitInfo] = defaultdict(
            lambda: {"hg_list": [], "best_bitscore": 0.0}
        )
        with open(align_result_file) as oarf:
            for line in oarf:
                line = line.strip()
                ls = line.split("\t")
                kq, hg, pident, bs, qcovhsp, scovhsp = ls
                bs = float(bs)
                pident = float(pident)
                qcovhsp = float(qcovhsp)
                scovhsp = float(scovhsp)  
                if not skip_key_protein_validation:
                    if pident <= 99.0: continue # identity
                    if qcovhsp <= 99.0: continue # query coverage
                    if scovhsp <= 99.0: continue # subject coverage
                if bs > kq_top_hits[kq]["best_bitscore"]:
                    kq_top_hits[kq] = {"hg_list": [hg], "best_bitscore": bs}
                elif bs == kq_top_hits[kq]["best_bitscore"]:
                    kq_top_hits[kq]["hg_list"].append(hg)

        with open(key_protein_queries_fasta) as okpqf:
            for rec in SeqIO.parse(okpqf, "fasta"):
                try:
                    if not (len(kq_top_hits[rec.id]["hg_list"]) >= 1):
                        raise AssertionError()
                    for hg in kq_top_hits[rec.id]["hg_list"]:
                        key_hgs.add(hg)
                except Exception as e2:
                    msg = f"Found no mapping for key query protein {rec.id} to general proteins."
                    log_object.error(msg)
                    sys.stderr.write(msg + "\n")
                    log_object.error(traceback.format_exc())
                    sys.exit(1)

    except Exception as e:
        log_object.error(
            "Issues with mapping key proteins to general proteins or consensus sequence of ortholog groups."
        )
        sys.stderr.write(
            "Issues with mapping key proteins to general proteins or consensus sequence of ortholog groups.\n"
        )
        sys.stderr.write(str(e) + "\n")
        raise RuntimeError(traceback.format_exc())

    return key_hgs


def load_target_genome_info_from_pkl(
    target_genomes_pkl_dir,
    sample, 
    diamond_results,
    min_distinct_genes,
    log_object,
) -> Tuple[Dict[str, Any], Dict[str, Any], Set[str]]:
    """
    Description:
    This function loads target genome information from a pickle file.
    ********************************************************************************************************************
    Parameters:
    - sample: The target genome sample identifier.
    - target_genomes_pkl_dir: The directory where pickle files on target genome CDS location information are stored.
    - diamond_results: A dictionary of diamond results.
    - min_distinct_genes: The minimum number of distinct genes per scaffold to be considered for a gene cluster 
                          segment.
    - log_object: A logging object.
    ********************************************************************************************************************
    Return:
    - target_genome_info: A dictionary of dictionaries:
            - gene_locations: A multi - tiered dictionary with location information for CDS features.
            - scaffold_genes: A double dictionary with keys corresponding to samples and secondary keys to scaffolds and
                              values to CDS features found on them.
            - boundary_genes: A dictionary with keys corresponding to samples and values to sets of CDS features nearby
                              scaffold edges.
    ********************************************************************************************************************
    """
    target_genome_info: Tuple[Dict[str, Any], Dict[str, Any], Set[str]] = (defaultdict(dict), defaultdict(set), set())
    try:
        sample_info_pkl_file = target_genomes_pkl_dir + sample + ".pkl"
        genbank_info = None
        try:
            if not (os.path.isfile(sample_info_pkl_file)):
                raise AssertionError()
            with open(sample_info_pkl_file, "rb") as handle:
                genbank_info = pickle.load(handle)
        except Exception as e:
            msg = f"Could not find pickle file with CDS coordinate information for sample {sample}"
            sys.stderr.write(msg + "\n")
            log_object.error(msg)
            util.report_cause(fallback_traceback=True)
            sys.exit(1)

        gene_locs = genbank_info[0]
        scaff_genes = genbank_info[1]
        bound_genes = genbank_info[2]

        gene_locs_filt: Dict[str, Any] = {}
        scaff_genes_filt: Dict[str, Any] = defaultdict(set)
        bound_genes_filt: Set[str] = set([])
        for scaffold in scaff_genes:
            if (
                len(
                    scaff_genes[scaffold].intersection(
                        diamond_results[sample]['lts']
                    )
                )
                >= min_distinct_genes
            ):
                scaff_genes_filt[scaffold] = scaff_genes[scaffold]
                for lt in scaff_genes[scaffold]:
                    gene_locs_filt[lt] = gene_locs[lt]
                bound_genes_filt = bound_genes_filt.union(
                    scaff_genes[scaffold].intersection(bound_genes)
                )
        
        target_genome_info = (gene_locs_filt, scaff_genes_filt, bound_genes_filt)

    except Exception as e:
        log_object.error(f"Issues with loading target genome information for sample {sample} file.")
        sys.stderr.write(str(e) + "\n")
        sys.exit(1)

    return target_genome_info

GC_STATE = 0
BG_STATE = 1


class GeneClusterHMM(BaseHMM):
    """
    Two-state HMM (state 0 = gene cluster, state 1 = background) used by fai to delineate gene clusters.

    Each gene is a symbol: "background" if it lacks a hit to the query gene cluster, otherwise the ortholog group it
    hit. Every ortholog group gets the same emission probability, and the per-state emission values are used as given
    (they need not sum to 1 across symbols), matching the earlier pomegranate-based model. hmmlearn does the
    forward-backward decoding; only the emission log-likelihoods are supplied here.

    Plain numpy/hmmlearn objects, so the model pickles cleanly for multiprocessing.
    """

    def __init__(self, labels=("background",), emission_logprob=None, startprob=(0.5, 0.5),
                 transmat=((0.9, 0.1), (0.1, 0.9))):
        super().__init__(n_components=2, algorithm="map", params="", init_params="")
        self.labels = list(labels)
        self.emission_logprob = emission_logprob
        self.startprob = startprob
        self.transmat = transmat
        self.startprob_ = numpy.asarray(startprob, dtype=float)
        self.transmat_ = numpy.asarray(transmat, dtype=float)
        self._label_index = {label: i for i, label in enumerate(self.labels)}

    def _compute_log_likelihood(self, X):
        return self.emission_logprob[numpy.asarray(X, dtype=int)[:, 0]]

    def encode(self, hgs_ordered):
        """ Convert an ordered list of ortholog group labels (or "background") into the model's symbol matrix. """
        try:
            return numpy.array([[self._label_index[hg]] for hg in hgs_ordered], dtype=int)
        except KeyError as e:
            raise ValueError(f"Ortholog group {e} is not part of the gene-cluster HMM.") from None


def build_gene_cluster_hmm(hgs, gc_emission_prob_with_hit, gc_emission_prob_without_hit,
                           gc_to_gc_transition_prob, bg_to_bg_transition_prob) -> GeneClusterHMM:
    """
    Description:
    Build the two-state gene cluster / background HMM used by fai's HMM delineation mode.
    ********************************************************************************************************************
    Parameters:
    - hgs: The ortholog groups of the query gene cluster.
    - gc_emission_prob_with_hit: Emission probability of a gene with a hit, in the gene-cluster state.
    - gc_emission_prob_without_hit: Emission probability of a gene without a hit, in the gene-cluster state.
    - gc_to_gc_transition_prob: Probability of staying in the gene-cluster state.
    - bg_to_bg_transition_prob: Probability of staying in the background state.
    ********************************************************************************************************************
    Returns:
    - A GeneClusterHMM whose labels are ["background"] + hgs.
    ********************************************************************************************************************
    """
    for name, value in [("gc_emission_prob_with_hit", gc_emission_prob_with_hit),
                        ("gc_emission_prob_without_hit", gc_emission_prob_without_hit),
                        ("gc_to_gc_transition_prob", gc_to_gc_transition_prob),
                        ("bg_to_bg_transition_prob", bg_to_bg_transition_prob)]:
        if not 0.0 <= float(value) <= 1.0:
            raise ValueError(f"{name} must be between 0.0 and 1.0, got {value}.")

    labels = ["background"] + [hg for hg in hgs if hg != "background"]
    n_hgs = len(labels) - 1
    gc_probs = [gc_emission_prob_without_hit] + [gc_emission_prob_with_hit] * n_hgs
    bg_probs = [1.0 - gc_emission_prob_without_hit] + [1.0 - gc_emission_prob_with_hit] * n_hgs
    with numpy.errstate(divide="ignore"):
        emission_logprob = numpy.log(numpy.array([gc_probs, bg_probs], dtype=float).T)

    transmat = [
        [gc_to_gc_transition_prob, 1.0 - gc_to_gc_transition_prob],
        [1.0 - bg_to_bg_transition_prob, bg_to_bg_transition_prob],
    ]
    # Start probabilities are equal; the previous model's end probabilities were also equal (0.5/0.5), which only
    # scales every path by the same constant and so does not affect posterior decoding.
    return GeneClusterHMM(labels=labels, emission_logprob=emission_logprob, startprob=[0.5, 0.5], transmat=transmat)


def predict_gene_cluster_states(model: GeneClusterHMM, hgs_ordered) -> List[int]:
    """ Posterior (forward-backward) decoding of gene-cluster (0) vs background (1) state for each gene in order. """
    if len(hgs_ordered) == 0:
        return []
    return model.predict(model.encode(hgs_ordered)).tolist()


def identify_gc_instances(
    query_information,
    target_information,
    diamond_results_summary_dir,
    diamond_results,
    target_genomes_pkl_dir,
    work_dir,
    log_object,
    min_hits=5,
    min_key_hits=3,
    draft_mode=False,
    gc_to_gc_transition_prob=0.9,
    bg_to_bg_transition_prob=0.9,
    gc_emission_prob_with_hit=0.95,
    gc_emission_prob_without_hit=0.2,
    syntenic_correlation_threshold=0.8,
    max_int_genes_for_merge=0,
    kq_evalue_threshold=1e-20,
    flanking_context=1000,
    threads=1,
    block_size=SAMPLE_BLOCK_SIZE,
    gc_delineation_mode="GENE-CLUMPER",
    min_distinct_genes=3,
) -> None:
    """
    Description:
    This function sets up identification of homologous / orthologous instances of gene clusters in target genomes based
    on DIAMOND blastp results and location information of genes in the query and target genomes. It is void and produces
    two types of files in subdirectories of work_dir, one being an extracted instance of homologous gene - clusters from
    target genomes and the other a TSV with information for individual CDS genes on how well they are matching query
    sequence.
    ********************************************************************************************************************
    Parameters:
    - query_information: A multi-tiered dictionary with query sequence information:
            - protein_to_hg: Dictionary mapping individual query sequences to ortholog groups.
            - query_fasta: The path to a FASTA of query sequences.
        - comp_gene_info: Dictionary linking to query CDS location / syntenic information.
    - target_information: A dictionary where keys are target genome sample identifiers and values are paths to GenBank files.
    - diamond_results_summary_dir: Directory with per-block (zstd-compressed) summaries of DIAMOND blastp hits, as
                                   written by process_diamond_blastp_shards with the same block_size.
    - work_dir: The workspace where to write files to for the analysis / function.
    - log_object: A logging object.
    - min_hits: The minimum number of non-redundant or ortholog group query sequences needed to be hit for gene cluster
                detection.
    - min_key_hits: The minimum number of key non-redundant or ortholog group query sequences needed to be hit for gene
                    cluster detection.
    - draft_mode: Whether to use draft mode and regard cutoffs in aggregate across candidate gene cluster instances
                  found near scaffold edges.
    - gc_to_gc_transition_prob: The gene cluster to gene cluster transition probability to set for HMM.
    - bg_to_bg_transition_prob: The background to background transition probability to set for HMM.
    - gc_emission_prob_with_hit: The gene cluster emission probability for the gene cluster state when a DIAMOND hit is
                                 present at a CDS / ORF.
    - gc_emission_prob_without_hit: The gene cluster emission probability for the gene cluster state  when no DIAMOND
                                    hit is present at a CDS / ORF.
    - syntenic_correlation_threshold: The threshold for the absolute syntenic correlation between the candidate gene
                                                                      cluster homologous instance and the query gene cluster.
    - max_int_genes_for_merge: The maximum number of intermediate genes between those aligning to the query proteins to
                               be regarded as co-located enough to correspond to a potential single instance of the
                               gene cluster.
    - kq_evalue_threshold: The E-value threshold for "key" designated ortholog groups to be regarded as actually being
                           hit.
    - flanking_context: The number of basepairs to take around identified gene clusters when extracting the gene cluster
                        specific GenBank from the full target genome GenBank.
    - threads: The number of threads to use.
    - block_size: The maximum number of samples to consider at a time. This might be good to make an adjustable option
                  for those interested in applying fai to metagenomes, where setting it to a much lower value likely
                  makes sense.
    - gc_delineation_mode: The method to use for gene-cluster delineation, can either be "GENE-CLUMPER" or "HMM"
    - min_distinct_genes: The minimum number of distinct genes per scaffold to be considered for a gene-cluster 
                          segment.
    ********************************************************************************************************************
    """
    try:
        work_dir = os.path.abspath(work_dir) + '/'

        # Check platform
        import platform
        system_platform = platform.system()
        if system_platform == 'Linux':
            platform_type = 'linux'
        elif system_platform == 'Darwin':
            platform_type = 'macos'
            # Set environment variables to limit threading for macOS multiprocessing compatibility
            os.environ['OMP_NUM_THREADS'] = '1'  # Limit OpenMP threads
            os.environ['MKL_NUM_THREADS'] = '1'  # Limit MKL threads
            os.environ['OPENBLAS_NUM_THREADS'] = '1'  # Limit OpenBLAS threads
            os.environ['VECLIB_MAXIMUM_THREADS'] = '1'  # Limit Accelerate framework threads
            # Ensure multiprocessing uses spawn method on macOS
            multiprocessing.set_start_method('spawn', force=True)
        else:
            platform_type = 'other'
        
        msg = f"Detected platform: {system_platform} ({platform_type})\nWill use {'multiprocessing (fork)' if platform_type == 'linux' else 'multiprocessing (spawn)' if platform_type == 'macos' else 'single-threaded'} processing" 
        log_object.info(msg)
        sys.stdout.write(msg + '\n')

        gc_genbanks_dir = work_dir + "GeneCluster_Genbanks/"
        gc_info_dir = work_dir + "GeneCluster_Info/"
        util.setup_ready_directory([gc_genbanks_dir, gc_info_dir])

        # unpack information in dictionaries
        query_gene_info = query_information["comp_gene_info"]
        single_query_mode = query_information["single_query_mode"]
        lt_to_hg = query_information["protein_to_hg"]
        all_hgs = set(lt_to_hg.values()) if lt_to_hg else set()
        key_hgs = query_information["key_hgs"]

        model = None
        model_labels = None
        if gc_delineation_mode == "HMM":
            model = build_gene_cluster_hmm(
                sorted(all_hgs),
                gc_emission_prob_with_hit,
                gc_emission_prob_without_hit,
                gc_to_gc_transition_prob,
                bg_to_bg_transition_prob,
            )
            model_labels = model.labels

        gc_hmm_evalues_file = (
            work_dir + "GeneCluster_NewInstances_HMMEvalues.txt"
        )
        gc_list_file = work_dir + "GeneCluster_Homolog_Listings.txt"

        # open handle to file where expanded GCF listings will be written
        assert target_information is not None
        # blocks must match those used when writing the per-block DIAMOND hit summaries (process_diamond_blastp_shards)
        for block_index, block_samp_list in enumerate(sample_blocks(target_information.keys(), block_size)):
            block_samp_set = set(block_samp_list)
            
            # start process of finding new BGCs
            sample_lt_to_hg = defaultdict(dict)
            sample_hgs = defaultdict(set)
            sample_lt_to_evalue = defaultdict(dict)
            sample_lt_to_identity = defaultdict(dict)
            sample_lt_to_sqlratio = defaultdict(dict)
            sample_lt_to_bitscore = defaultdict(dict)

            sample_data = []
            with util.open_text_file(block_summary_file(diamond_results_summary_dir, block_index)) as drsf:
                # Locus Tag\tHomolog Group\tSample\tBest Bitscore\tE-value\tIdentity\tQuery Length\tSubject Length
                next(drsf, None)  # header
                for line in drsf:
                    line = line.strip()
                    lt, hg, sample, bitscore, evalue, identity, qlen, slen = line.split("\t")
                    sqlratio = float(slen) / float(qlen)
                    if sample in block_samp_set:
                        
                        sample_data.append([lt, hg, sample, float(bitscore), float(evalue), float(identity), float(sqlratio)])

            # previously was just selecting one at random since we already had selected hits based 
            # on best bitscore in process_diamond_results() - this is just if there is a tie - 
            # then use E-value to try to distinguish            
            already_accounted_lts = set([])
            for sdr in sorted(sample_data, key=itemgetter(4)):
                lt, hg, sample, bitscore, evalue, identity, sqlratio = sdr
                if lt in already_accounted_lts: continue
                sample_lt_to_hg[sample][lt] = hg
                sample_hgs[sample].add(hg)
                sample_lt_to_bitscore[sample][lt] = bitscore
                sample_lt_to_evalue[sample][lt] = evalue
                sample_lt_to_identity[sample][lt] = identity
                sample_lt_to_sqlratio[sample][lt] = sqlratio
                already_accounted_lts.add(lt)
            
            identify_gc_segments_input = []
            for sample in sample_hgs:
                # load sample gene 
                sample_gene_location, sample_scaffold_genes, sample_boundary_genes = load_target_genome_info_from_pkl(target_genomes_pkl_dir, sample, diamond_results, min_distinct_genes, log_object)
                
                # if len(sample_hgs[sample]) < 3: continue
                hgs_ordered_dict = {}
                lts_ordered_dict = {}
                for scaffold in sample_scaffold_genes:
                    # Type assertion to help the type checker
                    assert sample_gene_location is not None 
                    lts_with_start = []
                    for lt in sample_scaffold_genes[scaffold]:
                        try:
                            lts_with_start.append(
                                [lt, sample_gene_location[lt]["start"]]  # type: ignore
                            )
                        except KeyError:
                            log_object.error(f"KeyError: {lt} not found in {sample}")
                            util.report_cause()
                            sys.exit(1)
                        except Exception as e:
                            log_object.error(f"Error getting start for {lt} in {sample}")
                            log_object.error(e)
                            sys.exit(1)

                    hgs_ordered = []
                    lts_ordered = []
                    for lt_start in sorted(lts_with_start, key=itemgetter(1)):
                        lt, start = lt_start
                        lts_ordered.append(lt)
                        if lt in sample_lt_to_hg[sample]:
                            hgs_ordered.append(sample_lt_to_hg[sample][lt])
                        else:
                            hgs_ordered.append("background")
                    hgs_ordered_dict[scaffold] = hgs_ordered
                    lts_ordered_dict[scaffold] = lts_ordered

                identify_gc_segments_input.append(
                    [
                        sample,
                        min_hits,
                        min_key_hits,
                        key_hgs,
                        kq_evalue_threshold,
                        syntenic_correlation_threshold,
                        max_int_genes_for_merge,
                        flanking_context,
                        draft_mode,
                        gc_delineation_mode,
                        sample_gene_location,
                        sample_boundary_genes,
                        hgs_ordered_dict,
                        lts_ordered_dict,
                        target_information[sample],
                        sample_lt_to_bitscore[sample],
                        sample_lt_to_evalue[sample],
                        sample_lt_to_identity[sample],
                        sample_lt_to_sqlratio[sample],
                        gc_genbanks_dir,
                        gc_info_dir,
                        query_gene_info,
                        lt_to_hg,
                        model,
                        model_labels,
                        single_query_mode
                    ]
                )

            msg = (
                f"Delineating gene clusters for {len(identify_gc_segments_input)} genomes determined as candidates for potential carriage."
            )
            log_object.info(msg)
            sys.stdout.write(msg + '\n')

            # Use multiprocessing on Linux/macOS, otherwise single-threaded processing. A worker error stops the run
            # (WorkerError); only a failure of the pool itself falls back to single-threaded processing.
            run_serially = platform_type not in ['linux', 'macos']
            if not run_serially:
                try:
                    util.robust_multiprocess_executor(
                        worker_function=_identify_gc_instances_worker,
                        inputs=identify_gc_segments_input,
                        pool_size=threads,
                        error_strategy="report_and_stop", 
                        log_object=log_object,
                        description="gene cluster identification"
                    )
                    log_object.info("Successfully completed multiprocessing gene cluster identification")
                except util.WorkerError:
                    raise
                except Exception as e:
                    log_object.error(f"Multiprocessing failed: {e}")
                    log_object.info("Falling back to single-threaded processing")
                    run_serially = True
            if run_serially:
                util.robust_serial_executor(
                    worker_function=_identify_gc_instances_worker,
                    inputs=identify_gc_segments_input,
                    error_strategy="report_and_stop",
                    log_object=log_object,
                    description="gene cluster identification"
                )

        util.concat_files(
            sorted(glob.glob(os.path.join(gc_info_dir, "**", "*.bgcs.txt"), recursive=True)),
            gc_list_file,
            append=True,
        )
        util.concat_files(
            sorted(glob.glob(os.path.join(gc_info_dir, "**", "*.hg_evalues.txt"), recursive=True)),
            gc_hmm_evalues_file,
            append=True,
        )

    except Exception as e:
        msg = "Issues with managing algorithm to delineate gene cluster homolog segments."
        log_object.error(msg)
        sys.stderr.write(msg + "\n")
        util.report_cause(fallback_traceback=True)
        sys.exit(1)


def _identify_gc_instances_worker(input_args) -> Tuple[str, Optional[str]]:
    """ This function is the actual one which identifies gene clusters. """
    try:
        (
            sample,
            min_hits,
            min_key_hits,
            key_hgs,
            kq_evalue_threshold,
            syntenic_correlation_threshold,
            max_int_genes_for_merge,
            flanking_context,
            draft_mode,
            gc_delineation_mode,
            gene_locations,
            boundary_genes,
            hgs_ordered_dict,
            lts_ordered_dict,
            target_annotation_info,
            sample_lt_to_bitscore,
            sample_lt_to_evalue,
            sample_lt_to_identity,
            sample_lt_to_sqlratio,
            gc_genbanks_dir, 
            gc_info_dir, 
            query_gene_info, 
            lt_to_hg, 
            model, 
            model_labels, 
            single_query_mode
        ) = input_args
        bg_set = set(["background"])

        # Type assertions to help the type checker
        assert hgs_ordered_dict is not None
        assert lts_ordered_dict is not None
        assert gene_locations is not None
        assert target_annotation_info is not None
        assert sample_lt_to_identity is not None
        assert sample_lt_to_sqlratio is not None
        assert sample_lt_to_bitscore is not None
        assert sample_lt_to_evalue is not None
        assert boundary_genes is not None
        assert lt_to_hg is not None
        if gc_delineation_mode == "HMM":
            assert model_labels is not None
            assert model is not None

        if single_query_mode:
            with open(gc_info_dir + sample + ".bgcs.txt", "w") as gc_sample_listing_handle:
                with open(gc_info_dir + sample + ".hg_evalues.txt", "w") as gc_hg_evalue_handle:
                    sample_gc_id = 1
                    # regions are collected and extracted together, reading the genome GenBank only once
                    regions_to_extract = []
                    for scaffold in hgs_ordered_dict:
                        hgs_ordered = hgs_ordered_dict[scaffold]
                        lts_ordered = lts_ordered_dict[scaffold]
                        for i, hg in enumerate(hgs_ordered):
                            if hg != "background":
                                lt = lts_ordered[i]
                                gc_genbank_file = (
                                    gc_genbanks_dir
                                    + sample
                                    + "_fai-gene-cluster-"
                                    + str(sample_gc_id)
                                    + ".gbk"
                                )
                                sample_gc_id += 1

                                min_gc_pos = (
                                    gene_locations[lt]["start"] - flanking_context
                                )
                                max_gc_pos = (
                                    gene_locations[lt]["end"] + flanking_context
                                )

                                regions_to_extract.append(
                                    (gc_genbank_file, scaffold, min_gc_pos, max_gc_pos)
                                )
                                gc_sample_listing_handle.write(
                                    "\t".join([sample, gc_genbank_file]) + "\n"
                                )

                                identity = 0.0
                                sqlratio = 0.0
                                bitscore = 0.0
                                if lt in sample_lt_to_identity:
                                    identity = sample_lt_to_identity[lt]
                                if lt in sample_lt_to_sqlratio:
                                    sqlratio = sample_lt_to_sqlratio[lt]
                                if lt in sample_lt_to_bitscore:
                                    bitscore = sample_lt_to_bitscore[lt]
                                gc_hg_evalue_handle.write(
                                    "\t".join(
                                        [
                                            gc_genbank_file,
                                            sample,
                                            lt,
                                            hg,
                                            str(bitscore),
                                            str(identity),
                                            str(sqlratio),
                                            str(hg in key_hgs),
                                        ]
                                    )
                                    + "\n"
                                )

                    util.create_genbanks(target_annotation_info, regions_to_extract)

            return ('success', None)

        sample_gc_predictions = []
        if gc_delineation_mode == "GENE-CLUMPER":
            gcs_id = 1
            for scaffold in hgs_ordered_dict:
                hgs_ordered = hgs_ordered_dict[scaffold]
                lts_ordered = lts_ordered_dict[scaffold]
                tmp = [] # type: ignore
                dist_counter = 0
                hg_counter = 0
                last_hg_i = 0
                active_search_mode = False
                for i, hg in enumerate(hgs_ordered):
                    if hg != "background":
                        tmp.append(i)
                        dist_counter = 0
                        hg_counter += 1
                        last_hg_i = i
                        active_search_mode = True
                    elif (
                        dist_counter < max_int_genes_for_merge
                        and active_search_mode
                    ):
                        tmp.append(i)
                        dist_counter += 1
                    else:
                        if hg_counter >= 3:
                            gc_state_lts = [] # type: ignore
                            gc_state_hgs = [] # type: ignore

                            begin_i = min(tmp)
                            if last_hg_i == (len(hgs_ordered) - 1):
                                gc_state_lts = lts_ordered[begin_i:]
                                gc_state_hgs = hgs_ordered[begin_i:]
                            else:
                                gc_state_lts = lts_ordered[
                                    min(tmp) : last_hg_i + 1
                                ]
                                gc_state_hgs = hgs_ordered[
                                    min(tmp) : last_hg_i + 1
                                ]
                            boundary_lt_featured = False
                            features_key_hg = False
                            key_hgs_detected = set([])
                            if (
                                len(
                                    key_hgs.intersection(
                                        set(gc_state_hgs).difference(bg_set)
                                    )
                                )
                                > 0
                            ):
                                for j, lt in enumerate(gc_state_lts):
                                    curr_hg = gc_state_hgs[j]
                                    if (
                                        curr_hg in key_hgs
                                        and lt in sample_lt_to_evalue
                                        and sample_lt_to_evalue[lt]
                                        <= kq_evalue_threshold
                                    ):
                                        key_hgs_detected.add(curr_hg)
                                        features_key_hg = True
                            if (
                                len(
                                    boundary_genes.intersection(
                                        set(gc_state_lts).difference(bg_set)
                                    )
                                )
                                > 0
                            ):
                                boundary_lt_featured = True
                            sample_gc_predictions.append(
                                [
                                    gc_state_lts,
                                    gc_state_hgs,
                                    len(gc_state_lts),
                                    len(set(gc_state_hgs).difference(bg_set)),
                                    len(key_hgs_detected),
                                    scaffold,
                                    boundary_lt_featured,
                                    features_key_hg,
                                    gcs_id,
                                    key_hgs_detected,
                                ]
                            )
                            gcs_id += 1
                        tmp = [] # type: ignore
                        hg_counter = 0
                        dist_counter = 0
                        active_search_mode = False
                if hg_counter >= 3:
                    gc_state_lts = [] # type: ignore
                    gc_state_hgs = [] # type: ignore
                    begin_i = min(tmp)
                    if last_hg_i == (len(hgs_ordered) - 1):
                        gc_state_lts = lts_ordered[begin_i:]
                        gc_state_hgs = hgs_ordered[begin_i:]
                    else:
                        gc_state_lts = lts_ordered[min(tmp) : last_hg_i + 1]
                        gc_state_hgs = hgs_ordered[min(tmp) : last_hg_i + 1]
                    boundary_lt_featured = False
                    features_key_hg = False
                    key_hgs_detected = set([])
                    if (
                        len(
                            key_hgs.intersection(
                                set(gc_state_hgs).difference(bg_set)
                            )
                        )
                        > 0
                    ):
                        for j, lt in enumerate(gc_state_lts):
                            curr_hg = gc_state_hgs[j]
                            if (
                                curr_hg in key_hgs
                                and lt in sample_lt_to_evalue
                                and sample_lt_to_evalue[lt]
                                <= kq_evalue_threshold
                            ):
                                features_key_hg = True
                                key_hgs_detected.add(curr_hg)
                    if (
                        len(
                            boundary_genes.intersection(
                                set(gc_state_lts).difference(bg_set)
                            )
                        )
                        > 0
                    ):
                        boundary_lt_featured = True
                    sample_gc_predictions.append(
                        [
                            gc_state_lts,
                            gc_state_hgs,
                            len(gc_state_lts),
                            len(set(gc_state_hgs).difference(bg_set)),
                            len(key_hgs_detected),
                            scaffold,
                            boundary_lt_featured,
                            features_key_hg,
                            gcs_id,
                            key_hgs_detected,
                        ]
                    )
                    gcs_id += 1
        else:
            gcs_id = 1
            for scaffold in hgs_ordered_dict:
                hgs_ordered = hgs_ordered_dict[scaffold]
                lts_ordered = lts_ordered_dict[scaffold]

                hmm_predictions = predict_gene_cluster_states(model, hgs_ordered)

                int_bg_coords = set([])
                tmp = [] # type: ignore
                for i, hg_state in enumerate(hmm_predictions):
                    if hg_state == 1:
                        tmp.append(i)
                    else:
                        if len(tmp) <= max_int_genes_for_merge:
                            for lt_index in tmp:
                                int_bg_coords.add(lt_index)
                        tmp = []

                gc_state_lts = [] # type: ignore
                gc_state_hgs = [] # type: ignore
                for i, hg_state in enumerate(hmm_predictions):
                    lt = lts_ordered[i]
                    hg = hgs_ordered[i]
                    if hg_state == 0:
                        gc_state_lts.append(lt)
                        gc_state_hgs.append(hg)
                    if hg_state == 1 or i == (len(hmm_predictions) - 1):
                        if len(set(gc_state_hgs).difference(bg_set)) >= 3:
                            boundary_lt_featured = False
                            features_key_hg = False
                            key_hgs_detected = set([])
                            if (
                                len(
                                    key_hgs.intersection(
                                        set(gc_state_hgs).difference(bg_set)
                                    )
                                )
                                > 0
                            ):
                                for j, lt in enumerate(gc_state_lts):
                                    curr_hg = gc_state_hgs[j]
                                    if (
                                        curr_hg in key_hgs
                                        and lt in sample_lt_to_evalue
                                        and sample_lt_to_evalue[lt]
                                        <= kq_evalue_threshold
                                    ):
                                        features_key_hg = True
                                        key_hgs_detected.add(curr_hg)
                            if (
                                len(
                                    boundary_genes.intersection(
                                        set(gc_state_lts).difference(bg_set)
                                    )
                                )
                                > 0
                            ):
                                boundary_lt_featured = True
                            sample_gc_predictions.append(
                                [
                                    gc_state_lts,
                                    gc_state_hgs,
                                    len(gc_state_lts),
                                    len(set(gc_state_hgs).difference(bg_set)),
                                    len(key_hgs_detected),
                                    scaffold,
                                    boundary_lt_featured,
                                    features_key_hg,
                                    gcs_id,
                                    key_hgs_detected,
                                ]
                            )
                            gcs_id += 1
                        gc_state_lts = [] # type: ignore
                        gc_state_hgs = [] # type: ignore

                gc_state_lts = [] # type: ignore
                gc_state_hgs = [] # type: ignore
                gc_states = [] # type: ignore
                for i, hg_state in enumerate(hmm_predictions):
                    lt = lts_ordered[i]
                    hg = hgs_ordered[i]
                    if hg_state == 0 or i in int_bg_coords:
                        gc_state_lts.append(lt)
                        gc_state_hgs.append(hg)
                        gc_states.append(hg_state)
                    elif hg_state == 1:
                        if len(set(gc_state_hgs).difference(bg_set)) >= 3:
                            if "1" in set(gc_states):
                                gc_state_lts = [] # type: ignore
                                gc_state_hgs = [] # type: ignore
                                gc_states = [] # type: ignore
                                continue
                            boundary_lt_featured = False
                            features_key_hg = False
                            key_hgs_detected = set([])
                            if (
                                len(
                                    key_hgs.intersection(
                                        set(gc_state_hgs).difference(bg_set)
                                    )
                                )
                                > 0
                            ):
                                for j, lt in enumerate(gc_state_lts):
                                    curr_hg = gc_state_hgs[j]
                                    if (
                                        curr_hg in key_hgs
                                        and lt in sample_lt_to_evalue
                                        and sample_lt_to_evalue[lt]
                                        <= kq_evalue_threshold
                                    ):
                                        features_key_hg = True
                                        key_hgs_detected.add(curr_hg)
                            if (
                                len(
                                    boundary_genes.intersection(
                                        set(gc_state_lts).difference(bg_set)
                                    )
                                )
                                > 0
                            ):
                                boundary_lt_featured = True
                            sample_gc_predictions.append(
                                [
                                    gc_state_lts,
                                    gc_state_hgs,
                                    len(gc_state_lts),
                                    len(set(gc_state_hgs).difference(bg_set)),
                                    len(key_hgs_detected),
                                    scaffold,
                                    boundary_lt_featured,
                                    features_key_hg,
                                    gcs_id,
                                    key_hgs_detected,
                                ]
                            )
                            gcs_id += 1
                        gc_state_lts = [] # type: ignore
                        gc_state_hgs = [] # type: ignore
                        gc_states = [] # type: ignore
                    if i == (len(hmm_predictions) - 1):
                        if len(set(gc_state_hgs).difference(bg_set)) >= 3:
                            if "1" in set(gc_states):
                                gc_state_lts = [] # type: ignore
                                gc_state_hgs = [] # type: ignore
                                gc_states = [] # type: ignore
                                continue
                            boundary_lt_featured = False
                            features_key_hg = False
                            key_hgs_detected = set([])
                            if (
                                len(
                                    key_hgs.intersection(
                                        set(gc_state_hgs).difference(bg_set)
                                    )
                                )
                                > 0
                            ):
                                for j, lt in enumerate(gc_state_lts):
                                    curr_hg = gc_state_hgs[j]
                                    if (
                                        curr_hg in key_hgs
                                        and lt in sample_lt_to_evalue
                                        and sample_lt_to_evalue[lt]
                                        <= kq_evalue_threshold
                                    ):
                                        features_key_hg = True
                                        key_hgs_detected.add(curr_hg)
                            if (
                                len(
                                    boundary_genes.intersection(
                                        set(gc_state_lts).difference(bg_set)
                                    )
                                )
                                > 0
                            ):
                                boundary_lt_featured = True
                            sample_gc_predictions.append(
                                [
                                    gc_state_lts,
                                    gc_state_hgs,
                                    len(gc_state_lts),
                                    len(set(gc_state_hgs).difference(bg_set)),
                                    len(key_hgs_detected),
                                    scaffold,
                                    boundary_lt_featured,
                                    features_key_hg,
                                    gcs_id,
                                    key_hgs_detected,
                                ]
                            )
                            gcs_id += 1
                        gc_state_lts = [] # type: ignore
                        gc_state_hgs = [] # type: ignore
                        gc_states = [] # type: ignore

        if len(sample_gc_predictions) == 0:
            return ('success', None)

        sorted_sample_gc_predictions = [
            x
            for x in sorted(sample_gc_predictions, key=itemgetter(3), reverse=True)
        ]

        cumulative_edge_hgs = set([])
        cumulative_edge_key_hgs = set([])
        visited_scaffolds_with_edge_gc_segment = set([])
        sample_gc_predictions_filtered = []
        sample_edge_gc_predictions_filtered = []

        for gc_segment in sorted_sample_gc_predictions:
            if (gc_segment[3] >= min_hits and gc_segment[4] >= min_key_hits) or (
                gc_segment[3] >= 3
                and gc_segment[6]
                and not gc_segment[5] in visited_scaffolds_with_edge_gc_segment
            ):
                # code to determine whether syntenically, the considered segment aligns with what is expected.
                # (skipped if input mode was 3)
                input_mode_3 = False
                if query_gene_info is None:
                    input_mode_3 = True

                best_corr = "irrelevant"
                if not input_mode_3 and syntenic_correlation_threshold > 0.0:
                    best_corr = None
                    segment_hg_order = []
                    segment_hg_direction = []
                    gc_hg_orders = defaultdict(list)
                    gc_hg_directions = defaultdict(list)

                    copy_count_of_hgs_in_segment = defaultdict(int)
                    for hg in gc_segment[1]:
                        copy_count_of_hgs_in_segment[hg] += 1

                    for gi, g in enumerate(gc_segment[0]):
                        hg = gc_segment[1][gi]
                        if copy_count_of_hgs_in_segment[hg] != 1:
                            continue
                        gene_midpoint = (
                            gene_locations[g]["start"]
                            + gene_locations[g]["end"]
                        ) / 2.0
                        segment_hg_order.append(gene_midpoint)
                        segment_hg_direction.append(
                            gene_locations[g]["direction"]
                        )

                        for gc in query_gene_info: # type: ignore
                            g_matching = []
                            for g in query_gene_info[gc]: # type: ignore
                                if g in lt_to_hg:
                                    hg_of_g = lt_to_hg[g]
                                    if hg_of_g == hg:
                                        g_matching.append(g)
                            if len(g_matching) == 1:
                                gc_gene_midpoint = (
                                    query_gene_info[gc][g_matching[0]]["start"] # type: ignore
                                    + query_gene_info[gc][g_matching[0]]["end"] # type: ignore
                                ) / 2.0
                                gc_hg_orders[gc].append(gc_gene_midpoint)
                                gc_hg_directions[gc].append(
                                    query_gene_info[gc][g_matching[0]]["direction"] # type: ignore
                                )
                            else:
                                gc_hg_orders[gc].append(None)
                                gc_hg_directions[gc].append(None)

                    best_corr = None
                    for gc in query_gene_info: # type: ignore
                        try:
                            if not (len(segment_hg_order) == len(gc_hg_orders[gc])):
                                raise AssertionError()
                            list1_same_dir = []
                            list2_same_dir = []
                            list1_comp_dir = []
                            list2_comp_dir = []
                            for iter, hgval1 in enumerate(segment_hg_order):
                                hgdir1 = segment_hg_direction[iter]
                                hgval2 = gc_hg_orders[gc][iter]
                                hgdir2 = gc_hg_directions[gc][iter]
                                if hgval1 is None or hgval2 is None:
                                    continue
                                if hgdir1 is None or hgdir2 is None:
                                    continue
                                if hgdir1 == hgdir2:
                                    list1_same_dir.append(hgval1)
                                    list2_same_dir.append(hgval2)
                                else:
                                    list1_comp_dir.append(hgval1)
                                    list2_comp_dir.append(hgval2)

                            if len(list1_same_dir) >= 3:
                                corr, pval = pearsonr(
                                    list1_same_dir, list2_same_dir
                                )
                                corr = abs(corr)  # type: ignore
                                if (pval < 0.1) and (  # type: ignore
                                    (best_corr and best_corr < corr)
                                    or (not best_corr)
                                ):
                                    best_corr = corr
                            if len(list1_comp_dir) >= 3:
                                corr, pval = pearsonr(
                                    list1_comp_dir, list2_comp_dir
                                )
                                corr = abs(corr)  # type: ignore
                                if (pval < 0.1) and (  # type: ignore
                                    (best_corr and best_corr < corr)
                                    or (not best_corr)
                                ):
                                    best_corr = corr
                        except Exception as e:
                            util.report_cause(fallback_traceback=True)
                            sys.exit(1)
                            # pass
                if best_corr != "irrelevant" and (
                    best_corr is None or best_corr < syntenic_correlation_threshold
                ):
                    continue
                if gc_segment[3] >= min_hits and gc_segment[4] >= min_key_hits:
                    sample_gc_predictions_filtered.append(gc_segment + [best_corr])
                    if gc_segment[6]:
                        cumulative_edge_hgs = cumulative_edge_hgs.union(
                            set(gc_segment[1])
                        )
                        cumulative_edge_key_hgs = cumulative_edge_key_hgs.union(
                            set(gc_segment[9])
                        )
                        visited_scaffolds_with_edge_gc_segment.add(gc_segment[5])
                elif (
                    gc_segment[3] >= 3
                    and gc_segment[6]
                    and not gc_segment[5] in visited_scaffolds_with_edge_gc_segment
                ):
                    sample_edge_gc_predictions_filtered.append(
                        gc_segment + [best_corr]
                    )
                    visited_scaffolds_with_edge_gc_segment.add(gc_segment[5])
                    cumulative_edge_hgs = cumulative_edge_hgs.union(
                        set(gc_segment[1])
                    )
                    cumulative_edge_key_hgs = cumulative_edge_key_hgs.union(
                        set(gc_segment[9])
                    )

        if len(sample_edge_gc_predictions_filtered) >= 1 and draft_mode:
            if (
                len(cumulative_edge_hgs) >= min_hits
                and len(cumulative_edge_key_hgs) >= min_key_hits
            ):
                sample_gc_predictions_filtered += (
                    sample_edge_gc_predictions_filtered
                )

        # dereplicate nested segments!
        redundant_gcs = set([])
        for i, gcs1 in enumerate(
            sorted(sample_gc_predictions_filtered, key=itemgetter(2), reverse=True)
        ):
            for j, gcs2 in enumerate(
                sorted(
                    sample_gc_predictions_filtered, key=itemgetter(2), reverse=True
                )
            ):
                if i < j and len(set(gcs1[0]).intersection(set(gcs2[0]))) > 0:
                    redundant_gcs.add(gcs2[8])

        dereplicated_sample_gc_predictions_filtered = []
        for gcs in sample_gc_predictions_filtered:
            if gcs[8] in redundant_gcs:
                continue
            dereplicated_sample_gc_predictions_filtered.append(gcs)

        with open(gc_info_dir + sample + ".bgcs.txt", "w") as gc_sample_listing_handle:
            gc_hg_evalue_handle = open(gc_info_dir + sample + ".hg_evalues.txt", "w")
            gc_corr_info_handle = open(gc_info_dir + sample + ".corr_info.txt", "w")

            sample_gc_id = 1
            # regions are collected and extracted together, reading the genome GenBank only once
            regions_to_extract = []
            for gc_segment in dereplicated_sample_gc_predictions_filtered:
                gc_genbank_file = (
                    gc_genbanks_dir
                    + sample
                    + "_fai-gene-cluster-"
                    + str(sample_gc_id)
                    + ".gbk"
                )
                sample_gc_id += 1

                gc_segment_scaff = gc_segment[5]
                min_gc_pos = (
                    min([gene_locations[g]["start"] for g in gc_segment[0]])
                    - flanking_context
                )
                max_gc_pos = (
                    max([gene_locations[g]["end"] for g in gc_segment[0]])
                    + flanking_context
                )

                regions_to_extract.append(
                    (gc_genbank_file, gc_segment_scaff, min_gc_pos, max_gc_pos)
                )
                gc_sample_listing_handle.write(
                    "\t".join([sample, gc_genbank_file]) + "\n"
                )
                gc_corr_info_handle.write(
                    "\t".join([sample, gc_genbank_file, str(gc_segment[10])]) + "\n"
                )

                for i, lt in enumerate(gc_segment[0]):
                    hg = gc_segment[1][i]
                    identity = 0.0
                    sqlratio = 0.0
                    bitscore = 0.0
                    if lt in sample_lt_to_identity:
                        identity = sample_lt_to_identity[lt]
                    if lt in sample_lt_to_sqlratio:
                        sqlratio = sample_lt_to_sqlratio[lt]
                    if lt in sample_lt_to_bitscore:
                        bitscore = sample_lt_to_bitscore[lt]
                    gc_hg_evalue_handle.write(
                        "\t".join(
                            [
                                gc_genbank_file,
                                sample,
                                lt,
                                hg,
                                str(bitscore),
                                str(identity),
                                str(sqlratio),
                                str(hg in key_hgs),
                            ]
                        )
                        + "\n"
                    )

            gc_hg_evalue_handle.close()
            gc_corr_info_handle.close()

        util.create_genbanks(target_annotation_info, regions_to_extract)

        return ('success', None)
    except Exception as e:
        error_msg = f"An issue occurred with identifying gene cluster instances for sample {input_args[0]}: {str(e)}"
        return ('error', error_msg)


def filter_paralogous_segments_and_concatenate_into_multi_record_genbanks(
    hmm_work_dir, homologous_gbk_dir, single_query_mode, log_object
) -> None:
    """
    Description:
    This function allows for filtering paralogous gene cluster segments identified in target genomes where the segments
    have overlap in the query non-redundant / ortholog groups they are mapping to. Gene cluster instances are considered
    paralogous if at least two overlapping query hits are found. Gene clusters retained are prioritized based on
    cumulative bitscores to query genes.
    ********************************************************************************************************************
    Parameters:
    - hmm_work_dir: The workspace where identifyGCInstances() saved extracted GenBanks for gene cluster instances from \
                    target genomes along with homology information to the query gene cluster.
    - homologous_gbk_dir: The directory where to save the final extracted gene cluster files.
    - single_query_mode: Whether the query is a single protein.
    - log_object: A logging object.
    ********************************************************************************************************************
    """
    try:
        hmm_work_dir = os.path.abspath(hmm_work_dir) + '/'
        gbk_info_dir = hmm_work_dir + "GeneCluster_Info/"
        gbk_filt_dir = hmm_work_dir + "GeneCluster_Filtered_Segments/"
        util.setup_ready_directory([gbk_filt_dir])

        sample_gc_segs = defaultdict(set)
        sample_gcs_segs_to_filter = defaultdict(set)
        sample_gcs_hgs = defaultdict(lambda: defaultdict(set))
        sample_gcs_hg_bs = defaultdict(lambda: defaultdict(float))
        for gcs_info in os.listdir(gbk_info_dir):
            if not gcs_info.endswith(".hg_evalues.txt"):
                continue
            sample = gcs_info.split(".hg_evalues.txt")[0]
            gcs_info_file = gbk_info_dir + gcs_info
            with open(gcs_info_file) as ogif:
                for line in ogif:
                    gc_seg_gbk, _, lt, hg, bits, iden, sqlr, is_key_hg = (
                        line.split("\t")
                    )
                    sample_gc_segs[sample].add(gc_seg_gbk)
                    if hg != "background":
                        sample_gcs_hgs[sample][gc_seg_gbk].add(hg)
                        sample_gcs_hg_bs[sample][gc_seg_gbk] += float(bits)

        # greedy removal of paralogs by prioritization of higher aggregate bitscore gene - cluster segments
        for sample in sample_gcs_hgs:
            gbk_filt_file = gbk_filt_dir + sample + ".txt"
            with open(gbk_filt_file, "w") as gbk_filt_handle:
                for i, gcs1 in enumerate(
                    sorted(
                        sample_gcs_hg_bs[sample], key=itemgetter(1), reverse=True
                    )
                ):
                    for j, gcs2 in enumerate(
                        sorted(
                            sample_gcs_hg_bs[sample],
                            key=itemgetter(1),
                            reverse=True,
                        )
                    ):
                        if i >= j:
                            continue
                        gcs1_hg = sample_gcs_hgs[sample][gcs1]
                        gcs2_hg = sample_gcs_hgs[sample][gcs2]
                        # consider segments paralogous if more than 2 reference proteins / ortholog groups are overlapping
                        # suggesting paralogy beyond fragmentation that might have split a gene in two.
                        # This is not the case if single_query_mode is True.
                        intersection_hgs = gcs1_hg.intersection(gcs2_hg)
                        if len(intersection_hgs) >= 2 or (
                            single_query_mode and len(intersection_hgs) >= 1
                        ):
                            sample_gcs_segs_to_filter[sample].add(gcs2)
                            gbk_filt_handle.write(str(gcs2) + "\n")

        for sample in sample_gc_segs:
            retained_gcs = set([])
            for gcs in sample_gc_segs[sample]:
                if not gcs in sample_gcs_segs_to_filter[sample]:
                    retained_gcs.add(gcs)
            if len(retained_gcs) == 1:
                shutil.move(list(retained_gcs)[0], homologous_gbk_dir)
            elif len(retained_gcs) > 1:
                with open(
                    homologous_gbk_dir + sample + ".gbk", "a+"
                ) as sample_gbk_handle:
                    for gcs_gbk in retained_gcs:
                        with open(gcs_gbk) as oggbk:
                            for line in oggbk:
                                sample_gbk_handle.write(line)

    except Exception as e:
        log_object.error("Issues resolving paralogous gene segments.")
        sys.stderr.write("Issues resolving paralogous gene segments.\n")
        util.report_cause(fallback_traceback=True)
        sys.exit(1)


def plot_overviews(
    target_annotation_info,
    hmm_work_dir,
    protein_to_hg,
    plot_work_dir,
    plot_result_pdf,
    plot_naming_pdf,
    log_object,
    height=10,
    width=15,
) -> None:
    """
    Description:
    This function serves to create plots from fai analysis showing the similarity of homologous / orthologous
    gene cluster instances identified to the query gene cluster proteins.
    ********************************************************************************************************************
    Parameters:
    - target_annotation_info: Dictionary linking sample names to their respective full genome GenBank files.
    - hmm_work_dir: The workspace where identifyGCInstances() saved extracted GenBanks for gene cluster instances from \
                            target genomes along with homology information to the query gene cluster.
    - protein_to_hg: A dictionary linking query proteins to non-redundant clusters / ortholog groups.
    - plot_work_dir: The workspace where to write intermediate files used for creating plots.
    - plot_result_pdf: A final plot (PDF) showcasing homology of target genome gene clusters to the query \
                       gene cluster(s).
    - plot_naming_pdf: A plot (PDF) based table which shows the mapping of query protein identifiers to ortholog group \
                       or non-redudnant cluster identifiers.
    - log_object: A logging object.
    - height: The height of the plot in inches.
    - width: The width of the plot in inches.
    ********************************************************************************************************************
    """
    try:
        hmm_work_dir = os.path.abspath(hmm_work_dir) + '/'

        que_info_table_file = plot_work_dir + "Name_Mapping.txt"
        with open(que_info_table_file, "w") as que_info_table_handle:
            que_info_table_handle.write("NR_ID\tCDS_Locus_Tags\n")
            hg_prots = defaultdict(set)
            for p in protein_to_hg:
                hg_prots[protein_to_hg[p]].add(p)
            for hg in hg_prots:
                que_info_table_handle.write(
                    str(hg) + "\t"
                    + ", ".join(
                        sorted(
                            list([x.split("|")[-1].strip() for x in hg_prots[hg]])
                        )
                    )
                    + "\n"
                )

        gbk_info_dir = hmm_work_dir + "GeneCluster_Info/"
        gbk_filt_dir = hmm_work_dir + "GeneCluster_Filtered_Segments/"

        relevant_samples = set([])
        relevant_lts = set([])
        if os.path.isdir(gbk_info_dir):
            for gcs_info in os.listdir(gbk_info_dir):
                if (
                    not gcs_info.endswith(".hg_evalues.txt")
                    or os.path.getsize(gbk_info_dir + gcs_info) == 0
                ):
                    continue
                sample = gcs_info.split(".hg_evalues.txt")[0]
                relevant_samples.add(sample)
                with open(gbk_info_dir + gcs_info) as ogif:
                    for line in ogif:
                        line = line.strip()
                        ls = line.split("\t")
                        relevant_lts.add(ls[2])
        lt_dists_to_scaff_starts = defaultdict(dict)
        lt_dists_to_scaff_ends = defaultdict(dict)
        lt_starts = defaultdict(dict)
        for sample in target_annotation_info:
            if not sample in relevant_samples:
                continue
    
            with util.open_text_file(target_annotation_info[sample]) as osgbk:
                for rec in SeqIO.parse(osgbk, "genbank"):
                    scaff_len = len(str(rec.seq))
                    scaff_is_relevant = False
                    for feature in rec.features:
                        if feature.type != "CDS":
                            continue
                        lt = None
                        try:
                            lt = feature.qualifiers.get("locus_tag")[0]
                        except Exception as e:
                            try:
                                lt = feature.qualifiers.get("protein_id")[0]
                            except Exception as e:
                                sys.stderr.write(
                                    f"The GenBank {target_annotation_info[sample]}, cataloging a homologous instance to the query gene cluster had at least one CDS without either a locus_tag or protein_id feature.\n"
                                )
                                util.report_cause()
                                sys.exit(1)
                        if lt in relevant_lts:
                            scaff_is_relevant = True
                            break
                    if not scaff_is_relevant:
                        continue
                    for feature in rec.features:
                        if feature.type != "CDS":
                            continue
                        lt = None
                        try:
                            lt = feature.qualifiers.get("locus_tag")[0]
                        except Exception as e:
                            try:
                                lt = feature.qualifiers.get("protein_id")[0]
                            except Exception as e:
                                sys.stderr.write(
                                    f"The GenBank {target_annotation_info[sample]}, cataloging a homologous instance to the query gene cluster had at least one CDS without either a locus_tag or protein_id feature.\n"
                                )
                                util.report_cause()
                                sys.exit(1)

                        start, end, direction, all_coords = (
                            util.process_location_string(str(feature.location)) # type: ignore
                        )

                        lt_dists_to_scaff_starts[sample][lt] = start
                        lt_dists_to_scaff_ends[sample][lt] = scaff_len - end
                        lt_starts[sample][lt] = start

        sample_filt_gcs = defaultdict(set)
        if os.path.isdir(gbk_filt_dir):
            for gfi in os.listdir(gbk_filt_dir):
                sample = gfi.split(".txt")[0]
                with open(gbk_filt_dir + gfi) as ogfi:
                    for gcs in ogfi:
                        sample_filt_gcs[sample].add(gcs.strip())

        plot_input_file = plot_work_dir + "Bar_Plot_Input.txt"
        with open(plot_input_file, "w") as plot_input_handle:
            plot_input_handle.write(
                "\t".join(
                    [
                        "sample",
                        "segment_title",
                        "gene",
                        "key",
                        "gene_order",
                        "identity",
                        "sql_ratio",
                    ]
                )
                + "\n"
            )
            for gcs_info in os.listdir(gbk_info_dir):
                if not gcs_info.endswith(".hg_evalues.txt"):
                    continue
                sample = gcs_info.split(".hg_evalues.txt")[0]
                gcs_info_file = gbk_info_dir + gcs_info

                dists_to_scaff_start = defaultdict(list)
                dists_to_scaff_end = defaultdict(list)
                filt_status = "Retained"
                with open(gcs_info_file) as ogif:
                    for line in ogif:
                        gc_seg_gbk, _, lt, hg, bits, iden, sqlr, is_key_hg = (
                            line.strip().split("\t")
                        )
                        dists_to_scaff_start[gc_seg_gbk].append(
                            lt_dists_to_scaff_starts[sample][lt]
                        )
                        dists_to_scaff_end[gc_seg_gbk].append(
                            lt_dists_to_scaff_ends[sample][lt]
                        )

                seg_to_name: Dict[str, Any] = {}
                for gc_seg_gbk in dists_to_scaff_start:
                    gcs_id = (
                        gc_seg_gbk.split("/")[-1]
                        .split("fai-gene-cluster-")[1]
                        .split(".gbk")[0]
                    )
                    filt_status = "Retained"
                    if gc_seg_gbk in sample_filt_gcs[sample]:
                        filt_status = "Filtered"
                    gcs_name = (
                        sample
                        + f" Segment_{gcs_id} ("
                        + filt_status
                        + "); "
                        + str(min(dists_to_scaff_start[gc_seg_gbk]))
                        + " bp to Scaffold Start; "
                        + str(min(dists_to_scaff_end[gc_seg_gbk]))
                        + " bp to Scaffold End"
                    )
                    seg_to_name[gc_seg_gbk] = gcs_name

                hg_iters = defaultdict(int)
                with open(gcs_info_file) as ogif:
                    for line in ogif:
                        gc_seg_gbk, _, lt, hg, bits, iden, sqlr, is_key_hg = (
                            line.strip().split("\t")
                        )
                        if hg == "background":
                            hg = "nov_aux"
                        hg_iters[hg] += 1
                        hg_id = hg + "|" + str(hg_iters[hg])
                        plot_input_handle.write(
                            "\t".join(
                                [
                                    sample,
                                    seg_to_name[gc_seg_gbk],
                                    hg_id,
                                    is_key_hg,
                                    str(lt_starts[sample][lt]),
                                    str(iden),
                                    str(sqlr),
                                ]
                            )
                            + "\n"
                        )


        from zol import plotting
        plotting.plot_table_pdf(que_info_table_file, plot_naming_pdf)
        plotting.plot_segments(plot_input_file, plot_result_pdf, height=height, width=width)

    except Exception as e:
        msg = "Issues with plotting overviews of homologous gene cluster segments identified."
        log_object.error(msg)
        sys.stderr.write(msg + "\n")
        util.report_cause(fallback_traceback=True)
        sys.exit(1)


def plot_tree_heatmap(
    homologous_gbk_dir,
    hmm_work_dir,
    species_tree,
    plot_phylo_dir,
    plot_result_pdf,
    log_object,
    height=10,
    width=10,
) -> None:
    """
    Description:
    This function serves to create a phyloheatmap plot when / if a species tree is provided to fai. Takes
    the best bitscore available for a query protein / ortholog group per sample in a locus deemed homologous.
    ********************************************************************************************************************
    Parameters:
    - homologous_gbk_dir: Results directory with homologous gene cluster GenBank files.
    - hmm_work_dir: Directory with itermediate files and information from HMM search for homologous gene cluster
                    instances.
    - species_tree: Species tree in Newick format.
    - plot_phylo_dir: Directory with intermediate files for plotting phyloheatmap.
    - plot_pdf: The final PDF result of the phyloheatmap.
    - log_object: A logging object.
    - height: The height of the plot in inches.
    - width: The width of the plot in inches.
    ********************************************************************************************************************
    """
    try:
        hmm_work_dir = os.path.abspath(hmm_work_dir) + '/'

        heatmap_info_file = plot_phylo_dir + "Heatmap_Track.txt"
        heatmap_rows = []

        sample_final_lts = defaultdict(set)
        for gbk in os.listdir(homologous_gbk_dir):
            sample = ".gbk".join(gbk.split(".gbk")[:-1]).split(
                "_fai-gene-cluster-"
            )[0]
            with open(homologous_gbk_dir + gbk) as ogf:
                for rec in SeqIO.parse(ogf, "genbank"):
                    for feature in rec.features:
                        if feature.type != "CDS":
                            continue
                        lt = None
                        try:
                            lt = feature.qualifiers.get("locus_tag")[0]
                        except Exception as e:
                            sys.stderr.write(
                                f"The GenBank {homologous_gbk_dir + gbk}, cataloging a homologous instance to the query gene cluster had at least one CDS without a locus_tag feature.\n" # type: ignore
                            )
                            util.report_cause()
                            sys.exit(1)
                        sample_final_lts[sample].add(lt)
        if len(sample_final_lts) == 0:
            log_object.warning(
                "Unable to generate phylogenetic heatmap because no gene cluster instances were detected in target genomes."
            )
            sys.stderr.write(
                "Warning: unable to generate phylogenetic heatmap because no gene cluster instances were detected in target genomes.\n"
            )
            return

        t = Tree(species_tree)
        all_samples_in_tree = set([])
        samples_accounted = set([])
        for node in t.traverse("postorder"): # type: ignore
            if node.is_leaf() and node.name != "":
                all_samples_in_tree.add(node.name)

        if (
            len(all_samples_in_tree.intersection(set(sample_final_lts.keys())))
            == 0
        ):
            msg = "Unable to generate phylogenetic heatmap because species tree provided doesn't match target genomes searched against."
            log_object.warning(msg)
            sys.stderr.write(msg + "\n")
            return

        gbk_info_dir = hmm_work_dir + "GeneCluster_Info/"
        example_og = None
        if os.path.isdir(gbk_info_dir):
            for gcs_info in os.listdir(gbk_info_dir):
                if (
                    not gcs_info.endswith(".hg_evalues.txt")
                    or os.path.getsize(gbk_info_dir + gcs_info) == 0
                ):
                    continue
                sample = gcs_info.split(".hg_evalues.txt")[0]
                if not sample in all_samples_in_tree:
                    msg = f"Sample {sample} not in phylogeny or has an unexpected name alteration."
                    log_object.warning(msg)
                    sys.stderr.write(msg + "\n")
                    continue
                samples_accounted.add(sample)
                with open(gbk_info_dir + gcs_info) as ogif:
                    for line in ogif:
                        line = line.strip()
                        ls = line.split("\t")
                        lt = ls[2]
                        if not lt in sample_final_lts[sample]:
                            continue
                        og = ls[3]
                        if og == "background":
                            continue
                        example_og = og
                        bitscore = float(ls[4])
                        heatmap_rows.append(f"{sample}\t{og}\t{bitscore}\n")

        if example_og is None:
            log_object.warning(
                "Unable to generate phylogenetic heatmap because no gene cluster instances were detected in target genomes."
            )
            sys.stderr.write(
                "Warning: unable to generate phylogenetic heatmap because no gene cluster instances were detected in target genomes.\n"
            )
            return

        for sample in sorted(all_samples_in_tree):
            if not sample in samples_accounted:
                heatmap_rows.append(sample + f"\t{example_og}\tNA\n")

        with open(heatmap_info_file, "w") as heatmap_info_file_handle:
            heatmap_info_file_handle.write("label\tquery_prot_id\tbitscore\n")
            heatmap_info_file_handle.writelines(heatmap_rows)

        from zol import plotting
        plotting.plot_phylo_heatmap(species_tree, heatmap_info_file, plot_result_pdf, height=height, width=width)

    except Exception as e:
        log_object.error("Issues with creating phyloheatmap.")
        sys.stderr.write("Issues with creating phyloheatmap.\n")
        util.report_cause(fallback_traceback=True)
        sys.exit(1)


def create_overview_spreadsheet_and_tiny_aai_plot(
    hmm_work_dir,
    protein_to_hg,
    key_hgs,
    spreadsheet_result_file,
    spreadsheet_result_tsv_dir,
    homologous_gbk_dir,
    plot_work_dir,
    tiny_aai_plot_file,
    log_object,
) -> None:
    """
    Description:
    This function serves to create an overview spreadsheet detailing the candidate homologous gene cluster segments
    identified. It will also generate a tiny AAI vs. shared genes plot.
    ********************************************************************************************************************
    Parameters:
    - hmm_work_dir: The workspace where identifyGCInstances() saved extracted GenBanks for gene cluster instances from \
                            target genomes along with homology information to the query gene cluster.
    - protein_to_hg: A dictionary linking query proteins to non-redundant clusters / ortholog groups.
    - key_hgs: Key homolog groups as defined by user.
    - spreadsheet_result_file: The path to the spreadsheet XLSX file with the overview info.
    - spreadsheet_result_tsv_dir: Directory where to write the TSV files that ultimately go into the final spreadsheet.
    - plot_work_dir: The workspace where to write intermediate files used for creating plots.
    - tiny_aai_plot_file: A final plot (PDF) showcasing average amino acid identity of homologous gene clusters identified \
                          to proportion of gene clusters represented.
    - log_object: A logging object.
    ********************************************************************************************************************
    """
    try:
        hmm_work_dir = os.path.abspath(hmm_work_dir) + '/'
        gbk_info_dir = hmm_work_dir + "GeneCluster_Info/"

        hg_prots = defaultdict(set)
        hg_list = []
        hg_headers = []
        for p in protein_to_hg:
            hg_prots[protein_to_hg[p]].add(p)
            hg_list.append(protein_to_hg[p])
            hg = protein_to_hg[p]
            hg_name = hg
            if hg in key_hgs:
                hg_name = "KEY-PROTEIN_" + protein_to_hg[p]
            hg_headers.append(hg_name + "-PID")
            hg_headers.append(hg_name + "-SQR")

        gcgbk_to_corr = defaultdict(lambda: float("nan"))
        for f in os.listdir(gbk_info_dir):
            info_file = gbk_info_dir + f
            if f.endswith(".corr_info.txt") and os.path.getsize(info_file) > 0:
                with open(info_file) as oif:
                    for line in oif:
                        line = line.strip()
                        sample, gc_genbank_file, corr_value = line.split("\t")
                        if corr_value == "irrelevant":
                            gcgbk_to_corr[gc_genbank_file] = float("nan")
                        else:
                            gcgbk_to_corr[gc_genbank_file] = round(
                                float(corr_value), 3
                            )

        tiny_aai_plot_input = plot_work_dir + "Tiny_AAI_Plot_Data.txt"
        tiny_aai_plot_handle = open(tiny_aai_plot_input, "w")
        tiny_aai_plot_handle.write("Sample\tAAI\tProp_Genes_Found\tMean_Syntenic_Correlation\n")

        tot_tsv_file = spreadsheet_result_tsv_dir + "total_gcs.tsv"
        igc_tsv_file = spreadsheet_result_tsv_dir + "individual_gcs.tsv"

        tot_tsv_handle = open(tot_tsv_file, "w")
        igc_tsv_handle = open(igc_tsv_file, "w")

        data_tot_header = [
            "sample",
            "aggregate-bitscore",
            "aai-to-query",
            "mean-sequence-to-query-ratio",
            "proportion-query-genes-found",
            "avg-syntenic-correlation",
            "number-background-genes",
            "number-gene-clusters",
            "copy-counts",
        ] + hg_headers
        data_igc_header = [
            "sample",
            "gene-cluster-id",
            "scaffold",
            "start-coordinate",
            "end-coordinate",
            "aggregate-bitscore",
            "aai-to-query",
            "mean-sequence-to-query-ratio",
            "proportion-query-genes-found",
            "syntenic-correlation",
            "number-background-genes",
            "copy-counts",
        ] + hg_headers

        tot_tsv_handle.write("\t".join(data_tot_header) + "\n")
        igc_tsv_handle.write("\t".join(data_igc_header) + "\n")

        tot_row_count = 0
        igc_row_count = 0

        igc_aggregate_bitscores = []
        aggregate_bitscores = []
        for f in os.listdir(gbk_info_dir):
            info_file = gbk_info_dir + f
            if (
                f.endswith(".hg_evalues.txt")
                and os.path.getsize(info_file) > 0
            ):
                tot_hg_best_hits = defaultdict(lambda: [0.0, [], []])
                tot_hg_hits_copy_counts = defaultdict(int)
                tot_bg_genes = 0
                igc_hg_best_hits = defaultdict(
                    lambda: defaultdict(lambda: [0.0, [], []])
                )
                igc_hg_hits_copy_counts = defaultdict(lambda: defaultdict(int))
                igc_bg_genes = defaultdict(int)
                tot_corr_values = []
                sample = None
                igcs = set([])
                with open(info_file) as oif:
                    for line in oif:
                        line = line.strip("\n")
                        (
                            igc,
                            sample,
                            lt,
                            hg,
                            bitscore,
                            identity,
                            sqlratio,
                            is_key_hg,
                        ) = line.split("\t")
                        bitscore = float(bitscore)
                        igcs.add(igc)
                        tot_corr_values.append(gcgbk_to_corr[igc])
                        if hg == "background":
                            tot_bg_genes += 1
                            igc_bg_genes[igc] += 1
                        else:
                            tot_hg_hits_copy_counts[hg] += 1
                            igc_hg_hits_copy_counts[igc][hg] += 1

                            if tot_hg_best_hits[hg][0] < bitscore: # type: ignore
                                identity = float(identity)
                                sqlratio = float(sqlratio)
                                tot_hg_best_hits[hg] = [
                                    bitscore,
                                    [identity],
                                    [sqlratio],
                                ]
                            elif tot_hg_best_hits[hg][0] == bitscore:
                                identity = float(identity)
                                sqlratio = float(sqlratio)
                                tot_hg_best_hits[hg][1].append(identity) # type: ignore
                                tot_hg_best_hits[hg][2].append(sqlratio) # type: ignore

                            if igc_hg_best_hits[igc][hg][0] < bitscore: # type: ignore
                                identity = float(identity)
                                sqlratio = float(sqlratio)
                                igc_hg_best_hits[igc][hg] = [
                                    bitscore,
                                    [identity],
                                    [sqlratio],
                                ]
                            elif igc_hg_best_hits[igc][hg][0] == bitscore:
                                identity = float(identity)
                                sqlratio = float(sqlratio)
                                igc_hg_best_hits[igc][hg][1].append(identity) # type: ignore
                                igc_hg_best_hits[igc][hg][2].append(sqlratio) # type: ignore

                copy_counts = []
                hg_pid_info = []
                hg_sql_info = []
                hg_pid_sql_pair = []
                aggregate_bitscore = 0.0
                for hg in hg_list:
                    max_pid = 0.0
                    max_sql = 0.0
                    aggregate_bitscore += tot_hg_best_hits[hg][0] # type: ignore
                    for i, pid in enumerate(tot_hg_best_hits[hg][1]): # type: ignore
                        if pid >= max_pid:
                            max_pid = pid
                            max_sql = tot_hg_best_hits[hg][2][i] # type: ignore
                    hg_pid_info.append(max_pid)
                    hg_sql_info.append(max_sql)
                    hg_pid_sql_pair.append(
                        str(round(max_pid, 3)) + "\t" + str(round(max_sql, 3))
                    )
                    copy_counts.append(str(tot_hg_hits_copy_counts[hg]))

                mean_aai = statistics.mean([x for x in hg_pid_info if x > 0.0])
                mean_sql = statistics.mean([x for x in hg_sql_info if x > 0.0])
                prop_hg_found = len(
                    [x for x in hg_pid_info if x > 0.0]
                ) / float(len(hg_list))
                copy_count_string = ",".join(copy_counts)
                avg_corr_value = statistics.mean(tot_corr_values)
                number_gcs = len(igcs)

                tot_row = [
                    str(x)
                    for x in [
                        sample,
                        round(aggregate_bitscore, 3),
                        round(mean_aai, 3),
                        round(mean_sql, 3),
                        round(prop_hg_found, 3),
                        round(avg_corr_value, 3),
                        tot_bg_genes,
                        number_gcs,
                        copy_count_string,
                    ]
                ] + hg_pid_sql_pair
                aggregate_bitscores.append(aggregate_bitscore)
                tot_tsv_handle.write("\t".join(tot_row) + "\n")
                tiny_aai_plot_handle.write(
                    "\t".join(
                        [
                            str(x)
                            for x in [
                                sample,
                                mean_aai,
                                prop_hg_found,
                                avg_corr_value,
                            ]
                        ]
                    )
                    + "\n"
                )
                tot_row_count += 1

                for igc in igcs:
                    igc_copy_counts = []
                    igc_hg_pid_info = []
                    igc_hg_sql_info = []
                    igc_hg_pid_sql_pair = []
                    igc_aggregate_bitscore = 0.0
                    for hg in hg_list:
                        max_pid = 0.0
                        max_sql = 0.0
                        igc_aggregate_bitscore += igc_hg_best_hits[igc][hg][0] # type: ignore
                        for i, pid in enumerate(igc_hg_best_hits[igc][hg][1]): # type: ignore
                            if pid >= max_pid:
                                max_pid = pid
                                max_sql = igc_hg_best_hits[igc][hg][2][i] # type: ignore
                        igc_hg_pid_info.append(max_pid)
                        igc_hg_sql_info.append(max_sql)
                        igc_hg_pid_sql_pair.append(
                            str(round(max_pid, 3))
                            + "\t"
                            + str(round(max_sql, 3))
                        )
                        igc_copy_counts.append(
                            str(igc_hg_hits_copy_counts[igc][hg])
                        )

                    igc_mean_aai = statistics.mean(
                        [x for x in igc_hg_pid_info if x > 0.0]
                    )
                    igc_mean_sql = statistics.mean(
                        [x for x in igc_hg_sql_info if x > 0.0]
                    )
                    igc_prop_hg_found = len(
                        [x for x in igc_hg_pid_info if x > 0.0]
                    ) / float(len(hg_list))
                    copy_count_string = ",".join(igc_copy_counts)

                    igc_path = igc
                    if not os.path.isfile(igc):
                        igc_path = homologous_gbk_dir + igc.split("/")[-1]
                        if not (os.path.isfile(igc_path)):
                            sys.stderr.write(
                                f"Warning: Gene cluster GenBank file not found / confirmed for: {igc}\n"
                            )
                            log_object.warning(
                                f"Gene cluster GenBank file not found / confirmed for: {igc}"
                            )

                    # Extract scaffold and coordinates from GenBank file
                    scaffold_name = "NA"
                    start_coord = "NA"
                    end_coord = "NA"
                    try:
                        with open(igc_path) as igc_handle:
                            for rec in SeqIO.parse(igc_handle, "genbank"):
                                # Get scaffold from record ID
                                scaffold_name = rec.id
                                # Parse coordinates from comment field in annotations
                                comment = None
                                if hasattr(rec, 'annotations') and 'comment' in rec.annotations:
                                    comment = rec.annotations['comment']
                                
                                if comment:
                                    # Expected format: "Extracted region: scaffold:start..end"
                                    if "Extracted region:" in comment:
                                        for line in comment.split('\n'):
                                            if "Extracted region:" in line:
                                                coord_part = line.split("Extracted region:")[1].strip()
                                                if ":" in coord_part and ".." in coord_part:
                                                    parts = coord_part.split(":")
                                                    if len(parts) >= 2:
                                                        coord_range = parts[-1]
                                                        if ".." in coord_range:
                                                            start_str, end_str = coord_range.split("..")
                                                            # Convert to integers if possible
                                                            try:
                                                                start_coord = int(start_str.strip())
                                                                end_coord = int(end_str.strip())
                                                            except ValueError:
                                                                # Keep as "NA" if conversion fails
                                                                pass
                                                break
                                break  # Only process first record
                    except Exception as e:
                        sys.stderr.write(
                            f"Warning: Could not extract coordinates from GenBank file {igc_path}: {str(e)}\n"
                        )

                    # Build row with proper typing - coordinates as numbers when available, strings when NA
                    igc_row = [
                        sample,
                        igc_path,
                        scaffold_name,
                        start_coord,  # Will be int or "NA"
                        end_coord,    # Will be int or "NA"
                        round(igc_aggregate_bitscore, 3),
                        round(igc_mean_aai, 3),
                        round(igc_mean_sql, 3),
                        igc_prop_hg_found,
                        gcgbk_to_corr[igc],
                        igc_bg_genes[igc],
                        copy_count_string,
                    ]
                    # Convert to strings for TSV output
                    igc_row_str = [str(x) for x in igc_row] + igc_hg_pid_sql_pair
                    igc_aggregate_bitscores.append(igc_aggregate_bitscore)
                    igc_tsv_handle.write("\t".join(igc_row_str) + "\n")
                    igc_row_count += 1

        tot_tsv_handle.close()
        igc_tsv_handle.close()
        tiny_aai_plot_handle.close()

        # Generate Excel spreadsheet
        writer = pd.ExcelWriter(spreadsheet_result_file, engine="xlsxwriter")
        workbook = writer.book
        dd_sheet = workbook.add_worksheet("Data Dictionary")
        dd_sheet.write(
            0,
            0,
            'Data Dictionary describing columns of "Overview" spreadsheets can be found below and on fai\'s Wiki page at:',
        )
        dd_sheet.write(
            1,
            0,
            "https://github.com/Kalan-Lab/zol/wiki/2.-more-info-on-fai#explanation-of-report",
        )

        warn_format = workbook.add_format(
            {
                "text_wrap": False,
                "border": 1,
                'font_color': '#FF0000',
                'bold': True
            }
        )

        dd_sheet.write(
            2,
            0,
            "WARNING: If DIAMOND linclust was used to collapse redundancy in the prepTG database (default behavior since v1.6.1), sequence similarity stats reported for some hits are proxied based on their representative proteins.",
            warn_format
        )
        dd_sheet.write(
            3,
            0,
            "To add detail to your methods descriptions, check out additional citations for dependency software at:",
        )
        dd_sheet.write(
            4,
            0,
            "https://github.com/Kalan-Lab/zol/wiki/6.-dependencies",
        )

        wrap_format = workbook.add_format(
            {
                "text_wrap": True,
                "valign": "vcenter",
                "align": "center",
                "border": 1,
            }
        )
        header_format = workbook.add_format(
            {
                "bold": True,
                "text_wrap": True,
                "valign": "top",
                "fg_color": "#FFFFFF",
                "border": 1,
            }
        )

        data_dict_fai = data_dictionary.fai_dd()
        data_dict_fai_df = util.load_table_in_panda_data_frame_from_string(
            data_dict_fai
        )
        worksheet_dd = writer.sheets["Data Dictionary"]
        worksheet_dd.set_column(1, 3, 50)

        for col_num, value in enumerate(data_dict_fai_df.columns.values): # type: ignore
            worksheet_dd.write(6, col_num + 1, value, header_format)

        colnames = ["Column", "Description", "Notes"]
        for index, row in data_dict_fai_df.iterrows(): # type: ignore
            row_ind = int(index) + 7 # type: ignore
            format = wrap_format
            for col_ind in range(0, 3):
                col_name = colnames[col_ind]
                worksheet_dd.write(row_ind, col_ind + 1, row[col_name], format)

        numeric_columns = set(
            [
                "aggregate-bitscore",
                "aai-to-query",
                "mean-sequence-to-query-ratio",
                "proportion-query-genes-found",
                "syntenic-correlation",
                "avg-syntenic-correlation",
                "number-background-genes",
                "number-gene-clusters",
                "start-coordinate",
                "end-coordinate",
            ]
            + hg_headers
        )

        tot_results_df = util.load_table_in_panda_data_frame(
            tot_tsv_file, numeric_columns
        )
        tot_results_df.to_excel( # type: ignore
            writer, sheet_name="Genome Wide - Report", index=False, na_rep="NA"
        )

        igc_results_df = util.load_table_in_panda_data_frame(
            igc_tsv_file, numeric_columns
        )
        igc_results_df.to_excel( # type: ignore
            writer,
            sheet_name="Gene Cluster Instance - Report",
            index=False,
            na_rep="NA",
        )

        tot_worksheet = writer.sheets["Genome Wide - Report"]
        igc_worksheet = writer.sheets["Gene Cluster Instance - Report"]

        tot_worksheet.conditional_format(
            "A1:HA1",
            {
                "type": "cell",
                "criteria": "!=",
                "value": "NA",
                "format": header_format,
            },
        )
        igc_worksheet.conditional_format(
            "A1:HA1",
            {
                "type": "cell",
                "criteria": "!=",
                "value": "NA",
                "format": header_format,
            },
        )

        tot_row_count += 1
        igc_row_count += 1
        if len(aggregate_bitscores) > 0:
            # aggregate bitscore coloring
            tot_worksheet.conditional_format(
                "B2:B" + str(tot_row_count),
                {
                    "type": "2_color_scale",
                    "min_color": "#d9d9d9",
                    "max_color": "#949494",
                    "min_type": "num",
                    "max_type": "num",
                    "min_value": min(aggregate_bitscores),
                    "max_value": max(aggregate_bitscores),
                },
            )
            # aggregate-bitscore is in column F for igc (shifted by 3 due to scaffold + coordinates)
            igc_worksheet.conditional_format(
                "F2:F" + str(igc_row_count),
                {
                    "type": "2_color_scale",
                    "min_color": "#d9d9d9",
                    "max_color": "#949494",
                    "min_type": "num",
                    "max_type": "num",
                    "min_value": min(igc_aggregate_bitscores),
                    "max_value": max(igc_aggregate_bitscores),
                },
            )

            # mean aai
            tot_worksheet.conditional_format(
                "C2:C" + str(tot_row_count),
                {
                    "type": "2_color_scale",
                    "min_color": "#d8eaf0",
                    "max_color": "#83c1d4",
                    "min_value": 0.0,
                    "max_value": 100.0,
                    "min_type": "num",
                    "max_type": "num",
                },
            )
            # aai-to-query is in column G for igc (shifted by 3)
            igc_worksheet.conditional_format(
                "G2:G" + str(igc_row_count),
                {
                    "type": "2_color_scale",
                    "min_color": "#d8eaf0",
                    "max_color": "#83c1d4",
                    "min_value": 0.0,
                    "max_value": 100.0,
                    "min_type": "num",
                    "max_type": "num",
                },
            )

            # mean sequence-to-query ratio
            tot_worksheet.conditional_format(
                "D2:D" + str(tot_row_count),
                {
                    "type": "2_color_scale",
                    "min_color": "#dbd5e8",
                    "max_color": "#afa1cf",
                    "min_value": 0.0,
                    "max_value": 1.0,
                    "min_type": "num",
                    "max_type": "num",
                },
            )
            # mean-sequence-to-query-ratio is in column H for igc (shifted by 3)
            igc_worksheet.conditional_format(
                "H2:H" + str(igc_row_count),
                {
                    "type": "2_color_scale",
                    "min_color": "#dbd5e8",
                    "max_color": "#afa1cf",
                    "min_value": 0.0,
                    "max_value": 1.0,
                    "min_type": "num",
                    "max_type": "num",
                },
            )

            # proportion query genes found
            tot_worksheet.conditional_format(
                "E2:E" + str(tot_row_count),
                {
                    "type": "2_color_scale",
                    "min_color": "#fce8ee",
                    "max_color": "#eb8da9",
                    "min_value": 0.0,
                    "max_value": 1.0,
                    "min_type": "num",
                    "max_type": "num",
                },
            )
            # proportion-query-genes-found is in column I for igc (shifted by 3)
            igc_worksheet.conditional_format(
                "I2:I" + str(igc_row_count),
                {
                    "type": "2_color_scale",
                    "min_color": "#fce8ee",
                    "max_color": "#eb8da9",
                    "min_value": 0.0,
                    "max_value": 1.0,
                    "min_type": "num",
                    "max_type": "num",
                },
            )

            # syntenic-correlation
            tot_worksheet.conditional_format(
                "F2:F" + str(tot_row_count),
                {
                    "type": "2_color_scale",
                    "min_color": "#fcf6eb",
                    "max_color": "#d1c0a1",
                    "min_value": 0.0,
                    "max_value": 1.0,
                    "min_type": "num",
                    "max_type": "num",
                },
            )
            # syntenic-correlation is in column J for igc (shifted by 3)
            igc_worksheet.conditional_format(
                "J2:J" + str(igc_row_count),
                {
                    "type": "2_color_scale",
                    "min_color": "#fcf6eb",
                    "max_color": "#d1c0a1",
                    "min_value": 0.0,
                    "max_value": 1.0,
                    "min_type": "num",
                    "max_type": "num",
                },
            )

            # hg_index for igc needs to account for 3 extra columns (scaffold, start-coordinate, end-coordinate)
            hg_index = 9
            hg_index_igc = 12
            for hg in hg_list:
                columnid_pid = util.determine_column_name_based_on_index(hg_index)
                columnid_sql = util.determine_column_name_based_on_index(
                    hg_index + 1
                )
                columnid_pid_igc = util.determine_column_name_based_on_index(hg_index_igc)
                columnid_sql_igc = util.determine_column_name_based_on_index(
                    hg_index_igc + 1
                )
                hg_index += 2
                hg_index_igc += 2

                # percent identity
                tot_worksheet.conditional_format(
                    columnid_pid + "2:" + columnid_pid + str(tot_row_count), # type: ignore
                    {
                        "type": "2_color_scale",
                        "min_color": "#d8eaf0",
                        "max_color": "#83c1d4",
                        "min_value": 0.0,
                        "max_value": 100.0,
                        "min_type": "num",
                        "max_type": "num",
                    },
                )
                igc_worksheet.conditional_format(
                    columnid_pid_igc + "2:" + columnid_pid_igc + str(igc_row_count), # type: ignore
                    {
                        "type": "2_color_scale",
                        "min_color": "#d8eaf0",
                        "max_color": "#83c1d4",
                        "min_value": 0.0,
                        "max_value": 100.0,
                        "min_type": "num",
                        "max_type": "num",
                    },
                )

                # sequence-to-query ratio
                tot_worksheet.conditional_format(
                    columnid_sql + "2:" + columnid_sql + str(tot_row_count), # type: ignore
                    {
                        "type": "2_color_scale",
                        "min_color": "#dbd5e8",
                        "max_color": "#afa1cf",
                        "min_value": 0.0,
                        "max_value": 1.0,
                        "min_type": "num",
                        "max_type": "num",
                    },
                )
                igc_worksheet.conditional_format(
                    columnid_sql_igc + "2:" + columnid_sql_igc + str(igc_row_count), # type: ignore
                    {
                        "type": "2_color_scale",
                        "min_color": "#dbd5e8",
                        "max_color": "#afa1cf",
                        "min_value": 0.0,
                        "max_value": 1.0,
                        "min_type": "num",
                        "max_type": "num",
                    },
                )

        else:
            tot_worksheet.write(
                1,
                0,
                "No hits to query gene cluster in target genomes found at requested thresholds! Try loosening parameters perhaps.",
            )
            igc_worksheet.write(
                1,
                0,
                "No hits to query gene cluster in target genomes found at requested thresholds! Try loosening parameters perhaps.",
            )

        # Freeze the first row of both sheets
        tot_worksheet.freeze_panes(1, 0)
        igc_worksheet.freeze_panes(1, 0)

        # close workbook
        workbook.close()

        # plot tiny AAI figure

        from zol import plotting
        plotting.plot_tiny_aai(tiny_aai_plot_input, tiny_aai_plot_file)

    except Exception:
        msg = "Issues with producing spreadsheet / AAI plot overviews of homologous gene cluster segments identified."
        log_object.error(msg)
        sys.stderr.write(msg + "\n")
        util.report_cause(fallback_traceback=True)
        sys.exit(1)

