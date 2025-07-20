from Bio.SeqFeature import FeatureLocation
from collections import defaultdict
from ete3 import Tree
from operator import itemgetter
import asyncio
import decimal
import gzip
import itertools
import logging
import math
import multiprocessing
import os
from typing import Dict, List, Optional, Union, Any, Tuple, Generator, TypedDict
import pickle
import resource
import shutil
import statistics
import subprocess
import sys
import traceback

import pandas as pd
import pandas as pd
import warnings
import warnings
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from scipy import stats
import aiofile
import aiohttp
import copy
import importlib.metadata
os.environ['KMP_WARNINGS'] = 'off'
import numpy as np
import pandas as pd
import pyhmmer
import tqdm
from zol import fai, zol


def assess_multiprocessing_mode() -> str:
    """
    Assess the current multiprocessing start method.
    
    Returns:
        str: The current multiprocessing start method ('spawn', 'fork', 'forkserver', or 'unknown')
    """
    try:
        import multiprocessing
        start_method = multiprocessing.get_start_method()
        return start_method
    except Exception as e:
        # Fallback method to infer from the platform
        try:
            import platform
            if platform.system() == 'Windows':
                return 'spawn'
            elif platform.system() == 'Darwin':  # macOS
                return 'fork'
            else:  # Linux and others
                return 'fork'
        except:
            return "unknown"


def set_multiprocessing_start_method(method: str = "fork") -> bool:
    """
    Set the multiprocessing start method.
    
    Args:
        method: The start method to use ('spawn', 'fork', 'forkserver')
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        import multiprocessing
        multiprocessing.set_start_method(method, force=True)
        return True
    except Exception as e:
        return False


def log_multiprocessing_info(log_object) -> None:
    """
    Log multiprocessing information for debugging purposes.
    
    Args:
        log_object: A logging object
    """
    current_method = assess_multiprocessing_mode()
    log_object.info(f"Current multiprocessing start method: {current_method}")
    
    import platform
    log_object.info(f"Platform: {platform.system()}")
    log_object.info(f"Python version: {platform.python_version()}")


# TypedDict definitions for better type safety
class BestHitInfo(TypedDict):
    hg_list: List[str]
    best_bitscore: float
    identity_list: List[float]
    sql_ratio_list: List[float]

class SampleBestHitInfo(TypedDict):
    bitscore: float
    identity: float
    sql_ratio: float
    locus_tag: str

# TypedDict definitions for process_diamond_for_g_cto_ribo_ratio function
class GciQueryTopHitInfo(TypedDict):
    hits: List[str]
    best_bitscore: float
    identity_list: List[float]

class GciTopHitInfo(TypedDict):
    query: str
    bitscore: float
    max_hit: str
    max_identity: float

class GciDataInfo(TypedDict):
    gci: str
    ribo_aai: float
    gc_aai: float
    ribo_prot_prop: float
    gc_prot_prop: float

try:
    package_name = "zol"
    package_version = str(importlib.metadata.version(package_name))
except importlib.metadata.PackageNotFoundError:
    package_version = "NA"

valid_alleles = set(["A", "C", "G", "T"])


def process_location_string(location_string) -> None:
    """
    Description:
    Function to process a location string from a GenBank file.
    ********************************************************************************************************************
    Parameters:
    - location_string: The location string to process.
    ********************************************************************************************************************
    Returns:
    - A list of items: [start, end, direction, all_coords]
    ********************************************************************************************************************
    """
    try:
        all_starts = []
        all_ends = []
        all_coords = []
        direction = None
        if "order" in str(location_string):
            all_directions = [] # type: ignore
            for exon_coord in location_string[6:-1].split(", "):
                start = (
                    min(
                        [
                            int(x.strip(">").strip("<"))
                            for x in exon_coord[1:].split("]")[0].split(":")
                        ]
                    )
                    + 1
                )
                end = max(
                    [
                        int(x.strip(">").strip("<"))
                        for x in exon_coord[1:].split("]")[0].split(":")
                    ]
                )
                direction = exon_coord.split("(")[1].split(")")[0]
                all_starts.append(start)
                all_ends.append(end)
                all_directions.append(direction)
                all_coords.append([start, end, direction])
            start = min(all_starts)
            end = max(all_ends)
            assert len(set(all_directions)) == 1
            direction = all_directions[0]
        elif "join" in location_string:
            all_directions = []
            for exon_coord in location_string[5:-1].split(", "):
                start = (
                    min(
                        [
                            int(x.strip(">").strip("<"))
                            for x in exon_coord[1:].split("]")[0].split(":")
                        ]
                    )
                    + 1
                )
                end = max(
                    [
                        int(x.strip(">").strip("<"))
                        for x in exon_coord[1:].split("]")[0].split(":")
                    ]
                )
                direction = exon_coord.split("(")[1].split(")")[0]
                all_starts.append(start)
                all_ends.append(end)
                all_directions.append(direction)
                all_coords.append([start, end, direction])
            start = min(all_starts)
            end = max(all_ends)
            assert len(set(all_directions)) == 1
            direction = all_directions[0]
        elif "{" not in location_string:
            start = (
                min(
                    [
                        int(x.strip(">").strip("<"))
                        for x in location_string[1:].split("]")[0].split(":")
                    ]
                )
                + 1
            )
            end = max(
                [
                    int(x.strip(">").strip("<"))
                    for x in location_string[1:].split("]")[0].split(":")
                ]
            )
            direction = location_string.split("(")[1].split(")")[0]
            all_starts.append(start)
            all_ends.append(end)
            all_coords.append([start, end, direction])
        else:
            msg = \
     'Error: There appears to be a location operator that is neither "join" nor "order". This is currently not supported in zol.'
            sys.stderr.write(msg + "\n")
            raise RuntimeError(
                f"Error processing location string {location_string}"
            )

        return [start, end, direction, all_coords] # type: ignore
    except Exception as e:
        sys.stderr.write(traceback.format_exc())
        raise RuntimeError(
            f"Error processing location string {location_string}"
        )

def _download_files(urls, resdir):
    """
    Download files from the given URLs and save them to the specified directory.
    Note, this function was taken from:
    https://gist.github.com/darwing1210/c9ff8e3af8ba832e38e6e6e347d9047a
    ********************************************************************
    Parameters:
    - urls: List of URLs to download.
    - resdir: Directory to save the downloaded files.
    ********************************************************************
    """
    os.makedirs(resdir, exist_ok=True)
    sema = asyncio.BoundedSemaphore(5)

    async def fetch_file(session, url):
        fname = url.split("/")[-1]
        try:
            async with sema:
                async with session.get(url) as resp:
                    try:
                        assert resp.status == 200
                    except Exception as e:
                        sys.stderr.write(f'Issue downloading {url}\n')
                    data = await resp.read()

            async with aiofile.async_open(
                os.path.join(resdir, fname), "wb"
            ) as outfile:
                await outfile.write(data)
        except Exception as e:
            sys.stderr.write(f'Issue downloading {url}\n')
    async def main():
        async with aiohttp.ClientSession() as session:
            tasks = [fetch_file(session, url) for url in urls]
            await asyncio.gather(*tasks)

    loop = asyncio.get_event_loop()
    loop.run_until_complete(main())
    loop.close()

def download_gtdb_genomes(
    taxa_name,
    gtdb_release,
    gtdb_dir,
    gtdb_listing_file,
    taxon_listing_file,
    log_object,
    sanity_check=False,
    automated_download=False,
) -> None:
    """
    Download GTDB genomes from NCBI Genbank.
    **********************************************************
    Parameters:
            - log_object: Logger object for logging messages.
            - taxa_name: Taxa name to search for in GTDB.
            :param gtdb_release: GTDB release version.
            :param genomes: List of genomes to download.
            :param outdir: Output directory for downloaded genomes.
            :param sanity_check: Boolean flag for sanity check.
    """

    msg = "Using axel to download GTDB listing."
    sys.stdout.write(msg + "\n")
    log_object.info(msg)
    axel_cmd = [
        "axel",
        f"https://github.com/raufs/gtdb_gca_to_taxa_mappings/raw/main/GTDB_{gtdb_release}_Information_with_Genome_URLs.txt.gz",
        "-o",
        gtdb_listing_file,
    ]
    run_cmd_via_subprocess(axel_cmd, log_object, check_files=[gtdb_listing_file])

    msg = (
        f"Beginning by assessing which genomic assemblies are available for the taxa {taxa_name} in GTDB {gtdb_release}"
    )
    sys.stdout.write(msg + "\n")
    log_object.info(msg)

    with open(taxon_listing_file, "w") as select_genome_listing_handle:
        genome_url_paths = [] 
        url_file_to_polished_name: Dict[str, Any] = {}
        with gzip.open(gtdb_listing_file, "rt") as ogtdb:
            for i, line in enumerate(ogtdb):
                line = line.strip("\n")
                if i == 0:
                    select_genome_listing_handle.write(line + "\n")
                    continue
                ls = line.split("\t")
                gca, gtdb_genus, gtdb_species, genome_url, version_match = (
                    line.split("\t")
                )
                if gca == "none":
                    continue
                if genome_url == "NA":
                    continue
                if len(taxa_name.split()) == 1:
                    if gtdb_genus == taxa_name:
                        select_genome_listing_handle.write(line + "\n")
                        genome_url_paths.append(genome_url)
                        filename = genome_url.split("/")[-1]
                        url_file_to_polished_name[filename] = (
                            "_".join(gtdb_species.split())
                            + f"_{gca}.fasta.gz"
                        )
                elif len(taxa_name.split()) == 2:
                    if gtdb_species == taxa_name:
                        select_genome_listing_handle.write(line + "\n")
                        genome_url_paths.append(genome_url)
                        filename = genome_url.split("/")[-1]
                        url_file_to_polished_name[filename] = (
                            "_".join(gtdb_species.split())
                            + f"_{gca}.fasta.gz"
                        )

    

    genome_count = len(genome_url_paths)
    if genome_count == 0:
        msg = \
     "Error: no genomes found to belong the genus or species specified in GTDB."
        sys.stderr.write(msg + "\n")
        log_object.info(msg)
        sys.exit(1)
    else:
        if not automated_download:
            try:
                response = input(
                    f"Will be downloading {genome_count} genomes for the taxon {taxa_name}. Note, each\ngenome in compressed FASTA formatting is approximately 1 - 2 MB.\nIf downloading thousands of genomes, this can lead to significant disk space being\nused. Do you wish to continue? (yes / no): "
                )
                if response.lower() != "yes":
                    os.system(
                        "User does NOT want to download genomes for taxa from NCBI. Exiting ..."
                    )
                    sys.exit(1)
            except Exception as e:
                msg = \
     "Error: user did not respond to download genomes for taxa from NCBI. Exiting ..."
                sys.stderr.write(msg + "\n")
                log_object.info(msg)
                sys.exit(1)
        try:
            _download_files(genome_url_paths, gtdb_dir)
        except Exception as e:
            msg = "Error downloading genomes from NCBI. Exiting ..."
            sys.stderr.write(msg + "\n")
            log_object.info(msg)
            sys.exit(1)

        final_genome_count = 0
        for f in os.listdir(gtdb_dir):
            genome_file = gtdb_dir + f
            if not os.path.isfile(genome_file):
                continue
            with gzip.open(genome_file, 'r') as fh:
                try:
                    fh.read(1)
                except gzip.BadGzipFile:
                    msg = f'Warning: genome {genome_file}, file with .gz suffix does not actually appear to be gzipped, skipping.'
                    log_object.warning(msg)
                    sys.stderr.write(msg + "\n")
                    os.remove(genome_file)
                    continue
            if sanity_check:
                try:
                    assert is_fasta(genome_file)
                    final_genome_count += 1
                except AssertionError as e:
                    msg = (
                        f"Warning: genome {f} is not a valid FASTA file. Removing ..."
                    )
                    try:
                        os.remove(genome_file)
                    except Exception as e:
                        pass
                    sys.stderr.write(msg + "\n")
                    log_object.info(msg)
            else:
                final_genome_count += 1
            gca = "_".join(f.split("_")[:2])
            polished_filename = url_file_to_polished_name[f]
            renamed_gfile = gtdb_dir + polished_filename
            os.rename(genome_file, renamed_gfile)

        msg = (
            f'Was able to download {final_genome_count} of {genome_count} genomes belonging to taxa {taxa_name} in GTDB {gtdb_release}.'
        )
        sys.stdout.write(msg + "\n")
        log_object.info(msg)


def print_progress_bar(
    iteration,
    total,
    prefix="",
    suffix="",
    decimals=1,
    length=100,
    fill="█",
    print_end="\r",
) -> None:
    """
        Function from: https://stackoverflow.com/questions/3173320/text-progress-bar-in-terminal-with-block-characters
    Call in a loop to create terminal progress bar
        @params:
        iteration   - Required  : current iteration (Int)
                total       - Required  : total iterations (Int)
        prefix      - Optional  : prefix string (Str)
        suffix      - Optional  : suffix string (Str)
        decimals    - Optional  : positive number of decimals in percent complete (Int)
        length      - Optional  : character length of bar (Int)
        fill        - Optional  : bar fill character (Str)
        printEnd    - Optional  : end character (e.g. "\r", "\r\n") (Str)
    """

    percent = ("{0:." + str(decimals) + "f}").format(
        100 * (iteration / float(total))
    )
    filled_length = int(length * iteration // total)
    bar = fill * filled_length + "-" * (length - filled_length)
    sys.stderr.write(f"\r{prefix} |{bar}| {percent}% {suffix}")
    # Print New Line on Complete
    if iteration == total:
        sys.stderr.write("\n")


def process_diamond_linclust_cluster_file(diamond_linclust_cluster_file, pickle_file, log_object) -> None:
    """
    Description:
    This function parses the DIAMOND linclust cluster file and pickles the dictionary of cluster information.
    ********************************************************************************************************************
    Parameters:
    - diamond_linclust_cluster_file: The path to the DIAMOND linclust cluster file (or the FASTA file used to create the cluster file - to 
                          accomodate for databases created prior to v1.6.0).
    - pickle_file: The path to the pickle file to store the dictionary of cluster information.
    - log_object: A logging object.
    ********************************************************************************************************************
    """
    try:
        rep_prot_to_nonreps = {}
        cluster_counts = {}
        multi_prot_cluster_reps = set([])

        with open(diamond_linclust_cluster_file) as occf:
            for line in occf:
                cluster_id, _ = line.strip().split('\t')
                if cluster_id not in cluster_counts:
                    cluster_counts[cluster_id] = 0
                cluster_counts[cluster_id] += 1

        for cluster_id, count in cluster_counts.items():
            if count > 1:
                multi_prot_cluster_reps.add(cluster_id)

        del cluster_counts

        if len(multi_prot_cluster_reps) == 0:
            with open(pickle_file, "wb") as handle:
                pickle.dump(rep_prot_to_nonreps, handle)
        else:
            rep_prot_to_nonreps_list = {}
            with open(diamond_linclust_cluster_file) as occf:
                for line in occf:
                    cluster_id, cluster_member = line.strip().split('\t')
                    if cluster_id in multi_prot_cluster_reps:
                        if cluster_id not in rep_prot_to_nonreps_list:
                            rep_prot_to_nonreps_list[cluster_id] = [cluster_member]
                        else:
                            rep_prot_to_nonreps_list[cluster_id].append(cluster_member)

            for cluster_id in rep_prot_to_nonreps_list:
                rep_prot_to_nonreps[cluster_id] = tuple(rep_prot_to_nonreps_list[cluster_id])

            with open(pickle_file, "wb") as handle:
                pickle.dump(rep_prot_to_nonreps, handle)

    except Exception as e:
        msg = f"Issues parsing DIAMOND linclust cluster file:\n{diamond_linclust_cluster_file}.\n"
        sys.stderr.write(msg + "\n")
        log_object.info(msg)
        sys.exit(1)



def memory_limit(mem) -> None:
    """
    Description:
    Experimental function to limit memory.
    ********************************************************************************************************************
    Parameters:
    - mem: The memory limit in GB.
    ********************************************************************************************************************
    """
    max_virtual_memory = mem * 1000000000
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (max_virtual_memory, hard))
    print(resource.getrlimit(resource.RLIMIT_AS))


def run_cmd_via_subprocess(
    cmd,
    log_object,
    check_files=[],
    check_directories=[],
    stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
) -> None:
    """
    Description:
    Function to run some command via subprocess.
    ********************************************************************************************************************
    Parameters:
    - cmd: The command as a list.
    - log_object: A logging object.
    - check_files: Files to check the existence of assuming successful run of the command.
    - check_directories: Directories to check the existence of assuming successful run of the command.
    - stdout: Where to have subprocess direct standard output.
    - stderr: Where to have subprocess direct standard errorr.
    ********************************************************************************************************************
    """
    log_object.info(f"Running {' '.join(cmd)}")
    try:
        subprocess.call(
            " ".join(cmd),
            shell=True,
            stdout=stdout,
            stderr=stderr,
            executable="/bin/bash",
        )
        for cf in check_files:
            assert os.path.isfile(cf)
        for cd in check_directories:
            assert os.path.isdir(cd)
        log_object.info(f"Successfully ran: {' '.join(cmd)}")
    except Exception as e:
        log_object.error(f"Had an issue running: {' '.join(cmd)}")
        log_object.error(traceback.format_exc())
        raise RuntimeError(f"Had an issue running: {' '.join(cmd)}")


def clean_up_sample_name(original_name) -> None:
    """
    Description:
    Function to clean up sample names for troublesome characters that makes unix based file creation tricky.
    ********************************************************************************************************************
    Parameters:
    - original_name: The original name of the sample.
    ********************************************************************************************************************
    """
    return (
        original_name.replace("#", "")
        .replace("*", "_")
        .replace(":", "_")
        .replace(";", "_")
        .replace(" ", "_")
        .replace(":", "_")
        .replace("|", "_")
        .replace('"', "_")
        .replace("'", "_")
        .replace("=", "_")
        .replace("-", "_")
        .replace("(", "")
        .replace(")", "")
        .replace("/", "")
        .replace("\\", "")
        .replace("[", "")
        .replace("]", "")
        .replace(",", "")
    )


def read_in_annotation_files_for_expanded_sample_set(
    expansion_listing_file, full_dir, log_object=None
) -> Dict[str, Dict[str, str]]:
    """
    Description:
    Function to read in GenBank paths from expansion listing file and load into dictionary with keys corresponding to
    sample IDs.
    ********************************************************************************************************************
    Parameters:
    - expansion_listing_file: A tab-delimited file with two columns: (1) sample ID (2) GenBank file name. \
    - full_dir: The path to where target genome GenBanks are stored.
    - log_object: A logging object.
    ********************************************************************************************************************
    Returns:
    - sample_annotation_data: A dictionary of dictionaries with primary keys as sample names and the secondary key as
                              "genbank" with final values being paths to the corresponding GenBank file for a sample
                              target genome.
    ********************************************************************************************************************
    """
    sample_annotation_data = defaultdict(dict)
    try:
        with open(expansion_listing_file) as oalf:
            for line in oalf:
                line = line.strip()
                sample, genbank = line.split("\t")
                genbank = full_dir + genbank.split("/")[-1]
                sample = clean_up_sample_name(sample)
                try:
                    assert os.path.isfile(genbank) and os.path.isfile(genbank)
                    sample_annotation_data[sample]["genbank"] = genbank
                except Exception as e:
                    if log_object:
                        log_object.warning(
                            f"Ignoring sample {sample}, because at least one of two annotation files does not seem to exist."
                        )
                    else:
                        sys.stderr.write(
                            f"Ignoring sample {sample}, because at least one of two annotation files does not seem to exist.\n"
                        )
        assert len(sample_annotation_data) >= 1
        return sample_annotation_data
    except Exception as e:
        if log_object:
            log_object.error(
                "Input file listing the location of annotation files for samples leads to incorrect paths or something else went wrong with processing of it. Exiting now ..."
            )
            log_object.error(traceback.format_exc()) # type: ignore
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def extract_scaffold_from_gzipped_genbank(filename, scaffold_name) -> Optional[SeqRecord]:
    """
    Description:
    Note, this is an AI generated function. 
    Extracts a specific scaffold record from a gzipped GenBank file by its name.
    ********************************************************************************************************************
    Parameters:
    - filename: Path to the gzipped GenBank file (e.g., "genome.gb.gz").
    - scaffold_name: The name (ID) of the scaffold record to extract.
    ********************************************************************************************************************
    Returns:
        Bio.SeqRecord.SeqRecord or None: The SeqRecord object if found, otherwise None.
    """
    try:
        if filename.endswith(".gz"):
            with gzip.open(filename, "rt") as handle:  # Open the gzipped file in text mode
                for record in SeqIO.parse(handle, "genbank"):
                    if record.id == scaffold_name:
                        return record
        else:
            with open(filename) as handle:  # Open the gzipped file in text mode
                for record in SeqIO.parse(handle, "genbank"):
                    if record.id == scaffold_name:
                        return record
        return None  # Scaffold not found
    except Exception as e:
        print(f"Error: {e}")
        return None


def create_genbank(
    full_genbank_file, new_genbank_file, scaffold, start_coord, end_coord
) -> None:
    """
    Description:
    Function to prune full genome-sized GenBank for only features in BGC of interest.
    ********************************************************************************************************************
    Parameters:
    - full_genbank_file: GenBank file for full genome.
    - new_genbank_file: Path to gene cluster specific GenBank to be created.
    - scaffold: Scaffold identifier.
    - start_coord: Start coordinate.
    - end_coord: End coordinate.
    ********************************************************************************************************************
    """
    try:
        with open(new_genbank_file, "w") as ngf_handle:
            pruned_coords = set(range(start_coord, end_coord + 1))
            rec = extract_scaffold_from_gzipped_genbank(full_genbank_file, scaffold)
            if rec is None:
                raise RuntimeError(f"Scaffold {scaffold} not found in {full_genbank_file}")
            original_seq = str(rec.seq) # type: ignore
            filtered_seq = ""
            start_coord = max(start_coord, 1)
            if end_coord >= len(original_seq):
                filtered_seq = original_seq[start_coord - 1 :]
            else:
                filtered_seq = original_seq[start_coord - 1 : end_coord]

            new_seq_object = Seq(filtered_seq)
            updated_rec = copy.deepcopy(rec)
            updated_rec.seq = new_seq_object # type: ignore

            updated_features = []
            for feature in rec.features: # type: ignore
                start, end, direction, all_coords = process_location_string(
                    str(feature.location) # type: ignore
                )

                feature_coords = set(range(start, end + 1))
                edgy_feat = "False"
                if (start <= 2000) or (end + 1 >= len(original_seq) - 2000):
                    edgy_feat = "True"
                part_of_cds_hanging = False
                if len(feature_coords.intersection(pruned_coords)) > 0:
                    fls = []
                    for sc, ec, dc in all_coords:
                        exon_coord = set(range(sc, ec + 1))
                        if len(exon_coord.intersection(pruned_coords)) == 0:
                            part_of_cds_hanging = True
                            continue
                        updated_start = sc - start_coord + 1
                        updated_end = ec - start_coord + 1
                        if ec > end_coord:
                            # note overlapping genes in prokaryotes are possible so avoid proteins that overlap
                            # with boundary proteins found by the HMM.
                            if feature.type == "CDS":
                                part_of_cds_hanging = True
                                continue
                            else:
                                updated_end = (
                                    end_coord - start_coord + 1
                                )  # ; flag1 = True
                        if sc < start_coord:
                            if feature.type == "CDS":
                                part_of_cds_hanging = True
                                continue
                            else:
                                updated_start = 1  # ; flag2 = True
                        strand = 1
                        if dc == "-":
                            strand = -1
                        fls.append(
                            FeatureLocation(
                                updated_start - 1, updated_end, strand=strand
                            )
                        )
                    if len(fls) > 0 and not part_of_cds_hanging:
                        updated_location = fls[0]
                        if len(fls) > 1:
                            updated_location = sum(fls)
                        feature.location = updated_location # type: ignore
                        feature.qualifiers["near_scaffold_edge"] = edgy_feat
                        updated_features.append(feature)
            updated_rec.features = updated_features # type: ignore
            SeqIO.write(updated_rec, ngf_handle, "genbank") # type: ignore

    except Exception as e:
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def multi_process(input) -> None:
    """
    Description:
    This is a generalizable function to be used with multiprocessing to parallelize list of commands. Inputs should
    correspond to space separated command (as list), with last item in list corresponding to a logging object handle for \
    logging progress.
    ********************************************************************************************************************
    Parameters:
    - input: A list corresponding to a command to run with the last item in the list corresponding to a logging object
             for the function.
    ********************************************************************************************************************
    """
    input_cmd = input[:-1]
    log_object = input[-1]
    log_object.info(f"Running the following command: {' '.join(input_cmd)}")
    try:
        subprocess.call(
            " ".join(input_cmd),
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            executable="/bin/bash",
        )
        log_object.info(f"Successfully ran: {' '.join(input_cmd)}")
    except Exception as e:
        log_object.error(f"Had an issue running: {' '.join(input_cmd)}")
        sys.stderr.write(f"Had an issue running: {' '.join(input_cmd)}")
        log_object.error(e)
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def setup_ready_directory(directories, delete_if_exist=False) -> None:
    """
    Description:
    This is a generalizable function to create directories.
    ********************************************************************************************************************
    Parameters:
    - dictionaries: A list of paths to directories to create or recreate (after removing). \
    ********************************************************************************************************************
    """
    try:
        assert type(directories) is list
        for d in directories:
            if os.path.isdir(d):
                if delete_if_exist:
                    response = input(
                        f"The directory {d}\nalready exists, so it will be deleted it and recreated. Do you wish to proceed? (yes / no): "
                    )
                    if response.lower() != "yes":
                        sys.stderr.write(
                            f"Deletion of directory {d} was not requested!\n"
                        )
                        return
                    os.system(f"rm -rf {d}")
                    sys.stderr.write(
                        f"Warning: directory {d} was deleted and recreated!\n"
                    )
                    os.system(f"mkdir {d}")
                else:
                    sys.stderr.write(
                        f"Warning: directory {d} exists! Overwriting\n"
                    )
            else:
                os.system(f"mkdir {d}")
    except Exception as e:
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def is_fasta(fasta) -> bool:
    """
    Description:
    Function to validate if a file is actually a FASTA file.
    ********************************************************************************************************************
    Parameters:
    - fasta: A file that should be in FASTA format.
    ********************************************************************************************************************
    Returns:
    - True or False statement depending on whether file is in FASTA format.
    ********************************************************************************************************************
    """
    warnings.filterwarnings("ignore")
    try:
        recs = 0
        if fasta.endswith(".gz"):
            with gzip.open(fasta, "rt") as ogf:
                for rec in SeqIO.parse(ogf, "fasta"):
                    recs += 1
                    break
        else:
            with open(fasta) as of:
                for rec in SeqIO.parse(of, "fasta"):
                    recs += 1
                    break
        if recs > 0:
            return True
        else:
            return False
    except Exception as e:
        return False


def is_genbank(gbk, check_for_cds=False) -> bool:
    """
    Description:
    Function to validate if a file is actually a GenBank file.
    ********************************************************************************************************************
    Parameters:
    - gbk: A file that should be in GenBank format.
    - check_for_cds: Whether to also check that the GenBank contains CDS features.
    ********************************************************************************************************************
    Returns:
    - True or False statement depending on whether file is in GenBank format.
    ********************************************************************************************************************
    """
    warnings.filterwarnings("ignore")
    try:
        recs = 0
        cds_flag = False
        assert (
            gbk.endswith(".gbk")
            or gbk.endswith(".gbff")
            or gbk.endswith(".gbk.gz")
            or gbk.endswith(".gbff.gz")
            or gbk.endswith(".gb")
            or gbk.endswith(".gb.gz")
            or gbk.endswith("genbank")
            or gbk.endswith(".genbank.gz")
        )
        if gbk.endswith(".gz"):
            with gzip.open(gbk, "rt") as ogf:
                for rec in SeqIO.parse(ogf, "genbank"):
                    if check_for_cds:
                        for feature in rec.features:
                            if feature.type == "CDS":
                                cds_flag = True
                    if not check_for_cds or cds_flag:
                        recs += 1
                        break
        else:
            with open(gbk) as ogf:
                for rec in SeqIO.parse(ogf, "genbank"):
                    if check_for_cds:
                        for feature in rec.features:
                            if feature.type == "CDS":
                                cds_flag = True
                    if not check_for_cds or cds_flag:
                        recs += 1
                        break
        if recs > 0:
            return True
        else:
            return False
    except Exception as e:
        return False


def check_valid_genbank(
    gbk,
    quality_assessment=False,
    draft_assessment=False,
    use_either_lt_or_pi=False,
    require_translation=False,
) -> bool:
    """
    Description:
    Function to check whether gene cluster GenBanks provided to zol as input meets the criteria requested by the user.
    ********************************************************************************************************************
    Parameters:
    - gbk: The path to the GenBank file.
    - quality_assessment: Whether to check that most bases in the gene cluster are non-ambiguous.
    - draft_assessment: Whether to check that gene cluster does not lie on edge of the scaffold.
    - use_either_lt_or_pi: Whether protein_id is acceptable to use if locus_tag unavailable (currently for usage in fai,
                           not zol).
    ********************************************************************************************************************
    Returns:
    - True or False statement depending on whether file meets criteria for inclusion in zol analysis.
    ********************************************************************************************************************
    """
    try:
        number_of_cds = 0
        lt_has_comma = False
        lt_is_none = False
        seqs = ""
        recs = 0
        edgy_cds = False
        with open(gbk) as ogbk:
            for rec in SeqIO.parse(ogbk, "genbank"):
                for feature in rec.features:
                    if feature.type == "CDS":
                        number_of_cds += 1
                        lt = None
                        pi = None
                        try:
                            lt = feature.qualifiers.get("locus_tag")[0]
                        except Exception as e:
                            pass
                        try:
                            pi = feature.qualifiers.get("protein_id")[0]
                        except Exception as e:
                            pass
                        if use_either_lt_or_pi:
                            if lt == None and pi != None:
                                lt = pi
                        if lt == None:
                            lt_is_none = True
                        if lt and "," in lt:
                            lt_has_comma = True

                        if require_translation:
                            try:
                                protein_seq = feature.qualifiers.get(
                                    "translation"
                                )[0]
                            except Exception as e:
                                return False

                        try:
                            edgy_cds = (
                                feature.qualifiers.get("near_scaffold_edge")[0]
                                == "True"
                            )
                            # print(edgy_cds)
                        except Exception as e:
                            pass
                seqs += str(rec.seq)
                recs += 1
        prop_missing = sum(
            [1 for bp in seqs if not bp in set(["A", "C", "G", "T"])]
        ) / len(seqs)
        if number_of_cds > 0 and not lt_is_none and not lt_has_comma:
            if (quality_assessment and prop_missing >= 0.1) or (
                draft_assessment and (recs > 1 or edgy_cds)
            ):
                return False
            else:
                return True
        else:
            return False
    except Exception as e:
        return False


def convert_genbank_to_cds_prots_fasta(
    gbk, protein, log_object, use_either_lt_or_pi=False
) -> None:
    """
    Description:
    This function extracts protein sequences for CDS features from a GenBank.
    ********************************************************************************************************************
    Parameters:
    - gbk: The path to the GenBank file.
    - protein: The path to the file to which to write the protein sequences to in FASTA format.
    - log_object: A logging object.
    - use_either_lt_or_pi: Whether protein_id is acceptable to use if locus_tag unavailable.
    ********************************************************************************************************************
    """
    try:
        with open(protein, "w") as prot_handle:
            with open(gbk) as ogbk:
                for rec in SeqIO.parse(ogbk, "genbank"):
                    for feature in rec.features:
                        if feature.type != "CDS":
                            continue
                        lt = None
                        pi = None
                        try:
                            lt = feature.qualifiers.get("locus_tag")[0]
                        except Exception as e:
                            pass
                        try:
                            pi = feature.qualifiers.get("protein_id")[0]
                        except Exception as e:
                            pass
                        if use_either_lt_or_pi:
                            if lt == None and pi != None:
                                lt = pi
                        prot_seq = feature.qualifiers.get("translation")[0]
                        prot_handle.write(f">{lt}\n" + str(prot_seq) + "\n")
        
    except Exception as e:
        sys.stderr.write(
            "Difficulties in processing input GenBank and converting to protein fasta. Please make sure translation and locus_tag fields are available for GenBank!\n"
        )
        log_object.error( #
            "Difficulties in processing input GenBank and converting to protein fasta. Please make sure translation and locus_tag fields are available for GenBank!"
        )
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def parse_genbank_for_cds_proteins_and_dna(
    gbk, log_object, allow_edge_cds=True, feature_type="CDS"
) -> List[Dict[str, Any]]:
    """
    Description:
    This function parses GenBank for CDS protein and nucleotide sequences.
    ********************************************************************************************************************
    Parameters:
    - gbk: Path to the GenBank file.
    - log_object: A logging object.
    - allow_edge_cds: Whether to regard CDS features near scaffold edges.
    - feature_type: The focal type of feature. Default is CDS.
    ********************************************************************************************************************
    Returns:
    - A list which can be expanded to the following dictionaries:
            - proteins: A dictionary mapping locus tag identifiers to protein sequences.
            - nucleotides: A dictionary mapping locus tag identifiers to nucleotide sequences.
            - upstream_regions: A dictionary mapping locus tag identifiers to upstream region sequences.
    ********************************************************************************************************************
    """
    try:
        proteins: Dict[str, Any] = {}
        nucleotides: Dict[str, Any] = {}
        upstream_regions: Dict[str, Any] = {}
        with open(gbk) as ogbk:
            for rec in SeqIO.parse(ogbk, "genbank"):
                full_sequence = str(rec.seq).upper()
                for feature in rec.features:
                    if feature.type != feature_type:
                        continue
                    lt = feature.qualifiers.get("locus_tag")[0]

                    start, end, direction, all_coords = process_location_string(
                        str(feature.location) # type: ignore
                    )

                    min_sc = start
                    max_ec = end
                    nucl_seq = ""
                    for sc, ec, dc in sorted(
                        all_coords, key=itemgetter(0), reverse=False
                    ):
                        if ec >= len(full_sequence):
                            nucl_seq += full_sequence[sc - 1 :]
                        else:
                            nucl_seq += full_sequence[sc - 1 : ec]

                    upstream_region = None
                    if direction == "-":
                        nucl_seq = str(Seq(nucl_seq).reverse_complement())
                        if ec + 100 >= len(full_sequence):
                            upstream_region = str(
                                Seq(
                                    full_sequence[max_ec : max_ec + 100]
                                ).reverse_complement()
                            )
                    else:
                        if sc - 100 >= 0:
                            upstream_region = str(
                                Seq(full_sequence[min_sc - 101 : min_sc - 1])
                            )

                    final_prot_seq = None
                    final_nucl_seq = None
                    final_upstream_region = None
                    edgy_cds = False

                    try:
                        final_upstream_region = feature.qualifiers.get(
                            "paf_upstream"
                        )[0].replace(" ", "")
                    except Exception as e:
                        final_upstream_region = upstream_region

                    try:
                        final_nucl_seq = feature.qualifiers.get(
                            "paf_nucl_seq"
                        )[0].replace(" ", "")
                    except Exception as e:
                        final_nucl_seq = nucl_seq
                    try:
                        final_prot_seq = feature.qualifiers.get("translation")[
                            0
                        ]
                    except Exception as e:
                        final_prot_seq = str(Seq(nucl_seq).translate())
                    try:
                        edgy_cds = (
                            feature.qualifiers.get("near_scaffold_edge")[0]
                            == "True"
                        )
                    except Exception as e:
                        edgy_cds = False
                    if allow_edge_cds or not edgy_cds:
                        proteins[lt] = final_prot_seq
                        nucleotides[lt] = final_nucl_seq
                        if upstream_region != None:
                            upstream_regions[lt] = final_upstream_region

        return [proteins, nucleotides, upstream_regions]
    except Exception as e:
        sys.stderr.write(f"Issues with parsing the GenBank {gbk}\n")
        log_object.error(f"Issues with parsing the GenBank {gbk}")
        sys.stderr.write(str(e) + "\n")
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def get_version() -> str:
    """
    Description:
    Parses the version of the zol suite from the setup.py program.
    ********************************************************************************************************************
    """
    return package_version 


def default_to_regular(d) -> Any:
    """
    Description:
    Function to convert defaultdict of defaultdict to dict of dict
    Taken from Martijn Pieters response in StackOverflow:
    https://stackoverflow.com/questions/26496831/how-to-convert-defaultdict-of-defaultdicts-of-defaultdicts-to-dict-of-dicts-o
    ********************************************************************************************************************
    Parameters:
    - d: Defaultdict dictionary.
    ********************************************************************************************************************
    Returns:
    - A regular old dictionary.
    ******************************************************************************************************************
    """
    if isinstance(d, defaultdict):
        d = {k: default_to_regular(v) for k, v in d.items()}
    return d


def create_logger_object(log_file) -> logging.Logger:
    """
    Description:
    This function creates a logging object.
    ********************************************************************************************************************
    Parameters:
    - log_file: Path to file to which to write logging.
    ********************************************************************************************************************
    Returns:
    - logger: A logging object.
    ********************************************************************************************************************
    """

    logger = logging.getLogger("task_logger")
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        "%(asctime)s - %(levelname)s - %(message)s", f"%Y-%m-%d %H:%M:%S"
    )
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    return logger


def close_logger_object(log_object) -> None:
    """
    Description:
    This function closes a logging object.
    ********************************************************************************************************************
    Parameters:
    - log_object: A logging object.
    ********************************************************************************************************************
    """

    handlers = log_object.handlers[:]
    for handler in handlers:
        handler.close()
        log_object.removeHandler(handler)


def log_parameters_to_file(parameter_file, parameter_names, parameter_values) -> None:
    """
    Description:
    This function serves to create a parameters input file for major programs, e.g. fai and zol.
    ********************************************************************************************************************
    Parameters:
    - parameter_file: The path to the file where to write parameter information to. Will overwrite each time.
    - parameter_names: A list containing the parameter names.
    - parameter_values: A list in the same order as parameter_names which contains the respective arguments provided.
    ********************************************************************************************************************
    """
    with open(parameter_file, "w") as parameter_handle:
        for i, pv in enumerate(parameter_values):
            pn = parameter_names[i]
            parameter_handle.write(pn + ": " + str(pv) + "\n")
        
def is_integer(x) -> bool:
    """
    Description:
    This function checks whether the input variable corresponds to an integer.
    ********************************************************************************************************************
    Parameters:
    - x: Input variable.
    ********************************************************************************************************************
    Returns:
    - True or False statement depending on whether input variable is an integer.
    ********************************************************************************************************************
    """
    try:
        x = int(x)
        return True
    except Exception as e:
        return False


def is_numeric(x) -> bool:
    """
    Description:
    This function checks whether the input variable is numeric.
    ********************************************************************************************************************
    Parameters:
    - x: Input variable.
    ********************************************************************************************************************
    Returns:
    - True or False statement depending on whether input variable is numeric.
    ********************************************************************************************************************
    """
    try:
        x = float(x)
        return True
    except Exception as e:
        return False


def cast_to_numeric(x) -> Any:
    """
    Description:
    This function attempts to cast a variable into a float. A special exception is whether "< 3 segregating sites!" is
    the value of the variable, which will simply be retained as a string.
    ********************************************************************************************************************
    Parameters:
    - x: Input variable.
    ********************************************************************************************************************
    Returns:
    - A float casting of the variable's value if numeric or "nan" if not.
    ********************************************************************************************************************
    """
    try:
        if x == "< 3 segregating sites!":
            return x
        else:
            x = float(x)
            return x
    except Exception as e:
        return float("nan")


def check_core_homolog_groups_exist(ortho_matrix_file) -> bool:
    """
    Description:
    This function checks whether at least one core ortholog group exists within an ortholog group by sample matrix.
    ********************************************************************************************************************
    Parameters:
    - orthogroup_matrix_file: The ortholog group vs sample matrix file, where cells correspond to locus tag identifiers.
    ********************************************************************************************************************
    Returns:
    - True or False depending on whether a core ortholog group is found.
    ********************************************************************************************************************
    """

    try:
        core_hgs = set([])
        with open(ortho_matrix_file) as omf:
            for i, line in enumerate(omf):
                if i == 0:
                    continue
                line = line.rstrip("\n")
                ls = line.split("\t")
                hg = ls[0]
                sample_count = 0
                for lts in ls[1:]:
                    if not lts.strip() == "":
                        sample_count += 1
                if sample_count / float(len(ls[1:])) == 1.0:
                    core_hgs.add(hg)
        assert len(core_hgs) != 0
        return True
    except Exception as e:
        return False


def process_genomes_using_miniprot(
    reference_proteome,
    sample_genomes,
    additional_miniprot_outdir,
    additional_proteomes_directory,
    additional_genbanks_directory,
    log_object,
    threads=1,
    locus_tag_length=3,
) -> None:
    """
    Description:
    This function oversees processing of input genomes to create proteome and GenBank files using miniprot.
    ********************************************************************************************************************
    Parameters:
    - reference_proteome: The reference proteome (in FASTA format) to use to map CDS features onto target sample \
                          genomes.
    - sample_genomes: A dictionary mapping sample identifiers to the path of their genomes in FASTA format.
    - additional_miniprot_outdir: Workspace where miniprot (intermediate) results should be written to directly. \
    - additional_proteomes_directory: Directory where final proteome files (in FASTA format) for target genomes will be \
                                      saved.
    - additional_genbanks_directory: Directory where final GenBank files for target genomes will be saved.
    - log_object: A logging object.
    - threads: The number of threads to use.
    - locus_tag_length: The length of the locus tags to generate.
    ********************************************************************************************************************
    """
    try:
        alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        possible_locustags = sorted(
            list(
                [
                    "".join(list(x))
                    for x in list(
                        itertools.product(alphabet, repeat=locus_tag_length)
                    )
                ]
            )
        )

        miniprot_cmds = []
        for i, sample in enumerate(sorted(sample_genomes)):
            sample_assembly = sample_genomes[sample]
            sample_locus_tag = "".join(list(possible_locustags[i]))

            sample_mp_db = additional_miniprot_outdir + sample + ".mpi"
            sample_mp_gff = additional_miniprot_outdir + sample + ".gff"
            sample_mp_gbk = additional_genbanks_directory + sample + ".gbk"
            sample_mp_faa = additional_proteomes_directory + sample + ".faa"

            miniprot_index_cmd = [
                "miniprot",
                "-t1",
                "-d",
                sample_mp_db,
                sample_assembly,
            ]
            miniprot_run_cmd = [
                "miniprot",
                "--trans",
                "--gff",
                "-I",
                "-t1",
                sample_mp_db,
                reference_proteome,
                ">",
                sample_mp_gff,
            ]
            miniprot_process_cmd = [
                "convertMiniprotGffToGbkAndProt.py",
                "-g",
                sample_mp_gff,
                "-f",
                sample_assembly,
                "-l",
                sample_locus_tag,
                "-og",
                sample_mp_gbk,
                "-op",
                sample_mp_faa,
            ]
            miniprot_cmds.append(
                miniprot_index_cmd
                + [";"]
                + miniprot_run_cmd
                + [";"]
                + miniprot_process_cmd
                + [log_object]
            )

        msg = f"Running miniprot for {len(miniprot_cmds)} genomes"
        log_object.info(msg)
        sys.stdout.write(msg + "\n")

        p = multiprocessing.Pool(threads)
        for _ in tqdm.tqdm(
            p.imap_unordered(multi_process, miniprot_cmds),
            total=len(miniprot_cmds),
        ):
            pass
        p.close()

        gzip_cmds = []
        for sample in sample_genomes:
            try:
                sample_mp_gbk = additional_genbanks_directory + sample + ".gbk"
                gzip_cmds.append(['gzip', sample_mp_gbk, log_object])
                sample_mp_faa = (
                    additional_proteomes_directory + sample + ".faa"
                )
                assert os.path.isfile(sample_mp_gbk) and os.path.isfile(sample_mp_faa)
            except Exception as e:
                sys.stderr.write(
                    f"Unable to validate successful genbank / predicted - proteome creation for sample {sample}"
                )
                sys.stderr.write(traceback.format_exc())
                sys.exit(1)
   
        p = multiprocessing.Pool(threads)
        msg = f"Gzipping {len(gzip_cmds)} genome-wide GenBank files"
        log_object.info(msg)
        sys.stdout.write(msg + "\n")
        for _ in tqdm.tqdm(
            p.imap_unordered(multi_process, gzip_cmds),
            total=len(gzip_cmds),
        ):
            pass
        p.close()
   
    except Exception as e:
        log_object.error(
            "Problem with creating commands for running miniprot or convertMiniprotGffToGbkAndProt.py. Exiting now ..."
        )
        log_object.error(traceback.format_exc())
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)

def process_genomes_using_prodigal(
    sample_genomes,
    prodigal_outdir,
    prodigal_proteomes,
    prodigal_genbanks,
    log_object,
    threads=1,
    locus_tag_length=3,
    gene_calling_method="pyrodigal",
    meta_mode=False,
    avoid_locus_tags=set([]),
):
    """
    Description:
    This function oversees processing of input genomes to create proteome and GenBank files using p(y)rodigal. \
    ********************************************************************************************************************
    Parameters:
    - sample_genomes: A dictionary mapping sample identifiers to the path of their genomes in FASTA format.
    - prodigal_outdir: Workspace where prodigal (intermediate) results should be written to directly. \
    - prodigal_proteomes: Directory where final proteome files (in FASTA format) for target genomes will be saved. \
    - prodigal_genbanks_directory: Directory where final GenBank files for target genomes will be saved.
    - log_object: A logging object.
    - threads: The number of threads to use.
    - locus_tag_length: The length of the locus tags to generate.
    - gene_calling_method: Whether to use pyrodigal (default), prodigal, or prodigal - gv. \
    - meta_mode: Whether to run pyrodigal/prodigal in metagenomics mode.
    - avoid_locus_tags: Whether to avoid using certain locus tags.
    ********************************************************************************************************************
    """
    try:
        alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        possible_locustags = sorted(
            list(
                set(
                    [
                        "".join(list(x))
                        for x in list(
                            itertools.product(
                                alphabet, repeat=locus_tag_length
                            )
                        )
                    ]
                ).difference(avoid_locus_tags)
            )
        )

        prodigal_cmds = []
        for i, sample in enumerate(sorted(sample_genomes)):
            sample_assembly = sample_genomes[sample]
            sample_locus_tag = "".join(list(possible_locustags[i]))

            prodigal_cmd = [
                "runProdigalAndMakeProperGenbank.py",
                "-i",
                sample_assembly,
                "-s",
                sample,
                "-gcm",
                gene_calling_method,
                "-l",
                sample_locus_tag,
                "-o",
                prodigal_outdir,
            ]
            if meta_mode:
                prodigal_cmd += ["-m"]
            prodigal_cmds.append(prodigal_cmd + [log_object])

        msg = f"Running {gene_calling_method} for {len(prodigal_cmds)} genomes"
        log_object.info(msg)
        sys.stdout.write(msg + "\n")

        p = multiprocessing.Pool(threads)
        for _ in tqdm.tqdm(
            p.imap_unordered(multi_process, prodigal_cmds),
            total=len(prodigal_cmds),
        ):
            pass
        p.close()

        gzip_cmds = []
        for sample in sample_genomes:
            try:
                assert os.path.isfile(prodigal_outdir + sample + ".faa") and os.path.isfile(prodigal_outdir + sample + ".gbk")
                os.system(
                    f"mv {prodigal_outdir + sample + '.gbk'} {prodigal_genbanks}"
                )
                gzip_cmds.append(['gzip', prodigal_genbanks + sample + '.gbk', log_object])
                os.system(
                    f"mv {prodigal_outdir + sample + '.faa'} {prodigal_proteomes}"
                )
            except Exception as e:
                sys.stderr.write(
                    f"Unable to validate successful genbank/predicted - proteome creation for sample {sample}\n"
                )
                sys.stderr.write(traceback.format_exc())
                sys.exit(1)
        
        p = multiprocessing.Pool(threads)
        msg = f"Gzipping {len(gzip_cmds)} genome-wide GenBank files"
        log_object.info(msg)
        sys.stdout.write(msg + "\n")
        for _ in tqdm.tqdm(
            p.imap_unordered(multi_process, gzip_cmds),
            total=len(gzip_cmds),
        ):
            pass
        p.close()

    except Exception as e:
        log_object.error(
            "Problem with creating commands for running prodigal via script runProdigalAndMakeProperGenbank.py. Exiting now ..."
        )
        log_object.error(traceback.format_exc())
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def process_genomes_as_genbanks(
    sample_genomes,
    proteomes_directory,
    genbanks_directory,
    gene_name_mapping_outdir,
    error_directory,
    error_summary_file,
    log_object,
    threads=1,
    locus_tag_length=3,
    avoid_locus_tags=set([]),
    rename_locus_tags=False,
    error_no_lt=False,
    error_no_translation=False,
):
    """
    Description:
    This function oversees processing of input genomes as GenBanks with CDS features already available.
    ********************************************************************************************************************
    Parameters:
    - sample_genomes: A dictionary mapping sample identifiers to the path of their genomes in GenBank format with CDS
                      features available.
    - proteomes_directory: Directory where final proteome files (in FASTA format) for target genomes will be saved. \
    - genbanks_directory: Directory where final GenBank files for target genomes will be saved.
    - gene_name_mapping_outdir: Directory where mapping files for original locus tags to new locus tags will be saved.
    - error_directory: Directory where error files should be written. Should already be created!
    - error_summary_file: File where summary of errors should be written.
    - log_object: A logging object.
    - threads: The number of threads to use.
    - locus_tag_length: The length of the locus tags to generate.
    - avoid_locus_tags: Whether to avoid using certain locus tags.
    - rename_locus_tags: Whether to rename locus tags.
    ********************************************************************************************************************
    Returns:
    - sample_genomes_updated: Dictionary mapping sample names to paths of final / processed sample GenBanks.
    ********************************************************************************************************************
    """

    sample_genomes_updated: Dict[str, Any] = {}
    process_cmds = []
    try:
        alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        possible_locustags = sorted(
            list(
                set(
                    [
                        "".join(list(x))
                        for x in list(
                            itertools.product(
                                alphabet, repeat=locus_tag_length
                            )
                        )
                    ]
                ).difference(avoid_locus_tags)
            )
        )

        for i, sample in enumerate(sorted(sample_genomes)):
            sample_locus_tag = possible_locustags[i]
            sample_genbank = sample_genomes[sample]
            process_cmd = [
                "processNCBIGenBank.py",
                "-i",
                sample_genbank,
                "-s",
                sample,
                "-g",
                genbanks_directory,
                "-p",
                proteomes_directory,
                "-n",
                gene_name_mapping_outdir,
                "-e",
                error_directory,
                "-l",
                sample_locus_tag,
            ]

            if error_no_translation:
                process_cmd += ["-ent"]
            if error_no_lt:
                process_cmd += ["-enl"]

            if rename_locus_tags:
                process_cmd += ["-r"]

            process_cmds.append(process_cmd + [log_object])

        msg = (
            f"Attempting to process/reformat {len(process_cmds)} genomes provided as GenBank files"
        )
        log_object.info(msg)
        sys.stdout.write(msg + "\n")

        p = multiprocessing.Pool(threads)
        for _ in tqdm.tqdm(
            p.imap_unordered(multi_process, process_cmds),
            total=len(process_cmds),
        ):
            pass
        p.close()

        error_summary_file_handle = open(error_summary_file, 'w')
        samples_with_noted_errors = set([])
        for f in os.listdir(error_directory):
            if f.endswith('.txt'):
                sample = '.'.join(f.split('.')[:-1])
                sample_err_msg = f'{sample}\t'
                error_found = False
                with open(error_directory + f, 'r') as error_file_handle:
                    for line in error_file_handle:
                        sample_err_msg += line.strip() + ' '
                        error_found = True
                if error_found:
                    samples_with_noted_errors.add(sample)
                    sample_err_msg += '\n'
                    error_summary_file_handle.write(sample_err_msg) # type: ignore

        successfully_processed = 0
        gzip_cmds = []
        for sample in sample_genomes:
            faa_file = proteomes_directory + sample + ".faa"
            gbk_file = genbanks_directory + sample + ".gbk"
            map_file = gene_name_mapping_outdir + sample + ".txt"
            try:
                assert (
                    os.path.isfile(faa_file)
                    and os.path.isfile(gbk_file)
                    and os.path.isfile(map_file)
                )
                assert (
                    os.path.getsize(faa_file) > 0
                    and os.path.getsize(gbk_file) > 0
                    and os.path.getsize(map_file) > 0
                )
                sample_genomes_updated[sample] = (
                    genbanks_directory + sample + ".gbk"
                )
                gzip_cmds.append(["gzip", gbk_file, log_object])
                successfully_processed += 1
            except AssertionError:
                if os.path.isfile(faa_file):
                    os.system("rm -f " + faa_file)
                if os.path.isfile(gbk_file):
                    os.system("rm -f " + gbk_file)
                if os.path.isfile(map_file):
                    os.system("rm -f " + map_file)
                msg = f"Unable to validate successful genbank reformatting / predicted - proteome creation for sample {sample}\n"
                sys.stderr.write(msg)
                if not sample in samples_with_noted_errors:
                    error_summary_file_handle.write(sample + '\t' + msg) 
                
        error_summary_file_handle.close()

        msg = f"Successfully processed {successfully_processed} genomes of {len(process_cmds)} genomes attempted to be processed!\n"
        sys.stdout.write(msg + '\n')
        log_object.info(msg)

        msg = f"Genomes with errors are reported in the following file:\n{error_summary_file}"
        sys.stdout.write(msg + '\n')
        log_object.info(msg)

        p = multiprocessing.Pool(threads)
        msg = f"Gzipping {len(gzip_cmds)} genome-wide GenBank files"
        log_object.info(msg)
        sys.stdout.write(msg + "\n")
        for _ in tqdm.tqdm(
            p.imap_unordered(multi_process, gzip_cmds),
            total=len(gzip_cmds),
        ):
            pass
        p.close()        

        sys.stdout.write(
            f"Successfully processed {successfully_processed} genomes!\n"
        )

    except Exception as e:
        log_object.error(
            "Problem with processing existing Genbanks to (re)create genbanks/proteomes. Exiting now ..." \
        )
        log_object.error(traceback.format_exc())
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)
    return sample_genomes_updated


def determine_genome_format(inputs) -> None:
    """
    Description:
    This function determines whether a target sample genome is provided in GenBank or FASTA format.
    ********************************************************************************************************************
    Parameters:
    - input a list which can be expanded to the following:
            - sample: The sample identifier / name.
            - genome_File: The path to the genome file.
            - format_assess_dir: The directory where the genome type information for the sample will be written.
            - log_object: A logging object.
    ********************************************************************************************************************
    """

    sample, genome_file, format_assess_dir, log_object = inputs
    try:
        gtype = "unknown"
        if is_fasta(genome_file):
            gtype = "fasta"
        if is_genbank(genome_file, check_for_cds=True):
            if gtype == "fasta":
                gtype = "unknown"
            else:
                gtype = "genbank"
        with open(format_assess_dir + sample + ".txt", "w") as sample_res_handle:
            sample_res_handle.write(f"{sample}\t{gtype}\n")
        
    except Exception as e:
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def parse_sample_genomes(
    genome_listing_file,
    format_assess_dir,
    format_predictions_file,
    log_object,
    threads=1,
) -> Tuple[Dict[str, Any], str]:
    """
    Description:
    This function parses the input sample target genomes and determines whether they are all provided in the same format
    and whether everything aligns with expectations.
    ********************************************************************************************************************
    Parameters:
    - genome_listing_file: A tab separated file with two columns: (1) sample name, (2) path to genome file. \
    - format_assess_dir: The directory / workspace where genome format information will be saved.
    - format_predictions_file: The file where to concatenate genome format information.
    - log_object: A logging object.
    - threads: The number of threads to use.
    ********************************************************************************************************************
    Returns:
    - sample_genomes: A dictionary which maps sample names to genome file paths (note, unknown format files will be
                      dropped).
    - format_prediction: The format prediction for genome files.
    ********************************************************************************************************************
    """
    try:
        sample_genomes: Dict[str, Any] = {}
        assess_inputs = []
        with open(genome_listing_file) as oglf:
            for line in oglf:
                line = line.strip()
                ls = line.split("\t")
                sample, genome_file = ls
                assess_inputs.append(
                    [sample, genome_file, format_assess_dir, log_object]
                )
                try:
                    assert os.path.isfile(genome_file)
                except Exception as e:
                    log_object.warning(
                        f"Problem with finding genome file {genome_file} for sample {sample}, skipping"
                    )
                    continue
                if sample in sample_genomes:
                    log_object.warning(
                        f"Skipping genome {genome_file} for sample {sample} because a genome file was already provided for this sample"
                    )
                    continue
                sample_genomes[sample] = genome_file

        msg = (
            f"Attempting to determine format (FASTA or GenBank) for {len(assess_inputs)} genomes"
        )
        log_object.info(msg)
        sys.stdout.write(msg + "\n")

        p = multiprocessing.Pool(threads)
        for _ in tqdm.tqdm(
            p.imap_unordered(determine_genome_format, assess_inputs),
            total=len(assess_inputs),
        ):
            pass
        p.close()

        os.system(
            f"find {format_assess_dir} -maxdepth 1 -type f | xargs cat >> {format_predictions_file}"
        )

        format_prediction = "mixed"
        gtypes = set([])
        with open(format_predictions_file) as ofpf:
            for line in ofpf:
                line = line.strip()
                sample, gtype = line.split("\t")
                if gtype == "unknown":
                    sys.stderr.write(
                        f"unsure about format for genome {sample_genomes[sample]} for sample {sample}, skipping inclusion...\n"
                    )
                    log_object.warning(
                        f"unsure about format for genome {sample_genomes[sample]} for sample {sample}, skipping inclusion..."
                    )
                    del sample_genomes[sample]
                else:
                    gtypes.add(gtype)

        if len(gtypes) == 1:
            format_prediction = list(gtypes)[0]

        return [sample_genomes, format_prediction] # type: ignore

    except Exception as e:
        log_object.error(
            "Problem with creating commands for running Prodigal. Exiting now ..."
        )
        log_object.error(traceback.format_exc())
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def filter_records_near_scaffold_edge(
    gbk, filt_genbank_file, log_object, quality_assessment=False
) -> None:
    """
    Description:
    This function filters specific records in a GenBank if they are near a scaffold edge.
    ********************************************************************************************************************
    Parameters:
    - gbk: The GenBank file.
    - filt_genbank_file: The filtered GenBank file.
    - log_object: A logging object.
    - quality_assessment: Whether to perform quality assessment and drop the GenBank if >10% of nucleotides are
                          ambiguous.
    ********************************************************************************************************************
    """
    try:
        number_of_cds = 0
        seqs = ""
        recs = 0
        recs_with_edgy_cds = set([])
        with open(gbk) as ogbk:
            for rec_it, rec in enumerate(SeqIO.parse(ogbk, "genbank")):
                edgy_cds = False
                for feature in rec.features:
                    if feature.type == "CDS":
                        number_of_cds += 1
                        try:
                            if (
                                feature.qualifiers.get("near_scaffold_edge")[0]
                                == "True"
                            ):
                                recs_with_edgy_cds.add(rec_it)
                                edgy_cds = True
                        except Exception as e:
                            pass
                if not edgy_cds:
                    seqs += str(rec.seq)
                recs += 1

        if len(seqs) == 0:
            return
        prop_missing = sum([1 for bp in seqs if not bp in set(["A", "C", "G", "T"])]) / len(seqs)
        recs_without_edgy_cds = recs - len(recs_with_edgy_cds)
        if (
            number_of_cds > 0
            and (prop_missing <= 0.1 or not quality_assessment)
            and recs_without_edgy_cds > 0
        ):
            with open(filt_genbank_file, "w") as out_handle:
                with open(gbk) as ogbk:
                    for rec_it, rec in enumerate(SeqIO.parse(ogbk, "genbank")):
                        if rec_it in recs_with_edgy_cds:
                            continue
                        SeqIO.write(rec, out_handle, "genbank")
            

    except Exception as e:
        sys.stderr.write(
            f"Issue parsing GenBank {gbk} and CDS locus tag renaming.\n"
        )
        log_object.error(
            f"Issue parsing GenBank {gbk} and CDS locus tag renaming."
        )
        sys.stderr.write(str(e) + "\n")
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def rename_cds_locus_tag(
    gbk,
    lt,
    rn_genbank_file,
    log_object,
    quality_assessment=False,
    draft_assessment=False,
) -> None:
    """
    Description:
    This function renames or creates locus tags in a GenBank.
    ********************************************************************************************************************
    Parameters:
    - gbk: The GenBank file.
    - lt: The new locus tag prefix.
    - rn_genbank_file: The path to the GenBank file to be created with the new locus tag names.
    - log_object: A logging object.
    - quality_assessment: Whether to perform quality assessment and drop the GenBank if >10% of nucleotides are
                          ambiguous.
    - draft_assessment: Whether to perform draft quality assessment and drop the GenBank if there are no records not
                        nearby scaffold edges.
    ********************************************************************************************************************
    """
    try:
        number_of_cds = 0
        seqs = ""
        recs = 0
        recs_with_edgy_cds = set([])
        with open(gbk) as ogbk:
            for rec_it, rec in enumerate(SeqIO.parse(ogbk, "genbank")):
                edgy_cds = False
                for feature in rec.features:
                    if feature.type == "CDS":
                        number_of_cds += 1
                        try:
                            if (
                                feature.qualifiers.get("near_scaffold_edge")[0]
                                == "True"
                            ):
                                recs_with_edgy_cds.add(rec_it)
                                edgy_cds = True
                        except Exception as e:
                            pass
                if not edgy_cds:
                    seqs += str(rec.seq)
                recs += 1
        if len(seqs) == 0:
            return
        prop_missing = sum([1 for bp in seqs if not bp in set(["A", "C", "G", "T"])]) / len(seqs)
        recs_without_edgy_cds = recs - len(recs_with_edgy_cds)
        if (
            number_of_cds > 0
            and (prop_missing <= 0.1 or not quality_assessment)
            and (recs_without_edgy_cds > 0 or not draft_assessment)
        ):
            with open(rn_genbank_file, "w") as out_handle:
                locus_tag_iterator = 1
                with open(gbk) as ogbk:
                    for rec_it, rec in enumerate(SeqIO.parse(ogbk, "genbank")):
                        if rec_it in recs_with_edgy_cds and draft_assessment:
                            continue
                        for feature in rec.features:
                            if feature.type != "CDS":
                                continue
                            new_locus_tag = lt + "_"
                            if locus_tag_iterator < 10:
                                new_locus_tag += "00000" + str(locus_tag_iterator)
                            elif locus_tag_iterator < 100:
                                new_locus_tag += "0000" + str(locus_tag_iterator)
                            elif locus_tag_iterator < 1000:
                                new_locus_tag += "000" + str(locus_tag_iterator)
                            elif locus_tag_iterator < 10000:
                                new_locus_tag += "00" + str(locus_tag_iterator)
                            elif locus_tag_iterator < 100000:
                                new_locus_tag += "0" + str(locus_tag_iterator)
                            else:
                                new_locus_tag += str(locus_tag_iterator)
                            feature.qualifiers["locus_tag"] = new_locus_tag
                            locus_tag_iterator += 1
                        SeqIO.write(rec, out_handle, "genbank")
            
    except Exception as e:
        sys.stderr.write(
            f"Issue parsing GenBank {gbk} and CDS locus tag renaming.\n"
        )
        log_object.error(
            f"Issue parsing GenBank {gbk} and CDS locus tag renaming."
        )
        sys.stderr.write(str(e) + "\n")
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def parse_gbk(
    gbk, prefix, log_object, use_either_lt_or_pi=False, feature_type="CDS"
) -> Dict[str, Any]:
    """
    Description:
    This function parses CDS coordinate information from a GenBank.
    ********************************************************************************************************************
    Parameters:
    - gbk: The GenBank file.
    - prefix: The prefix to append to locus tags (often gene cluster name) to make them use. \
    - log_object: A logging object.
    - use_either_lt_or_pi: Use protein_id qualifier if locus_tag is unavailable for CDS feature.
    ********************************************************************************************************************
    Returns:
    - gc_gene_locations: A dictionary for CDS locations where keys correspond to "prefix|locus_tag" and the values are
                         another dictionary with the keys scaffold, start position, end position.
    ********************************************************************************************************************
    """
    try:
        gc_gene_locations: Dict[str, Any] = {}
        with open(gbk) as ogbk:
            for rec in SeqIO.parse(ogbk, "genbank"):
                for feature in rec.features:
                    if feature.type != feature_type:
                        continue
                    lt = None
                    pi = None
                    try:
                        lt = feature.qualifiers.get("locus_tag")[0]
                    except Exception as e:
                        pass
                    try:
                        pi = feature.qualifiers.get("protein_id")[0]
                    except Exception as e:
                        pass
                    if use_either_lt_or_pi:
                        if lt == None and pi != None:
                            lt = pi

                    start, end, direction, all_coords = process_location_string(
                        str(feature.location) # type: ignore
                    )

                    location = {
                        "scaffold": rec.id,
                        "start": start,
                        "end": end,
                        "direction": direction,
                    }
                    gc_gene_locations[prefix + "|" + lt] = location
        return gc_gene_locations
    except Exception as e:
        sys.stderr.write(f"Issue parsing GenBank {gbk}\n")
        log_object.error(f"Issue parsing GenBank {gbk}") # type: ignore
        sys.stderr.write(str(e) + "\n")
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def determine_possible_lts() -> List:
    """
    Description:
    This function creates a sorted list of possible locus tag prefices of length 4 each.
    ********************************************************************************************************************
    Returns:
    - possible_locustags: A sorted list of possible locus tag prefices of length 4 each.
    ********************************************************************************************************************
    """
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    possible_locustags = sorted(
        list(
            ["".join(list(lt)) for lt in itertools.product(alphabet, repeat=4)]
        )
    )
    return possible_locustags


def gather_annotation_from_dict_for_homolo_group(hg, db, annot_dict) -> str:
    """
    Description:
    This function formats the annotation information for the final zol TSV and XLSX for an ortholog group for a
    particular database.
    ********************************************************************************************************************
    Parameters:
    - hg: The ortholog group identifier.
    - db: The database identifier.
    - annot_dict: The dictionary with annotation information.
    ********************************************************************************************************************
    Returns:
    - A string with the annotation and the E - value in parentheses or "NA" if an annotation is not available.
    ********************************************************************************************************************
    """
    try:
        assert db in annot_dict
        annot_set_filt = set(
            [x for x in annot_dict[db][hg][0] if x.strip() != ""]
        )
        assert len(annot_set_filt) > 0
        return (
            f"{'; '.join(sorted(annot_set_filt))} ({max(annot_dict[db][hg][1])})"
        )
    except Exception as e:
        return "NA"


def gather_value_from_dict_for_homolog_group(hg, info_dict) -> str:
    """
    Description:
    This function formats / gathers information for an ortholog group from a dictionary.
    ********************************************************************************************************************
    Parameters:
    - hg: The ortholog group identifier.
    - info_dict: The dictionary with information to be extracted for the ortholog group.
    ********************************************************************************************************************
    Returns:
    - A string with the information for the ortholog group or "NA" if information is not available.
    ********************************************************************************************************************
    """
    try:
        return info_dict[hg]
    except Exception as e:
        return "NA"


def load_table_in_panda_data_frame_from_string(input_string) -> pd.DataFrame:
    """
    Description:
    This function processes a string and stores it as a pandas dataframe.
    ********************************************************************************************************************
    Parameters:
    - input_file: The input TSV file, with first row corresponding to the header.
    - numeric_columns: Set of column names which should have numeric data.
    ********************************************************************************************************************
    Returns:
    - panda_df: A pandas DataFrame object reprsentation of the input TSV file.
    ********************************************************************************************************************
    """
    panda_df = None
    try:
        data = []
        for line in input_string.split("\n"):
            if line.strip() == "":
                continue
            line = line.strip("\n")
            ls = line.split("\t")
            data.append(ls)

        panda_dict: Dict[str, Any] = {}
        for ls in zip(*data):
            key = ls[0]
            cast_vals = ls[1:]
            panda_dict[key] = cast_vals
        panda_df = pd.DataFrame(panda_dict)

    except Exception as e:
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)
    return panda_df


def load_table_in_panda_data_frame(input_file, numeric_columns) -> pd.DataFrame:
    """
    Description:
    This function formats reads a TSV file and stores it as a pandas dataframe.
    ********************************************************************************************************************
    Parameters:
    - input_file: The input TSV file, with first row corresponding to the header.
    - numeric_columns: Set of column names which should have numeric data.
    ********************************************************************************************************************
    Returns:
    - panda_df: A pandas DataFrame object reprsentation of the input TSV file.
    ********************************************************************************************************************
    """
    panda_df = None
    try:
        data = []
        with open(input_file) as oif:
            for line in oif:
                line = line.strip("\n")
                ls = line.split("\t")
                data.append(ls)

        panda_dict: Dict[str, Any] = {}
        for ls in zip(*data):
            key = ls[0]
            cast_vals = ls[1:] # type: ignore
            if key in numeric_columns:
                cast_vals = []
                for val in ls[1:]:
                    cast_vals.append(cast_to_numeric(val))
            panda_dict[key] = cast_vals
        panda_df = pd.DataFrame(panda_dict)

    except Exception as e:
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)
    return panda_df


def chunks(lst, n) -> Generator[List[Any], None, None]:
    """
        Description:
        Function to yield successive n - sized chunks from lst.
        Solution taken from: https://stackoverflow.com/questions/312443/how-do-you-split-a-list-into-evenly-sized-chunks
        ********************************************************************************************************************
        Parameters:
        - lst: A list.
        - n: The chunk size.
        ********************************************************************************************************************
        Yields:
        - chunks of size n.
        ********************************************************************************************************************
    """
    for i in range(0, len(lst), n):
        yield lst[i : i + n]


def parse_genbank_and_find_boundary_genes(inputs) -> None:
    """
    Description:
    Function to parse a full genome GenBank and writes a dictionary of genes per scaffold, gene to scaffold, and a
    set of CDS features which lie on the boundary of scaffolds (within 2000 bp).
    ********************************************************************************************************************
    Parameters:
    - input: A list which can be expanded to the following:
            - sample: The sample / target genome identifier.
            - sample_genbank: The sample / target genome GenBank file.
            - pkl_result_file: The path to the pickle file to write with the sample / target genome information mentioned
                               in the description.
    ********************************************************************************************************************
    """

    distance_to_scaffold_boundary = 2000
    gene_location: Dict[str, Any] = {}
    scaffold_genes = defaultdict(set)
    boundary_genes = set([])
    gene_id_to_order = defaultdict(dict)
    gene_order_to_id = defaultdict(dict)

    sample, sample_genbank, pkl_result_file = inputs
    osg = None
    if sample_genbank.endswith(".gz"):
        osg = gzip.open(sample_genbank, "rt")
    else:
        osg = open(sample_genbank)

    for rec in SeqIO.parse(osg, "genbank"):
        scaffold = rec.id
        scaffold_length = len(str(rec.seq))
        boundary_ranges = set(
            range(1, distance_to_scaffold_boundary + 1)
        ).union(
            set(
                range(
                    scaffold_length - distance_to_scaffold_boundary,
                    scaffold_length + 1,
                )
            )
        )
        gene_starts = []
        for feature in rec.features:
            if not feature.type == "CDS":
                continue
            locus_tag = feature.qualifiers.get("locus_tag")[0]

            start, end, direction, all_coords = process_location_string(
                str(feature.location) # type: ignore
            )

            gene_location[locus_tag] = {
                "scaffold": scaffold,
                "start": start,
                "end": end,
                "direction": direction,
            }
            scaffold_genes[scaffold].add(locus_tag)

            gene_range = set(range(start, end + 1))
            if len(gene_range.intersection(boundary_ranges)) > 0:
                boundary_genes.add(locus_tag)

            gene_starts.append([locus_tag, start])

        for i, g in enumerate(sorted(gene_starts, key=itemgetter(1))):
            gene_id_to_order[scaffold][g[0]] = i
            gene_order_to_id[scaffold][i] = g[0]
    

    sample_data = [
        gene_location,
        dict(scaffold_genes),
        boundary_genes,
        dict(gene_id_to_order),
        dict(gene_order_to_id),
    ]

    with open(pkl_result_file, "wb") as pickle_file:
        pickle.dump(sample_data, pickle_file, protocol=pickle.HIGHEST_PROTOCOL)


def convert_genome_genbank_to_fasta(inputs) -> None:
    input_gbk, output_fasta = inputs
    if input_gbk.endswith(".gz"):
        with open(output_fasta, "w") as output_fasta_handle:
            with gzip.open(input_gbk, "rt") as oig:
                for rec in SeqIO.parse(oig, "genbank"):
                    output_fasta_handle.write(f">{rec.id}\n{str(rec.seq)}\n")
    else:
        with open(output_fasta, "w") as output_fasta_handle:
            with open(input_gbk) as oig:
                for rec in SeqIO.parse(oig, "genbank"):
                    output_fasta_handle.write(f">{rec.id}\n{str(rec.seq)}\n")
                
def create_nj_tree(
    additional_genbanks_directory,
    species_tree,
    workspace_dir,
    log_object,
    threads=1,
) -> None:
    """
    Description:
    Function to create a species tree using ANI estimates + a neighbor joining approach.
    ********************************************************************************************************************
    Parameters:
    - additional_genbanks_directory: Directory with full genomes in GenBank format.
    - workspace_dir: Workspace directory.
    - log_object: A logging object.
    - threads: The number of threads to use.
    ********************************************************************************************************************
    """

    try:
        workspace_dir = os.path.abspath(workspace_dir) + "/"
        tmp_genome_fasta_dir = workspace_dir + "Genome_FASTAs/"
        setup_ready_directory([tmp_genome_fasta_dir])

        all_genomes_listing_file = workspace_dir + "Genome_FASTAs_Listing.txt"
        with open(all_genomes_listing_file, "w") as all_genomes_listing_handle:
            conversion_inputs = []
            all_samples = set([])
            for f in os.listdir(additional_genbanks_directory):
                if not f.endswith(".gbk.gz"):
                    continue
                input_gbk = additional_genbanks_directory + f
                sample = ".gbk.gz".join(f.split(".gbk.gz")[:-1])
                output_fasta = tmp_genome_fasta_dir + (sample + ".fasta") 
                all_samples.add(sample)
                conversion_inputs.append([input_gbk, output_fasta])
                all_genomes_listing_handle.write(f"{output_fasta}\n")        

        # parallelize conversion
        msg = f"Converting genome GenBank to FASTA for {len(conversion_inputs)} genomes"
        log_object.info(msg)
        sys.stdout.write(msg + "\n")

        p = multiprocessing.Pool(threads)
        for _ in tqdm.tqdm(
            p.imap_unordered(convert_genome_genbank_to_fasta, conversion_inputs),
            total=len(conversion_inputs),
        ):
            pass
        p.close()

        # run skani triangle
        skani_result_file = workspace_dir + "Skani_Triangle_Edge_Output.txt"
        skani_triangle_cmd = [
            "skani",
            "triangle",
            "-E",
            "-l",
            all_genomes_listing_file,
            "-t",
            str(threads),
            "-o",
            skani_result_file,
        ]
        try:
            subprocess.call(
                " ".join(skani_triangle_cmd),
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                executable="/bin/bash",
            )
            assert os.path.isfile(skani_result_file)
            log_object.info(
                f"Successfully ran: {' '.join(skani_triangle_cmd)}"
            )
        except Exception as e:
            log_object.error(
                f"Had an issue with running skani: {' '.join(skani_triangle_cmd)}"
            )
            sys.stderr.write(
                f"Had an issue with running skani: {' '.join(skani_triangle_cmd)}\n"
            )
            log_object.error(e)
            sys.stderr.write(traceback.format_exc())
            sys.exit(1)

        shutil.rmtree(tmp_genome_fasta_dir)

        dist_matrix_file = workspace_dir + "Skani_Based_Distance_Matrix.txt"
        with open(dist_matrix_file, "w") as dist_matrix_handle:
            stos_dists = defaultdict(lambda: defaultdict(lambda: 1.0))
            with open(skani_result_file) as osrf:
                for i, line in enumerate(osrf):
                    if i == 0:
                        continue
                    line = line.strip()
                    ls = line.split("\t")
                    s1 = ".fasta".join(ls[0].split("/")[-1].split(".fasta")[:-1])
                    s2 = ".fasta".join(ls[1].split("/")[-1].split(".fasta")[:-1])
                    if s1 != s2:
                        ani = float(ls[2]) / 100.0
                        dist_ani = 1.0 - ani
                        stos_dists[s1][s2] = dist_ani
                        stos_dists[s2][s1] = dist_ani
                    else:
                        stos_dists[s1][s2] = 0.0

            dist_matrix_handle.write(
                "sample\t" + "\t".join(sorted(all_samples)) + "\n"
            )
            for s1 in sorted(all_samples):
                printlist = [s1]
                for s2 in sorted(all_samples):
                    printlist.append(str(stos_dists[s1][s2]))
                dist_matrix_handle.write("\t".join(printlist) + "\n")
            

        rscript_path = workspace_dir + "generateNjTree.R"
        unrooted_tree_file = workspace_dir + "Unrooted_Species_Tree.nwk"
        generate_nj_tree(
            rscript_path, dist_matrix_file, unrooted_tree_file, log_object
        )
        nj_tree_cmd = ["Rscript", rscript_path]
        try:
            subprocess.call(
                " ".join(nj_tree_cmd),
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                executable="/bin/bash",
            )
            assert os.path.isfile(unrooted_tree_file)
            log_object.info(f"Successfully ran: {' '.join(nj_tree_cmd)}")
        except Exception as e:
            log_object.error(
                f"Had an issue with generating neigbhor - joining tree: {' '.join(nj_tree_cmd)}" \
            )
            sys.stderr.write(
                f"Had an issue with generating neighbor - joining tree: {' '.join(nj_tree_cmd)}\n"
            )
            log_object.error(e)
            sys.stderr.write(traceback.format_exc())
            sys.exit(1)

        # Midpoint the tree using ete3
        t = Tree(unrooted_tree_file)
        R = t.get_midpoint_outgroup()
        t.set_outgroup(R)
        t.write(outfile=species_tree, format=1)

    except Exception as e:
        sys.stderr.write("Issues with creating species tree.\n")
        sys.stderr.write(traceback.format_exc())
        log_object.error("Issues with creating species tree.")
        log_object.error(traceback.format_exc())
        sys.exit(1)


def determine_column_name_based_on_index(index) -> str:
    """
    Function to determine spreadsheet column name for a given index
    """
    # offset at 0
    num_to_char: Dict[str, Any] = {}
    alphabet = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    alphabet_first_spot = [""] + alphabet
    for i, c in enumerate(alphabet):
        num_to_char[str(i)] = c
    level = math.floor(index / 26)
    remainder = index % 26
    columname = alphabet_first_spot[level]
    columname += alphabet[remainder]
    return columname


def clean_up(clean_up_dirs_and_files, log_object) -> None:
    """
    Basic function to clean - up disk heavy files / directories that are intermediate.
    """
    try:
        for df in clean_up_dirs_and_files:
            if os.path.isfile(df):
                log_object.warning(f"Deleting the file {df}")
                os.system(f"rm -fi {df}")
            elif os.path.isdir(df):
                log_object.warning(f"Deleting the file {df}")
                shutil.rmtree(df)
            else:
                log_object.error(f"Couldn't find {df} to delete!")
    except Exception as e:
        sys.stderr.write("Issues with cleaning up files / directories.\n")
        sys.stderr.write(traceback.format_exc())
        log_object.error("Issues with cleaning up files / directories.\n")
        log_object.error(traceback.format_exc())
        sys.exit(1)

def diamond_blast_and_get_best_hits(
    cluster_id,
    query_protein_fasta,
    key_protein_fasta,
    target_genomes_db,
    workspace_dir,
    log_object,
    identity_cutoff=40.0,
    coverage_cutoff=70.0,
    evalue_cutoff=1e-3,
    blastp_mode="very - sensitive",
    prop_key_prots_needed=0.0,
    threads=1,
) -> None:
    """
    Description:

    ********************************************************************************************************************
    Parameters:
    - cluster_id: Identifier of the query gene cluster (e.g. single BGC or phage) in consideration. \
    - query_protein_fasta: FASTA file with all protein sequences for the query gene cluster.
    - key_protein_fasta: FASTA file with all the key (subset) protein sequences for the query gene cluster. \
    - target_genomes_db: prepTG results directory of target genomes.
    - workspace_dir: Workspace directory.
    - log_object: A logging object.
    - identity_cutoff: Minimum percent identity for homologs to query gene cluster proteins.
    - coverage_cutoff: Minimum query coverage for homologs to query gene cluster proteins.
    - evalue_cutoff: Maximum E - value for homologs to query gene cluster proteins.
    - blastp_mode: Sensitivity mode for DIAMOND blastp.
    - prop_key_prots_needed: The proportion of key proteins needed for the query gene cluster to be deemed present in a
                             target genome.
    - threads: The number of threads to use.
    ********************************************************************************************************************
    """
    try:
        genome_wide_tsv_result_file = workspace_dir + "total_gcs.tsv"
        target_genome_dmnd_db = target_genomes_db + "Target_Genomes_DB.dmnd"

        # run diamond (use's fai function: runDiamondBlastp)
        diamond_results_file = workspace_dir + "DIAMOND_Results.txt"
        fai.run_diamond_blastp(
            target_genome_dmnd_db,
            query_protein_fasta,
            workspace_dir,
            log_object,
            diamond_sensitivity=blastp_mode,
            evalue_cutoff=evalue_cutoff,
            compute_query_coverage=True,
            threads=threads,
        )

        all_hgs = set([])
        key_hgs = set([])

        with open(query_protein_fasta) as oqpf:
            for rec in SeqIO.parse(oqpf, "fasta"):
                all_hgs.add(rec.id)

        with open(key_protein_fasta) as okpf:
            for rec in SeqIO.parse(okpf, "fasta"):
                key_hgs.add(rec.id)

        # parse results (based closely on fai function: processDiamondBlastp)
        # Use TypedDict for better type safety
        best_hit_per_lt: Dict[str, Dict[str, BestHitInfo]] = defaultdict(
            lambda: defaultdict(
                lambda: {
                    "hg_list": [],
                    "best_bitscore": 0.0,
                    "identity_list": [],
                    "sql_ratio_list": []
                }
            )
        )
        with open(diamond_results_file) as orf:
            for line in orf:
                line = line.strip()
                ls = line.split()
                hg = ls[0]
                sample = ls[1].split("|")[0]
                lt = ls[1].split("|")[1]
                identity = float(ls[2])
                qcovhsp = float(ls[7])
                if qcovhsp < coverage_cutoff or identity < identity_cutoff:
                    continue
                bitscore = float(ls[4])
                qlen = float(ls[5])
                slen = float(ls[6])
                sql_ratio = float(slen) / float(qlen)
                
                # Type assertion to ensure we're working with the correct types
                current_best = best_hit_per_lt[sample][lt]
                if bitscore > current_best["best_bitscore"]:
                    current_best["best_bitscore"] = bitscore
                    current_best["hg_list"] = [hg]
                    current_best["identity_list"] = [identity]
                    current_best["sql_ratio_list"] = [sql_ratio]
                elif bitscore == current_best["best_bitscore"]:
                    current_best["hg_list"].append(hg)
                    current_best["identity_list"].append(identity)
                    current_best["sql_ratio_list"].append(sql_ratio)

        sample_best_hit_per_hg: Dict[str, Dict[str, List[SampleBestHitInfo]]] = defaultdict(lambda: defaultdict(list))
        for sample in best_hit_per_lt:
            for lt in best_hit_per_lt[sample]:
                current_best = best_hit_per_lt[sample][lt]
                for hgi, hg in enumerate(current_best["hg_list"]):
                    sample_best_hit_per_hg[sample][hg].append({
                        "bitscore": current_best["best_bitscore"],
                        "identity": current_best["identity_list"][hgi],
                        "sql_ratio": current_best["sql_ratio_list"][hgi],
                        "locus_tag": lt
                    })

        with open(genome_wide_tsv_result_file, "w") as outf_handle:
            outf_handle.write("genome\tblank1\taai\tblank2\tshared_gene_prop\n")
            for sample in sample_best_hit_per_hg:
                top_identities: List[float] = []
                hgs_found = set([])
                key_hgs_found = set([])
                for hg in sample_best_hit_per_hg[sample]:
                    for i, hginfo in enumerate(
                        sorted(
                            sample_best_hit_per_hg[sample][hg],
                            key=lambda x: (x["bitscore"], x["identity"], x["sql_ratio"]),
                            reverse=False,
                        )
                    ):
                        if i == 0:
                            top_identities.append(hginfo["identity"])
                            hgs_found.add(hginfo["locus_tag"])
                            if hg in key_hgs:
                                key_hgs_found.add(hg)

                key_hgs_found_prop = float(len(key_hgs_found)) / len(key_hgs)
                if key_hgs_found_prop < prop_key_prots_needed:
                    continue
                aai = statistics.mean(top_identities)
                sgp = float(len(hgs_found)) / len(all_hgs)
                printrow = [sample, "NA", str(aai), "NA", str(sgp)]
                outf_handle.write("\t".join(printrow) + "\n")
        

    except Exception as e:
        sys.stderr.write(
            "Issues with running simple BLASTp for abon, atpoc, or apos.\n"
        )
        sys.stderr.write(traceback.format_exc())
        log_object.error(
            "Issues with running simple BLASTp for abon, atpoc, or apos.\n"
        )
        log_object.error(traceback.format_exc())
        sys.exit(1)


def diamond_blast(inputs) -> None:
    og_prot_file, og_prot_dmnd_db, og_blast_file, log_object = inputs

    makedb_cmd = [
        "diamond",
        "makedb",
        "--ignore - warnings",
        "--in",
        og_prot_file,
        "-d",
        og_prot_dmnd_db,
        "--threads",
        "1"
    ]
    search_cmd = [
        "diamond",
        "blastp",
        "--ignore - warnings",
        "-p",
        "1",
        "-d",
        og_prot_dmnd_db,
        "-q",
        og_prot_file,
        "-o",
        og_blast_file,
    ]
    try:
        subprocess.call(
            " ".join(makedb_cmd),
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            executable="/bin/bash",
        )
        assert os.path.exists(og_prot_dmnd_db)
        log_object.info(f"Successfully ran: {' '.join(makedb_cmd)}")
    except Exception as e:
        log_object.error(
            f"Had an issue running DIAMOND makedb: {' '.join(makedb_cmd)}"
        )
        sys.stderr.write(
            f"Had an issue running DIAMOND makedb: {' '.join(makedb_cmd)}\n"
        )
        log_object.error(e)
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)

    try:
        subprocess.call(
            " ".join(search_cmd),
            shell=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            executable="/bin/bash",
        )
        assert os.path.exists(og_blast_file)
        log_object.info(f"Successfully ran: {' '.join(search_cmd)}")
    except Exception as e:
        log_object.error(
            f"Had an issue running DIAMOND blastp: {' '.join(search_cmd)}"
        )
        sys.stderr.write(
            f"Had an issue running DIAMOND blastp: {' '.join(search_cmd)}\n"
        )
        log_object.error(e)
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def determine_fai_param_recommendations(
    genbanks, ortho_matrix_file, hg_prot_dir, outdir, log_object, threads=1
) -> None:
    """
    Description:
    Function to determine parameter recommendations for running fai from known instances based on zol orthology
    inference.
    ********************************************************************************************************************
    Parameters:
    - genbanks: list of gene cluster GenBank files.
    - ortho_matrix_file: zol orthology inference.
    - hg_prot_dir: directory of ortholog group protein sequences, each ortholog group is its own FASTA file.
    - outdir: workspace directory
    - log_object: A logging object.
    ********************************************************************************************************************
    """
    try:
        outdir = os.path.abspath(outdir) + "/"
        lt_to_og: Dict[str, Any] = {}
        og_gbks = defaultdict(set)
        og_gbk_lts = defaultdict(lambda: defaultdict(set))
        og_lts = defaultdict(set)
        gbks = []
        with open(ortho_matrix_file) as oomf:
            for i, line in enumerate(oomf):
                line = line.strip("\n")
                ls = line.split("\t")
                if i == 0:
                    gbks = ls[1:]
                else:
                    og = ls[0]
                    for j, lts in enumerate(ls[1:]):
                        gbk = gbks[j]
                        for lt in lts.split(", "):
                            if lt.strip() != "":
                                lt = "|".join(lt.split("|")[1:])
                                lt_to_og[lt] = og
                                og_lts[og].add(lt)
                                og_gbk_lts[og][gbk].add(lt)
                                og_gbks[og].add(gbk)

        self_blast_dir = outdir + "Self_Blast_OG_Proteins/"
        setup_ready_directory([self_blast_dir])

        diamond_self_blasting_inputs = []
        for f in os.listdir(hg_prot_dir):
            og = f.split(".faa")[0]
            og_prot_file = hg_prot_dir + f
            og_blast_file = self_blast_dir + og + ".txt"
            og_prot_dmnd_db = hg_prot_dir + f
            diamond_self_blasting_inputs.append(
                [og_prot_file, og_prot_dmnd_db, og_blast_file, log_object]
            )

        msg = f"Running reflexive DIAMOND BLASTp for {len(diamond_self_blasting_inputs)} ortholog groups"
        log_object.info(msg)
        sys.stdout.write(msg + "\n")

        p = multiprocessing.Pool(threads)
        for _ in tqdm.tqdm(
            p.imap_unordered(diamond_blast, diamond_self_blasting_inputs),
            total=len(diamond_self_blasting_inputs),
        ):
            pass
        p.close()

        og_min_eval_params = outdir + "OG_Information.txt"
        with open(og_min_eval_params, "w") as omep_handle:
            maximum_evalues = []
            near_core_ogs_maximum_evalues = []
            omep_handle.write(
                "\t".join(
                    [
                        "og",
                        "maximum_evalue",
                        "conservation",
                        "is_near_core",
                        "is_single_copy",
                    ]
                )
                + "\n"
            )
            near_core_ogs = set([])
            for f in os.listdir(self_blast_dir):
                self_blast_file = self_blast_dir + f    
                og = f.split(".txt")[0]
                maximum_evalue = -1.0
                with open(self_blast_file) as osbf:
                    for line in osbf:
                        line = line.strip()
                        ls = line.split("\t")
                        evalue = float(ls[-2])
                        if evalue > maximum_evalue:
                            maximum_evalue = evalue
                maximum_evalues.append(maximum_evalue)

                og_gbk_count = 0
                singlecopy = True
                for gbk in og_gbk_lts[og]:
                    og_gbk_count += 1
                    lt_count = 0
                    for lt in og_gbk_lts[og][gbk]:
                        lt_count += 1
                    if lt_count > 1:
                        singlecopy = False

                conservation = og_gbk_count / float(len(genbanks))
                nearcore = False
                if conservation >= 0.8:
                    nearcore = True
                    near_core_ogs.add(og)
                    near_core_ogs_maximum_evalues.append(maximum_evalue)

                omep_handle.write(
                    f"{og}\t{maximum_evalue}\t{conservation}\t{nearcore}\t{singlecopy}\n"
                )
        

        og_list = defaultdict(list)
        cds_counts = []
        prop_genes_nc = []
        gbk_og_counts = defaultdict(lambda: defaultdict(int))
        for gbk in genbanks:
            cds_count = 0
            nc_cds_count = 0
            with open(gbk) as ogbk:
                for i, rec in enumerate(SeqIO.parse(ogbk, "genbank")):
                    for feature in rec.features:
                        if feature.type != "CDS":
                            continue
                        lt = feature.qualifiers.get("locus_tag")[0]
                        loc_str = str(feature.location)
                        start, end, direction, all_coords = process_location_string(loc_str) # type: ignore
                        og = lt_to_og[lt]
                        gbk_og_counts[gbk][og] += 1
                        nc = False
                        if og in near_core_ogs:
                            nc = True
                            nc_cds_count += 1
                        og_list[f"{gbk}|{str(i)}"].append([lt, og, start, nc])
                        cds_count += 1
            prop_genes_nc.append(nc_cds_count / float(cds_count))
            cds_counts.append(cds_count)

        cds_between_ncs = []
        for rec in og_list:
            last_nc = None
            for i, lt_info in enumerate(
                sorted(og_list[rec], key=itemgetter(2))
            ):
                is_nc = lt_info[3]
                if is_nc:
                    if last_nc != None:
                        cbn = i - last_nc
                        cds_between_ncs.append(cbn)
                    last_nc = i

        gbk_scores: Dict[str, Any] = {}
        for i, gbk1 in enumerate(genbanks):
            go1cs = gbk_og_counts[gbk1]
            g1ogs = set(go1cs.keys())
            sum_jaccard = 0.0
            for j, gbk2 in enumerate(genbanks):
                go2cs = gbk_og_counts[gbk2]
                g2ogs = set(go2cs.keys())
                ogs_intersect = g1ogs.intersection(g2ogs)
                ogs_union = g1ogs.union(g2ogs)
                ogs_jaccard = len(ogs_intersect) / float(len(ogs_union))
                sum_jaccard += ogs_jaccard
            gbk_scores[gbk1] = sum_jaccard

        ref_gbk = None
        for i, gbk in enumerate(
            sorted(gbk_scores.items(), key=itemgetter(1), reverse=True)
        ):
            if i == 0:
                ref_gbk = gbk[0]

        # extract near - core OG sequences from ref GBK
        near_core_prots_faa_file = (
            outdir + "NearCore_Proteins_from_Representative.faa"
        )
        with open(near_core_prots_faa_file, "w") as ncpff_handle:
            ref_cds_nc = 0
            ref_cds = 0
            if ref_gbk is not None:
                with open(ref_gbk) as orgf:
                    for rec in SeqIO.parse(orgf, "genbank"):
                        for feature in rec.features:
                            if not feature.type == "CDS":
                                continue
                            lt = feature.qualifiers.get("locus_tag")[0]
                            seq = feature.qualifiers.get("translation")[0]
                            ref_cds += 1
                            og = lt_to_og[lt]
                            if og in near_core_ogs:
                                ref_cds_nc += 1
                                ncpff_handle.write(f">{lt}\n{seq}\n")
        

        prop_ref_cds_nc = "NA"
        if float(ref_cds) > 0:
            prop_ref_cds_nc = ref_cds_nc / float(ref_cds)

        max_distance_between_ncs = "NA"
        if len(cds_between_ncs) > 0:
            max_distance_between_ncs = max(cds_between_ncs)
        median_cds_count = "NA"
        if len(cds_counts) > 0:
            median_cds_count = statistics.median(cds_counts)
        maximum_of_maximum_evalues = "NA"
        if len(maximum_evalues) > 0:
            maximum_of_maximum_evalues = max(maximum_evalues)
        maximum_of_near_core_ogs_maximum_evalues = "NA"
        if len(near_core_ogs_maximum_evalues) > 0:
            maximum_of_near_core_ogs_maximum_evalues = max(
                near_core_ogs_maximum_evalues
            )
        median_prop_cds_nc = "NA"
        if len(prop_genes_nc) > 0:
            median_prop_cds_nc = statistics.median(prop_genes_nc)

        parameter_recommendations_file = (
            outdir + "Parameter_Recommendations_for_fai.txt"
        )
        with open(parameter_recommendations_file, "w") as prf_handle:
            prf_handle.write(
                "=============================================================\n"
            )
            prf_handle.write(
                "Recommendations for running fai to find additional instances of gene cluster:\n"
            )
            prf_handle.write(
                "-------------------------------------------------------------\n"
            )
            prf_handle.write(
                "Note, this functionality assumes that the known instances of the gene cluster\nare representative of the gene cluster / taxonomic diversity you will be searching.\n"
            )
            prf_handle.write(
                "=============================================================\n"
            )
            prf_handle.write("General statistics:\n")
            prf_handle.write(
                "=============================================================\n"
            )
            prf_handle.write(
                f"Maximum of maximum E-values observed for any OG\t{maximum_of_maximum_evalues}\n"
            )
            prf_handle.write(
                f"Maximum of near-core OG  E-values observed:\t{maximum_of_near_core_ogs_maximum_evalues}\n"
            )
            prf_handle.write(
                f"Maximum distance between near-core OGs:\t{max_distance_between_ncs}\n"
            )
            prf_handle.write(f"Median CDS count:\t{median_cds_count}\n")
            prf_handle.write(
                f"Median proportion of CDS which are near-core (conserved in 80 percent of gene clusters):\t{median_prop_cds_nc}\n"
            )
            prf_handle.write(
                f"Best representative query gene cluster instance to use:\t{ref_gbk}\n"
            )
            prf_handle.write(
                "=============================================================\n"
            )
            prf_handle.write(
                "Parameter recommendations - threads set to 4 by default\n"
            )
            prf_handle.write(
                "please provide the path to the prepTG database yourself!\n"
            )
            prf_handle.write(
                "=============================================================\n"
            )
            prf_handle.write(
                "Lenient / Sensitive Recommendations for Exploratory Analysis:\n"
            )
            fai_cmd = [
                "fai",
                "--threads",
                "4",
                "--output-dir",
                "fai_Search_Results/",
                "--draft-mode",
            ]
            if maximum_of_maximum_evalues != "NA":
                fai_cmd += [
                    "--evalue-cutoff",
                    str(max(maximum_of_maximum_evalues, 1e-10)),
                ]
            else:
                fai_cmd += ["--evalue-cutoff", "1e-10"]

            if prop_ref_cds_nc != "NA":
                fai_cmd += ["--min-prop", str(max(prop_ref_cds_nc - 0.25, 0.1))]
            else:
                fai_cmd += ["--min-prop", "0.1"]

            fai_cmd += ["--syntenic-correlation-threshold", "0.0"]
            if max_distance_between_ncs != "NA":
                fai_cmd += [
                    "--max-genes-disconnect",
                    str(max_distance_between_ncs + 3),
                ]
            else:
                fai_cmd += ["--max-genes-disconnect", "10"]
            prf_handle.write(" ".join(fai_cmd) + "\n")
            prf_handle.write(
                "-------------------------------------------------------------\n"
            )

            prf_handle.write("Strict / Specific Recommendations:\n")
            fai_cmd = [
                "fai",
                "--threads",
                "4",
                "--output - dir",
                "fai_Search_Results/",
                "--draft - mode",
                "--filter - paralogs",
            ]

            if maximum_of_maximum_evalues != "NA":
                fai_cmd += ["--evalue - cutoff", str(maximum_of_maximum_evalues)]

            if prop_ref_cds_nc != "NA":
                fai_cmd += ["--min - prop", str(max(prop_ref_cds_nc, 0.25))]

            if max_distance_between_ncs != "NA":
                fai_cmd += [
                    "--max-genes-disconnect",
                    str(max_distance_between_ncs),
                ]

            fai_cmd += [
                "key-protein-queries",
                near_core_prots_faa_file,
                "--key-protein-min-prop",
                "0.5",
            ]
            if maximum_of_near_core_ogs_maximum_evalues != "NA":
                fai_cmd += [
                    "--key-protein-evalue-cutoff",
                    str(maximum_of_near_core_ogs_maximum_evalues),
                ]

            prf_handle.write(" ".join(fai_cmd) + "\n")
            prf_handle.write(
                "-------------------------------------------------------------\n"
            )
        

        os.system(f"cat {parameter_recommendations_file}")

    except Exception as e:
        sys.stderr.write(
            "Issue with determining parameter recommendations for running fai based on quick zol analysis of known gene cluster instances.\n"
        )
        sys.stderr.write(traceback.format_exc())
        log_object.error(
            "Issue with determining parameter recommendations for running fai based on quick zol analysis of known gene cluster instances.\n"
        )
        log_object.error(traceback.format_exc())
        sys.exit(1)


def run_pyhmmer_for_ribo_prots(
    best_tg_gbk_file, tg_query_prots_file, ribo_norm_dir, log_object, threads=1
) -> None:
    """
    Description:
    Annotate ribosomal proteins from Hug et al. 2016 (using HMMs as provided in GToTree by Lee 2019) in a reference genome. \
    ********************************************************************************************************************
    Parameters:
    - best_tg_gbk_file: The full genome GenBank file for the reference / select genome.
    - tg_query_prots_file: The FASTA file to write ribosomal proteins identified to (note will append to file). \
    - ribo_norm_dir: Output workspace to write intermediate files to.
    - log_object: A logging object.
    - threads: The number of threads to use [Default is 1].
    ********************************************************************************************************************
    """
    try:
        rp_db_file = None
        try:
            zol_data_directory = str(os.getenv("ZOL_DATA_PATH")).strip()
            db_locations = None
            if zol_data_directory != "None":
                try:
                    zol_data_directory = (
                        os.path.abspath(zol_data_directory) + "/"
                    )
                    db_locations = (
                        zol_data_directory + "database_location_paths.txt"
                    )
                except Exception as e:
                    pass

            if db_locations == None or not os.path.exists(db_locations):
                sys.stderr.write(
                    "Warning: databases do not appear to be setup or setup properly - so unable to annotate!\n"
                )

            if db_locations is not None:
                with open(db_locations) as odl:
                    for line in odl:
                        line = line.strip()
                        if len(line.split("\t")) != 4:
                            continue
                        name, _, db_file, _ = line.split("\t")
                        if name == "riboprots":
                            rp_db_file = db_file

            assert rp_db_file != None and os.path.exists(rp_db_file)
        except Exception as e:
            msg = \
     "Issues validating that the ribosomal proteins file is available. Downloading from Zenodo to output directory."
            log_object.warning(msg)
            sys.stderr.write(msg + "\n")

            curr_path = os.path.abspath(os.getcwd()) + "/"
            try:
                download_path = "/".join((db_locations.split("/")[:-1])) + "/" # type: ignore
                download_links = [
                    "https://zenodo.org/records/13858489/files/Universal-Hug-et-al.hmm?download=1"
                ]
            except Exception as e:
                msg = "Unable to download ribosomal proteins file from Zenodo."
                log_object.error(msg)
                sys.stderr.write(msg + "\n")
                sys.exit(1)

            # Download
            os.chdir(download_path)
            rp_db_file = download_path + "Universal-Hug-et-al.hmm"
            try:
                for dl in download_links:
                    axel_download_dbs_cmd = [
                        "axel",
                        "-a",
                        "-n",
                        str(threads),
                        dl,
                    ]
                    os.system(" ".join(axel_download_dbs_cmd))
                    assert os.path.exists(rp_db_file)
                    os.system(" ".join(["hmmpress", rp_db_file]))
            except Exception as e:
                sys.stderr.write("Error occurred during downloading!\n")
                log_object.error(
                    "Error occurred during downloading with axel.\n"
                )
                sys.exit(1)
            os.chdir(curr_path)

        hmm_lengths: Dict[str, Any] = {}
        z = 0
        try:
            with pyhmmer.plan7.HMMFile(rp_db_file) as hmm_file:
                for hmm in hmm_file:
                    hmm_lengths[hmm.name.decode()] = len(hmm.consensus) # type: ignore
                    z += 1
        except Exception as e:
            raise RuntimeError("Problem getting HMM consensus lengths!")

        best_tg_faa_file = ribo_norm_dir + "Reference_Genome_All_Proteins.faa"
        with open(best_tg_faa_file, "w") as best_tg_faa_handle:
            try:
                if best_tg_gbk_file.endswith('.gz'):
                    with gzip.open(best_tg_gbk_file, 'rt') as obtgf:
                        for rec in SeqIO.parse(obtgf, 'genbank'):
                            for feat in rec.features:
                                if feat.type == 'CDS':
                                    lt = feat.qualifiers.get('locus_tag')[0]
                                    prot_seq = feat.qualifiers.get('translation')[0]
                                    best_tg_faa_handle.write(f">{lt}\n{prot_seq}\n")
                else: 
                    with open(best_tg_gbk_file) as obtgf:   
                        for rec in SeqIO.parse(obtgf, "genbank"):
                            for feat in rec.features:
                                if feat.type == "CDS":
                                    lt = feat.qualifiers.get("locus_tag")[0]
                                    prot_seq = feat.qualifiers.get("translation")[0]
                                    best_tg_faa_handle.write(f">{lt}\n{prot_seq}\n")
            except Exception as e:
                raise RuntimeError("Problem processing full genome GenBank file.")
        

        alphabet = pyhmmer.easel.Alphabet.amino()
        sequences = []
        with pyhmmer.easel.SequenceFile(
            best_tg_faa_file, digital=True, alphabet=alphabet
        ) as seq_file:
            sequences = list(seq_file)

        reference_ribo_prots = set([])
        with pyhmmer.plan7.HMMFile(rp_db_file) as hmm_file:
            for hits in pyhmmer.hmmsearch(
                hmm_file,
                sequences,
                bit_cutoffs="trusted",
                Z=int(z),
                cpus=threads,
            ):
                for hit in hits:
                    # solution for calcualting coverage taken from pcamargo's answer in a pyhmmer ticket on Github: https://github.com / althonos / pyhmmer / issues / 27
                    n_aligned_positions = len(
                        hit.best_domain.alignment.hmm_sequence
                    ) - hit.best_domain.alignment.hmm_sequence.count(".")
                    hmm_coverage = (
                        n_aligned_positions
                        / hmm_lengths[hit.best_domain.alignment.hmm_name.decode()]
                    )
                    if hmm_coverage >= 0.25:
                        reference_ribo_prots.add(hit.name.decode())

        with open(tg_query_prots_file, "a+") as tg_query_prots_handle:
            with open(best_tg_faa_file) as obtff:
                for rec in SeqIO.parse(obtff, "fasta"):
                    if rec.id in reference_ribo_prots:
                        tg_query_prots_handle.write(
                            f">Ribosomal_Protein|{rec.id}\n{str(rec.seq)}\n"
                        )
        

    except Exception as e:
        raise RuntimeError(
            "Problem with running pyhmmer for finding ribosomal proteins in reference genome!"
        )


def run_pyhmmer_for_vo_gfor_salt(inputs) -> None:
    """
    Description:
    Annotate a single sample's predicted proteome file for VOG HMMs.
    ********************************************************************************************************************
    Parameters:
    - inputs: A list of length 5:
            - db_file: HMM database.
            - z: Size of database, used for E-value computation.
            - protein_faa: sample's proteome FASTA file.
            - annotation_result_file: Path to output file where to write annotation information.
            - threads: number of threads to use for search.
    ********************************************************************************************************************
    """
    db_file, z, protein_faa, annotation_result_file, threads = inputs
    try:
        hmm_lengths: Dict[str, Any] = {}
        try:
            with pyhmmer.plan7.HMMFile(db_file) as hmm_file:
                for hmm in hmm_file:
                    hmm_lengths[hmm.name.decode()] = len(hmm.consensus) # type: ignore
        except Exception as e:
            raise RuntimeError("Problem getting HMM consensus lengths!")

        alphabet = pyhmmer.easel.Alphabet.amino()
        sequences = []
        with pyhmmer.easel.SequenceFile(
            protein_faa, digital=True, alphabet=alphabet
        ) as seq_file:
            sequences = list(seq_file)

        with open(annotation_result_file, "w") as outf:
            with pyhmmer.plan7.HMMFile(db_file) as hmm_file:
                for hits in pyhmmer.hmmsearch(
                    hmm_file, sequences, Z=int(z), cpus=threads
                ):
                    for hit in hits:
                        # solution for calcualting coverage taken from pcamargo's answer in a pyhmmer ticket on Github: https://github.com/althonos/pyhmmer/issues/27
                        n_aligned_positions = len(
                            hit.best_domain.alignment.hmm_sequence
                        ) - hit.best_domain.alignment.hmm_sequence.count(".")
                        hmm_coverage = (
                            n_aligned_positions
                            / hmm_lengths[hit.best_domain.alignment.hmm_name.decode()] # type: ignore
                        )
                        outf.write(
                            "\t".join(
                                [
                                    hits.query.name.decode(),
                                    "NA",
                                    hit.name.decode(),
                                    "NA",
                                    str(hit.evalue),
                                    str(hit.score),
                                    str(hmm_coverage),
                                ]
                            )
                            + "\n"
                        )
        
    except Exception as e:
        raise RuntimeError(
            "Problem running pyhmmer for annotating MGEs in target genomes!"
        )


def annotate_mges(inputs) -> None:
    """
    Description:
    Annotate MGEs for a single sample's predicted proteome file.
    ********************************************************************************************************************
    Parameters:
    - inputs: A list of length 6:
            - sample: sample / genome name.
            - faa_file: path to proteome for sample / genome.
            - vog_annot_file: path to the pyhmmer results for VOG annotations (will be written to - should not exist). \
            - mobsuite_annot_file: path to the DIAMOND blastp results for MOB-suite annotations (will be written to - should not exist). \
            - is_annot_file: path to the DIAMOND blastp results for ISfinder annotations (will be written to - should not exist). \
            - log_object: a logging object.
    ********************************************************************************************************************
    """
    sample = "NA"
    try:
        (
            sample,
            faa_file,
            vog_annot_file,
            mobsuite_annot_file,
            is_annot_file,
            log_object,
        ) = inputs

        zol_data_directory = str(os.getenv("ZOL_DATA_PATH")).strip()
        db_locations = None
        conda_setup_success = None
        if zol_data_directory != "None":
            try:
                zol_data_directory = os.path.abspath(zol_data_directory) + "/"
                db_locations = (
                    zol_data_directory + "database_location_paths.txt"
                )
            except Exception as e:
                pass

        if db_locations == None or not os.path.exists(db_locations):
            sys.stderr.write(
                "Warning: databases do not appear to be setup or setup properly - so unable to annotate!\n"
            )

        if db_locations is not None:
            with open(db_locations) as odl: # type: ignore
                for line in odl:
                    line = line.strip()
                    if len(line.split("\t")) != 4:
                        continue
                    name, _, db_file, z = line.split("\t")
                    if name == "vog":
                        vog_pyhmmer_input = [
                            db_file,
                            z,
                            faa_file,
                            vog_annot_file,
                            1,
                        ]
                        run_pyhmmer_for_vo_gfor_salt(vog_pyhmmer_input)
                    elif name == "mobsuite" or name == "isfinder":
                        annot_result_file = is_annot_file
                        if name == "mobsuite":
                            annot_result_file = mobsuite_annot_file
                        search_cmd = [
                            "diamond",
                            "blastp",
                            "--ignore-warnings",
                            "-fast",
                            "--outfmt",
                            "6",
                            "qseqid",
                            "sseqid",
                            "pident",
                            "evalue",
                            "qcovhsp",
                            "-p",
                            "1",
                            "-d",
                            db_file,
                            "-q",
                            faa_file,
                            "-o",
                            annot_result_file,
                        ]
                        run_cmd_via_subprocess(
                            search_cmd, log_object, check_files=[annot_result_file]
                        )
    except Exception as e:
        sys.stderr.write(
            f"Issues with MGE annotation commands for sample {sample}.\n"
        )
        log_object.error(
            f"Issues with MGE annotation commands for sample {sample}."
        )
        sys.stderr.write(traceback.format_exc())
        log_object.error(traceback.format_exc())
        sys.exit(1)


def process_mge_annotations(inputs) -> None:
    """
    Description:
    Function to process MGE annotation results from DIAMOND blastp and pyhmmer searching for a single sample as well as
    corresponding genomic GenBank file to determine summarized information for
    ********************************************************************************************************************
    Parameters:
    - inputs: A list of length 6:
            - sample: sample / genome name.
            - gbk_file: path to GenBank for sample / genome.
            - summary_file: The output file for the sample with information on MGE locations.
            - vog_annot_file: path to the pyhmmer results for VOG annotations (should already exist). \
            - mobsuite_annot_file: path to the DIAMOND blastp results for MOB-suite annotations (should already exist). \
            - is_annot_file: path to the DIAMOND blastp results for ISfinder annotations (should already exist). \
            - log_object: a logging object.
    ********************************************************************************************************************
    """
    try:
        # CURRENTLY HARDCODED!
        max_dmnd_annotation_evalue = 1e-5
        max_hmm_annotation_evalue = 1e-5
        min_percent_identity = 40.0  # only used for diamond
        min_coverage = 70.0
        (
            sample,
            gbk_file,
            summary_file,
            vog_annot_file,
            mobsuite_annot_file,
            is_annot_file,
            log_object,
        ) = inputs

        try:
            assert (
                os.path.isfile(vog_annot_file)
                and os.path.isfile(mobsuite_annot_file)
                and os.path.isfile(is_annot_file)
            )
        except Exception as e:
            sys.stderr.write(
                f"Issue validating the existence of one or more of the annotation files for sample: {sample}\n"
            )
            log_object.error(
                f"Issue validating the existence of one or more of the annotation files for sample: {sample}"
            )
            # TODO: Consider making this an end - program error
            return

        vog_hits = set([])
        with open(vog_annot_file) as ovaf:
            for line in ovaf:
                line = line.rstrip("\n")
                if line.startswith("#"):
                    continue
                ls = line.split("\t")
                _, _, query, _, evalue, score, coverage = ls
                evalue = decimal.Decimal(evalue)
                coverage = float(coverage) * 100.0
                if (
                    evalue <= max_hmm_annotation_evalue
                    and coverage >= min_coverage
                ):
                    vog_hits.add(query)

        mobsuite_hits = set([])
        with open(mobsuite_annot_file) as omaf:
            for line in omaf:
                line = line.strip()
                query, _, pident, evalue, qcovhsp = line.split("\t")
                pident = float(pident)
                evalue = decimal.Decimal(evalue)
                qcovhsp = float(pident)
                if (
                    qcovhsp >= min_coverage
                    and evalue <= max_dmnd_annotation_evalue
                    and pident >= min_percent_identity
                ):
                    mobsuite_hits.add(query)

        isfinder_hits = set([])
        with open(is_annot_file) as oiaf:
            for line in oiaf:
                line = line.strip()
                query, _, pident, evalue, qcovhsp = line.split("\t")
                pident = float(pident)
                evalue = decimal.Decimal(evalue)
                qcovhsp = float(pident)
                if (
                    qcovhsp >= min_coverage
                    and evalue <= max_dmnd_annotation_evalue
                    and pident >= min_percent_identity
                ):
                    isfinder_hits.add(query)

        with open(summary_file, "w") as summary_handle:
            summary_handle.write(
                "\t".join(
                    [
                        "sample",
                        "scaffold",
                        "total_cds",
                        "mobsuite_cds_hits",
                        "mobsuite_cds_prop",
                        "vog_cds_hits",
                        "vog_cds_prop",
                        "is_cds_hits",
                        "is_cds_prop",
                        "is_cds_boundaries",
                    ]
                )
                + "\n"
            )
            with gzip.open(gbk_file, "rt") as ogbf:
                for rec in SeqIO.parse(ogbf, "genbank"):
                    scaffold_id = rec.id
                    scaffold_is_boundary_coords = set([])
                    total_cds = 0
                    mhits = 0
                    vhits = 0
                    ihits = 0
                    for feat in rec.features:
                        if feat.type != "CDS":
                            continue
                        lt = feat.qualifiers.get("locus_tag")[0]
                        total_cds += 1
                        if lt in isfinder_hits:
                            ihits += 1
                            start = (
                                min(
                                    [
                                        int(x.strip(">").strip("<"))
                                        for x in str(feat.location)[1:]
                                        .split("]")[0]
                                        .split(":")
                                    ]
                                )
                                + 1
                            )
                            end = max(
                                [
                                    int(x.strip(">").strip("<"))
                                    for x in str(feat.location)[1:]
                                    .split("]")[0]
                                    .split(":")
                                ]
                            )
                            scaffold_is_boundary_coords.add(start)
                            scaffold_is_boundary_coords.add(end)
                        if lt in mobsuite_hits:
                            mhits += 1
                        if lt in vog_hits:
                            vhits += 1
                    if total_cds > 0:
                        summary_handle.write(
                            "\t".join(
                                [
                                    str(x)
                                    for x in [
                                        sample,
                                        scaffold_id,
                                        total_cds,
                                        mhits,
                                        round(mhits / total_cds, 3),
                                        vhits,
                                        round(vhits / total_cds, 3),
                                        ihits,
                                        round(ihits / total_cds, 3),
                                    ]
                                ]
                                + [
                                    ", ".join(
                                        [
                                            str(y)
                                            for y in scaffold_is_boundary_coords
                                        ]
                                    )
                                ]
                            )
                            + "\n"
                        )
        

    except Exception as e:
        sys.stderr.write(
            f"Issues with MGE annotation processing for sample {sample}.\n"
        )
        log_object.error(
            f"Issues with MGE annotation processing for sample {sample}."
        )
        sys.stderr.write(traceback.format_exc())
        log_object.error(traceback.format_exc())
        sys.exit(1)


def process_diamond_for_gc_to_ribo_ratio(
    diamond_results_file,
    rep_prot_to_nonreps,
    query_faa,
    fai_ind_gc_file,
    gc_to_ribo_aai_stat_file,
    log_object,
) -> None:
    """
    Description:
    Function to process results from DIAMOND blastping ribosomal + focal gene cluster proteins from reference genome
    to the remainder of target genomes to determine the Beta-RD statistic.
    ********************************************************************************************************************
    Parameters:
    - diamond_results_file: DIAMOND blastp results searching for query gene cluster and ribosomal proteins from reference
                            genome against all target genomes.
    - rep_prot_to_nonreps: a dictionary mapping representative protein IDs to non-representative protein IDs.   
    - query_faa: the query FASTA file containing all ribosomal + gene cluster protein seqeunces.
    - fai_ind_gc_file: fai individual gene cluster instances tsv.
    - gc_to_ribo_aai_stat_file: Path to output file to write statistics per gene cluster instance.
    - log_object: a logging object.
    ********************************************************************************************************************
    """
    try:
        gc_prots = defaultdict(set)
        genome_gcis = defaultdict(set)
        lt_to_gci: Dict[str, Any] = {}
        with open(fai_ind_gc_file) as ofigf:
            for i, line in enumerate(ofigf):
                if i == 0:
                    continue
                line = line.strip()
                ls = line.split("\t")
                gc_gbk = ls[1]
                gc_gbk_file = gc_gbk.split("/")[-1]
                gci = ".gbk".join(gc_gbk_file.split(".gbk")[:-1])
                genome = ".gbk".join(gc_gbk_file.split(".gbk")[:-1])
                if "_fai-gene-cluster" in genome:
                    genome = genome.split("_fai-gene-cluster")[0]
                genome_gcis[genome].add(gci)
                if gc_gbk.endswith('.gz'):
                    with gzip.open(gc_gbk, 'rt') as ogg:
                        for rec in SeqIO.parse(ogg, "genbank"):
                            for feat in rec.features:
                                if feat.type == "CDS":
                                    lt = feat.qualifiers.get("locus_tag")[0]
                                    gc_prots[genome].add(lt)
                                    lt_to_gci[genome + "|" + lt] = gci
                else:
                    with open(gc_gbk) as ogg:
                        for rec in SeqIO.parse(ogg, "genbank"):
                            for feat in rec.features:
                                if feat.type == "CDS":
                                    lt = feat.qualifiers.get("locus_tag")[0]
                                    gc_prots[genome].add(lt)
                                    lt_to_gci[genome + "|" + lt] = gci
        # Use TypedDict for better type safety
        gci_query_top_hits: Dict[str, Dict[str, GciQueryTopHitInfo]] = defaultdict(
            lambda: defaultdict(
                lambda: {
                    "hits": [],
                    "best_bitscore": 0.0,
                    "identity_list": []
                }
            )
        )
        with open(diamond_results_file) as odrf:
            for line in odrf:
                line = line.strip()
                ls = line.split("\t")
                query = ls[0]
                query_class = query.split("|")[0]
                rep_hit_id = ls[1]
                all_hits = [rep_hit_id]
                if rep_prot_to_nonreps != None:
                    if rep_hit_id in rep_prot_to_nonreps:
                        all_hits = rep_prot_to_nonreps[rep_hit_id]
                    else:
                        all_hits = [rep_hit_id]

                for hit in all_hits:
                    genome, lt = hit.split("|")
                    if (
                        query_class == "Gene_Cluster_Protein"
                        and not lt in gc_prots[genome]
                    ):
                        continue
                    identity = float(ls[2])
                    bitscore = float(ls[4])
                    qcovhsp = float(ls[7])
                    if qcovhsp <= 25.0:
                        continue
                    if query_class == "Ribosomal_Protein":
                        gcis = genome_gcis[genome]
                    else:
                        gcis = set([lt_to_gci[hit]])
                    for gci in gcis:
                        current_hit = gci_query_top_hits[gci][query]
                        if bitscore > current_hit["best_bitscore"]:
                            current_hit["hits"] = [hit]
                            current_hit["best_bitscore"] = bitscore
                            current_hit["identity_list"] = [identity]
                        elif bitscore == current_hit["best_bitscore"]:
                            current_hit["hits"].append(hit)
                            current_hit["identity_list"].append(identity)

        all_ribo_prots = set([])
        all_gc_prots = set([])
        with open(query_faa) as oqf:
            for rec in SeqIO.parse(oqf, "fasta"):
                if rec.id.startswith("Gene_Cluster_Protein|"):
                    all_gc_prots.add(rec.id)
                else:
                    all_ribo_prots.add(rec.id)

        gci_top_hits: Dict[str, List[GciTopHitInfo]] = defaultdict(list)
        for gci in gci_query_top_hits:
            for query in gci_query_top_hits[gci]:
                max_identity = 0.0
                max_hit = None
                current_hit = gci_query_top_hits[gci][query]
                for i, hit in enumerate(current_hit["hits"]):
                    if current_hit["identity_list"][i] >= max_identity:
                        max_identity = current_hit["identity_list"][i]
                        max_hit = hit
                gci_top_hits[gci].append({
                    "query": query,
                    "bitscore": current_hit["best_bitscore"],
                    "max_hit": max_hit, # type: ignore
                    "max_identity": max_identity
                })

        gc_aais = []
        rp_aais = []
        data = []
        for gci in gci_top_hits:
            ribo_identities = []
            gc_identities = []
            accounted_lts = set([])
            hits = gci_top_hits[gci]
            for query_hit_data in sorted(
                hits, key=lambda x: x["bitscore"], reverse=True
            ):
                query = query_hit_data["query"]
                max_ident = query_hit_data["max_identity"]
                max_hit = query_hit_data["max_hit"]
                if max_hit is not None and max_hit not in accounted_lts:
                    if query.split("|")[0] == "Ribosomal_Protein":
                        ribo_identities.append(max_ident)
                    else:
                        gc_identities.append(max_ident)
                if max_hit is not None:
                    accounted_lts.add(max_hit)
            if len(ribo_identities) > 0 and len(gc_identities) > 0:
                ribo_aai = round(statistics.mean(ribo_identities), 3)
                gc_aai = round(statistics.mean(gc_identities), 3)
                ribo_prot_prop = round(
                    len(ribo_identities) / len(all_ribo_prots), 3
                )
                gc_prot_prop = round(len(gc_identities) / len(all_gc_prots), 3)
                gc_aais.append(gc_aai)
                rp_aais.append(ribo_aai)
                data.append(
                    [gci, ribo_aai, gc_aai, ribo_prot_prop, gc_prot_prop]
                )

        slope, intercept, r, p, se = stats.linregress(rp_aais, gc_aais)

        with open(gc_to_ribo_aai_stat_file, "w") as gc_to_ribo_aai_stat_handle:
            gc_to_ribo_aai_stat_handle.write(
                "\t".join(
                    [
                        "gc_instance",
                        "distance_to_linear_regression",
                        "ribo_aai",
                        "gc_aai",
                        "prop_ribo_prots",
                        "prop_gc_prots",
                    ]
                )
                + "\n"
            )
            for rd in data:
                gci, ribo_aai, gc_aai, ribo_prot_prop, gc_prot_prop = rd
                expected_gc_aai = ribo_aai * slope + intercept
                diff_obs_to_exp = gc_aai - expected_gc_aai
                row_data = [
                    gci,
                    diff_obs_to_exp,
                    ribo_aai,
                    gc_aai,
                    ribo_prot_prop,
                    gc_prot_prop,
                ]
                gc_to_ribo_aai_stat_handle.write(
                    "\t".join([str(x) for x in row_data]) + "\n"
                )
        
    except Exception as e:
        sys.stderr.write(
            "Issues with processing DIAMOND results for aligning ribosomal and gene cluster proteins from the reference to the remainder of the target genomes."
        )
        log_object.error(
            "Issues with processing DIAMOND results for aligning ribosomal and gene cluster proteins from the reference to the remainder of the target genomes."
        )
        sys.stderr.write(traceback.format_exc())
        log_object.error(traceback.format_exc())
        sys.exit(1)


def make_gc_vs_ribo_prot_aai_scatterplot(
    rscript_file,
    gc_to_ribo_aai_stat_file,
    gc_ribo_aai_plot_pdf_file,
    log_object,
) -> None:
    """
    Description:
    Function to plot the gene cluster to ribosomal protein AAI relationship as a scatterplot.
    ********************************************************************************************************************
    Parameters:
    - gc_to_ribo_aai_stat_file: The gene cluster and ribosomal proteins AAI relationship file.
    - gc_ribo_aai_plot_pdf_file: The path to the resulting PDF with the scatterplot showing the GC to ribo proteins AAI
                                                             relationship.
    - log_object: a logging object.
    ********************************************************************************************************************
    """
    try:
        generate_salt_gc_vs_ribo_aai_plot( # type: ignore
            rscript_file,
            gc_to_ribo_aai_stat_file,
            gc_ribo_aai_plot_pdf_file,
            log_object,
        )
        plot_cmd = ["Rscript", rscript_file]
        try:
            subprocess.call(
                " ".join(plot_cmd),
                shell=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                executable="/bin/bash",
            )
            assert os.path.isfile(gc_ribo_aai_plot_pdf_file)
            log_object.info(f"Successfully ran: {' '.join(plot_cmd)}")
        except Exception as e:
            log_object.error(
                f"Had an issue running R based plotting - potentially because of R setup issues in conda: {' '.join(plot_cmd)}" \
            )
            sys.stderr.write(
                f"Had an issue running R based plotting - potentially because of R setup issues in conda: {' '.join(plot_cmd)}\n"
            )
            log_object.error(e)
            sys.stderr.write(traceback.format_exc())
            sys.exit(1)
    except Exception as e:
        sys.stderr.write(
            "Issues with creating a scatterplot of gene cluster vs. ribosomal protein AAIs.\n"
        )
        log_object.error(
            "Issues with creating a scatterplot of gene cluster vs. ribosomal protein AAIs."
        )
        sys.stderr.write(traceback.format_exc())
        log_object.error(traceback.format_exc())
        sys.exit(1)


def consolidate_salty_spreadsheet(
    fai_ind_gc_file,
    genome_gbks,
    codoff_result_dir,
    mge_annot_result_file,
    gc_to_ribo_aai_stat_file,
    result_tsv_file,
    result_xlsx_file,
    log_object,
) -> None:
    """
    Description:
    Consolidate results from all three analyses to inform on lateral transfer into a table and create an auto - colored
    spreadsheet XLSX file.
    ********************************************************************************************************************
    Parameters:
    - fai_ind_gc_file: fai individual gene cluster instances tsv.
    - genome_gbks: A dictionary mapping a genome id / name to a full genome GenBank file.
    - codoff_result_dir: results directory from running codoff for each gene cluster against their respective background
                                     genomes.
    - mge_annot_result_file: overview file with MGE annotations for scaffolds from target genomes.
    - gc_to_ribo_aai_stat_file: Path to output file to write statistics per gene cluster instance.
    - result_tsv_file: The path to the output file to write the consolidated table to in TSV format.
    - result_xlsx_file: The path to the output file to write the consolidated table to in XLSX format.
    - log_object: a logging object.
    ********************************************************************************************************************
    """
    try:

        header = [
            "gene cluster (GC) instance",
            "GC gbk path",
            "genome",
            "scaffold",
            "scaffold length (bp)",
            "scaffold CDS count",
            "GC CDS count",
            "codoff empirical P-value",
            "GC AAI observed - expectation",
            "GC AAI between genome and reference genome",
            "ribosomal protein AAI between genome and reference genome",
            "distance to IS-associated element",
            "scaffold CDS proportion IS-associated elements",
            "scaffold CDS proportion VOGs",
            "scaffold CDS proportion plasmid - associated",
        ]

        gc_codoff_pvals = defaultdict(lambda: "NA")
        for f in os.listdir(codoff_result_dir):
            gc = ".txt".join(f.split(".txt")[:-1])
            gc_codoff_file = codoff_result_dir + f
            with open(gc_codoff_file) as ogcf:
                for line in ogcf:
                    line = line.strip()
                    ls = line.split("\t")
                    if ls[0] == "Empirical P-value":
                        gc_codoff_pvals[gc] = ls[1]

        genome_scaffold_annot_info = defaultdict(
            lambda: defaultdict(lambda: ["NA"] * 8)
        )
        if mge_annot_result_file != None and os.path.isfile(mge_annot_result_file):
            with open(mge_annot_result_file) as omarf:
                for i, line in enumerate(omarf):
                    if i == 0:
                        continue
                    line = line.strip("\n")
                    (
                        sample,
                        scaffold,
                        total_cds,
                        mobsuite_cds_hits,
                        mobsuite_cds_prop,
                        vog_cds_hits,
                        vog_cds_prop,
                        is_cds_hits,
                        is_cds_prop,
                        is_cds_boundaries,
                    ) = line.split("\t")
                    genome_scaffold_annot_info[sample][scaffold] = [
                        total_cds,
                        mobsuite_cds_hits,
                        mobsuite_cds_prop,
                        vog_cds_hits,
                        vog_cds_prop,
                        is_cds_hits,
                        is_cds_prop,
                        is_cds_boundaries,
                    ]

        gc_aai_stats = defaultdict(lambda: ["NA"] * 5)
        with open(gc_to_ribo_aai_stat_file) as ogtrasf:
            for i, line in enumerate(ogtrasf):
                if i == 0:
                    continue
                line = line.strip("\n")
                (
                    gc,
                    ratio_statistic,
                    ribo_aai,
                    gc_aai,
                    prop_ribo_prots,
                    prop_gc_prots,
                ) = line.split("\t")
                gc_aai_stats[gc] = [
                    ratio_statistic,
                    ribo_aai,
                    gc_aai,
                    prop_ribo_prots,
                    prop_gc_prots,
                ]

        num_rows = 0
        with open(result_tsv_file, "w") as tsv_outf:
            tsv_outf.write("\t".join(header) + "\n")
            with open(fai_ind_gc_file) as ofigf:
                for i, line in enumerate(ofigf):
                    if i == 0:
                        continue
                    line = line.strip()
                    ls = line.split("\t")
                    gc_gbk = ls[1]
                    gc_gbk_file = gc_gbk.split("/")[-1]
                    gc = ".gbk".join(gc_gbk_file.split(".gbk")[:-1])
                    genome = ".gbk".join(gc_gbk_file.split(".gbk")[:-1])
                    if "_fai-gene-cluster" in genome:
                        genome = genome.split("_fai-gene-cluster")[0]

                    scaffold = "NA"
                    gc_cds_count = 0

                    cds_coords: Dict[str, Any] = {}
                    scaffold_lens: Dict[str, Any] = {}
                    if genome_gbks[genome].endswith('.gz'):
                        with gzip.open(genome_gbks[genome], 'rt') as ogg:
                            for rec in SeqIO.parse(ogg, "genbank"):
                                scaffold_lens[rec.id] = len(str(rec.seq))
                                for feat in rec.features:
                                    if feat.type == "CDS":
                                        lt = feat.qualifiers.get("locus_tag")[0]
                                        loc_str = str(feat.location)
                                        start, end, direction, all_coords = process_location_string(loc_str) # type: ignore
                                        cds_coords[lt] = [start, end]
                    else:
                        with open(genome_gbks[genome]) as ogg:
                            for rec in SeqIO.parse(ogg, "genbank"):
                                scaffold_lens[rec.id] = len(str(rec.seq))
                                for feat in rec.features:
                                    if feat.type == "CDS":
                                        lt = feat.qualifiers.get("locus_tag")[0]
                                        loc_str = str(feat.location)
                                        start, end, direction, all_coords = process_location_string(loc_str) # type: ignore
                                        cds_coords[lt] = [start, end]

                    min_gc_coord = 1e100
                    max_gc_coord = 0
                    with open(gc_gbk) as ogg:
                        for rec in SeqIO.parse(ogg, "genbank"):
                            scaffold = rec.id
                            for feat in rec.features:
                                if feat.type == "CDS":
                                    gc_cds_count += 1
                                    lt = feat.qualifiers.get("locus_tag")[0]
                                    cds_start = cds_coords[lt][0]
                                    cds_end = cds_coords[lt][1]
                                    if cds_start < min_gc_coord:
                                        min_gc_coord = cds_start
                                    if cds_end > max_gc_coord:
                                        max_gc_coord = cds_end

                    (
                        scaffold_cds_count,
                        mobsuite_cds_hits,
                        mobsuite_cds_prop,
                        vog_cds_hits,
                        vog_cds_prop,
                        is_cds_hits,
                        is_cds_prop,
                        is_cds_boundaries,
                    ) = genome_scaffold_annot_info[genome][scaffold]
                    (
                        ratio_statistic,
                        ribo_aai,
                        gc_aai,
                        prop_ribo_prots,
                        prop_gc_prots,
                    ) = gc_aai_stats[gc]
                    gc_codoff_pval = gc_codoff_pvals[gc]

                    min_is_distance = 1e100
                    for is_coord in is_cds_boundaries.split(", "):
                        if is_numeric(is_coord):
                            is_coord = int(is_coord)
                            dist = None
                            if (
                                is_coord >= min_gc_coord
                                and is_coord <= max_gc_coord
                            ):
                                dist = 0
                            else:
                                dist = min(
                                    [
                                        abs(is_coord - min_gc_coord),
                                        abs(is_coord - max_gc_coord),
                                    ]
                                )
                            if dist < min_is_distance:
                                min_is_distance = dist
                    if min_is_distance == 1e100:
                        min_is_distance = "NA"

                    scaffold_length = scaffold_lens[scaffold]

                    row = [
                        gc,
                        gc_gbk,
                        genome,
                        scaffold,
                        scaffold_length,
                        scaffold_cds_count,
                        gc_cds_count,
                        gc_codoff_pval,
                        ratio_statistic,
                        gc_aai,
                        ribo_aai,
                        min_is_distance,
                        is_cds_prop,
                        vog_cds_prop,
                        mobsuite_cds_prop,
                    ]
                    tsv_outf.write("\t".join([str(x) for x in row]) + "\n")
                    num_rows += 1
        

        # Generate Excel spreadsheet
        writer = pd.ExcelWriter(result_xlsx_file, engine="xlsxwriter")
        workbook = writer.book
        dd_sheet = workbook.add_worksheet("Data Dictionary")
        dd_sheet.write(
            0,
            0,
            'Data Dictionary describing columns of "SALT Results" spreadsheet can be found on zol\'s Wiki page at:',
        )
        dd_sheet.write(
            1,
            0,
            "https://github.com/Kalan-Lab/zol/wiki/5.4-horizontal-or-lateral-transfer-assessment-of-gene-clusters-using-salt",
        )

        na_format = workbook.add_format(
            {"font_color": "#a6a6a6", "bg_color": "#FFFFFF", "italic": True}
        )
        header_format = workbook.add_format(
            {
                "bold": True,
                "text_wrap": True,
                "valign": "top",
                "fg_color": "#D7E4BC",
                "border": 1,
            }
        )

        numeric_columns = {
            "scaffold length (bp)",
            "scaffold CDS count",
            "total scaffold CDS count",
            "GC CDS count",
            "codoff empirical P-value",
            "GC AAI observed - expectation",
            "distance to IS-associated element",
            "scaffold CDS proportion IS-associated elements",
            "scaffold CDS proportion VOGs",
            "scaffold CDS proportion plasmid-associated",
            "GC AAI between genome and reference genome",
            "ribosomal protein AAI between genome and reference genome",
        }

        results_df = load_table_in_panda_data_frame(
            result_tsv_file, numeric_columns
        )
        results_df.to_excel(
            writer, sheet_name="SALT Results", index=False, na_rep="NA"
        )
        worksheet = writer.sheets["SALT Results"]
        worksheet.conditional_format(
            "A2:BA" + str(num_rows + 1),
            {
                "type": "cell",
                "criteria": "==",
                "value": '"NA"',
                "format": na_format,
            },
        )
        worksheet.conditional_format(
            "A1:BA1",
            {
                "type": "cell",
                "criteria": "!=",
                "value": "NA",
                "format": header_format,
            },
        )

        # codoff p - value
        worksheet.conditional_format(
            "H2:H" + str(num_rows + 1),
            {
                "type": "3_color_scale",
                "min_color": "#f07878",
                "mid_color": "#f7ca81",
                "max_color": "#f7de99",
                "min_value": 0.0,
                "mid_value": 0.05,
                "max_value": 1.0,
                "min_type": "num",
                "mid_type": "num",
                "max_type": "num",
            },
        )

        # diff
        worksheet.conditional_format(
            "I2:I" + str(num_rows + 1),
            {
                "type": "3_color_scale",
                "min_color": "#917967",
                "mid_color": "#FFFFFF",
                "max_color": "#ba8dc9",
                "min_value": -100.0,
                "mid_value": 0.0,
                "max_value": 100.0,
                "min_type": "num",
                "mid_type": "num",
                "max_type": "num",
            },
        )

        # AAI columns
        worksheet.conditional_format(
            "J2:J" + str(num_rows + 1),
            {
                "type": "2_color_scale",
                "min_color": "#dfedf5",
                "max_color": "#8fa9b8",
                "min_value": 0.0,
                "max_value": 100.0,
                "min_type": "num",
                "max_type": "num",
            },
        )
        worksheet.conditional_format(
            "K2:K" + str(num_rows + 1),
            {
                "type": "2_color_scale",
                "min_color": "#dfedf5",
                "max_color": "#8fa9b8",
                "min_value": 0.0,
                "max_value": 100.0,
                "min_type": "num",
                "max_type": "num",
            },
        )

        # dist to IS / transposon
        worksheet.conditional_format(
            "L2:L" + str(num_rows + 1),
            {
                "type": "2_color_scale",
                "min_color": "#ed87ad",
                "max_color": "#f5c6d8",
                "min_value": 0.0,
                "max_value": 10000.0,
                "min_type": "num",
                "max_type": "num",
            },
        )

        # IS / VOG / plasmid
        worksheet.conditional_format(
            "M2:M" + str(num_rows + 1),
            {
                "type": "2_color_scale",
                "min_color": "#98ad93",
                "max_color": "#78b56b",
                "min_value": 0.0,
                "max_value": 1.0,
                "min_type": "num",
                "max_type": "num",
            },
        )
        worksheet.conditional_format(
            "N2:N" + str(num_rows + 1),
            {
                "type": "2_color_scale",
                "min_color": "#98ad93",
                "max_color": "#78b56b",
                "min_value": 0.0,
                "max_value": 1.0,
                "min_type": "num",
                "max_type": "num",
            },
        )
        worksheet.conditional_format(
            "O2:O" + str(num_rows + 1),
            {
                "type": "2_color_scale",
                "min_color": "#98ad93",
                "max_color": "#78b56b",
                "min_value": 0.0,
                "max_value": 1.0,
                "min_type": "num",
                "max_type": "num",
            },
        )

        workbook.close()

    except Exception as e:
        sys.stderr.write(
            "Issues with creating the final salt XLSX spreadsheet.\n"
        )
        log_object.error(
            "Issues with creating the final salt XLSX spreadsheet."
        )
        sys.stderr.write(traceback.format_exc())
        log_object.error(traceback.format_exc())
        sys.exit(1)


# R functions
def cluster_heatmap_r(
    medlen_data_file,
    heatmap_data_file,
    pdf_file,
    height,
    width,
    rscript_path,
    log_object,
) -> None:
    try:

        with open(rscript_path, "w") as rph:
            rph.write("library(ggplot2)\n")
            rph.write("library(cowplot)\n\n")

            rph.write('medlen.data_file <- "' + medlen_data_file + '"\n')
            rph.write('heatmap.data_file <- "' + heatmap_data_file + '"\n')
            rph.write('pdf_file <- "' + pdf_file + '"\n')
            rph.write("height <- as.numeric(" + str(height) + ")\n")
            rph.write("width <- as.numeric(" + str(width) + ")\n\n")

            rph.write(
                'medlen.data <- read.table(medlen.data_file, header=T, sep="\\t")\n'
            )
            rph.write(
                'heatmap.data <- read.table(heatmap.data_file, header=T, sep="\\t")\n\n'
            )

            rph.write(
                'gg_ml <- ggplot(medlen.data, aes(x = \
        reorder(og, og_order), y = med_length)) + theme_classic() + xlab("") + \n'
            )
            rph.write(
                'ylab("Median Length\n(kbp)") + geom_bar(stat= \
        "identity", fill="black") + \n'
            )
            rph.write(
                "theme(axis.title.x= \
        element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())\n\n"
            )

            rph.write(
                "gg_hm <- ggplot(heatmap.data, aes(x = \
        reorder(og, og_order), y = genbank, fill=as.factor(og_presence), label=copy_count)) + \n"
            )
            rph.write(
                'theme_classic() + xlab("ortholog group IDs in Consensus Order") + ylab("") + geom_tile(color= \
        "white", show.legend=F) + geom_text() + \n'
            )
            rph.write(
                "theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + \n"
            )
            rph.write('scale_fill_manual(values=c("#FFFFFF", "#889cbd"))\n\n')

            rph.write("pdf(pdf_file, height=height, width=width)\n")
            rph.write(
                'print(plot_grid(gg_ml, gg_hm, ncol= \
        1, axis="l", align="v", rel_heights=c(1, 4)))\n'
            )
            rph.write("dev.off()\n")
            
    except Exception as e:
        sys.stderr.write("Issues with creating Rscript clusterHeatmap.R.\n")
        log_object.error("Issues with creating Rscript clusterHeatmap.R.")
        sys.stderr.write(traceback.format_exc())
        log_object.error(traceback.format_exc())
        sys.exit(1)


def phylo_heatmap_r(
    phylo_tree_file,
    heatmap_data_file,
    pdf_file,
    height,
    width,
    rscript_path,
    log_object,
) -> None:
    try:
        with open(rscript_path, "w") as rph:
            rph.write("library(ggtree)\n")
            rph.write("library(ggplot2)\n")
            rph.write("library(ape)\n")
            rph.write("library(dplyr)\n")
            rph.write("library(aplot)\n\n")

            rph.write('phylo.tree_file <- "' + phylo_tree_file + '"\n')
            rph.write('heatmap.data_file <- "' + heatmap_data_file + '"\n')
            rph.write('pdf_file <- "' + pdf_file + '"\n')
            rph.write("pdf.height <- as.numeric(" + str(height) + ")\n")
            rph.write("pdf.width <- as.numeric(" + str(width) + ")\n\n")

            rph.write("phylo.tree <- read.tree(phylo.tree_file)\n")
            rph.write(
                'heatmap.data <- read.table(heatmap.data_file, header=T, sep="\\t")\n\n'
            )

            rph.write("pdf(pdf_file, height=pdf.height, width=pdf.width)\n")
            rph.write("gg_tr <- ggtree(phylo.tree)\n")
            rph.write(
                "gg_hm <- ggplot(heatmap.data, aes(x= \
        query_prot_id, y=label, fill=bitscore)) + \n"
            )
            rph.write(
                'theme_classic() + scale_fill_gradient(low= \
        "grey", high="black", na.value="white") + \n'
            )
            rph.write(
                'xlab("Query proteins / Homolog groups") + ylab("") + geom_tile(color= \
        "white", show.legend=F) + \n'
            )
            rph.write(
                "theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))\n"
            )
            rph.write("gg_hm %>% insert_left(gg_tr, width=0.4)\n")
            rph.write("dev.off()\n")
            
    except Exception as e:
        sys.stderr.write("Issues with creating Rscript phyloHeatmap.R.\n")
        log_object.error("Issues with creating Rscript phyloHeatmap.R.")
        sys.stderr.write(traceback.format_exc())
        log_object.error(traceback.format_exc())
        sys.exit(1)


def plot_segments_r(
    que_info_file,
    seq_data_file,
    pdf_file,
    pdf_file2,
    height,
    width,
    rscript_path,
    log_object,
) -> None:
    try:
        with open(rscript_path, "w") as rph:
            rph.write("library(ggplot2)\n")
            rph.write("library(gridExtra)\n\n")

            rph.write('que.info.file <- "' + que_info_file + '"\n')
            rph.write('seg.data.file <- "' + seq_data_file + '"\n')
            rph.write('pdf_file <- "' + pdf_file + '"\n')
            rph.write('pdf_file2 <- "' + pdf_file2 + '"\n')
            rph.write("height <- as.numeric(" + str(height) + ")\n")
            rph.write("width <- as.numeric(" + str(width) + ")\n\n")

            rph.write(
                'naming.data <- read.table(que.info.file, header=T, sep="\\t")\n'
            )
            rph.write("pdf(pdf_file2, height=30, width=10)\n")
            rph.write("grid.table(naming.data)\n")
            rph.write("dev.off()\n\n")

            rph.write(
                'segment.data <- read.table(seg.data.file, header=T, sep="\\t")\n'
            )
            rph.write('colors <- c("#000000", "#FFFFFF")\n')
            rph.write('names(colors) <- c("True", "False")\n\n')

            rph.write("samples <- unique(segment.data$sample)\n")
            rph.write("pdf(pdf_file, height=height, width=width)\n")
            rph.write("for (s in samples) {\n")
            rph.write(
                "sample.segment.data <- segment.data[segment.data$sample==s,]\n"
            )
            rph.write("print(unique(sample.segment.data$segment_title))\n")
            rph.write(
                "g<-ggplot(sample.segment.data, aes(x= \
        reorder(gene, gene_order), y=sql_ratio, fill=identity, color=key)) + \n"
            )
            rph.write(
                'geom_bar(stat= \
        "identity", position="dodge", size=1.5) + geom_hline(yintercept=1.0, color="blue", linetype=2) + \n'
            )
            rph.write(
                'ylab("CDS to Query Length Ratio") + facet_grid(.~segment_title, space= \
        "free", scales="free",  labeller = label_wrap_gen(width = 50, multi_line = TRUE)) + \n'
            )
            rph.write(
                'xlab("CDS Classifications") + theme_bw() + theme(legend.position= \
        "bottom", axis.text.x = element_text(angle = 45, vjust=1, hjust=1)) + \n'
            )
            rph.write(
                'scale_color_manual(values= \
        colors) + scale_fill_gradient(limits=c(0.0, 100.0), breaks=c(0.0, 50.0, 100.0), low="grey", high="red")\n'
            )
            rph.write("print(g)\n")
            rph.write("}\n")
            rph.write("dev.off()\n\n")

        
    except Exception as e:
        sys.stderr.write("Issues with creating Rscript plotSegments.R.\n")
        log_object.error("Issues with creating Rscript plotSegments.R.")
        sys.stderr.write(traceback.format_exc())
        log_object.error(traceback.format_exc())
        sys.exit(1)


def plot_tiny_aair(info_file, pdf_file, rscript_path, log_object) -> None:
    try:
        with open(rscript_path, "w") as rph:
            rph.write("library(ggplot2)\n")
            rph.write("library(gridExtra)\n")

            rph.write('info.file <- "' + info_file + '"\n')
            rph.write('pdf_file <- "' + pdf_file + '"\n\n')

            rph.write('info.dat <- read.table(info.file, header=T, sep="\\t")\n\n')

            rph.write("pdf(pdf_file, height=10, width=10)\n")
            rph.write(
                "ggplot(info.dat, aes(x= \
        AAI, y=Prop_Genes_Found, color=Mean_Syntenic_Correlation)) + \n"
            )
            rph.write(
                'geom_point(alpha= \
        0.7) + theme_bw() + scale_color_gradient(low="#e6ffbd", high="#1754b0") + \n'
            )
            rph.write(
                'guides(color=guide_legend("Syntenic Correlation to Query")) + \n'
            )
            rph.write(
                'xlab("Average amino acid identity") + ylab("Proportion of query proteins with match")\n' \
            )
            rph.write("dev.off()\n")

        
    except Exception as e:
        sys.stderr.write("Issues with creating Rscript plotTinyAAI.R.\n")
        log_object.error("Issues with creating Rscript plotTinyAAI.R.")
        sys.stderr.write(traceback.format_exc())
        log_object.error(traceback.format_exc())
        sys.exit(1)


def generate_syntenic_visual_r(input_file, pdf_file, height, width, rscript_path) -> None: 
    try:
        with open(rscript_path, "w") as rph:
            rph.write("library(ggplot2)\n")
            rph.write("library(gggenes)\n")

            rph.write('input.file <- "' + input_file + '"\n')
            rph.write('pdf.file <- "' + pdf_file + '"\n')
            rph.write("height <- as.numeric(" + str(height) + ")\n")
            rph.write("width <- as.numeric(" + str(width) + ")\n\n")

            rph.write(
                'input.data <- read.table(file=input.file, sep="\\t", header=T)\n\n'
            )

            rph.write("pdf(pdf.file, height=height, width=width)\n")
            rph.write(
                'ggplot(input.data, aes(xmin= \
        Start, xmax = End, y = "", forward = Direction, label=SC)) + \n'
            )
            rph.write("geom_gene_arrow(aes(fill=Metric)) + theme_classic() + \n")
            rph.write(
                'scale_fill_gradient2(low= \
        "#e05c74", mid="#f2f2f2", high="#2087b3", breaks =c(-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0),\n'
            )
            rph.write(
                'labels= \
        c("", "-2", "", "0", "", "2", ""), limits=c(-3, 3), na.value="grey50", guide="colourbar", aesthetics="fill") + \n'
            )
            rph.write(
                'geom_gene_label(align= \
        "centre", min.size=5) + theme(legend.position="bottom")\n'
            )
            rph.write("dev.off()\n")

        
    except Exception as e:
        sys.stderr.write(
            f"Issues with creating / running Rscript {rscript_path}.\n"
        )
        sys.stderr.write(traceback.format_exc())
        sys.exit(1)


def generate_salt_gc_vs_ribo_aai_plot(
    rscript_path, input_data_file, pdf_file, log_object
) -> None:
    try:
        with open(rscript_path, "w") as rph:
            rph.write("library(ggplot2)\n\n")

            rph.write('input.data_file <- "' + input_data_file + '"\n')
            rph.write('pdf_file <- "' + pdf_file + '"\n\n')

            rph.write('dat <- read.table(input.data_file, header=T, sep="\\t")\n')
            rph.write("pdf(pdf_file, height=5, width=5)\n")
            rph.write(
                "ggplot(dat, aes(x=ribo_aai, y=gc_aai)) + geom_point(alpha=0.7) +\n"
            )
            rph.write(
                'geom_smooth(method="lm", formula= y~x, linetype=2, color="red") +\n'
            )
            rph.write(
                'geom_abline(slope= \
        1, yintercept=0, linetype=2, color="grey") + theme_bw() +\n'
            )
            rph.write('xlab("Ribosomal Protein AAI") + ylab("Gene Cluster AAI")\n')
            rph.write("dev.off()\n")
    except Exception as e:
        sys.stderr.write(
            "Issues with creating Rscript generateSaltGCVsRiboAAIPlot.R.\n"
        )
        log_object.error(
            "Issues with creating Rscript generateSaltGCVsRiboAAIPlot.R.\n"
        )
        sys.stderr.write(traceback.format_exc())
        log_object.error(traceback.format_exc())
        sys.exit(1)


def generate_nj_tree(rscript_path, input_dist_file, output_tree_file, log_object) -> None: 
    try:
        with open(rscript_path, "w") as rph:
            rph.write("library(ape)\n\n")
            rph.write(f'input.dist.file <- "{input_dist_file}"\n')
            rph.write(f'output.nwk.file <- "{output_tree_file}"\n\n')
            rph.write(
                'dat <- read.table(input.dist.file, header=T, sep="\\t", row.names=1)\n'
            )
            rph.write("d <- as.dist(as.matrix(dat))\n")
            rph.write("njt <- nj(d)\n")
            rph.write("write.tree(njt, file=output.nwk.file)\n")

    except Exception as e:
        sys.stderr.write("Issues with creating Rscript generateNjTree.R.\n")
        log_object.error("Issues with creating Rscript generateNjTree.R.")
        sys.stderr.write(traceback.format_exc())
        log_object.error(traceback.format_exc())
        sys.exit(1)

def create_fake_diamond_linclust_file(fasta_file: str, diamond_linclust_cluster_file: str, log_object) -> None:
    """
    Description:
    Creates a fake DIAMOND linclust cluster file where each protein is its own representative.
    This is used when DIAMOND linclust clustering is skipped but downstream analysis expects
    a cluster file format.
    ********************************************************************************************************************
    Parameters:
    - diamond_linclust_cluster_file: The path to the DIAMOND linclust cluster file to create.
    - log_object: A logging object.
    ********************************************************************************************************************
    Returns:
    - None
    ********************************************************************************************************************
    """
    try:

        # Check if the FASTA file exists
        if not os.path.isfile(fasta_file):
            msg = f"FASTA file not found: {fasta_file}"
            log_object.error(msg)
            sys.stderr.write(msg + '\n')
            sys.exit(1)
        
        msg = f"Creating fake DIAMOND linclust cluster file from: {fasta_file}"
        log_object.info(msg)
        sys.stdout.write(msg + '\n')
        
        # Create the fake cluster file
        with open(diamond_linclust_cluster_file, 'w') as cluster_handle:
            cluster_id = 0
            with open(fasta_file) as fasta_handle:
                for rec in SeqIO.parse(fasta_handle, 'fasta'):
                    cluster_handle.write(f"{rec.id}\t{rec.id}\n")
        
        log_object.info(f"Successfully created fake DIAMOND linclust cluster file: {diamond_linclust_cluster_file}")
        
    except Exception as e:
        msg = f"Error creating fake DIAMOND linclust cluster file: {e}"
        log_object.error(msg)
        sys.stderr.write(msg + '\n')
        sys.exit(1)
