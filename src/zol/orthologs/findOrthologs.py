#! /usr/bin/env python3

"""
Program: findOrthologs.py
Author: Rauf Salamzade
Kalan Lab
UW Madison, Department of Medical Microbiology and Immunology
"""

# BSD 3 - Clause License
#
# Copyright (c) 2023-2025, Kalan - Lab
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
#    list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and / or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from collections import defaultdict
from operator import attrgetter, itemgetter
import argparse
import os
from typing import Dict, List, Optional, Union, Any, Tuple
import shutil
import sys
import traceback
from Bio import SeqIO
from zol import util

zol_exec_directory = str(os.getenv("ZOL_EXEC_PATH")).strip()
split_diamond_results_prog = None
rbh_prog = None
if zol_exec_directory != "None":
    try:
        zol_exec_directory = os.path.abspath(zol_exec_directory) + "/"
        rbh_prog = zol_exec_directory + "runRBH"
        split_diamond_results_prog = zol_exec_directory + "splitDiamondResults"
    except Exception as e:
        pass
if (
    rbh_prog == None
    or split_diamond_results_prog == None
    or not os.path.isfile(rbh_prog)
    or not os.path.isfile(split_diamond_results_prog)
):
    sys.stderr.write(
        "Issues in setup of the zol-suite (in findOrthologs.py) - please describe your installation process and post an issue on GitHub!\n" \
    )
    sys.exit(1)


def create_parser() -> argparse.Namespace:
    """Parse arguments"""
    parser = argparse.ArgumentParser(
        description="""
	Program: findOrthologs.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

	Performs reflexive alignment of proteins across different protein-sets / genomes and 
    determines orthologs, co-orthologs, and in-paralogs.
	""",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "-p",
        "--proteome_dir",
        help="Path to directory with genomes. Should end with .faa,\n"
             ".faa, or .fasta",
        required=True,
    )
    parser.add_argument(
        "-o", "--output_dir", help="Output directory.", required=True
    )
    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        help="E-value cutoff for determining orthologs.",
        required=False,
        default=0.001,
    )
    parser.add_argument(
        "-i",
        "--identity",
        type=float,
        help="Percent identity cutoff for determining orthologs.",
        required=False,
        default=30.0,
    )
    parser.add_argument(
        "-q",
        "--coverage",
        type=float,
        help="Bi-directional coverage cutoff for determining orthologs.",
        required=False,
        default=50.0,
    )
    parser.add_argument(
        "-mi",
        "--mcl-inflation",
        type=float,
        help= "The inflation parameter value for MCL. If -1 will\n"
              "use single-linkage clustering.",
        required=False,
        default=1.5,
    )
    parser.add_argument(
        "-dco",
        "--dc-orthogroup",
        action="store_true",
        help="Infer ortholog groups using diamond cluster.",
        required=False,
        default=False,
    )
    parser.add_argument(
        "-dcp",
        "--dc-params",
        help= "Parameters for performing diamond cluster based ortholog group\n"
              "clustering if requested via --dc-orthogroup.\n"
              "[Default is \"--approx-id 50 --mutual-cover 25\"].",
        required=False,
        default="--approx-id 50 --mutual-cover 25"
    )
    parser.add_argument(
        "-c",
        "--threads",
        type=int,
        help="Maximum number of threads to use. Default is 1.",
        default=1,
        required=False,
    )
    parser.add_argument(
        "-m",
        "--memory",
        type=int,
        help="Maximum memory in GB to use. Default is 4.",
        default=4,
        required=False,
    )
    args = parser.parse_args()
    return args


def refactor_proteomes_function(inputs) -> Tuple[str, Optional[str]]:
    (
        sample_name,
        prot_file,
        updated_prot_file,
        original_naming_file,
    ) = inputs
    try:
        with open(updated_prot_file, "w") as updated_prot_handle:
            original_naming_handle = open(original_naming_file, "w")
            with open(prot_file) as opf:
                for i, rec in enumerate(SeqIO.parse(opf, "fasta")):
                    original_naming_handle.write(
                        sample_name
                        + f"\t{prot_file}\t"
                        + "Protein_"
                        + str(i)
                        + "\t"
                        + rec.id + "\n"
                    )
                    updated_prot_handle.write(
                        f">{sample_name}|"
                        + "Protein_"
                        + str(i)
                        + "\n"
                        + str(rec.seq)
                        + "\n"
                    )
            original_naming_handle.close()
        return ('success', None)
    except Exception as e:
        error_msg = f"Error with refactoring proteome file {prot_file}: {str(e)}"
        return ('error', error_msg)


def one_vs_all_parse(inputs) -> Tuple[str, Optional[str]]:
    (
        sample,
        samp_algn_res,
        samp_forth_res,
        identity_cutoff,
        coverage_cutoff,
    ) = inputs
    try:
        rbh_cmd = [
            rbh_prog,
            samp_algn_res,
            str(identity_cutoff),
            str(coverage_cutoff),
            sample,
            ">",
            samp_forth_res,
        ]
        util.run_cmd_via_subprocess(rbh_cmd, check_files = [samp_forth_res])

        util.remove_file(samp_algn_res)
        return ('success', None)
    except Exception as e:
        error_msg = f"Error in one-vs-all parsing for sample {sample}: {str(e)}"
        return ('error', error_msg)

def create_final_results(
    concat_result_file,
    proteome_listing_file,
    original_naming_dir,
    result_tab_file,
    result_mat_file,
    log_object,
) -> None:
    try:
        samples = []
        with open(proteome_listing_file) as oplf:
            for line in oplf:
                line = line.strip()
                ls = line.split("\t")
                samples.append(".".join(ls[1].split("/")[-1].split(".")[:-1]))

        mapping: Dict[str, Any] = {}
        for f in os.listdir(original_naming_dir):
            name_file = original_naming_dir + f
            with open(name_file) as onf:
                for line in onf:
                    line = line.strip()
                    rn_sample, og_sample, rn_protein, og_protein = line.split(
                        "\t"
                    )
                    og_sample = ".".join(
                        og_sample.split("/")[-1].split(".")[:-1]
                    )
                    mapping[rn_sample + "|" + rn_protein] = tuple(
                        [og_sample, og_protein]
                    )

        with open(result_tab_file, "w") as result_tab_handle:
            with open(result_mat_file, "w") as result_mat_handle:
                sorted_samples = sorted(samples)
                result_mat_handle.write("Sample\t" + "\t".join(sorted_samples) + "\n")
                with open(concat_result_file) as ocrf:
                    for line in ocrf:
                        line = line.strip()
                        ls = line.split("\t")
                        og = ls[0]
                        prots = ls[1].split(" ")
                        samp_prots = defaultdict(set)
                        for p in sorted(prots):
                            s, op = mapping[p]
                            samp_prots[s].add(op)
                            result_tab_handle.write(f"{og}\t{s}\t{op}\n")
                        printlist = [og]
                        for s in sorted_samples:
                            printlist.append(", ".join(sorted(samp_prots[s])))
                        result_mat_handle.write("\t".join(printlist) + "\n")
        
    except Exception as e:
        log_object.error(e)
        sys.stderr.write(traceback.format_exc())
        sys.stderr.write(str(e))
        sys.exit(1)


def create_final_results_diamond_cluster(
    diamond_cluster_file, result_tab_file, result_mat_file, log_object
) -> None:
    try:
        cluster_id = None
        clust_proteins = defaultdict(lambda: defaultdict(set))
        clust_protein_counts = defaultdict(int)
        protein_to_clust: Dict[str, Any] = {}
        samples = set([])
        # Parse diamond cluster output
        cluster_representatives = {}  # Maps cluster rep to cluster members
        with open(diamond_cluster_file) as occf:
            for line in occf:
                line = line.strip()
                if not line:
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    continue
                rep_id = parts[0]
                mem_id = parts[1]
                
                if rep_id not in cluster_representatives:
                    cluster_representatives[rep_id] = set()
                cluster_representatives[rep_id].add(mem_id)
        
        # Process clusters
        for cluster_id, members in cluster_representatives.items():
            for lt in members:
                sample = lt.split("|")[0]
                samples.add(sample)
                protein_to_clust[lt] = cluster_id
                clust_proteins[cluster_id][sample].add(lt)
                clust_protein_counts[cluster_id] += 1

        sorted_samples = sorted(list(samples))

        result_tab_handle = open(result_tab_file, "w")
        result_mat_handle = open(result_mat_file, "w")
        result_mat_handle.write("Sample\t" + "\t".join(sorted_samples) + "\n")

        for ci, c in enumerate(
            sorted(
                clust_protein_counts.items(), key=itemgetter(1), reverse=True
            )
        ):
            og_id = "OG_" + str(ci + 1)
            mat_row = [og_id]
            for s in sorted_samples:
                samp_cluster_lts = []
                for l in clust_proteins[c[0]][s]:
                    samp_cluster_lts.append(l)
                    result_tab_handle.write(f"{og_id}\t{s}\t{l}\n")
                mat_row.append(", ".join(samp_cluster_lts))
            result_mat_handle.write("\t".join(mat_row) + "\n")
        
        result_tab_handle.close()
        result_mat_handle.close()
    except Exception as e:
        log_object.error(e)
        sys.stderr.write(traceback.format_exc())
        sys.stderr.write(str(e))
        sys.exit(1)


def find_orthologs() -> None:
    """
    Void function which runs primary workflow for program.
    """

    """
	PARSE INPUTS
	"""
    myargs = create_parser()

    proteome_dir = os.path.abspath(myargs.proteome_dir) + "/"
    outdir = os.path.abspath(myargs.output_dir) + "/"
    threads = myargs.threads
    coverage_cutoff = myargs.coverage
    identity_cutoff = myargs.identity
    evalue_cutoff = myargs.evalue
    mcl_inflation = myargs.mcl_inflation
    diamond_cluster_orthogroup_flag = myargs.dc_orthogroup
    diamond_cluster_params = myargs.dc_params
    memory = myargs.memory

    if not os.path.isdir(outdir):
        os.system(f"mkdir {outdir}")
    else:
        sys.stderr.write(
            "Note, output directory exists already, will perform steps that were not successfully completed!\n"
        )

    # create logging object
    log_file = outdir + "Progress.log"
    log_object = util.create_logger_object(log_file)

    version = util.get_version()
    sys.stdout.write(f"Running version: {version}\n")
    log_object.info(f"Running version: {version}")

    parameters_file = outdir + "Command_Issued.txt"
    sys.stdout.write(
        f"Appending command issued for future records to: {parameters_file}\n"
    )
    sys.stdout.write(f"Logging more details at: {log_file}\n")
    log_object.info("\nNEW RUN!!!\n**************************************")
    log_object.info(f"Running version {version}")
    log_object.info(
        f"Appending command issued for future records to: {parameters_file}"
    )

    with open(parameters_file, "a+") as parameters_handle:
        parameters_handle.write(" ".join(sys.argv) + "\n") 
    
    result_tab_file = outdir + "Orthogroups.Listing.tsv"
    result_mat_file = outdir + "Orthogroups.tsv"

    checkpoint_dir = outdir + "Checkpoint_Files/"
    if not os.path.isdir(checkpoint_dir):
        os.system(f"mkdir {checkpoint_dir}")

    if not diamond_cluster_orthogroup_flag:
        local_proteome_dir = outdir + "Proteomes/"
        if not os.path.isdir(local_proteome_dir):
            os.system(f"mkdir {local_proteome_dir}")

        proteome_name_dir = outdir + "Original_Naming_of_Proteomes/"
        if not os.path.isdir(proteome_name_dir):
            os.system(f"mkdir {proteome_name_dir}")

        msg = "--------------------\nStep 0\n--------------------\nProcessing proteomes"
        sys.stdout.write(msg + "\n")
        log_object.info(msg)

        proteome_listing_file = outdir + "Proteome_Listing.txt"
        step0_checkpoint_file = checkpoint_dir + "Step0_Checkpoint.txt"
        if not os.path.isfile(step0_checkpoint_file):
            try:
                refactor_proteomes = []
                with open(proteome_listing_file, "w") as proteome_listing_handle:
                    proteome_counter = 0
                    assert os.path.isdir(proteome_dir)
                    for f in sorted(os.listdir(proteome_dir)):
                        suffix = f.split(".")[-1]
                        if not suffix in set(["fa", "faa", "fasta"]):
                            continue
                        proteome_counter += 1
                        prot_file = proteome_dir + f
                        if not util.is_fasta(prot_file):
                            log_object.warning(
                                f"File {prot_file} not in valid FASTA format and ignored."
                            )
                            continue
                        sample_name = "Proteome_" + str(proteome_counter)
                        updated_prot_file = (
                            local_proteome_dir + (sample_name + ".faa")
                        )
                        original_naming_file = (
                            proteome_name_dir + (sample_name + ".txt")
                        )
                        refactor_proteomes.append(
                            [
                                sample_name,
                                prot_file,
                                updated_prot_file,
                                original_naming_file,
                            ]
                        )
                        proteome_listing_handle.write(f"{sample_name}\t{prot_file}\t{updated_prot_file}\n")

                # Use robust error handling for proteome refactoring
                result_summary = util.robust_multiprocess_executor(
                    worker_function=refactor_proteomes_function,
                    inputs=refactor_proteomes,
                    pool_size=threads,
                    error_strategy="report_and_continue",  # Continue even if some proteomes fail
                    log_object=log_object,
                    description="proteome refactoring"
                )

                success_prop = result_summary['success_count'] / result_summary['total_processed']
                if success_prop != 1.0:
                    msg = f"Issues with proteome refactoring for at least one proteome. Exiting now ..."
                    sys.stderr.write(msg + '\n')
                    log_object.error(msg)
                    sys.exit(1)

            except Exception as e:
                log_object.error(e)
                log_object.error(
                    "Difficulties with parsing directory of proteomes."
                )
                sys.exit(1)
            
            os.system(f"touch {step0_checkpoint_file}")

        """
		START WORKFLOW
		"""

        # Step 1: Concatenate all proteins into single multi-FASTA file
        msg = "--------------------\nStep 1\n--------------------\nConcatenating proteins into multi-FASTA"
        sys.stdout.write(msg + "\n")
        log_object.info(msg)
        concat_faa = outdir + "All_Proteins.faa"
        if not os.path.isfile(concat_faa):
            with open(concat_faa, "w") as cf_handle:
                for f in os.listdir(local_proteome_dir):
                    with open(local_proteome_dir + f) as olf:
                        for rec in SeqIO.parse(olf, "fasta"):
                            cf_handle.write(f">{rec.id}\n{rec.seq}\n")

        ########### Perform ortholog grouping using reflexive DIAMOND alignment + InParanoid-like algorithm

        # Step 2: Perform reflexive alignment via Diamond
        msg = "--------------------\nStep 2\n--------------------\nRunning reflexive alignment using Diamond\n"
        sys.stdout.write(msg + "\n")
        log_object.info(msg)
        step2_checkpoint_file = checkpoint_dir + "Step2_Checkpoint.txt"

        align_res_dir = outdir + "Alignments/"
        forth_res_dir = outdir + "Ortholog_Inference_Results/"

        sample_prots = defaultdict(set)
        sample_inputs = set([])
        sample_listing_file = outdir + "sample_listing.txt"
        with open(sample_listing_file, "w") as sample_listing_handle:
            visited = set([])
            try:
                with open(concat_faa) as ocf:
                    for rec in SeqIO.parse(ocf, "fasta"):
                        sample = rec.id.split("|")[0]
                        sample_prots[sample].add(rec.id)
                        samp_algn_res = align_res_dir + (sample + ".out")
                        samp_forth_res = forth_res_dir + (sample + ".abc-like")
                        sample_inputs.add(
                            tuple([sample, samp_algn_res, samp_forth_res])
                        )
                        if not sample in visited:
                            sample_listing_handle.write(f"{sample}\t{align_res_dir + sample}.out\n")
                        visited.add(sample)
            except Exception as e:
                log_object.error(e)
                raise RuntimeError(e)

        if not os.path.isfile(step2_checkpoint_file):
            os.system(f"mkdir {align_res_dir}")
            db_path = outdir + "All_Proteins.dmnd"
            diamond_makedb_cmd = [
                "diamond",
                "makedb",
                "--ignore-warnings",
                "--in",
                concat_faa,
                "--db",
                db_path,
                "--threads",
                str(threads)
            ]

            util.run_cmd_via_subprocess(diamond_makedb_cmd, check_files = [db_path])

            alignment_result_file = outdir + "Reflexive_Alignment.out"
            diamond_blastp_cmd = ["diamond", "blastp", "--ignore-warnings", "-d", db_path,
                                  "-q", concat_faa, "-o", alignment_result_file, "--outfmt",
                                  "6", "qseqid", "sseqid", "pident", "length", "mismatch",
                                  "gapopen", "qstart", "qend", "sstart", "send", "evalue",
                                  "bitscore", "qcovhsp", "-k0", "--very-sensitive", "-e",
                                  str(evalue_cutoff), "--masking", "0", "-p", str(threads)]
            util.run_cmd_via_subprocess(diamond_blastp_cmd, check_files = [alignment_result_file])

            split_diamond_cmd = [
                split_diamond_results_prog,
                alignment_result_file,
                sample_listing_file,
            ]
            
            util.run_cmd_via_subprocess(split_diamond_cmd, check_files = [alignment_result_file])
            os.system(f"touch {step2_checkpoint_file}")

        msg = "--------------------\nStep 3\n--------------------\nRunning RBH.\n"
        sys.stdout.write(msg + "\n")
        log_object.info(msg)
        scaled_find_orthos_result_file = outdir + "All.normalized.abc-like"
        step3_checkpoint_file = checkpoint_dir + "Step3_Checkpoint.txt"
        if not os.path.isfile(step3_checkpoint_file):
            os.system(f"mkdir {forth_res_dir}")

            find_orthos_result_file = outdir + "All.abc-like"

            parallelize_inputs = []
            for inp_tup in sample_inputs:
                parallelize_inputs.append(
                    list(inp_tup)
                    + [identity_cutoff, coverage_cutoff]
                )

            # Use robust error handling for one-vs-all parsing
            result_summary = util.robust_multiprocess_executor(
                worker_function=one_vs_all_parse,
                inputs=parallelize_inputs,
                pool_size=threads,
                error_strategy="report_and_stop",  # Continue even if some comparisons fail
                log_object=log_object,
                description="one-vs-all ortholog parsing"
            )

            success_prop = result_summary['success_count'] / result_summary['total_processed']
            if success_prop != 1.0:
                msg = f"Issues with one-vs-all ortholog parsing for at least one sample. Exiting now ..."
                sys.stderr.write(msg + '\n')
                log_object.error(msg)
                sys.exit(1)

            os.system(
                f"cat {forth_res_dir}* > {find_orthos_result_file}"
            )

            qh_pairs_accounted = set([])
            qshs_rbh = defaultdict(int)
            qshs_rbh_sum = defaultdict(float)
            with open(scaled_find_orthos_result_file, "w") as scaled_find_orthos_result_handle:
                with open(find_orthos_result_file) as osforf:
                    for line in osforf:
                        line = line.strip()
                        query, query_samp, hit, hit_samp, bs = line.split("\t")
                        qh_pair = tuple(sorted([query, hit]))
                        if qh_pair in qh_pairs_accounted:
                            continue
                        qh_pairs_accounted.add(qh_pair)
                        sample_pair = tuple(sorted([query_samp, hit_samp]))
                        qshs_rbh[sample_pair] += 1
                        qshs_rbh_sum[sample_pair] += float(bs)

                qshs_rbh_avg: Dict[str, Any] = {}
                for sp in qshs_rbh:
                    qshs_rbh_avg[sp] = qshs_rbh_sum[sp] / qshs_rbh[sp]

                qh_pairs_accounted = set([])
                with open(find_orthos_result_file) as osforf:
                    for line in osforf:
                        line = line.strip()
                        query, query_samp, hit, hit_samp, bs = line.split("\t")
                        qh_pair = tuple(sorted([query, hit]))
                        if qh_pair in qh_pairs_accounted:
                            continue
                        qh_pairs_accounted.add(qh_pair)
                        if query_samp == hit_samp and mcl_inflation > 0:
                            scaled_find_orthos_result_handle.write(
                                query
                                + f"\t{hit}\t"
                                + str(float(bs) / 100.0)
                                + "\n"
                            )
                        else:
                            sample_pair = tuple(sorted([query_samp, hit_samp]))
                            if mcl_inflation > 0:
                                scaled_find_orthos_result_handle.write(f"{query}\t{hit}\t{float(bs) / qshs_rbh_avg[sample_pair]}\n") # type: ignore
                            elif mcl_inflation == -1:
                                scaled_find_orthos_result_handle.write(f"{query}\t{hit}\n")
            

        msg = "--------------------\nStep 4\n--------------------\nRun MCL or single-linkage clustering to determine OrthoGroups.\n"
        sys.stdout.write(msg + "\n")
        log_object.info(msg)
        tmp_result_file = outdir + "Tmp_Results.txt"
        step4_checkpoint_file = checkpoint_dir + "Step4_Checkpoint.txt"
        if not os.path.isfile(step4_checkpoint_file):
            cluster_result_file = outdir + "MCL_Cluster_Results.txt"
            if mcl_inflation > 0:
                find_clusters_cmd = [
                    "mcl",
                    scaled_find_orthos_result_file,
                    "--abc",
                    "-I",
                    str(mcl_inflation),
                    "-te",
                    str(threads),
                    "-o",
                    cluster_result_file,
                ]
            elif mcl_inflation == -1:
                find_clusters_cmd = [
                    "slclust",
                    "<",
                    scaled_find_orthos_result_file,
                    ">",
                    cluster_result_file,
                ]
            else:
                msg = "MCL inflation parameter value is not valid. It is less than 0 but not -1 to indicate single-linkage clustering requested."
                sys.stderr.write(msg + "\n")
                log_object.error(msg)
                sys.exit(1)

            util.run_cmd_via_subprocess(find_clusters_cmd, check_files = [cluster_result_file])

            with open(tmp_result_file, "w") as result_handle:
                subclust_id = 1
                clustered = set([])
                with open(cluster_result_file) as omrf:
                    for line in sorted([l for l in omrf.readlines()]):
                        line = line.strip()
                        ls = line.split()
                        clustered = clustered.union(set(ls))
                        result_handle.write(f"OG_{subclust_id}\t{' '.join(ls)}\n")
                        subclust_id += 1

                with open(concat_faa) as ocf:
                    fasta_dict = SeqIO.to_dict(SeqIO.parse(ocf, "fasta"))
                    for rec in sorted(fasta_dict.values(), key=attrgetter("id")):
                        if not rec.id in clustered:
                            result_handle.write(f"OG_{subclust_id}\t{rec.id}\n")
                            subclust_id += 1
            
            util.remove_file(alignment_result_file)
            util.remove_file(cluster_result_file)
            shutil.rmtree(forth_res_dir)
            shutil.rmtree(align_res_dir)

            os.system(f"touch {step4_checkpoint_file}")

        step5_checkpoint_file = checkpoint_dir + "Step5_Checkpoint.txt"
        msg = "--------------------\nStep 5\n--------------------\nCreate final result files!\n"
        sys.stdout.write(msg + "\n")
        log_object.info(msg)

        if not os.path.isfile(step5_checkpoint_file):
            create_final_results(
                tmp_result_file,
                proteome_listing_file,
                proteome_name_dir,
                result_tab_file,
                result_mat_file,
                log_object,
            )
            os.system(f"touch {step5_checkpoint_file}")

    else:
        """
        START WORKFLOW
        """

        # Step 1: Concatenate all proteins into single multi - FASTA file
        msg = "--------------------\nStep 1\n--------------------\nConcatenating proteins into multi-FASTA\n"
        sys.stdout.write(msg + "\n")
        log_object.info(msg)

        concat_faa = outdir + "All_Proteins.faa"
        if not os.path.isfile(concat_faa):
            with open(concat_faa, "w") as cf_handle:
                for f in os.listdir(proteome_dir):
                    with open(proteome_dir + f) as olf:
                        for rec in SeqIO.parse(olf, "fasta"):
                            cf_handle.write(f">{rec.id}\n{rec.seq}\n")
                
        ########### Perform protein clustering using DIAMOND cluster

        # Step 2: Run DIAMOND cluster
        msg = "--------------------\nStep 2\n--------------------\nRunning DIAMOND cluster to determine protein clusters\n"
        sys.stdout.write(msg + "\n")
        log_object.info(msg)

        step2_checkpoint_file = checkpoint_dir + "DIAMOND_CLUSTER_Step2_Checkpoint.txt"

        diamond_cluster_file = outdir + "diamond_clusters.tsv"
        diamond_refined_cluster_file = outdir + "diamond_refined_clusters.tsv"
        if not os.path.isfile(step2_checkpoint_file):
            util.diamond_cluster_and_recluster(concat_faa, diamond_cluster_file, diamond_refined_cluster_file,
                                               log_object, threads=threads, diamond_params=diamond_cluster_params,
                                               mem=memory)
            os.system(f"touch {step2_checkpoint_file}")

        msg = "--------------------\nStep 3\n--------------------\nCreate final result files!\n"
        sys.stdout.write(msg + "\n")
        log_object.info(msg)
        step3_checkpoint_file = checkpoint_dir + "DIAMOND_CLUSTER_Step3_Checkpoint.txt"

        if not os.path.isfile(step3_checkpoint_file):
            create_final_results_diamond_cluster(
                diamond_refined_cluster_file, result_tab_file, result_mat_file, log_object
            )
            os.system(f"touch {step3_checkpoint_file}")

    # DONE!
    msg = f"--------------------\nDONE!\n--------------------\nOrthogroup by sample matrix can be found at: {result_mat_file}\n"
    sys.stdout.write(msg + "\n")
    log_object.info(msg)


if __name__ == "__main__":
    find_orthologs()
