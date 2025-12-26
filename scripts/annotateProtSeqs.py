#!/usr/bin/env python3

"""
Program: annotateProtSeqs.py
Author: Rauf Salamzade
Kalan Lab
UW Madison, Department of Medical Microbiology and Immunology
"""

# BSD 3-Clause License
#
# Copyright (c) 2023-2025, Kalan-Lab
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
#    and/or other materials provided with the distribution.
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
 # BSD 3-Clause License

import argparse
import os
from typing import Dict, List, Optional, Union, Any, Tuple
import sys

from Bio import SeqIO

from zol import util, zol


def create_parser() -> argparse.Namespace:
    """Parse arguments"""
    parser = argparse.ArgumentParser(
        description="""
	Program: annotateProtSeqs.py
	Author: Rauf Salamzade
	Affiliation: Kalan Lab, UW Madison, Department of Medical Microbiology and Immunology

    	Simple script to annotate a protein FASTA using zol's standard approach taken for consensus sequences for ortholog groups.
	""",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--protein_faa",
        type=str,
        help="Path to FASTA of protein sequences.",
        required=True,
    )
    parser.add_argument(
        "-o", "--output_dir", help="Path to output directory.", required=True
    )
    parser.add_argument(
        "-c",
        "--threads",
        type=int,
        help="The number of threads to use [Default is 1].",
        required=False,
        default=1,
    )
    args = parser.parse_args()

    return args


def annotate_prot_seqs() -> None:
    """
    PARSE INPUTS
    """
    myargs = create_parser()

    protein_faa = myargs.protein_faa # type: ignore
    outdir = os.path.abspath(myargs.output_dir) + "/" # type: ignore
    threads = myargs.threads # type: ignore

    try:
        assert util.is_fasta(protein_faa)
        protein_faa = os.path.abspath(protein_faa)
    except Exception as e:
        sys.stderr.write(f"Could not validate {protein_faa} as a FASTA file")

    util.setup_ready_directory([outdir])

    # create logging object
    log_file = outdir + "Progress.log"
    log_object = util.create_logger_object(log_file)

    # perform annotate consensus sequences
    annotations = zol.annotate_consensus_sequences(
        protein_faa, outdir, log_object, threads=threads
    )

    header = [
        "Protein ID",
        "KO Annotation (E-value)",
        "PGAP Annotation (E-value)",
        "PaperBLAST Annotation (E-value)",
        "CARD Annotation (E-value)",
        "IS Finder (E-value)",
        "MIBiG Annotation (E-value)",
        "VOG Annotation (E-value)",
        "VFDB Annotation (E-value)",
        "Pfam Domains",
        "Protein Sequence",
    ]

    seqs: Dict[str, Any] = {}
    with open(protein_faa) as opf:
        for rec in SeqIO.parse(opf, "fasta"):
            seqs[rec.id] = str(rec.seq)

    annotation_tsv_file = outdir + "Annotation_Results.tsv"
    with open(annotation_tsv_file, "w") as frt_handle:
        frt_handle.write("\t".join(header) + "\n")

        for prot in seqs:
            ko_annot = util.gather_annotation_from_dict_for_homolo_group(
                prot, "ko", annotations
            )
            pgap_annot = util.gather_annotation_from_dict_for_homolo_group(
                prot, "pgap", annotations
            )
            pb_annot = util.gather_annotation_from_dict_for_homolo_group(
                prot, "paperblast", annotations
            )
            card_annot = util.gather_annotation_from_dict_for_homolo_group(
                prot, "card", annotations
            )
            isf_annot = util.gather_annotation_from_dict_for_homolo_group(
                prot, "tn_is", annotations
            )
            mibig_annot = util.gather_annotation_from_dict_for_homolo_group(
                prot, "mibig", annotations
            )
            vog_annot = util.gather_annotation_from_dict_for_homolo_group(
                prot, "vog", annotations
            )
            vfdb_annot = util.gather_annotation_from_dict_for_homolo_group(
                prot, "vfdb", annotations
            )
            pfam_annots = "NA"
            if "pfam" in annotations and prot in annotations["pfam"]: # type: ignore
                pfam_annots = "; ".join(annotations["pfam"][prot][0]) # type: ignore
            prot_seq = seqs[prot]
            row = [
                prot,
                ko_annot,
                pgap_annot,
                pb_annot,
                card_annot,
                isf_annot,
                mibig_annot,
                vog_annot,
                vfdb_annot,
                pfam_annots,
                prot_seq,
            ]
            row = [str(x) for x in row]
            frt_handle.write("\t".join(row) + "\n")
    

    # close logging object and exit
    log_object.info( # type: ignore
        f"******************\nannotateProtSeqs.py finished!\n******************\nConsolidated TSV can be found at: {annotation_tsv_file}"
    )
    util.close_logger_object(log_object)
    sys.exit(0)


if __name__ == "__main__":
    annotate_prot_seqs()
