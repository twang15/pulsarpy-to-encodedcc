#!/usr/bin/env python3

"""
Fetches the given experiment from the Portal in JSON format, and serializes it to tab-delimited
format.
"""

import argparse
import os
import sys

import encode_utils.connection as euc
from encode_utils.parent_argparser import dcc_login_parser

# Check that Python3 is being used
v = sys.version_info
if v < (3, 3):
    raise Exception("Requires Python 3.3 or greater.")

EXP_TAB = "experiments.txt"
EXP_HEADER = [
  "name",
  "#accession",
  "upstream_identifier",
  "description",
  "#target_name",
  "document_ids",
  "submitter_comments",
  "notes",
]

FILE_TAB = "files.txt"
FILE_HEADER = [
  "name",
  "#accession",
  "#alias",
  "#platform",
  "#submitted_file_name",
  "#rep_uid",
  "#rep_accession",
  "#run_type",
  "#controlled_by",
  "#paired_with",
  "#paired_end",
  "#read_count",
  "#read_length",
  "#barcode",
  "#machine",
  "#lane",
]

REP_TAB = "replicates.txt"
REP_HEADER = [
  "name",
  "upstream_identifier",
  "chipseq_experiment_id",
  "biosample_id",
  "#biosample_alias",
  "biological_replicate_number",
  "technical_replicate_number",
  "antibody_id",
  "#antibody_accession",
  "submitter_comments",
  "notes",
]
BIO_TAB = "biosamples.txt"
BIO_HEADER = [
  "name",
  "#accession",
  "upstream_identifier",
  "part_of_id",
  "nih_institutional_certification",
  "pooled_from_biosample_ids",
  "treatment_ids",
  "document_ids",
  "biosample_type_id",
  "biosample_term_name_id",
  "vendor_id",
  "vendor_product_identifier",
  "lot_identifier",
  "donor_id",
  "passage_number",
  "date_biosample_taken",
  "submitter_comments",
  "notes",
]

LIB_TAB = "libraries.txt"
LIB_HEADER = [
  "name",
  "#accession",
  "upstream_identifier",
  "biosample_id",
  "nucleic_acid_term_id",
  "strand_specific",
  "document_ids",
  "size_range",
  "treatments",
  "vendor_id",
  "vendor_product_identifir",
  "lot_identifier",
  "library_fragmentation_method_id",
  "sequencing_library_prep_kit_id",
  "paired_end",
  "barcode_id",
  "paired_barcode_id",
  "submitter_comments",
  "notes",
]
GM_TAB = "gms.txt"
GM_HEADER = [
  "name",
  "accession",
  "upstream_identifier",
  "biosample_id ",
  "description",
  "document_ids",
  "category",
  "purpose",
  "#method",
  "#guide_rna_sequences",
  "#introduced_tags",
  "#reagents",
  "#characterizations",
  "crispr_construct_ids",
  "donor_construct_id",
  "notes",
]

CONN = "" # connection object to ENCODE Portal

def file_is_empty(name):
    return not os.stat(name).st_size

def portal_ids_to_aliases(ids):
    """
    Given a list of identifiers from the Portal, gets the first alias of each record specified by
    the identifier and returns the result in a list. If a particular record doesn't have any
    aliases, then the original identifier provided is used in place.
    """
    global  CONN
    res = []
    for i in ids:
        rec = CONN.get(i)
        aliases = rec.get("aliases", [])
        if not aliases:
            res.append(i)
        else:
            res.append(aliases[0])
    return res
        


def get_parser():
    parser = argparse.ArgumentParser(
        description = __doc__,
        parents=[dcc_login_parser],
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-e", "--exp", required=True, help="An identifier for an experiment record.")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory.")
    return parser

def main():
    global CONN
    parser = get_parser()
    args = parser.parse_args()
    exp_id = args.exp
    outdir = args.outdir
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    dcc_mode = args.dcc_mode

    if dcc_mode:
        CONN = euc.Connection(dcc_mode)
    else:
        # Default dcc_mode taken from environment variable DCC_MODE.
        CONN = euc.Connection()

    exp = CONN.get(exp_id)

    # Open output file handles
    # experiments file
    exp_file = os.path.join(outdir, EXP_TAB)
    expfh = open(exp_file, "a")
    if file_is_empty(exp_file):
        expfh.write("\t".join(EXP_HEADER) + "\n")
    # files file
    file_file = os.path.join(outdir, FILE_TAB)
    ffh = open(file_file, "a")
    if file_is_empty(file_file):
        ffh.write("\t".join(FILE_HEADER) + "\n")
    # replicates file
    rep_file = os.path.join(outdir, REP_TAB)
    repfh = open(rep_file, "a")
    if file_is_empty(rep_file):
        repfh.write("\t".join(REP_HEADER) + "\n")
    # genetic modifications file
    gm_file = os.path.join(outdir, GM_TAB)
    gmfh = open(gm_file, "a")
    if file_is_empty(gm_file):
        gmfh.write("\t".join(GM_HEADER) + "\n")
    # biosamples file
    bio_file = os.path.join(outdir, BIO_TAB)
    biofh = open(bio_file, "a")
    if file_is_empty(bio_file):
        biofh.write("\t".join(BIO_HEADER) + "\n")
    # libraries file
    lib_file = os.path.join(outdir, LIB_TAB)
    libfh = open(lib_file, "a")
    if file_is_empty(lib_file):
        libfh.write("\t".join(LIB_HEADER) + "\n")

    expfh.write("\t") # emtpy for name field in Pulsar
    expfh.write(exp["accession"] + "\t")
    exp_alias = exp["aliases"][0]
    expfh.write(exp_alias + "\t")
    expfh.write(exp["description"] + "\t")
    expfh.write(exp["target"]["name"] + "\t")
    document_aliases = portal_ids_to_aliases(exp["documents"])
    expfh.write(",".join(document_aliases) + "\t")
    submitter_comments = exp.get("submitter_comment", "")
    expfh.write(submitter_comments + "\t")
    expfh.write("\t") # empty for notes field in Pulsar
    expfh.write("\n")
    # START FILE FILE
    fastq_files = CONN.get_fastqfiles_on_exp(exp_id)
    for f in fastq_files:
        ffh.write("\t") # empty for name field in Pulsar
        ffh.write(f["accession"] + "\t")
        ffh.write(f["aliases"][0] + "\t")
        platform = f["platform"]["aliases"][-1]
        ffh.write(platform + "\t")
        ffh.write(f.get("submitted_file_name", "") + "\t")
        rep = f["replicate"]
        ffh.write(rep.get("uuid", "") + "\t")
        ffh.write(rep["aliases"][0] + "\t")
        ffh.write(f.get("run_type", "") + "\t")
        controlled_by = f.get("controlled_by", [])
        ffh.write(",".join(controlled_by) + "\t")
        ffh.write(f.get("paired_with", "") + "\t")
        ffh.write(f.get("paired_end", "") + "\t")
        ffh.write(str(f.get("read_count", "")) + "\t")
        ffh.write(str(f.get("read_length", "")) + "\t")
        fc = f.get("flowcell_details", {})
        if fc:
            fc = fc[0]
        ffh.write(fc.get("barcode", "") + "\t")
        ffh.write(fc.get("machine", "") + "\t")
        ffh.write(str(fc.get("lane", "")) + "\t")
        ffh.write("\n")

    # START REPLICATE FILE
    reps = exp["replicates"]
    for i in reps:
        repfh.write("\t") # empty for name field in Pulsar
        repfh.write(i["aliases"][0] + "\t")
        repfh.write(exp_alias + "\t")
        lib = i["library"]
        bio = lib["biosample"]
        repfh.write("\t") # empty for biosample_id fkey field in Pulsar
        biosample_alias = bio["aliases"][0]
        repfh.write(biosample_alias + "\t")
        repfh.write(str(i["biological_replicate_number"]) + "\t")
        repfh.write(str(i["technical_replicate_number"]) + "\t")
        repfh.write("\t") # empty for antibody_id fkey field in Pulsar
        antibody = i.get("antibody", "")
        if antibody:
            repfh.write(antibody["accession"] + "\t")
        else:
            repfh.write("\t")
        repfh.write(i.get("submitter_comment", "") + "\t")
        repfh.write("\t") # empty for notes field in Pulsar
        repfh.write("\n")
        # START BIOSAMPLE FILE
        biofh.write("\t") # empty for name field in Pulsar
        biofh.write(bio["accession"] + "\t") 
        biosample_upstream_id = bio["aliases"][0]
        biofh.write(biosample_upstream_id + "\t")
        biofh.write(bio.get("part_of", "") + "\t")
        biofh.write(bio.get("nih_institutional_certification", "") + "\t")
        pooled_from = bio.get("pooled_from", [])
        biofh.write(",".join(pooled_from) + "\t")
        treatment_dicts = bio.get("treatments", {})
        treatment_uuids = [x["uuid"] for x in treatment_dicts]
        treatment_aliases = portal_ids_to_aliases(treatment_uuids)
        biofh.write(",".join(treatment_aliases) + "\t")
        document_aliases = portal_ids_to_aliases(bio.get("documents", []))
        biofh.write(",".join(document_aliases) + "\t")
        biofh.write(bio["biosample_type"] + "\t")
        biofh.write(bio["biosample_term_name"] + "\t")
        biofh.write(bio["source"]["name"] + "\t")
        biofh.write(bio.get("product_id", "") + "\t")
        biofh.write(bio.get("lot_id", "") + "\t")
        biofh.write(bio["donor"]["aliases"][0] + "\t")
        biofh.write(bio.get("passage_number", "") + "\t")
        date_taken = bio.get("culture_start_date", "")
        if not date_taken:
            date_taken = bio.get("date_obtained", "")
        biofh.write(date_taken + "\t")
        biofh.write(bio.get("submitter_comment", "") + "\t")
        biofh.write("\t") # empty for notes field in Pulsar
        biofh.write("\n")
        # update gm file
        for gm_id in bio.get("genetic_modifications", []):
            gm = CONN.get(gm_id)
            gmfh.write("\t") # empty for name field in Pulsar
            gmfh.write(gm["accession"] + "\t")
            gmfh.write(gm["aliases"][0] + "\t")
            gmfh.write(biosample_upstream_id + "\t")
            gmfh.write(gm.get("description", "") + "\t")
            document_aliases = portal_ids_to_aliases(gm.get("documents", []))
            gmfh.write(",".join(document_aliases) + "\t")
            gmfh.write(gm.get("category", "") + "\t")
            gmfh.write(gm.get("purpose", "") + "\t")
            gmfh.write(gm.get("method", "") + "\t")
            guide_seqs = gm.get("guide_rna_sequences", [])
            gmfh.write(",".join(guide_seqs) + "\t")
            tags = gm.get("introduced_tags", [])
            gmfh.write(str(tags) + "\t")
            reagents = gm.get("reagents")
            gmfh.write(str(reagents) + "\t")
            chars = gm.get("characterizations", [])
            gmfh.write(",".join(chars) + "\t")
            gmfh.write("\t") # empty for crispr_construct_ids in Pulsar
            gmfh.write("\t") # empty for donor_construct_id in Pulsar
            gmfh.write("\t") # empty for notes field in pulsar
            gmfh.write("\n")
        # START LIBRARY FILE
        libfh.write("\t") # empty for name field in Pulsar
        libfh.write(lib["accession"] + "\t")
        libfh.write(lib["aliases"][0] + "\t")
        libfh.write(biosample_upstream_id + "\t")
        libfh.write(lib["nucleic_acid_term_name"] + "\t")
        strand_specific = str(lib.get("strand_specificity", False))
        libfh.write(strand_specific + "\t")
        document_aliases = portal_ids_to_aliases(lib["documents"])
        libfh.write(",".join(document_aliases) + "\t")
        libfh.write(lib.get("size_range", "") + "\t")
        treatment_aliases = portal_ids_to_aliases(lib["treatments"])
        libfh.write(",".join(treatment_aliases) + "\t")
        libfh.write(lib.get("source", "") + "\t")
        libfh.write(lib.get("product_id", "") + "\t")
        libfh.write(lib.get("lot_id", "") + "\t")
        libfh.write(lib.get("fragmentation_method", "") + "\t")
        libfh.write("\t") # empty for sequencing_library_prep_kit_id field in Pulsar
        libfh.write("\t") # empty for paired_end field in Pulsar
        libfh.write("\t") # empty for barcode_id in Pulsar
        libfh.write("\t") # empty for paired_barcode_id in Pulsar
  
        libfh.write(lib.get("submitter_comment", "") + "\t")
        libfh.write("\t") # empty for notes field in pulsar
        libfh.write("\n")

    expfh.close()
    ffh.close()
    repfh.close()
    biofh.close()
    libfh.close()
    gmfh.close()

if __name__ == "__main__":
    main()
