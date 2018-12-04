#!/usr/bin/env python3

"""
Imports the biosamples from the experiments given by an ENCODE Portal search URL.
"""

import pdb
import json
import time
import argparse

import encode_utils.connection as euc
import pulsarpy.models as models 

SPECIES = "Homo sapiens"

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-i", "--infile", help="JSON formatted file containing the records.")
    group.add_argument("-u", "--url", help="The search URL.")
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    url = args.url
    infile = args.infile
    if url:
        conn = euc.Connection("prod")
        results = conn.search(url=url)
    else:
        fh = open(infile)
        results = json.load(fh)
    admin = models.User.find_by({"first_name": "Admin"})
    if not admin:
        raise Exception("Could not find the Admin user in the database, which is needed for associating with new records.")
    created = 0
    patched = 0
    total = 0
    for rec in results:
        patch = False
        total += 1
        organism = rec["organism"]["scientific_name"]
        if organism != SPECIES:
          continue
        payload = {}
        label = rec["label"]
        payload["name"] = label
        # Check if the target already exists in the database.
        pulsar_record = models.Target.find_by({"name": label})
        upstream = rec["@id"].strip("/").split("/")[-1]
        if pulsar_record and upstream != pulsar_record["upstream_identifier"]:
            patch = True
        elif pulsar_record:
            continue # Can add support for patch operation later. 
        payload["upstream_identifier"] = upstream
        payload["user_id"] = admin["id"]
        xrefs = rec["dbxref"]
        for ref in xrefs:
            tokens = ref.split(":")
            prefix, ref = ref.rsplit(":", 1)
            if prefix == "ENSEMBL":
                payload["ensembl"] = ref
            elif prefix == "UniProtKB":
                payload["uniprotkb"] = ref
            elif prefix == "RefSeq":
                payload["refseq"] = ref
        print("Creating {}".format(payload))
        if patch:
            target = models.Target(pulsar_record["id"])
            target.patch(payload)
            patched += 1
            print("Patched: {}".format(patched))
        else:
            models.Target.post(payload)
            created += 1
            print("Created: {}".format(created))
        print("Total processed: {}".format(total))

if __name__ == "__main__":
    main()
