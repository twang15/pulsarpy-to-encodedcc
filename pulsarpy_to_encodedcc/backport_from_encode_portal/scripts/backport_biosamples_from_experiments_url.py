#!/usr/bin/env python3

"""
Imports the biosamples from the experiments given by an ENCODE Portal search URL.
"""

import argparse

import pulsarpy.models.backport_from_encode_portal.backport as bp
import encode_utils.connection as euc

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-u", "--url", required=True, help="The search URL.")
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    url = args.url
    conn = euc.Connection("prod")
    results = conn.search(url=url)
    for i in results:
        dcc_exp = conn.get(i["@id"])
        for rep in dcc_exp["replicates"]:
            biosample_id = rep["library"]["biosample"]["@id"]
            bp.biosample(rec_id=biosample_id)
