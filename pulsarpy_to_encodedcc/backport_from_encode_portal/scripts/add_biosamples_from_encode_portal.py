#!/usr/bin/env python3                                                                                 
# -*- coding: utf-8 -*-                                                                                
                                                                                                       
###                                                                                                    
# © 2018 The Board of Trustees of the Leland Stanford Junior University                                
# Nathaniel Watson                                                                                     
# nathankw@stanford.edu                                                                                
###

"""
Backports biosamples from the ENCODE Portal. Biosample identifiers on the ENCODE Portal can be 
provided in one of three ways:

    1. An input file with one ID per line. 
    2. On the command-line.
    3. Via an ENCODE Portal search URL.   
"""

import argparse
import pdb

import encode_utils.connection as euc
import pulsarpy_to_encodedcc.backport_from_encode_portal.backport as backport

ENC_CONN = euc.Connection("prod")
# Note that ex_url below has a double '%', where the second is used to escape the original. 
ex_url = "https://www.encodeproject.org/search/?type=Biosample&lab.title=Michael+Snyder%%2C+Stanford&award.rfa=ENCODE4&biosample_type=tissue"

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-u", "--url", help="""
    A search URL (wrapped in quotation marks to avoid shell expansion, to indicate which biosamples to import, i.e. {}.""".format(ex_url))
    group.add_argument("-i", "--infile", help="File containing biosample record IDs, one per line.")
    group.add_argument("-r", "--rec-ids", nargs="+", help="One or more biosample record identifiers separated by a space.")
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    url = args.url
    infile = args.infile
    rec_ids = args.rec_ids
    inputs = []
    if infile:
        with open(infile) as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                inputs.append(line)
    elif rec_ids:
        inputs = rec_ids
    else:
        res = ENC_CONN.search(url=url)
        for record in res: # i is a JSON object.
            # GET the object. Only part of the object is given in the search results. For example, 
            # the property 'genetic_modifications' isn't present. 
            inputs.append(record["@id"])

    for i in inputs:
        print("Backporting Biosample {}".format(i))
        backport.biosample(rec_id=i)
        
     
if __name__ == "__main__":
    main()
