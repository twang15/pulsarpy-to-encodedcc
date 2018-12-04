#!/usr/bin/env python3                                                                                 
# -*- coding: utf-8 -*-                                                                                
                                                                                                       
###                                                                                                    
# © 2018 The Board of Trustees of the Leland Stanford Junior University                                
# Nathaniel Watson                                                                                     
# nathankw@stanford.edu                                                                                
###

"""
Backports biosamples from the ENCODE Portal. A search URL for the ENCODE Portal must be given
in order to designate which biosamples to backport.  
"""

import argparse
import pdb

import pulsarpy.models
from encode_utils.connection import Connection
import encode_utils.profiles as eup

# Note that ex_url below has a double '%', where the second is used to escape the original. 
ex_url = "https://www.encodeproject.org/search/?type=Biosample&lab.title=Michael+Snyder%%2C+Stanford&award.rfa=ENCODE4&biosample_type=tissue"

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-u", "--url", required=True, help="""
    A search URL (wrapped in quotation marks to avoid shell expansion, to indicate which biosamples to import, i.e. {}.""".format(ex_url))
    return parser


def create_payload(rec_json):
    payload = {}
    payload["aliases"] = rec_json["aliases"]
    payload["description"] = rec_json["summary"]
    payload["submitter_comments"] = rec_json["submitter_comment"]

    payload["date_biosample_taken"] = rec_json["date_obtained"]
    payload["documents"] = rec_json["documents"]
    payload["donor"] = rec_json["donor"] 
    #payload["genetic_modifications"] = rec_json[
    payload["treatments"] = rec_json["documents"]

def main():
    biosample_profile = eup.Profile("biosample")
    parser = get_parser()
    args = parser.parse_args()
    url = args.url
    conn = Connection("prod")
    res = conn.search(url=url)
    for record in res: # i is a JSON object.
        # GET the object. Only part of the object is given in the search results. For example, 
        # the property 'genetic_modifications' isn't present. 
        rec_id = record["@id"]
        json_rec = conn.get(rec_id)
        json_rec = biosample_profile.filter_non_writable_props(json_rec, keep_identifying=True) 
        pdb.set_trace()
        payload = create_payload(json_rec)
        
     
if __name__ == "__main__":
    main()
