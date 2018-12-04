#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###Author
#Nathaniel Watson
#2018-10-22
#nathankw@stanford.edu
###

"""
Submits records from Pulsar to the ENCODE Portal. 
"""

import argparse
import inflection

from encode_utils.parent_argparser import dcc_login_parser    
import pulsarpy.dcc_submit as dcc_submit
from pulsarpy import models

MODEL_PROFILE = {}
MODEL_PROFILE["biosample"] = "Biosample"
MODEL_PROFILE["library"] = "Library"

def get_parser():
    parser = argparse.ArgumentParser(
        description=__doc__,
        parents=[dcc_login_parser],
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-p", "--profile-id", required=True, help="""                                  
        The ID of the ENCODE profile to submit to, i.e. use 'genetic_modification' for                            
        https://www.encodeproject.org/profiles/genetic_modification.json. Must be a value from the set
        {}""".format(list(MODEL_PROFILE.keys())))
    parser.add_argument("-i", "--infile", required=True, help="""
        The input file containing Pulsar record identifiers, one per row. The record identifiers should belong to
        a Pulsar model that is the equivalent of the DCC profile specified by --profile-id.  For example,
        if submitting Pulsar records to the ENCODE 'biosample' profile, the equivalent model in Pulsar is
        Biosample, thus your record identifiers in this file should thus be Biosample identifiers. Note that
        the record identifier to use must be either the record ID or record name if POSTING, or can additionally
        be the value of the Pulsar record's upstream_identifier attribute if PATCHING.""")
    parser.add_argument("--no-extend-arrays", action="store_true", help="""
        Only affects updating objects on the ENCODE Portal. By default, when updating an array
        attribute, the array will be extended with the provided values from the input file. However,
        including this command-line option means to first empty the array contents.""")
    parser.add_argument("--patch", action="store_true", help="""
        Presence of this option indicates to PATCH an existing DCC record rather than register a 
        new one.""")
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    no_extend_arrays = args.no_extend_arrays
    patch = args.patch
    dcc_mode = args.dcc_mode
    dcc_profile = args.profile_id
    infile = args.infile

    submit = dcc_submit.Submit(dcc_mode=dcc_mode, extend_arrays=not no_extend_arrays)
    post_method_name = "post_" + inflection.underscore(MODEL_PROFILE[dcc_profile])
    post_method = getattr(submit, post_method_name)

    fh = open(infile, 'r')
    rec_ids = []
    for line in fh:
        line = line.strip()
        if not line:
            continue
        rec_ids.append(line)
    fh.close()

    for i in rec_ids:
        post_method(rec_id=i, patch=patch)

if __name__ == "__main__":
    main()
