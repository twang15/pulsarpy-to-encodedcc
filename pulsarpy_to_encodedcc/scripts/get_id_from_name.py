#!/usr/bin/env python
# -*- coding: utf-8 -*-

###Author
#Nathaniel Watson
#2018-06-21
#nathankw@stanford.edu
###

"""
Given a file with a column of names of records of a given type of model in Pulsar, fetches the IDs.
The record IDs are written to a new file.
"""

import argparse
import pulsarpy.models

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-m", "--model", required=True, help="The name of the Rails model class of the records for which we need to get to IDs.")
    parser.add_argument("-i", "--infile", required=True, help="The input file containing record names, one per row.")
    parser.add_argument("-o", "--outfile", required=True, help="The output file containing a column of record IDs.  Each row corresponds to the same row in the input file.")
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    model = args.model
    model = getattr(pulsarpy.models, model)
    infile = args.infile
    outfile = args.outfile
    fh = open(infile, 'r')
    fout = open(outfile, 'w')
    for line in fh:
        name = line.strip()
        if not name:
            fout.write("\n")
            continue
        #rec = model.find_by({"id": name})
        rec = model.find_by({"name": name})
        #rec = model.find_by({"upstream_identifier": name})
        #fout.write(f'{rec["upstream_identifier"]} {name}\n')
        fout.write(f'{rec["id"]}\n')

    fh.close()
    fout.close()

if __name__ == "__main__":
    main()
