#!/usr/bin/env python

"""
Given one or more DCC experiment IDs, looks at all read2s that were submitted and updates each r2 file
object such that it's paired_with property points to the correct r1. This works by looking at the aliases
in the r2 file object to see if there is one with _R2_001 in it. If so, it sets paired_with to be
the same alias, but with that segment replace with _R1_001. Thus, this script is nice if submissions
went wrong with regard to the file pairings, and this is one way to fix that. 
"""

import argparse
import encode_utils.connection as euc
import re

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--infile", required=True, help="""
      The input file with a DCC experiment on each line.""")
    return parser

def main():
    conn = euc.Connection("prod")
    reg = re.compile("_R2_001")
    parser = get_parser()
    args = parser.parse_args()
    ids = []
    fh = open(args.infile)
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        ids.append(line)
    for i in ids:
        h = conn.get_fastqfile_replicate_hash(exp_id=i)
        for bio_rep in h:
            for tech_rep in h[bio_rep]:
                read_files = h[bio_rep][tech_rep].get(2)
                # read_files is a list of file objects 
                if not read_files:
                    continue
                for r in read_files:
                    aliases = r["aliases"]
                    for a in aliases:
                        match = reg.search(a)
                        if match:
                            paired_with_name = a.replace(reg.pattern, "_R1_001")
                            payload = {conn.ENCID_KEY: a}
                            payload["paired_with"] = paired_with_name
                            try:
                                conn.patch(payload=payload)
                            except Exception:
                                break
                            break
                    
if __name__ == "__main__":
    main()                
