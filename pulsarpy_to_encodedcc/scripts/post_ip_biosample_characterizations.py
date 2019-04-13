#!/usr/bin/env python3                                                                                 
# -*- coding: utf-8 -*-                                                                                
                                                                                                       
###Author                                                                                              
#Nathaniel Watson                                                                                      
#2019-04-12                                                                                            
#nathankw@stanford.edu                                                                                 
### 


"""
Submits Immunoblots in batch.
"""

import argparse

import pulsarpy_to_encodedcc.dcc_submit as d  

S = d.Submit("prod")

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-i", "--infile", required=True, help="""
      Tab-delimited input file where the first column contains Pulsar biosample IDs and the second contains
      Pulsar immunoblot_ids. Empty lines and lines starting with '#' are skipped.
    """)

    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    infile = args.infile
    
    items = []
    fh = open(infile)
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        items.append([int(x) for x in line.split("\t")])
     
    for i in items:
        b_id = i[0]
        ip_id = i[1]
        print("Biosample {}, Immunoblot {}".format(b_id, ip_id))
        S.post_ip_biosample_characterization(biosample_id=b_id, immunoblot_id=ip_id, patch=False)
        
if __name__ == "__main__":
   main() 
          
    
