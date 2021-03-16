import pulsarpy.models as models
import argparse

import pulsarpy_to_encodedcc.dcc_submit as d  

#S = d.Submit("prod")
S = d.Submit("dev")

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__)
#    parser.add_argument("-h", "--help", required=False, help="""
#      A single column input file contains Pulsar biosample upstream_IDs.
#      """)

    parser.add_argument("-i", "--infile", required=True, help="""
      A single column input file contains Pulsar biosample upstream_IDs.
      Empty lines and lines starting with '#' are skipped.
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
        items.append(line)

    model = getattr(models, "Atacseq")
    for item in items:
      #rec = model.find_by({"upstream_identifier": item})
      #rec_id = rec['id']
      rec_id = item
      
      # Post
      S.post_scAtacseq_exp(rec_id=rec_id, patch=False)
      
      # Patch
      #S.post_scAtacseq_exp(rec_id=rec_id, patch=False, patch_all=True)

    fh.close()

if __name__ == "__main__":
   main() 
