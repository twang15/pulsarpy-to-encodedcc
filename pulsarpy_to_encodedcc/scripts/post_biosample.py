import pulsarpy_to_encodedcc.dcc_submit as dcc_submit
import traceback
import argparse

from encode_utils.parent_argparser import dcc_login_parser

submitter = dcc_submit.Submit(dcc_mode='prod')
#chipSeq = dcc_submit.Submit(dcc_mode='dev')

def log_traceback(ex):
    tb_lines = traceback.format_exception(ex.__class__, ex, ex.__traceback__)
    for tb in tb_lines:
        print(tb)

def get_parser():
    parser = argparse.ArgumentParser(
            description = __doc__,
            parents=[dcc_login_parser],
            formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-i", "--infile", required=True,
            help="""
      each row corresponds to one upstream_identifier of one chipSeq experiment.
      """)

    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    infile = args.infile
  
    fh = open(infile)
    records = []
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
  
        records.append(line)
  
    fh.close()

    fout = open('out', 'w')

    for record in records:
        # patch
        #upstream_id = submitter.post_biosample(rec_id=record, patch=True)

        # post
        upstream_id = submitter.post_biosample(rec_id=record, patch=False)
        fout.write(f'{record} {upstream_id}\n')
    
    fout.close()

if __name__ == "__main__":
  main()
