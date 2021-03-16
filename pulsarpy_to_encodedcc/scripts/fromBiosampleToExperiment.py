import pulsarpy.models
import argparse

from encode_utils.parent_argparser import dcc_login_parser

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
  items = []
  for line in fh:
    line = line.strip()
    if not line or line.startswith("#"):
      continue

    items.append(line)

  fh.close()


  model = getattr(pulsarpy.models, "Biosample")
  fout = open("enc", "w")

  for item in items:
    rec = model.find_by({"id": item})

    if rec:
        fout.write(f"{item} {rec['crispr_modification_id']} {rec['part_of_id']} {rec['pcr_ids']} {rec['chipseq_experiment_ids']}\n")
    else:
      fout.write(f"{item}\n")

  fout.close()

if __name__ == "__main__":
  main()

