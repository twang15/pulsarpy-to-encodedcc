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


  model = getattr(pulsarpy.models, "ChipseqExperiment")
  #model = getattr(pulsarpy.models, "Atacseq")
  fout = open("bid", "w")

  for item in items:
    rec = model.find_by({"upstream_identifier": item})
    #rec = model.find_by({"id": item})

    for lid in rec['replicate_ids']:
      library_model = getattr(pulsarpy.models, "Library")
      library_rec = library_model.find_by({"id": lid})
      bid = library_rec["biosample_id"]
      biosample = pulsarpy.models.Biosample(bid)
      child_pcr = biosample.pcr_ids

      parent = biosample.part_of_id
      p = pulsarpy.models.Biosample(parent)
      parent_pcr = p.pcr_ids
      #bid = biosample.upstream_identifier

      #fout.write(f"{library_rec['upstream_identifier']} {bid} {rec['upstream_identifier']}\n")
      fout.write(f"{item} {bid} {child_pcr} {parent} {parent_pcr}\n")

  fout.close()

if __name__ == "__main__":
  main()

