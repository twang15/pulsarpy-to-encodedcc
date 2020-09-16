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

  #fh = open(infile)
  #items = []
  #for line in fh:
  #  line = line.strip()
  #  if not line or line.startswith("#"):
  #    continue

  #  items.append(line)

  #fh.close()


  model = getattr(pulsarpy.models, "ChipseqExperiment")
  #model = getattr(pulsarpy.models, "Atacseq")
  fout = open("bid", "w")
  subsetA = open(infile)
  As = []
  for line in subsetA:
      line = line.strip()
      As.append(line)

  subsetA.close()

  for chip in range(5, 362):
      rec = model.find_by({"id": chip})
      if rec and rec['upstream_identifier'] and rec['upstream_identifier'] not in As:
          fout.write(f"{rec['upstream_identifier']}\n")
  
  fout.close()

if __name__ == "__main__":
  main()

