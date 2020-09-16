import pulsarpy.models as models
import argparse

import pulsarpy_to_encodedcc.dcc_submit as d  

S = d.Submit("prod")

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
def post_biosample_ips(items, fout, forChildren=False):
    model = getattr(models, 'Biosample')
    child = 0

    for uid in items:
        if uid.startswith("ENC"):
          rec = model.find_by({"upstream_identifier": uid})
        else:
          rec = model.find_by({"id": uid})
        
        if not rec['chipseq_experiment_ids'] and rec['biosample_part_ids']:
          #items.extend([str(id) for id in rec['biosample_part_ids']])
          children = [str(id) for id in rec['biosample_part_ids']]
          success = post_biosample_ips(children, fout, True)

          if success:
            fout.write(uid + " Success\n")
            fout.flush()
            continue
          else:
            fout.write(uid + " Failed\n")

        b_id = rec['id']
        ip_id = rec['immunoblot_ids']
        #if not ip_id:
          #print("UPSTREAM: {}".format(uid))
          #fout.write(uid+"\n")
          #print("Biosample {}, Immunoblot {}".format(b_id, ip_id[0]))

        upstream_id = ""
        try:
          if b_id:
            if ip_id:
              print("Biosample {}, Immunoblot {}".format(b_id, ip_id[0]))
              upstream_id = S.post_ip_biosample_characterization(biosample_id=b_id, immunoblot_id=ip_id[0], patch=False)
            else:
              upstream_id = S.post_ip_biosample_characterization(biosample_id=b_id, immunoblot_id=None, patch=False)
        except ValueError as ve:
          if str(ve).startswith("Crispr") and not forChildren:
            fout.write(uid+" No GM\n")
            continue
          elif str(ve).startswith("IP") and not forChildren:
            fout.write(uid+" No IPs\n")
            continue
          elif str(ve).startswith("GelLane") and not forChildren:
            fout.write(f"{uid} GelLane didn't pass\n")
            continue
          elif str(ve).startswith("Replicate") and not forChildren:
            print (str(ve))
            fout.write(f"{uid} Replicate not used\n")
            continue
          else:
            raise ValueError("Unknown problem")
        except:
            raise ValueError("Unknown problem")
        
        if not upstream_id:
          fout.write(uid+"\n")
          fout.flush()
          raise ValueError("Unknown problem")
        elif forChildren:
          print ("post for Child biosamples {}".format(uid))
          child += 1
        else:
          fout.write(uid + " Success\n")
          fout.flush()

    if forChildren:
      if child > 0:
        return True
      else:
        return False
    else:
      return True

def main():
    parser = get_parser()
    args = parser.parse_args()
    infile = args.infile

    items = []
    fh = open(infile)
    fout = open("uid", "w")
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        items.append(line)

    post_biosample_ips(items, fout, False)
    fout.close()

if __name__ == "__main__":
   main() 
          
