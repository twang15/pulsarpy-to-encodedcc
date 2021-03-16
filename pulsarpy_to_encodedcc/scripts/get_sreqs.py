import pulsarpy.models as models
import argparse
from encode_utils.parent_argparser import dcc_login_parser
import os, glob
import sys
from dxpy.scripts.dx import main as dx_main
from dxpy.scripts.dx import set_args_list
import hashlib

args_list = ['download', 'project_name', '-r', '-f']
FASTQ_FOLDER = '/oak/stanford/scg/prj_ENCODE/SRES/'

def extract_chip_meta(accession, fout):
    '''
    accession: string, ENCODE accession ID for a chipseq experiment.
    fout: output file handler
    '''
    chip = models.ChipseqExperiment.find_by({'upstream_identifier': accession})
    
    # wild-type replicate: <= 1
    wt = chip['wild_type_control_id']
    wt_biosample = models.Biosample.find_by({'id': wt})

    count = 0
    if wt_biosample:
        count = len(wt_biosample['library_ids'])
    
    if count > 1:
        raise Exception('Biosample {} has zero or more than one libraries.'.format(wt))
    elif count == 1: 
        library_id = wt_biosample['library_ids'][0]
        extract_library_meta(accession, library_id, 'wild', fout)
    else:
        print(f"chip {accession} does not have wild-type control.\n")
    
    # biological replicates: >= 1
    count = 0
    if chip['replicate_ids']:
        count = len(chip['replicate_ids'])

    if count >= 1:
        for lib in chip['replicate_ids']:
            extract_library_meta(accession, lib, 'biological', fout)
    else:
        raise Exception('ChipseqExperiment {} does not have biological replicates.'.format(chip['upstream_identifier']))
    
    # input control replicate: <= 1
    control = chip['control_replicate_ids']
    count = 0
    if control:
        count = len(control)

    if count > 1:
        raise Exception('ChipseqExperiment {} has zero or more than one libraries.'.format(chip['upstream_identifier']))
    elif count == 1: 
        library_id = control[0]
        extract_library_meta(accession, library_id, 'control', fout)
    else:
        print(f"chip {accession} does not have input library control.\n")


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

def extract_library_meta(accession, library_id, replicate_type, fout):
    '''
    library_id: int, library id in pulsar
    replicate_type: one in {'wild', 'biological', 'control'}

    return: library_meta, dict
    '''
    library_meta = {'id': library_id}
    library_meta['type'] = replicate_type
    
    # library
    library = models.Library.find_by({'id': library_id})
    
    # sequencing request
    for sreq_id in library['sequencing_request_ids']:
        sreq = models.SequencingRequest.find_by({'id': sreq_id})

        # sequencing run and DX project
        if len(sreq['sequencing_run_ids']) != 1:
            raise Exception('Sequencing request {} has zero or more than one sequencing requests.'.format(sreq['id']))
        srun_id = sreq['sequencing_run_ids'][0]
        srun = models.SequencingRun.find_by({'id': srun_id})
        fout.write(f'{accession}\t{library_id}\t{sreq_id}\t{srun["data_storage_id"]}\n')

def sreqs_meta(fout):
    empty = []
    for sreq_id in range(1, 302):
        sreq = models.SequencingRequest.find_by({'id': sreq_id})

        if sreq:
            # sequencing run and DX project
            if len(sreq['sequencing_run_ids']) == 0:
                #raise Exception('Sequencing request {} has zero or more than one sequencing requests.'.format(sreq['id']))
                empty.extend([sreq_id])

            for srun_id in sreq['sequencing_run_ids']:
                srun = models.SequencingRun.find_by({'id': srun_id})
                fout.write(f'{sreq_id}\t{sreq["name"]}\t{srun["data_storage_id"]}\t{srun["name"]}\n')
                fout.flush()

    print (empty)

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

    fout = open('enc', 'w')
    sreqs_meta(fout)
#    for record in records:
#        extract_chip_meta(record, fout)
#        #fout.write(f'{record} {upstream_id}\n')

    fout.close()

if __name__ == "__main__":
  main()
