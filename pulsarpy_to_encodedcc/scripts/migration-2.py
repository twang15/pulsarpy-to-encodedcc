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

def set_download_project(project):
    args_list[1] = project + ":/"
    set_args_list(args_list)


def write_out(accession, fout, library, download=False):
    '''
    accession: string, ENCODE accession ID for a chipseq experiment.
    fout: output file handler
    library: dict, library meta data
    '''
    if download:
        fout.write(f"{accession}\t{library['type']}\t{library['id']}\t{library['barcode']}\
                 \t{library['sreq']}\t{library['dx_project']}\t{library['sres']}\t{library['read1_uri']}\
                 \t{library['read1_accession']}\t1\
                 \t{library['read1_file']}\t{library['read1_file_md5sum']}\n")
        fout.write(f"{accession}\t{library['type']}\t{library['id']}\t{library['barcode']}\
                 \t{library['sreq']}\t{library['dx_project']}\t{library['sres']}\t{library['read2_uri']}\
                 \t{library['read2_accession']}\t2\
                 \t{library['read2_file']}\t{library['read2_file_md5sum']}\n")
    else:
        fout.write(f"{accession}\t{library['type']}\t{library['id']}\t{library['barcode']}\
                 \t{library['sreq']}\t{library['dx_project']}\t{library['sres']}\t{library['read1_uri']}\
                 \t{library['read1_accession']}\t1\n")
        fout.write(f"{accession}\t{library['type']}\t{library['id']}\t{library['barcode']}\
                 \t{library['sreq']}\t{library['dx_project']}\t{library['sres']}\t{library['read2_uri']}\
                 \t{library['read2_accession']}\t2\n")
    
    fout.flush()

def download_fastq_file(library, read_num):
    '''
    Download fastq file from DNANexus.
    library: dict, library meta data
    read_num: string, read num for the current fastq file.
    '''

    project = library['dx_project']
    project_name = library['dx_project_name']
    path_to_sreq = FASTQ_FOLDER + str(library['sreq'])
    barcode = library['barcode']
    
    if os.path.isdir(path_to_sreq) and len(os.listdir(path_to_sreq)) != 0:
        print ("SReq dataset already downloaded")
        os.chdir(path_to_sreq)
    else:
        if not os.path.isdir(path_to_sreq):
            os.mkdir(path_to_sreq)
            #os.rmdir(path_to_sreq)
        
        os.chdir(path_to_sreq)
    
        # download the whole dataset for SREQ
        set_download_project(project_name)
        dx_main()

        # if download fails, clean the empty directory
        if len(os.listdir(path_to_sreq)) == 0:
            os.rmdir(path_to_sreq)

    # Check file existence
    # Assumption: there is only one fastq file for each read.
    files = glob.glob("**/*" + barcode + "*R" + read_num + "*.fastq.gz", recursive=True)

    # now look up the file's full path
    if len(files) == 1:
        return os.path.abspath(files[0])
    elif len(files) > 1:
        raise Exception(f'More than one fastq files exist for project {project} read {read_num}.')
        return None
    else:
        raise Exception(f'Fastq file for project {project} read {read_num} does not exist.')
        return None

       # # read 1
       # file_uri = library['read1_uri']
       # dx_file = dxpy.DXFile(dxid=file_uri, project=data_storage.project_identifier)
       # file_path = os.path.join(FASTQ_FOLDER, dx_file.name)

       # # read 2

def fastq_md5sum(fileName):
    '''
    fileName: string, the full path of the fastq file
    '''
    md5_hash = hashlib.md5()
    with open(fileName,"rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            md5_hash.update(byte_block)
        
        return md5_hash.hexdigest()

def download_and_compute_md5sum(library):
    read1_file = download_fastq_file(library, '1')
    read1_file_md5sum = fastq_md5sum(read1_file)

    read2_file = download_fastq_file(library, '2')
    read2_file_md5sum = fastq_md5sum(read2_file)

    return read1_file, read1_file_md5sum, read2_file, read2_file_md5sum

def download_and_compute_md5sum_tmp(library):
    return "read1_file", "read1_file_md5sum", "read2_file", "read2_file_md5sum"

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
    
    count = len(wt_biosample['library_ids'])
    if count > 1:
        raise Exception('Biosample {} has zero or more than one libraries.'.format(wt))
    elif count == 1:
        library_id = wt_biosample['library_ids'][0]
        extract_library_meta(accession, fout, library_id, 'wild')
    else:
        print(f"chip {accession} does not have wild-type control.\n")
    
    # biological replicates: >= 1
    count = 0
    if chip['replicate_ids']:
        count = len(chip['replicate_ids'])

    count = len(chip['replicate_ids'])
    if count >= 1:
        for library_id in chip['replicate_ids']:
            extract_library_meta(accession, fout, library_id, 'biological')
    else:
        raise Exception('ChipseqExperiment {} does not have biological replicates.'.format(chip['upstream_identifier']))
    
    # input control replicate: <= 1
    control = chip['control_replicate_ids']
    count = 0
    if control:
        count = len(control)

    count = len(control)
    if count > 1:
        raise Exception('ChipseqExperiment {} has zero or more than one libraries.'.format(chip['upstream_identifier']))
    elif count == 1: 
        library_id = control[0]
        extract_library_meta(accession, fout, library_id, 'control')
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

def extract_sreq_meta(library_id, sreq_id, replicate_type):
    '''
    library_id: int, library id in pulsar
    sreq_id: int, sequencing request id in pulsar
    replicate_type: one in {'wild', 'biological', 'control'}

    return: library_meta, dict
    '''
    library = models.Library.find_by({'id': library_id})
    library_meta = {'id': library_id}
    library_meta['type'] = replicate_type
    sreq = models.SequencingRequest.find_by({'id': sreq_id})
    
    # sequencing run and DX project
    if len(sreq['sequencing_run_ids']) != 1:
        raise Exception('Sequencing request {} has zero or more than one sequencing requests.'.format(sreq['id']))
    library_meta['sreq'] = sreq['id']

    srun = sreq['sequencing_run_ids'][0]
    srun = models.SequencingRun.find_by({'id': srun})
    ds = models.DataStorage.find_by({'id': srun['data_storage_id']})
    
    if not ds:
        raise Exception('Sequencing run {} has invalid DX project.'.format(srun['data_storage_id']))
    
    dx_project = ds['project_identifier']
    library_meta['dx_project'] = dx_project
    library_meta['dx_project_name'] = ds['name']
    
    # barcode
    barcode = library['barcode_id']
    barcode = models.Barcode.find_by({'id': barcode})
    barcode_sequence = barcode['sequence']
    library_meta['barcode'] = barcode_sequence
    
    # sequencing result
    if len(library['sequencing_result_ids']) == 0:
        raise Exception('Library {} has no sequencing result.'.format(library_id))
    elif len(library['sequencing_result_ids']) > 1:
        print('Library {} has more than one sequencing results.'.format(library_id))
        
        def find_sres(sres_list, library_id):
            '''
            sres_list: list, list of sres ids in pulsar of a srun
            library_id: int, library id in pulsar

            Return:
            sres_id: int, the sequencing result id in pulsar.
            '''
            for sres_id in sres_list:
                sres = models.SequencingResult.find_by({'id': sres_id})
                if sres['library_id'] == library_id:
                    return sres_id
        sres_id = find_sres(srun['sequencing_result_ids'], library_id)
    else:
        sres_id = library['sequencing_result_ids'][0]
    
    sres = models.SequencingResult.find_by({'id': sres_id})
    library_meta['sres'] = sres['id']

    # sequencing result ID
    sres_id = library['sequencing_result_ids'][0]
    
    # File DX ID
    # File ENCODE Accession
    
    # read 1
    read1_uri = sres['read1_uri']
    read1_accession = sres['read1_upstream_identifier']
    library_meta['read1_uri'] = read1_uri
    library_meta['read1_accession'] = read1_accession

    # read 2
    read2_uri = sres['read2_uri']
    read2_accession = sres['read2_upstream_identifier']
    library_meta['read2_uri'] = read2_uri
    library_meta['read2_accession'] = read2_accession
    
    f1, md5sum1, f2, md5sum2 = download_and_compute_md5sum(library_meta)
    library_meta['read1_file'] = f1
    library_meta['read1_file_md5sum'] = md5sum1
    library_meta['read2_file'] = f2
    library_meta['read2_file_md5sum'] = md5sum2

    return library_meta

def extract_library_meta(accession, fout, library_id, replicate_type):
    library = models.Library.find_by({'id': library_id})
    # sequencing request
    if len(library['sequencing_request_ids']) == 0:
        raise Exception('Library {} does not have sequencing requests.'.format(library_id))
    else:
        # one or more sequencing requests
        for sreq_id in library['sequencing_request_ids']:
            library_meta = extract_sreq_meta(library_id, sreq_id, replicate_type)
            write_out(accession, fout, library_meta, True)

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

    fout = open('enc', 'a')
    for record in records:
        extract_chip_meta(record, fout)
        #fout.write(f'{record} {upstream_id}\n')

    fout.close()

if __name__ == "__main__":
  main()
