import argparse
import os
import logging
import pulsarpy.models
import requests
import pulsarpy_to_encodedcc.dcc_submit as dcc_submit

submitter = dcc_submit.Submit(dcc_mode='prod')
#submitter = dcc_submit.Submit(dcc_mode='dev')

def get_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--no-extend-arrays", action="store_true", help="""
        Only affects updating objects on the ENCODE Portal. By default, when updating an array
        attribute, the array will be extended with the provided values from the input file. However,
        including this command-line option means to first empty the array contents.""")
    parser.add_argument("--patch", action="store_true", help="""
        Presence of this option indicates to PATCH an existing DCC record rather than register a 
        new one.""")
    
    parser.add_argument("--biosample", action="store_true", help="""
            Presence of this option indicates the input file is a list of biosample IDs.""")

    parser.add_argument("-i", "--infile", required=True, help="The input file containing record names, one per row.")
    parser.add_argument("-o", "--outfile", required=True, help="The output file containing a column of record IDs.  Each row corresponds to the same row in the input file.")
    return parser

def get_gel_lane_with_biosmaple_pcr(biosample_id, gel_id):

    '''
    Given a gel and a biosample, return the gel_lane for the biosample.
    If there are multiple gel lanes with the same biosample, then return
    all of them.

    Args:
      biosample_id: int, biosample id
      gel_id: int, gel id for its pcr

    Returns:
      `None` if the GelLane didn't pass. Otherwise, all related 
      `pulsarpy.models.GelLane` instances.

    Raises:
      `IpLaneException`: One of multiple issues that could be present as indicated by the error
       message, i.e., 

          * The Biosample doesn't have an associated Gel
          * There isn't a GelLane with the Biosample on it
    '''

    gel_lanes = []
    gel = pulsarpy.models.Gel(gel_id)
    for gel_lane_id in gel.gel_lane_ids:
        gel_lane = pulsarpy.models.GelLane(gel_lane_id)

        # could a gel_lane for pcr only be associated with a parent biosample?
        if biosample_id == gel_lane.biosample_id:
            gel_lanes.append(gel_lane)

    if not gel_lanes:
        biosample = pulsarpy.models.Biosample(biosample_id)
        children = biosample.biosample_part_ids
        for child in children:
            gel_lanes = get_gel_lane_with_biosmaple_pcr(child, gel_id)
            if gel_lanes:
                return gel_lanes

    return gel_lanes


def main():
    parser = get_parser()
    args = parser.parse_args()
    no_extend_arrays = args.no_extend_arrays
    patch = args.patch
    model = 'Pcr'
    pcr_model = getattr(pulsarpy.models, model)
    biosample_model = getattr(pulsarpy.models, "Biosample")

    # Genetic_modification that Pcr characterizes.
    # Crispr_model = getattr(pulsarpy.models, "CrisprModification")
    Crispr_model = pulsarpy.models.CrisprModification

    infile = args.infile
    outfile = args.outfile
    fh = open(infile, 'r')
    fout = open(outfile, 'w')
    for line in fh:
        if not line or line.startswith("#"):
            continue

        pcr_id = line.strip()
        if not pcr_id:
            fout.write("\n")
            continue
        pcr = pcr_model.find_by({"id": pcr_id})

        # A biosample has at most 1 pcr and a pcr can only be associated with one biosample.
        # Although a Pcr may characterize multiple biosamples of relations, such as parent/child, or siblings.
        biosample = biosample_model.find_by({'id': pcr['biosample_id']})
        crispr_id = biosample['crispr_modification_id']

        # required field https://www.encodeproject.org/profiles/genetic_modification_characterization
        crispr = Crispr_model(crispr_id)
        crispr_construct_id = crispr['crispr_construct_ids'][0]
        crispr_construct = pulsarpy.models.CrisprConstruct(crispr_construct_id)
        target_TF = pulsarpy.models.Target(crispr_construct['target_id'])['name']

        # cell line name
        btn = pulsarpy.models.BiosampleTermName(biosample['biosample_term_name_id']).name
        
        # use a brief caption
        caption = "genetic modification characterization for {} in cell line {}.".format(target_TF, btn)
        payload = {'caption': caption}
        
        # gel_lane
        gel = pulsarpy.models.Gel(pcr['gel_id'])
        gel_lanes = get_gel_lane_with_biosmaple_pcr(biosample['id'], gel.id)
        if not gel_lanes:
            raise IpLaneException("Could't find a GelLane that has Biosample {} on gel {}.".format(biosample_id, gel_id))

        # Submit payload
        if patch:  
            upstream_id = submitter.patch(payload, gel_lanes[0].upstream_identifier)
        else:
            upstream_id = submitter.post(payload=payload, dcc_profile="genetic_modification_characterization", pulsar_model=pulsarpy.models.GelLane, pulsar_rec_id=gel_lanes[0].id, repost=False)

            for gel_lane in gel_lanes[1:]:
                gel_lane.patch(payload={"upstream_identifier": upstream_id})
            

        fout.write(f"{pcr['id']} {upstream_id}\n")
    fh.close()
    fout.close()

if __name__ == "__main__":
    main()
