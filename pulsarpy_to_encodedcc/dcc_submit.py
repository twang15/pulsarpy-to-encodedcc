# -*- coding: utf-8 -*-

###
# © 2018 The Board of Trustees of the Leland Stanford Junior University
# Nathaniel Watson
# nathankw@stanford.edu
###

"""
Required environment variables
  1) Those that are required in the pulsarpy.models module for connecting to Pulsar LIMS:
     -PULSAR_API_URL
     -PULSAR_TOKEN
  2) Those that are required in the encode_utils.connection module to connect to the ENCODE Portal:
     -DCC_API_KEY
     -DCC_SECRET_KEY

Optional environment variables:
  1) DCC_MODE_ - Specifies which ENCODE Portal host to connect to. If not set, then must be provided
     when instantiating the Submit() class.

..  _DCC_MODE: https://encode-utils.readthedocs.io/en/latest/connection.html#encode_utils.connection.Connection.dcc_mode

"""

import base64
import logging
import os, glob
import re
import requests
import sys

import dxpy

import pulsarpy_to_encodedcc
from pulsarpy_to_encodedcc import FASTQ_FOLDER, log_error
from pulsarpy import models
import pulsarpy.utils
import encode_utils as eu
import encode_utils.aws_storage
import encode_utils.replicate
import encode_utils.connection as euc
import encode_utils.utils as euu
import pdb

error_logger = logging.getLogger(pulsarpy_to_encodedcc.ERROR_LOGGER_NAME)
# regex for finding one or more continuous spaces
space_reg = re.compile(r' +')


class IpLaneException(Exception):
    """
    Raised when posting an IP. The is a temporary class that'll be removed once the exceptions
    are handled property in the Submit.post_ip_lane method. 
    """

class ExpMissingReplicates(Exception):
    """
    Raised when trying to POST an experiment to the Portal (such as a control experiment) and
    there aren't any replicates (Biosample records) to attach to it.
    """

class MissingSequencingResult(Exception):
    pass

class MissingTargetUpstream(Exception):
    """
    Raised when submitting a record that tries to link to a DCC target, but the target record in Pulsar 
    doesn't have the upstream_identifier attribute set.
    """
    

class UpstreamNotSet(Exception):
    pass


class NoFastqFile(Exception):
    """
    Raised in Submit.post_fastq_file() when submitting either a R1 FASTQ file or a R2 FASTQ file,
    and the filepath isn't set in the corresponding SequencingResult record in Pulsar.
    """
    pass


#def dec(model_class):
#    """
#    A decorator that is to be used with the post_* methods defined in the Submit class defined below.
#    """
#
#    def wrapper(func):
#
#        def inner(self, rec_id, patch=False, *args, **kwargs):
#            """
#            Saves time by checking whether a record needs to be posted to the Portal before it's 
#            payload is constructed.  It need not be posted if in Pulsar it has the upstream_identifier
#            attribute set AND the 'patch' argument is False. The decorated method will thus only run
#            if upstream_identifier isn't set, or if the 'patch' argument is True. 
#            """
#            rec = model_class(rec_id) 
#            upstream = rec.get_upstream() 
#            if upstream and not patch:
#                # Then no need to post
#                return upstream
#            else:
#                self.func(rec_id=rec_id, patch=patch, *args, **kwargs)
#        return inner
#
#    return wrapper


class Submit():
    """
    Contains methods for submitting various types of objects in Pulsar to the ENCODE Portal. 
    """

    def __init__(self, dcc_mode=None, extend_arrays=True):
        if not dcc_mode:
            try:                                                                                    
                dcc_mode = os.environ["DCC_MODE"]                                                   
                print("Utilizing DCC_MODE environment variable.")                 
            except KeyError:                                                                        
                print("ERROR: You must supply the `dcc_mode` argument or set the environment variable DCC_MODE.")
                sys.exit(-1)                                                                        
        self.dcc_mode = dcc_mode
        self.ENC_CONN = euc.Connection(self.dcc_mode, submission=True)
        #: When patching, there is the option to extend array properties or overwrite their values.
        #: The default is to extend.
        self.extend_arrays = extend_arrays

    def filter_standard_attrs(self, payload):
        attrs = ["created_at", "id", "owner_id", "updated_at", "user_id"]
        for i in attrs:
            if i in payload:
                payload.pop(i)
        for i in payload:
            if i.startswith("_"):
                payload.pop(i)
        return payload

    def sanitize_prop_val(self, txt):
        """
        Replaces characters that can be problematic in property values on the ENCODE Portal. 
        For example, the '/' character in an alias is a problem since the alias is an identifying property
        that can be used in a URL to view the record. In this case, the '/' will be interpreted as a
        path separator. 

        Characters that get replaced: currently, just '/' with '-'. 

        Args:
            txt: `str`. The value to clean.

        Returns:
            `str`. The cleaned value that is submission acceptable. 
        """
        txt = txt.replace("/","-").strip()
        # Replace contiguous spaces with a single space
        return space_reg.sub(" ", txt)
        

    def get_vendor_id_from_encodeportal(self, pulsar_vendor_id):
        """
        Given a Pulsar Vendor record ID, returns the upstream identifier. 

        Raises:
            `UpstreamNotSet`: The Pulsar vendor.upstream_identifier attribute isn't set.
        """
        if not pulsar_vendor_id:
            return ""
        vendor = models.Vendor(pulsar_vendor_id)
        upstream = vendor.upstream_identifier
        if not upstream:
            msg = "Vendor {} with Pulsar ID {} does not have the upstream_identifier attribute set.".format(vendor.name, vendor.id)
            raise UpstreamNotSet(msg)
        return upstream
                  

    def patch(self, payload, upstream_id, dont_extend_arrays=False):
        """
        A wrapper over `encode_utils.connection.Connection.patch()`.

        Args:
            dont_extend_arrays: `bool`. Dynamic way to signal not to extend array property values. 
                If not True, then the boolean value of `self.extend_arrays` determines whether
                arrays are extended. 

        Returns:
            `dict`: The JSON response from the PATCH operation, or an empty dict if the record doesn't
                exist on the Portal.  See ``encode_utils.connection.Connection.patch()`` for more details.
        """
        payload[self.ENC_CONN.ENCID_KEY] = upstream_id
        dont_extend_arrays = True
        if dont_extend_arrays:
            extend = False
        else:
            extend = self.extend_arrays
        response_json = self.ENC_CONN.patch(payload=payload, extend_array_values=extend)
        if not response_json:
            raise Exception("Couldn't PATCH record on the Portal since it doesn't exist.")
        return upstream_id

    def post(self, payload, dcc_profile, pulsar_model, pulsar_rec_id, repost=False):
        """
        A wrapper over `encode_utils.connection.Connection.post()`. 

        First checks if the Pulsar record has an upstream_identifier set, and if set, returns
        it rather than attempting to re-post.

        Adds aliases to the payload being the record's record ID and name. 

        Sets the profile key in the payload.

        If the record is successfully posted to the prod ENCODE Portal, then sets the 
        upstream_identifier attribute in the Pulsar record.
    
        Args:
            payload: `dict`. The new record attributes to submit.
            dcc_profile: `str`. The name of the ENCODE Profile for this record, i.e. 'biosample',
                'genetic_modification'.
            pulsar_model: One of the defined subclasses of the ``models.Model`` class, i.e. 
                ``models.Model.Biosample``, which will be used to set the Pulsar record's 
                upstream_identifier attribute after a successful POST to the ENCODE Portal.
            pulsar_rec_id: `str`. The identifier of the Pulsar record to POST to the DCC.
        Returns:
            `str`: The upstream identifier for the new record on the ENCODE Portal, or the existing
            upstream identifier if the record already exists; see ``encode_utils.utils.get_record_id()``
            for more details. 
        """
        payload[self.ENC_CONN.PROFILE_KEY] = dcc_profile
        pulsar_rec = pulsar_model(pulsar_rec_id)
        upstream = pulsar_rec.upstream_identifier
        if upstream:
          # for other dcc_profile than biosample immunoblot (IP, "biosample_characterization"),
          # repost shall be False to aviod repost a record which already exists on ENCODE Portal.
          if not repost:
            # May sure that upstream exists. Could be that the upstream identifier belongs to a different
            # server, i.e. test or a demo, than the one we are currently connected to. 
            # No need to POST.
            if upstream.startswith("ENC"):
              # For sure this is a production accession and we should leave it alone.
                return upstream

            # Tao: not all upstream ids start with "ENC", e.g., Gel lanes has upstream ids
            #      like michael-snyder:pGL-156 (https://pulsar-encode.herokuapp.com/gels/38)
            # If it exists on ENCODE Portal, and we do not want to repost (because it is not
            # an IP), we return its upstream id immediately.
            exists_on_server = self.ENC_CONN.get(rec_ids=upstream)
            if exists_on_server:
              return upstream

          # Tao: if upstream id does not exist in pulsar or ENCODE Portal, 
          # or even if the upstream Id exists but we need a post for an
          # IP because child biosamples do not have their own IPs but share the same IP of their parent's,
          # we need a repost.
          else:
            assert(dcc_profile == "biosample_characterization"), "A repost only happens for biosample characterization (IP)"

        aliases = payload.get("aliases", [])
        abbrev_alias = pulsar_rec.abbrev_id()
        if abbrev_alias not in aliases:
            aliases.append(abbrev_alias)
        # Add value of 'name' property as an alias, if this property exists for the given model.
        try:
            name = self.sanitize_prop_val(pulsar_rec.name).replace('&', ' and ')
            if name:
                # Need to prepend the model abbreviation to the name since some names are the same
                # between models. For example, its common in Pulsar to have a Library named the
                # same as the Biosample it belongs to.
                alias_name = models.Model.PULSAR_LIMS_PREFIX + pulsar_rec.MODEL_ABBR + "-" + name
                if alias_name not in aliases:
                    aliases.append(alias_name)
        except KeyError:
            pass
        #payload["aliases"] = [aliases[0].replace(":", "-") + "-cs521-rep4"]
        payload["aliases"] = [alias.replace(":", "-") for alias in aliases]

        # make aliases unique since it is the key for ENCODE Portal to detect whether it is a duplicate.
        # Option 1: Myers approach, only submit one IP for both biosample
        #   chip-seq: https://www.encodeproject.org/experiments/ENCSR635VTR/
        #   biosample 1: ENCBS187ICM
        #   biosample 2: ENCBS739WGO
        # Option 2: Tao's approach, submit one IP for each biosample 
        #   chip-seq: https://www.encodeproject.org/experiments/ENCSR602QEJ/
        #   biosample 1: https://www.encodeproject.org/biosamples/ENCBS055BET/
        #   biosample 2: https://www.encodeproject.org/biosamples/ENCBS062KGP/
        # With the following discussion, we decided to take Option 1, by commenting out the following
        # two lines:
        #if repost:
        #  payload['aliases'] = [payload['aliases'][0] + "-" + payload['characterizes']]
        # Discussion:
        # Annika: We need to confirm one thing with you for IP submission. In case the IP is done on the 
        # parent, we submit it twice (one for each child), which is what we discussed before. 
        # However, the same IP shows up twice on the experiment then which is a little weird. 
        # Biosample 1: https://www.encodeproject.org/biosamples/ENCBS062KGP/ Biosample 
        # 2: https://www.encodeproject.org/biosamples/ENCBS055BET/ and corresponding 
        # experiment: https://www.encodeproject.org/experiments/ENCSR602QEJ/
        # Myers Mendehall does it differently: they also have the IP on the parent and 
        # submit it once for only one of the biosamples and the other one doesn't have an IP attached.
        #if repost:
        #  payload['aliases'] = [payload['aliases'][0] + "-" + payload['characterizes']]
    
        # `dict`. The POST response if the record didn't yet exist on the ENCODE Portal, or the
        # record JSON itself if it does already exist. Note that the dict. will be empty if the connection
        # object to the ENCODE Portal has the dry-run feature turned on.
        response_json = self.ENC_CONN.post(payload)
        upstream = euu.get_record_id(response_json)
        # Set value of the Pulsar record's upstream_identifier
        print("Setting the Pulsar record's upstream_identifier attribute to '{}'.".format(upstream))
        pulsar_rec.patch(payload={"upstream_identifier": upstream})
        print("upstream_identifier attribute set successfully.")
        return upstream

    def get_biosample_term_name_and_type(self, biosample):
        """
        Creates a dict. with the keys:
        
          biosample_term_name
          biosample_term_id
          biosample_type 

        Args:
            biosample: `pulsarpy.models.Biosample` instance.

        Returns:
            `dict`.
        """
        res = {} 
        btn = models.BiosampleTermName(biosample.biosample_term_name_id)
        res["biosample_term_name"] = btn.name
        res["biosample_term_id"] = btn.accession
        bty = models.BiosampleType(biosample.biosample_type_id)
        res["biosample_type"] = bty.name
        return res

    def get_exp_of_biosample(self, dcc_biosample_id):
        query_string =  "?searchTerm={}&type=Experiment".format(dcc_biosample_id)
        experiments = self.ENC_CONN.search(url=query_string)
        if len(experiments) > 1:
            raise Exception("Expected to find 1 experiment linked to biosample {}, instead found more: {}.".format(dcc_biosample_id, len(experiments)))
        return experiments[0]["accession"]

    def post_library_through_fastq(self, pulsar_library_id, dcc_exp_id, patch=False):
        """
        POSTS the Biosample, it's latest Library, and all SequencingResults for that Library.

        Args:
            pulsar_library_id: `int`. The ID of a Pulsar Library record.
            dcc_exp_id: `int`. The ID of the experiment record on the Portal to link the replicate to.
        """
        pulsar_library = models.Library(pulsar_library_id)
        
        # determine experiment type to delegate fastq submission APIs
        exp_type = "nonSC" # ChIP or bulkAtac
        if pulsar_library.atacseq_id:
            if pulsar_library.chipseq_experiment_id:
                raise Exception(f'A library must either be part of Atacseq or ChIP, but not both.')

            atacseq_exp = models.Atacseq(pulsar_library.atacseq_id)
            if atacseq_exp.single_cell and not atacseq_exp.snrna:
                if atacseq_exp.multiome:
                    exp_type = "scAtac-multiome"
                else:
                    exp_type = "scAtac-nonMultiome"
            elif not atacseq_exp.single_cell and atacseq_exp.snrna:
                if atacseq_exp.multiome == True:
                    exp_type = "snRNA-multiome"
                else:
                    exp_type = "snRNA-nonMultiome"

        # POST biosample record
        pulsar_biosample_id = pulsar_library.biosample_id
        biosample_upstream = self.post_biosample(pulsar_biosample_id, patch=patch, exp_type=exp_type)
        biosample = models.Biosample(pulsar_biosample_id)
        #if not dcc_exp_id:
        #    dcc_exp_id = self.get_exp_of_biosample(biosample_upstream)
        
        # POST library record
        library_upstream = self.post_library(rec_id=pulsar_library.id, exp_type=exp_type, patch=patch)
        
        # POST replicate record
        replicate_upstream = self.post_replicate(pulsar_library_id=pulsar_library.id, dcc_exp_id=dcc_exp_id, patch=patch)
        
        # POST file records for all sequencing results for the Library
        sres_ids = pulsar_library.sequencing_result_ids
        if not sres_ids:
            # TODO: shall we require existence of SRES for fastq submission?
            error_logger.error(f"Now trying post fastq file directly from SCG for Library {pulsar_library.id}")
            sreq = pulsar_library.sequencing_request_ids

            # the sequencing result does not exist either on DX or SCG.
            msg = "No SequencingResult for Library {} of Biosample {}, exiting.".format(pulsar_library.id, pulsar_biosample_id)
            error_logger.error(msg)
            raise MissingSequencingResult(msg)
        
        if exp_type == "nonSC":
            for i in sres_ids:
                # Single lane ChIP
                #print("ChIP-seq: single lane with 2 Reads.")
                #self.post_sres(pulsar_sres_id=i, enc_replicate_id=replicate_upstream, patch=patch)

                # Double lane ChIP
                print("ChIP-seq: double lanes with 2 Reads per lane.")
                self.post_sres_chip_seq(pulsar_sres_id=i, enc_replicate_id=replicate_upstream, patch=patch)
        elif exp_type == "snRNA-multiome":
            print("single nucleus RNA-seq multiome).")
            for i in sres_ids:
                self.post_sres_snRNA_seq(pulsar_sres_id=i, enc_replicate_id=replicate_upstream, patch=patch)
        elif exp_type == "snRNA-nonMultiome":
            print("single nucleus RNA-seq separated).")
            for i in sres_ids:
                self.post_sres_snRNA_seq(pulsar_sres_id=i, enc_replicate_id=replicate_upstream, patch=patch, multiome=False)
        elif exp_type == "scAtac-multiome":
            print("single cell Atac-seq, multiome.")
            for i in sres_ids:
                self.post_sres_scAtac_seq(pulsar_sres_id=i, enc_replicate_id=replicate_upstream, patch=patch)
        elif exp_type == "scAtac-nonMultiome":
            print("single cell Atac-seq, non-multiome (separated).")
            for i in sres_ids:
                self.post_sres_scAtac_seq(pulsar_sres_id=i, enc_replicate_id=replicate_upstream, patch=patch, multiome=False)
    
    def post_sres_snRNA_seq(self, pulsar_sres_id, enc_replicate_id, patch=False, multiome=True):
        """
        A wrapper over ``self.post_fastq_file()``.  Whereas ``self.post_fastq_file()`` only 
        uploads the FASTQ file for the given read number and lanes, this method calls ``self.post_fastq_file()``
        four times potentially

        For each of the two lanes (L001 and L002), there are 2 calls:
            one for R1 (read 1)
            one for R2 (read 2)
        """
        if multiome:
            exp_type = "snRNA-multiome"
        else:
            exp_type = "snRNA-nonMultiome"

        sres = models.SequencingResult(pulsar_sres_id)
        srun = models.SequencingRun(sres.sequencing_run_id)
        sreq = models.SequencingRequest(srun.sequencing_request_id)
        Lanes = ['L001', 'L002']

        for lane in Lanes:
            # R1: read 1
            self.post_fastq_file(pulsar_sres_id=sres.id, read_num=1, enc_replicate_id=enc_replicate_id, \
                        exp_type=exp_type, lane_num=lane, patch=patch)
            if not sreq.paired_end and sres.read2_uri:
                sres.patch({"paired_end": True})
    
            # R2: read 2
            if sreq.paired_end:
                # Submit read2
                self.post_fastq_file(pulsar_sres_id=sres.id, read_num=2, enc_replicate_id=enc_replicate_id, \
                        exp_type=exp_type, lane_num=lane, patch=patch)
    
    def post_sres_scAtac_seq(self, pulsar_sres_id, enc_replicate_id, patch=False, multiome=True):
        """
        A wrapper over ``self.post_fastq_file()``.  Whereas ``self.post_fastq_file()`` only 
        uploads the FASTQ file for the given read number and lanes, this method calls ``self.post_fastq_file()``
        six times potentially
        For each of the two lanes (L001 and L002), there are 3 calls:
            one for R1 (read 1)
            one for R3 (read 2)
            one for R2 (index read)
        """
        if multiome:
            exp_type = "scAtac-multiome"
        else:
            exp_type = "scAtac-nonMultiome"

        sres = models.SequencingResult(pulsar_sres_id)
        srun = models.SequencingRun(sres.sequencing_run_id)
        sreq = models.SequencingRequest(srun.sequencing_request_id)
        Lanes = ['L001', 'L002']

        for lane in Lanes:
            # R1: read 1
            self.post_fastq_file(pulsar_sres_id=sres.id, read_num=1, enc_replicate_id=enc_replicate_id, \
                        exp_type=exp_type, lane_num=lane, patch=patch)
            if not sreq.paired_end and sres.read2_uri:
                sres.patch({"paired_end": True})
    
            # R3: read 2
            if sreq.paired_end:
                # Submit read2
                self.post_fastq_file(pulsar_sres_id=sres.id, read_num=3, enc_replicate_id=enc_replicate_id, \
                        exp_type=exp_type, lane_num=lane, patch=patch)
    
            # R2: index read
            self.post_fastq_file(pulsar_sres_id=sres.id, read_num=2, enc_replicate_id=enc_replicate_id, \
                        exp_type=exp_type, lane_num=lane, patch=patch)
    
    def post_sres_chip_seq(self, pulsar_sres_id, enc_replicate_id, patch=False):
        """
        A wrapper over ``self.post_fastq_file()``.  Whereas ``self.post_fastq_file()`` only 
        uploads the FASTQ file for the given read number and lanes, this method calls ``self.post_fastq_file()``
        4 times potentially
        For each of the two lanes (L001 and L002), there are 2 calls:
            one for R1 (read 1)
            one for R2 (read 2)
        """
        sres = models.SequencingResult(pulsar_sres_id)
        srun = models.SequencingRun(sres.sequencing_run_id)
        sreq = models.SequencingRequest(srun.sequencing_request_id)
        Lanes = ['L001', 'L002']

        for lane in Lanes:
            # R1: read 1
            self.post_fastq_file(pulsar_sres_id=sres.id, read_num=1, enc_replicate_id=enc_replicate_id, \
                                 lane_num=lane, patch=patch)
            if not sreq.paired_end and sres.read2_uri:
                sres.patch({"paired_end": True})
    
            # R2: read 2
            if sreq.paired_end:
                # Submit read2
                self.post_fastq_file(pulsar_sres_id=sres.id, read_num=2, enc_replicate_id=enc_replicate_id, \
                                     lane_num=lane, patch=patch)
    
    def post_sres(self, pulsar_sres_id, enc_replicate_id, patch=False):
        """
        A wrapper over ``self.post_fastq_file()``.  Whereas ``self.post_fastq_file()`` only 
        uploads the FASTQ file for the given read number, this method calls ``self.post_fastq_file()``
        twice potentially, once for each FASTQ file in the Pulsar SequencingResult. Thus,
        if paired-end sequencing was done, ``self.post_fastq_file()`` will be called twice to upload
        the forward and reverse reads FASTQ files. 
        """
        sres = models.SequencingResult(pulsar_sres_id)
        srun = models.SequencingRun(sres.sequencing_run_id)
        sreq = models.SequencingRequest(srun.sequencing_request_id)
        self.post_fastq_file(pulsar_sres_id=sres.id, read_num=1, enc_replicate_id=enc_replicate_id, patch=patch)
        if not sreq.paired_end and sres.read2_uri:
            sres.patch({"paired_end": True})
        if sreq.paired_end:
            # Submit read2
            self.post_fastq_file(pulsar_sres_id=sres.id, read_num=2, enc_replicate_id=enc_replicate_id, patch=patch)

    def check_if_biosample_has_exp_on_portal(self, dcc_biosample_id):
        """
        Given a Portal biosample record ID, searches the Portal for associated experiment records.
        Any that are found are returned in a list. 

        Args:
            dcc_biosample_id: `str`. A biosample record identifier on the Portal.

        Returns:
            `list` of associated experiment records, where each is JSON-serialized. 

        Raises:
            `Exception`: The biosample is linked to more than one experiment.
        """
        if not dcc_biosample_id:
            return False
        exps = self.ENC_CONN.get_experiments_with_biosample(rec_id=dcc_biosample_id)
        if exps:
            # Should only exist on one experiment. exps is an array of >= 0 experiment records.
            if len(exps) > 1:
                accessions = []
                for i in exps:
                    accessions.append(i["accession"])
                msg = "Error: Biosample {} is associated to more than one Portal experiment: {}.".format(dcc_biosample_id, ", ".join(accessions))
                raise Exception(msg)
            return exps[0]
        return False
        
    def post_chipseq_ctl_exp(self, rec_id, wt_input=False, paired_input=False, exp_only=False, patch=False):
        """
        Creates a control experiment record on the ENCODE Portal for either the paired-input control
        biosample(s) or the wild-type input biosample on the Pulsar ChipseqExperiment.

        Args:
            rec_id: `int`. ID of a ChipseqExperiment record in Pulsar.
            wt_input: `bool`. True means to make a control experiment on the Portal for the wild-type
                input biosample on the Pulsar ChipseqExperiment. Note that either this or the
                `paired_input` parameter must be set to True and not both. 
            paired_input: `bool`. True means to make a control experiment on the Portal for the 
                paired-input control biosample(s) on the Pulsar ChipseqExperiment. Note that either 
                this or the `wild_type` parameter must be set to True and not both. 
            exp_only: `bool`. Only makes sense to use when the `patch` parameter is set to True.
                When `exp_only=True`, then don't PATCH Biosample records and everything downstream
                to the file records on the Portal (don't call `self.post_library_through_fastq()`).

        Returns:
            `str`: The ENCODE Portal accession of the control experiment. 

        Raises:
            `ValueError`: Both parameters `wt_input` and `paired_input` are set to False or True.
                Only one of them must be True. 
        """
        print(">>> IN post_chipseq_ctl_exp()")
        if (not wt_input and not paired_input) or (wt_input and paired_input):
            raise ValueError("Either the wt_input or the paired_input parameter must be set to True.")

        pulsar_exp = models.ChipseqExperiment(rec_id)
        input_ids = []
        if wt_input:
            experiment_type = "wild type"
            # Only 1 Wild Type input per experiment.
            if pulsar_exp.wild_type_control_id:
                input_ids.append(pulsar_exp.wild_type_control_id)
            else:
                raise ExpMissingReplicates("Chip-seq does not have wild-type control.")
        else:
            experiment_type = "paired-input"
            # Normally there will only be one paired_input control Biosample, but there could at times
            # be another. That happens when one of the reps fail, and another rep has to be made from a
            # different cell batch than the sibling rep on the experiment.
            input_ids.extend(
                models.Library(lib_id).biosample_id
                for lib_id in pulsar_exp.control_replicate_ids
            ) # Biosample records.
        if not input_ids:
            msg = "Can't submit {} control exp. for {}: no replicates.".format(experiment_type, pulsar_exp.abbrev_id())
            log_error(msg)
            raise ExpMissingReplicates(msg)
        inputs = [models.Biosample(x) for x in input_ids]
        dcc_exp = ""
        for i in inputs:
            dcc_exp = self.check_if_biosample_has_exp_on_portal(i.upstream_identifier)
            if dcc_exp:
                break
        payload = {}
        alias = ""
        for i in inputs:
            alias += i.abbrev_id()
        alias.rstrip()
        if wt_input:
            alias_prefix = "pWT-CTL_"
            payload["control_type"] = "wild type"
        else:
            alias_prefix = "pPI-CTL_"
            payload["control_type"] = "input library"
        alias = alias_prefix + alias 
        payload["aliases"] = [alias]
        payload.update(self.get_exp_core_payload_props(pulsar_exp_rec=pulsar_exp, assay_term_name="ChIP-seq"))
        payload["description"] = "ChIP-seq on human " + self.ENC_CONN.get(payload["biosample_ontology"])["term_name"]
  
        # Before POSTING experiment, check if it already exists on the Portal.
        # POST experiment. Don't use self.post() since there isn't a Pulsar model for a control experiment.
        # So, use encode-utils directly to POST.
        if patch:
            if not dcc_exp:
                msg = "Can't PATCH " + alias + " control experiment since it wasn't found on the Portal."
                raise Exception(msg)
            ctl_exp_accession = self.patch(payload=payload, upstream_id=dcc_exp["accession"])
        else:
            # post
            if not dcc_exp:
                payload[self.ENC_CONN.PROFILE_KEY] = "experiment"
                dcc_exp = self.ENC_CONN.post(payload=payload)
            ctl_exp_accession = dcc_exp["accession"]

        if (patch and not exp_only) or not patch:
            for b in inputs:
                for lib_id in b.library_ids:
                    self.post_library_through_fastq(pulsar_library_id=lib_id, dcc_exp_id=ctl_exp_accession, patch=patch)
        return ctl_exp_accession

    def post_snRNAseq_exp(self, rec_id, patch=False, patch_all=False):
        """
        Args:
            rec_id: `int`. ID of an single-cell RNA sequencing assay (experiment) record in Pulsar.
            patch: `bool`. True means to patch the DCC experiment record.
            patch_all: `bool`. True means to patch not just the experiment record, but its sub-entities
                also, i.e. biosamples, libraries, replicates, ... Setting this to True automatically
                sets `patch` to True as well.

        Returns:
        """
        if patch_all:
            patch = True
        pulsar_exp = models.Atacseq(rec_id)
        assert((not pulsar_exp.single_cell and pulsar_exp.snrna)), "Not a single-cell RNA sequencing assay."
        pulsar_exp_upstream = pulsar_exp.upstream_identifier
        payload = {}
        payload.update(self.get_exp_core_payload_props(pulsar_exp_rec=pulsar_exp, assay_term_name="single-cell RNA sequencing assay"))
        desc = pulsar_exp.description.strip()
        if desc:
            payload["description"] = desc
        # submit experiment
        if patch:
            dcc_exp_accession = self.patch(payload=payload, upstream_id=pulsar_exp_upstream)
        if patch_all or not patch:
            dcc_exp_accession = self.post(payload=payload, dcc_profile="experiment", pulsar_model=models.Atacseq, pulsar_rec_id=rec_id)
            self.post_experimental_reps(rec_id=rec_id, experiment_type="atac-seq", patch=patch)

        return dcc_exp_accession

    def post_scAtacseq_exp(self, rec_id, patch=False, patch_all=False):
        """
        Args:
            rec_id: `int`. ID of an AtacSeq experiment record in Pulsar. Should be a single-cell atac experiment. 
            patch: `bool`. True means to patch the DCC experiment record.
            patch_all: `bool`. True means to patch not just the experiment record, but its sub-entities
                also, i.e. biosamples, libraries, replicates, ... Setting this to True automatically
                sets `patch` to True as well.

        Returns:
        """
        if patch_all:
            patch = True
        pulsar_exp = models.Atacseq(rec_id)
        pulsar_exp_upstream = pulsar_exp.upstream_identifier
        payload = {}
        payload.update(self.get_exp_core_payload_props(pulsar_exp_rec=pulsar_exp, assay_term_name="single-nucleus ATAC-seq"))
        desc = pulsar_exp.description.strip()
        if desc:
            payload["description"] = desc
        # submit experiment
        if patch:
            dcc_exp_accession = self.patch(payload=payload, upstream_id=pulsar_exp_upstream)
        if patch_all or not patch:
            dcc_exp_accession = self.post(payload=payload, dcc_profile="experiment", pulsar_model=models.Atacseq, pulsar_rec_id=rec_id)
            self.post_experimental_reps(rec_id=rec_id, experiment_type="atac-seq", patch=patch)

        return dcc_exp_accession

    def post_bulk_atacseq_exp(self, rec_id, patch=False, patch_all=False):
        """
        Args:
            rec_id: `int`. ID of an AtacSeq experiment record in Pulsar. Should be a bulk and not
                a single-cell experiment. 
            patch: `bool`. True means to patch the DCC experiment record.
            patch_all: `bool`. True means to patch not just the experiment record, but its sub-entities
                also, i.e. biosamples, libraries, replicates, ... Setting this to True automatically
                sets `patch` to True as well.

        Returns:
        """
        if patch_all:
            patch = True
        pulsar_exp = models.Atacseq(rec_id)
        pulsar_exp_upstream = pulsar_exp.upstream_identifier
        payload = {}
        payload.update(self.get_exp_core_payload_props(pulsar_exp_rec=pulsar_exp, assay_term_name="ATAC-seq"))
        desc = pulsar_exp.description.strip()
        if desc:
            payload["description"] = desc
        # submit experiment
        if patch:
            dcc_exp_accession = self.patch(payload=payload, upstream_id=pulsar_exp_upstream)
        if patch_all or not patch:
            dcc_exp_accession = self.post(payload=payload, dcc_profile="experiment", pulsar_model=models.Atacseq, pulsar_rec_id=rec_id)
            self.post_experimental_reps(rec_id=rec_id, experiment_type="atac-seq", patch=patch)

        return dcc_exp_accession

    def post_chipseq_exp(self, rec_id, patch=False):
        """
        Args:
            rec_id: `int`. ID of a ChipseqExperiment record in Pulsar.

        Returns:
            `str`: The ENCODE Portal accession of the control experiment. 

        Raises:
            `ValueError`: Both parameters `wt_input` and `paired_input` are set to False or True.
                Only one of them must be True. 
        """
        pulsar_exp = models.ChipseqExperiment(rec_id)
        pulsar_exp_upstream = pulsar_exp.upstream_identifier
        payload = {}
        payload.update(self.get_exp_core_payload_props(pulsar_exp_rec=pulsar_exp, assay_term_name="ChIP-seq"))
        target = models.Target(pulsar_exp.target_id)
        target_upstream = target.upstream_identifier
        if not target_upstream:
            msg = "Target {} missing upstream identifier.".format(target.abbrev_id())
            log_error(msg)
            raise MissingTargetUpstream(msg)
        payload["target"] = target.upstream_identifier.replace('eGFP-', '').replace('3xFLAG-', '')
        #payload["description"] = pulsar_exp.description.strip()
        payload["description"] = target.upstream_identifier.rstrip('-human') + ' ChIP-seq on human ' + self.ENC_CONN.get(payload["biosample_ontology"])["term_name"]
        # submit experiment
        if patch:
            dcc_exp_accession = self.patch(payload=payload, upstream_id=pulsar_exp_upstream)
        else:
            dcc_exp_accession = self.post(payload=payload, dcc_profile="experiment", pulsar_model=models.ChipseqExperiment, pulsar_rec_id=rec_id)
            # Then POST WT-input and paired-input control experiments. The WT-input is shared across
            # multiple experiments from the same starting batch, so it's possible that it was POSTED
            # already during submission of a related experiment. 
            self.post_chipseq_control_experiments(rec_id=rec_id)
            # POST experimental biosampes
            self.post_experimental_reps(rec_id=rec_id, experiment_type="chip-seq")

        # Add-in/PATCH possible_controls property
        self.patch_chipseq_possible_controls(pulsar_exp.id)
        return dcc_exp_accession

    def patch_chipseq_possible_controls(self, pulsar_exp_id):
        possible_controls = self.get_chipseq_possible_controls(pulsar_exp_id)
        if not possible_controls:
            raise Exception("Either WT contorl or input library control or both should exist.")

        payload = {}
        payload["possible_controls"] = possible_controls
        exp = models.ChipseqExperiment(pulsar_exp_id)
        self.patch(payload=payload, upstream_id=exp.upstream_identifier, dont_extend_arrays=True)
        
    def get_chipseq_possible_controls(self, pulsar_exp_id):
        possible_controls = []
        exp = models.ChipseqExperiment(pulsar_exp_id)
        # Tao Wang: there may be not wild-type control
        if exp.wild_type_control_id:
            wt = models.Biosample(exp.wild_type_control_id)
            wt_upstream = wt.upstream_identifier
            wt_ctl_exp = self.check_if_biosample_has_exp_on_portal(wt_upstream)
            if not wt_ctl_exp:
                raise Exception("WT input {} on ChipseqExperiment {} doesn't have an upstream control experiment record.".format(wt.abbrev_id(), exp.abbrev_id()))
            possible_controls.append(wt_ctl_exp["accession"])
        
        # input library control may also missing. 
        # But either wild-type control or input library control should exist.
        pis = [
            models.Biosample(models.Library(lib_id).biosample_id)
            for lib_id in exp.control_replicate_ids
        ]
        for i in pis:
            upstream = i.upstream_identifier
            pi_ctl_exp = self.check_if_biosample_has_exp_on_portal(upstream)
            if not i and not pi_ctl_exp:
                raise Exception("Paired input {} on ChipseqExperiment {} doesn't have an upstream control experiment record.".format(i.abbrev_id(), exp.abbrev_id()))
            possible_controls.append(pi_ctl_exp["accession"])

        return list(set(possible_controls))
        
    def post_chipseq_control_experiments(self, rec_id):
        """
        Tao's note: following discussions w/ Annika, we accept one of them
        {WT-input, Paired-input, both}, while the previous restriction only
        accept 'both'.

        POSTS the WT input and the paired input controls that are associated to the indicated 
        ChipseqExperiment in Pulsar, turning each into an experiment record on the Portal.

        Args:
            rec_id: `int`. ID of a ChipseqExperiment record in Pulsar.
        """
        print(">>> IN post_chipseq_control_experiments()")
        controls = 0

        # First the WT-input:
        try:
            self.post_chipseq_ctl_exp(rec_id=rec_id, wt_input=True)
        except Exception as ex:
            print(ex)
            print("Wild-type control experiment is unavailable")
        else:
            controls += 1
            print("Wild-type control experiment is done")

        # Then the Paired-input, which is unique to this experiment.
        try:
            print ("posting input library control.")
            self.post_chipseq_ctl_exp(rec_id=rec_id, paired_input=True)
        except Exception as ex:
            print(ex)
            print("Input-library control experiment is unavailable.")
        else:
            print("Input-library control experiment is available.")
            controls += 1

        assert(controls >= 1), "Either Wild-type control or input-library control or both should be available"


    def post_experimental_reps(self, rec_id, experiment_type, patch=False):
        """
        POSTS the experimental replicates of a ChipseqExperiment or bulk Atacseq experiment object.

        Args:
            rec_id: `int`. ID of a ChipseqExperiment record in Pulsar.
            experiment_type: `str`. Either of chip-seq or atac-seq.
        """
        if experiment_type == "chip-seq":
            pulsar_exp = models.ChipseqExperiment(rec_id)
        elif experiment_type == "atac-seq":
            pulsar_exp = models.Atacseq(rec_id)
        else:
            raise Exception("Unknown experiment type '{}' passed to experiment_type parameter.".format(experiment_type))
        rep_ids = pulsar_exp.replicate_ids
        for i in rep_ids:
        #for i in [11817]:
            self.post_library_through_fastq(pulsar_library_id=i, dcc_exp_id=pulsar_exp.upstream_identifier, patch=patch)

    def get_exp_core_payload_props(self, pulsar_exp_rec, assay_term_name):
        """
        Args:
            pulsar_exp_rec: `str`. `pulsarpy.models` subclass being either ChipSeq or Atacseq.
            assay_term_name: `str`. Either 'ChIP-seq' or ATAC-seq.
        """
        payload = {}
        first_rep_library = models.Library(pulsar_exp_rec.replicate_ids[0])
        first_rep_biosample = models.Biosample(first_rep_library.biosample_id)
        # Add biosample_term_name, biosample_term_id, and biosample_type props
        btn = models.BiosampleTermName(first_rep_biosample.biosample_term_name_id)
        bty = models.BiosampleType(first_rep_biosample.biosample_type_id)
        payload["biosample_ontology"] = self.ENC_CONN.get_biosample_type(classification=bty.name, term_id=btn.accession)["@id"]
        payload["assay_term_name"] = assay_term_name
        payload["documents"] = self.post_documents(pulsar_exp_rec.document_ids)
        submitter_comment = pulsar_exp_rec.submitter_comments.strip()
        if submitter_comment:
            payload["submitter_comment"] = submitter_comment
        return payload
    
    def post_crispr_modification(self, rec_id, patch=False):
        rec = models.CrisprModification(rec_id)
        # CrisprConstruct(s)
        cc_ids = rec.crispr_construct_ids
        ccs = [models.CrisprConstruct(i) for i in cc_ids]
        # DonorConstruct
        dc = models.DonorConstruct(rec.donor_construct_id)
        dc_target = models.Target(dc.target_id)
        target_upstream = dc_target.upstream_identifier
        if not target_upstream:
            msg = "Target {} missing upstream identifier.".format(dc_target.abbrev_id())
            log_error(msg)
            raise MissingTargetUpstream(msg)

        payload = {}
        payload["category"] = rec.category # Required
        desc = rec.description.strip()
        if desc:
            payload["description"]
        payload["documents"] = self.post_documents(rec.document_ids)

        guide_seqs = list(c.guide_sequence for c in ccs)
        payload["guide_rna_sequences"] = guide_seqs

        if rec.category in ["insertion", "replacement"]:
            pass # The insert can be viewed in addgene. This doesn't look good to show on the Portal.
            #payload["introduced_sequence"] = dc.insert_sequence.upper()

        payload["method"] = "CRISPR"       # Required
        payload["modified_site_by_target_id"] = dc_target.upstream_identifier.replace('eGFP-', '').replace('3xFLAG-', '')
        payload["purpose"] = rec.purpose   # Required

        # Note that CrisprConstruct can also has_many construct_tags. Those are not part of the donor
        # insert though. 
        construct_tags = [models.ConstructTag(i) for i in dc.construct_tag_ids]
        construct_tag_names = [x.name for x in construct_tags]
        seen_tags = []
        introduced_tags = []
        for tag in construct_tag_names:
            if tag.startswith("eGFP"):
                # Pulsar has two eGFP tags that differ by the linker sequence:
                #    1) eGFP (MH170480)
                #    2) eGFP (MH170481)
                # The Portal, however, only has eGFP and it makes most sense to submit this as 
                # simply eGFP and mention the linker used elsewhere. 
                tag = "eGFP"
            if tag not in seen_tags:
                # Avoid potential for duplicate tags, which are not allowed of course on Portal.
                seen_tags.append(tag) 
                introduced_tags.append({"name": tag, "location": "C-terminal"})
        if not introduced_tags:
            # tags are required for modifications on the Portal.
            introduced_tags = [{"name": "eGFP", "location": "C-terminal"}]
        payload["introduced_tags"] = introduced_tags
        reagents = []
        for i in ccs + [dc]:
            addgene_id = getattr(i, "addgene_id")
            if addgene_id:
                r = {}
                r["source"] = "addgene"
                r["url"] = "http://www.addgene.org/" + addgene_id
                r["identifier"] = "addgene:" + addgene_id
                reagents.append(r)
        if reagents:
            payload["reagents"] = reagents
        # ex: ENCGM094ZOS

        if patch: 
            upstream_id = self.patch(payload=payload, upstream_id=rec.upstream_identifier)
        else:
            upstream_id = self.post(payload=payload, dcc_profile="genetic_modification", pulsar_model=models.CrisprModification, pulsar_rec_id=rec_id)
        return upstream_id
    
    def post_document(self, rec_id, patch=False):
        rec = models.Document(rec_id)
        payload = {}
        desc = rec.description.strip()
        if desc:
            payload["description"] = rec.description
        doc_type = models.DocumentType(rec.document_type_id)
        payload["document_type"] = doc_type.name
        content_type = rec.content_type
        # Create attachment for the attachment prop
        file_contents = rec.download()
        data = base64.b64encode(file_contents)
        temp_uri = str(data, "utf-8")
        href = "data:{mime_type};base64,{temp_uri}".format(mime_type=content_type, temp_uri=temp_uri)
        attachment = {}
        attachment["download"] = rec.name
        attachment["type"] = content_type 
        attachment["href"] = href
        payload["attachment"] = attachment
        if patch:
            upstream_id = self.patch(payload, rec.upstream_identifier)
        else:
            upstream_id = self.post(payload=payload, dcc_profile="document", pulsar_model=models.Document, pulsar_rec_id=rec_id)
        return upstream_id

    def post_documents(self, rec_ids, patch=False):
        upstreams = []
        for i in rec_ids:
            upstreams.append(self.post_document(rec_id=i, patch=patch))
        return upstreams

    def post_treatments(self, rec_ids, patch=False):
        upstreams = []
        for i in rec_ids:
            upstreams.append(self.post_treatment(rec_id=i, patch=patch))
        return upstreams

    def post_treatment(self, rec_id, patch=False):
        rec = models.Treatment(rec_id)
        payload = {}
        conc = rec.concentration
        if conc:
            payload["amount"] = conc
            conc_unit = models.Unit(rec.concentration_unit_id)
            payload["amount_units"] = conc_unit.name
        duration = rec.duration
        if duration:
            payload["duration"] = duration
            payload["duration_units"] = rec.duration_units
        temp = rec.temperature_celsius
        if temp:
            payload["temperature"] = temp
            payload["temperature_units"] = "Celsius"
        ttn = models.TreatmentTermName(rec.treatment_term_name_id)
        payload["treatment_term_id"] = ttn["accession"]
        payload["treatment_term_name"] = ttn["name"]
        payload["treatment_type"] = rec.treatment_type
        payload["documents"] = self.post_documents(rec.document_ids)
        # Submit
        if patch:
            upstream_id = self.patch(pyaload, rec.upstream_identifier)
        else:
            upstream_id = self.post(payload=payload, dcc_profile="treatment", pulsar_model=models.Treatment, pulsar_rec_id=rec_id)
        return upstream_id

    
    def post_vendor(self, rec_id, patch=False):
        """
        Vendors must be registered directly by the DCC personel. 
        """
        raise Exception("Vendors must be registered directly by the DCC personel.")

    def get_gel_lane_with_biosample(self, immunoblot_id, biosample_id, parent_biosample_id):
        """
        Given an Immunoblot record ID, and a Biosample record ID, returns the GelLane object with the 
        given Biosample. This method assumes that a Gel won't have more than one GelLane with the same
        Biosample. 

        Note that there should only be 1 Gel, even though the Rails Immunoblot model allows many -
        on the 'to fix list'. 

        Args:
            immunoblot_id: `int`. Immunoblot record ID.
            biosample_id: `int`. Biosample record ID.

        Returns:
           `None` if the GelLane didn't pass. Otherwise, a `pulsarpy.models.GelLane` instance.  
       
        Raises:
            `IpLaneException`: One of multiple issues that could be present as indicated by the error
              message, i.e. 

                * The Biosample doesn't have an associated Gel
                * There isn't a GelLane with the Biosample on it. 
        """
        ip = models.Immunoblot(immunoblot_id)
        if not ip.gel_ids:
            raise IpLaneException("IP {} for Biosample {} does not have a Gel.".format(ip.id, biosample_id))
        gel = models.Gel(ip.gel_ids[0])
        gl = "" # GelLane
        for gel_lane_id in gel.gel_lane_ids:
            gel_lane = models.GelLane(gel_lane_id)
            if biosample_id == gel_lane.biosample_id or parent_biosample_id == gel_lane.biosample_id:
                gl = gel_lane
        if not gl:
            raise IpLaneException("Could't find a GelLane that has Biosample {} on Immunoblot {}.".format(biosample_id, immunoblot_id))
        if not gl.attrs["pass"]:
            print("GelLane didn't pass for Biosample {}.".format(biosample_id))
            raise ValueError("GelLane didn't pass")
            #return None
        return gl

    def get_root_biosample_id(self, child_biosample_id):
      """
      To submit IP for a biosample, submission of its root ancester biosample's IP has to be done first.
      The root ancester is ALWAYS wild-type biosample and ALWAYS has its own IP to show background
      noise within itself. The root ancester is usually a biosample batch of a certain cell line.
      For example, cell line SK-N-SH (represented by SK-N-SH master B-2507) has several biosample 
      batches, such as SK-N-SH batch 0 B-2527, SK-N-SH batch 1 B-2531 and SK-N-SH batch 2 B-2535 among
      many other batches. Each of these batches is a unique root ancester for many child biosamples.
      The family tree looks like:
      
                             SK-N-SH master B-2507
                              /         |          \
                            /           |           \__________________\
                          /             |                               \
      SK-N-SH batch 1 B-2531 (root)   SK-N-SH batch 2 B-2535 (root)   SK-N-SH batch 0 B-2527 (root) ....
                /   |   \              /   |   \                          /  |  \
              /     |    \           /     |    \                        /   |   \
          B-2532  B-2533  ...
            / \
          /    \
       B-2534  ...

       So, all direct descedents of master B-2507 are biosample batches which are the root biosample for their
       respective biosample family tree.
      """
      root_biosample_id = None
      if not child_biosample_id:
        return root_biosample_id
      child_biosample = models.Biosample(child_biosample_id)
      parent_biosample_id = child_biosample.part_of_id
      
      while parent_biosample_id:
        root_biosample_id = child_biosample_id
        child_biosample_id = parent_biosample_id
        child_biosample = models.Biosample(child_biosample_id)
        parent_biosample_id = child_biosample.part_of_id

      return root_biosample_id


    def post_ip_biosample_characterization(self, immunoblot_id, biosample_id, patch=False):
        """
        Submits a Pulsar Immunoblot for a specific lane (biosample) on a Gel to the ENCODE biosample_characterization
        profile. Such an immunoblot is used to show whether the eGFP-tagged target (using CRISPR) 
        is expressed (has a band in the size range of the expected taget size).
        Only submit these after the ChipSeq experiment (and hence CrisprModification) has been submitted. Even though some Biosamples
        have a successful IP, they don't all need to be submitted. For example, in one case a Biosample
        was lost after a successful IP and hence couldn't do the crosslinking for ChIP later on.
        Another reason may be that we already have enough validated Biosamples to submit.

        This method makes the assumption that a given gel won't have more than one lane with the same
        Biosample.

        Returns:
            `None`: The Biosmaple isn't already registered on the Portal.
            `None`: The Biosample has an IP, but not one that passes (based on the GelLane.pass attribute)
            `None`: The Non-WT Biosample isn't yet registerd on the Portal
            `None`: The non-WT biosample that doesn't have a ChipSeq object.
            `int`: The ID of the created biosample_characterization record on the Portal. 
        """
        GEL_IMAGE_DIR = os.path.join(os.path.curdir, "gel_images")
        if not os.path.exists(GEL_IMAGE_DIR):
            os.mkdir(GEL_IMAGE_DIR)

        biosample = models.Biosample(biosample_id)
        if not biosample.upstream_identifier:
            # For now, don't sumbit until jadrian says otherwise. 
            print("Biosample missing upstream - skipping.")
            return None

        # At the early stage of ENCODE4, only parent biosamples gets biosample
        # characterization (IP), all its child biosamples do not have one as their
        # own but share their parent's IP.
        parent_biosample_id = None
        if not immunoblot_id:    
          parent_biosample_id = biosample.part_of_id
          parent_biosample = models.Biosample(parent_biosample_id)
          if parent_biosample['immunoblot_ids']:
            immunoblot_id = parent_biosample['immunoblot_ids'][0]
          else:
            print("Biosample and its parent do not have immunoblot - skipping")
            raise ValueError("IP does not exist on either parent or child biosample.")
            #return None

        ip = models.Immunoblot(immunoblot_id)
        gl = self.get_gel_lane_with_biosample(immunoblot_id=immunoblot_id, biosample_id=biosample_id, parent_biosample_id=parent_biosample_id)
        if not gl:
            return None

        payload = {}
        #if not biosample.wild_type:
        #    wt_biosample_id = self.get_root_biosample_id(biosample_id)
        #    assert(wt_biosample_id), "Each biosample must have a wild-type root biosample."
        #    ## Find WT parent that has an associated Immunoblot to use as control
        #    ## and set that Biosample's upstream_identifier as the value of the ENCODE property 
        #    ## biosample_characterization.wildtype_biosample. Currently, the WT biosample is determined
        #    ## soley by biosample_term_name.
        #    wt_biosample = models.Biosample(wt_biosample_id)
        #    wt_biosample_upstream = wt_biosample.upstream_identifier
        #    if not wt_biosample_upstream:
        #        print("POSTING WT parent Biosample.")
        #        wt_biosample_upstream = self.post_biosample(wt_biosample_id)
        #    # Then POST the Immunoblot linked to the WT Parent Biosample. Note that it's possible
        #    # but unlikely for a Biosample to be linked to multiple Immunoblots. In that case, the
        #    # first one will be submitted.
        #    if wt_biosample.immunoblot_ids:
        #      wt_ip_id = wt_biosample.immunoblot_ids[0]
        #    else:
        #      biosample_term_name = models.BiosampleTermName(biosample.biosample_term_name_id).name
        #      if biosample_term_name == "A549":
        #        wt_biosample_id = 2551
        #      elif biosample_term_name == "GM23338":
        #        wt_biosample_id = 2559
        #      elif biosample_term_name == "HepG2":
        #        wt_biosample_id = 2510
        #      elif biosample_term_name == "MCF-7":
        #        wt_biosample_id = 2515
        #      elif biosample_term_name == "SK-N-SH":
        #        wt_biosample_id = 11200
        #      else:
        #        msg = "Can't submit IP biosample_characterization for Biosample {} IP {} since the wild type biosample with its own Immunoblot can't be determined for biosample term name {}.".format(biosample.id, immunoblot_id, biosample_term_name)
        #        error_logger.error(msg)
        #        return None
        #      wt_biosample = models.Biosample(wt_biosample_id)
        #      wt_biosample_upstream = wt_biosample.upstream_identifier
        #      if not wt_biosample_upstream:
        #          print("POSTING WT parent Biosample.")
        #          wt_biosample_upstream = self.post_biosample(wt_biosample_id)
        #      wt_ip_id = wt_biosample.immunoblot_ids[0]
        #    print("POSTING WT parent Biosample's Immunoblot.")            
        #    self.post_ip_biosample_characterization(immunoblot_id=wt_ip_id, biosample_id=wt_biosample_id, patch=False)
        #    payload["wildtype_biosample"] = wt_biosample_upstream
        
        payload["characterization_method"] = "immunoblot"
        payload["characterizes"] = biosample.upstream_identifier
        payload["documents"] = self.post_documents(ip.document_ids)
        payload["review"] = {"lab": "richard-myers", "lane": gl.lane_number}
        # Process attachment property for gel image.
        # A Pulsar Gel object can have many GelImages (different exposure times), but Jess has indicated
        # to take the first one if multiple are present. 
        gel = models.Gel(gl.gel_id)
        if not gel.gel_image_ids:
            msg = "GelLane {} of Gel {} for Biosample {} is missing a GelImage.".format(gl.id, gel.id, biosample_id)
            error_logger.error(msg)
            raise IpLaneException(msg)
        gel_image = models.GelImage(sorted(gel.gel_image_ids)[0])
        # The image URI is expected to have public read permission.
        # Some paths store a // at the beginning to tell the browser to use the same protocol as it's currently
        # using (HTTP/HTTPS). In that case, just prefix it with 'https:'.
        image_uri = gel_image.image
        if image_uri.startswith("//"):
            image_uri = "https:" + image_uri
        image_basename = os.path.basename(image_uri)
        image_exists_locally = os.path.join(GEL_IMAGE_DIR, image_basename)
        if not os.path.exists(image_exists_locally):
            # Then download it
            stream = requests.get(image_uri, stream=True)
            fout = open(image_exists_locally, 'wb')
            for line in stream.iter_content(chunk_size=512):
                fout.write(line)
            fout.close()
        payload["attachment"] = {"path": image_exists_locally}
        # Caption
        btn = models.BiosampleTermName(biosample.biosample_term_name_id).name
        caption = "Immunoprecipitation was performed on nuclear extracts from biosample {}".format(biosample.upstream_identifier)
        if biosample.wild_type:
            caption += " ({} wild type)".format(btn)
        else:
            if not biosample.chipseq_experiment_ids:
                msg = "Biosample {} is not linked to any ChipSeq experiments.".format(biosample_id)
                error_logger.error(msg)
                raise ValueError("Replicate not used in Chip experiment.")
                #return None
                #raise IpLaneException(msg)
            elif len(biosample.chipseq_experiment_ids) > 1:
                msg = "Biosample {} is linked to more than 1 ChipSeq experiment. It is not known as to which one this IP relates.".format(biosample_id)
                error_logger.error(msg)
                raise IpLaneException(msg)
            chipseq_exp = models.ChipseqExperiment(biosample.chipseq_experiment_ids[0])
            if not chipseq_exp.upstream_identifier:
                msg = "ChipSeq experiment {} for Biosample {} needs to be submitted prior to submitting the IP biosample_characterization.".format(chipseq_exp.id, biosample_id)
                error_logger.error(msg)
                raise IpLaneException(msg)
            try:
              crispr_modification = models.CrisprModification(biosample.crispr_modification_id)
            except ValueError as ve:
              print("Cannot submit IP for Biosample {} without genetic modification in pulsar.".format(biosample_id))
              raise ValueError("Crispr genetic modification does not exist")
            if not crispr_modification.upstream_identifier:
                msg = "Biosample {} has a CrisprModification, but it isn't yet registered on the Portal.".format(biosample_id)
                error_logger.error(msg)
                raise IpLaneException(msg)
            crispr_construct = models.CrisprConstruct(crispr_modification.crispr_construct_ids[0])
            target = models.Target(crispr_construct.target_id)
            # Get biosample_replicate_number on experiment in Portal
            rep_hash = encode_utils.replicate.ExpReplicates(self.ENC_CONN, chipseq_exp.upstream_identifier).rep_hash
            try:
              brn = rep_hash[biosample.upstream_identifier]["brn"]
            except KeyError as ke:
              raise ValueError("Replicate not used in Chip experiment.")
            caption += " ({} eGFP-{} replicate {})".format(btn, target.name, brn)
        caption += " cells using anti-eGFP antibody."
        if btn in ['K562', 'GM12878', 'MCF-7']:
          caption += "The image shows western blot analysis of immunoprecipitation."
        elif btn in ['SK-N-SH', 'A549', 'HepG2', 'HEPG2', 'WTC-11']:
          caption += "The image shows western blot analysis of input"
          caption += " (lane 1),"
          if not biosample.wild_type:
              caption += " immunoprecipitate (lane 2), and mock immunoprecipitate using IgG (lane 3)."
          else:
              caption += " and immunoprecipitate (lane 2)."
        else:
            assert(False), "New biosample cell line name. Need to be in one of the above two btn list!"

        caption += " Molecular weight standard (Bio-Rad, cat. # 161-0374) contains 10"
        caption += " pre-stained recombinant proteins (250, 150, 100, 75, 50, 37, 25, 20, 15, and 10 kD)."
        
        if btn in ['K562', 'GM12878', 'MCF-7']:
          caption += "The target molecular weight is indicated."
        elif btn in ['SK-N-SH', 'A549', 'HepG2', 'HEPG2', 'WTC-11']:
          if not biosample.wild_type:
              caption += " The target molecular weight is {} kD as indicated with an arrow.".format(gl.expected_product_size)
          if gl.low_target_band_intensity:
              caption += " Lower size bands which may be due to potential degradation products are marked with asterisks."
        else:
          assert(False), "New biosample cell line name. Need to be in one of the above two btn list!"
       
        payload["caption"] = caption

        submitter_comment = ip.submitter_comments
        if submitter_comment:
            payload["submitter_comment"] = submitter_comment

        # Submit payload
        if patch:  
          upstream_id = self.patch(payload, gl.upstream_identifier)
        else:
          upstream_id = self.post(payload=payload, dcc_profile="biosample_characterization", pulsar_model=models.GelLane, pulsar_rec_id=gl.id, repost=False)
        return upstream_id

    def post_biosample(self, rec_id, patch=False, exp_type=None):
        rec = models.Biosample(rec_id)
        # The alias lab prefixes will be set in the encode_utils package if the DCC_LAB environment
        # variable is set.
        payload = {}
        # Add biosample_term_name, biosample_term_id, biosample_type props
        btn = models.BiosampleTermName(rec.biosample_term_name_id)
        bty = models.BiosampleType(rec.biosample_type_id)
        payload["biosample_ontology"] = self.ENC_CONN.get_biosample_type(classification=bty.name, term_id=btn.accession)["@id"]

        if exp_type in ["scAtac-multiome", "snRNA-multiome", "snRNA-nonMultiome"]:
            payload["subcellular_fraction_term_name"] = "nucleus"

        date_biosample_taken = rec.date_biosample_taken
        if date_biosample_taken:
            if bty.name == "tissue":
                payload["date_obtained"] = date_biosample_taken
            else:
                payload["culture_harvest_date"] = date_biosample_taken

        desc = rec.description.strip()
        if desc:
            payload["description"] = desc

        donor = models.Donor(rec.donor_id)
        donor_upstream = donor.get_upstream() 
        if not donor_upstream:
            raise Exception("Donor '{}' of biosample '{}' does not have its upstream set. Donors must be registered with the DCC directly.".format(donor.id, rec_id))
        payload["donor"] = donor_upstream

        lot_id = rec.lot_identifier
        if lot_id:
            payload["lot_id"] = lot_id

        nih_cert = rec.nih_institutional_certification
        assert(nih_cert), "Must have an NIH certificate!"
        if nih_cert:
            payload["nih_institutional_certification"] = nih_cert

        payload["organism"] = "human"

        passage_number = rec.passage_number
        if passage_number:
            payload["passage_number"] = passage_number

        starting_amount = rec.starting_amount
        if starting_amount:
            payload["starting_amount"] = starting_amount
            payload["starting_amount_units"] = rec["starting_amount_units"]
            #payload["starting_amount_units"] = models.Unit(rec.starting_amount_units_id).name

        submitter_comment = rec.submitter_comments
        if submitter_comment:
            payload["submitter_comment"] = submitter_comment

        preservation_method = rec.tissue_preservation_method
        if preservation_method:
            payload["preservation_method"] = preservation_method

        prod_id = rec.vendor_product_identifier
        if prod_id:
            payload["product_id"] = prod_id
    
        cm_id = rec.crispr_modification_id
        if cm_id:
            payload["genetic_modifications"] = [self.post_crispr_modification(cm_id)]
    
        payload["documents"] = self.post_documents(rec.document_ids)
    
        part_of_biosample_id = rec.part_of_id
        if part_of_biosample_id:
            part_of_biosample = models.Biosample(part_of_biosample_id)
            pob_upstream = part_of_biosample.get_upstream() 
            if not pob_upstream or not pob_upstream.startswith("ENCBS"):
                pob_upstream = self.post_biosample(part_of_biosample_id, exp_type=exp_type)
            payload["part_of"] = pob_upstream
    
        pooled_from_biosample_ids = rec.pooled_from_biosample_ids
        if pooled_from_biosample_ids:
            pooled_from_biosamples = [models.Biosample(p) for p in pooled_from_biosample_ids]
            payload["pooled_from"] = []
            for p in pooled_from_biosamples:
                p_upstream = p.get_upstream() 
                if not p_upstream:
                    p_upstream = self.post_biosample(p.id, exp_type=exp_type)
                payload["pooled_from"].append(p_upstream)
    
        if rec.vendor_id:
            payload["source"] = self.get_vendor_id_from_encodeportal(rec.vendor_id)
        else:
            payload["source"] = "michael-snyder"
    
        payload["treatments"] = self.post_treatments(rec.treatment_ids)
   
        if patch:  
            upstream_id = self.patch(payload, rec.upstream_identifier)
        else:
            upstream_id = self.post(payload=payload, dcc_profile="biosample", pulsar_model=models.Biosample, pulsar_rec_id=rec_id)
        return upstream_id

    def post_library(self, rec_id, patch=False, exp_type="nonSC"):
        """
        This method will check whether the biosample associated to this library is submitted. If it
        isn't, it will first submit the biosample. 
        """
        rec = models.Library(rec_id)
        payload = {}
        biosample = models.Biosample(rec.biosample_id)
        # If this Library record is a SingleCellSorting.library_prototype, then the Biosample it will
        # be linked to is the SingleCellSorting.sorting_biosample.
        payload["biosample"] = biosample.upstream_identifier
        payload["documents"] = self.post_documents(rec.document_ids)
        if exp_type in ['scAtac-multiome', 'scAtac-nonMultiome', "snRNA-multiome", "snRNA-nonMultiome"]:
            payload["construction_platform"] = "/platforms/NTR:0000452/"
        
        fragmentation_method_id = rec.library_fragmentation_method_id
        if fragmentation_method_id:
            fragmentation_method = models.LibraryFragmentationMethod(fragmentation_method_id)
            payload["fragmentation_methods"] = [fragmentation_method.name]
        if rec.lot_identifier:
            payload["lot_id"] = rec.lot_identifier
        payload["nucleic_acid_term_name"] = models.NucleicAcidTerm(rec.nucleic_acid_term_id).name
        
        if rec.vendor_product_identifier:
            payload["product_id"] = rec.vendor_product_identifier

        payload["size_range"] = rec.size_range
        #payload["strand_specificity"] = bool(rec.strand_specific)
        #if not rec.strand_specific:
        if rec.strand_specific or exp_type in ["snRNA-multiome", "snRNA-nonMultiome"]:
            payload["strand_specificity"] = "reverse"

        if rec.vendor_id:
            payload["source"] = self.get_vendor_id_from_encodeportal(rec.vendor_id)
        else:
            payload["source"] = "michael-snyder"
        ssc_id = rec.single_cell_sorting_id
        if ssc_id:
           barcode_details = self.get_barcode_details_for_ssc(ssc_id=ssc_id)
           payload["barcode_details"] = barcode_details

        # Submit payload
        if patch:  
            upstream_id = self.patch(payload, rec.upstream_identifier)
        else:
            upstream_id = self.post(payload=payload, dcc_profile="library", pulsar_model=models.Library, pulsar_rec_id=rec_id)
        return upstream_id

    def post_replicate(self, pulsar_library_id, dcc_exp_id, patch=False):
        """
        Submits a replicate record, linked to the specified library and experiment. 
        First, replicates on the experiment will be searched to see if a replicate already exists 
        for a specifc biosample and library combination, and if so then that repicate's JSON from 
        the ENCODE Portal is returned.

        If the associated experiment is ChIP-seq, and isn't a control experpiment, then the 
        replicate will be submitted with a link to antibody ENCAB728YTO (AB-9 in Pulsar), which is 
        the GFP-specific antibody used to pull down GFP-tagged TFs. 

        Args:
            pulsar_library_id: `int`. The ID of a Library record in Pulsar.
            dcc_exp_id: `int`. The ID of the experiment record on the Portal to link the replicate to.

        Returns:
            `str`: The replicate.uuid property value of the record on the ENCODE Portal.
        """
        # Required fields to submit to a replicate are:
        #  -biological_replicate_number
        #  -experiment
        #  -technical_replicate_number

        #dcc_lib = self.ENC_CONN.get(ignore404=False, rec_ids=dcc_library_id)       
        print(">>> In dcc_submit.post_replicate()")
        payload = {}
        lib = models.Library(pulsar_library_id)
        payload["library"] = lib.upstream_identifier
        biosample_id = lib.biosample_id
        biosample = models.Biosample(biosample_id)
        payload["aliases"] = [biosample.upstream_identifier + "-" + lib.upstream_identifier]
        
        payload["experiment"] = dcc_exp_id
        dcc_exp = self.ENC_CONN.get(rec_ids=dcc_exp_id)
        # Check if replicate already exists for this library
        exp_reps_instance = encode_utils.replicate.ExpReplicates(self.ENC_CONN, dcc_exp_id)
        rep_json = exp_reps_instance.get_rep(biosample_accession=biosample.upstream_identifier, library_accession=lib.upstream_identifier)
        brn = ""
        trn = ""
        if rep_json and not patch:
            return rep_json["uuid"]
        elif rep_json:
            brn = rep_json["biological_replicate_number"]
            trn = rep_json["technical_replicate_number"]
        else:
           # Then there isn't a replicate yet for this library, and maybe not even the biosample.
            if not biosample.upstream_identifier in exp_reps_instance.rep_hash:
                brn = exp_reps_instance.suggest_brn()
                trn = 1
            else:
                brn = exp_reps_instance.rep_hash[biosample.upstream_identifier]["brn"]
                trn = exp_reps_instance.suggest_trn(biosample.upstream_identifier)
        if lib.antibody_id:
            antibody = models.Antibody(lib.antibody_id) 
            if antibody.upstream_identifier:
                payload["antibody"] = antibody.upstream_identifier
        elif dcc_exp["assay_term_name"] == "ChIP-seq":
            # Only add antibody if not replicate on control experiment 
            if not dcc_exp.get("control_type"):
                payload["antibody"] = "ENCAB728YTO" #AB-9 in Pulsar
        #payload["aliases"] = 
        # Set biological_replicate_number and technical_replicate_number. For ssATAC-seq experiments,
        # these two attributes don't really make sense, but they are required to submit, so ...
            
        payload["biological_replicate_number"] = brn
        payload["technical_replicate_number"] = trn
        # Submit payload
        if patch:  
            upstream_id = self.patch(payload, rep_json["uuid"])
        else:
            # POST to ENCODE Portal. Don't use post() method defined here that is a wrapper over 
            # `encode_utils.connection.Connection.post`, since the wrapper only works if the record we
            # are POSTING has a corresponding record type on the Portal. Since Pulsar doesn't have a 
            # corresponding replicate model, we can't use the wrapper method. 
            payload[self.ENC_CONN.PROFILE_KEY] = "replicate"
            res_json = self.ENC_CONN.post(payload=payload)
            upstream_id = euu.get_record_id(res_json)
        return upstream_id

    def upload_fastq_files(self, dcc_biosample_id):
        pulsar_bio = models.Biosample(upstream_id=dcc_biosample_id)
        
        

    def post_fastq_file(self, pulsar_sres_id, read_num, enc_replicate_id, exp_type="nonSC", lane_num=None, patch=False):
        """
        Creates a file record on the ENCODE Portal. Checks the SequencingResult in Pulsar to see 
        where the file is stored. If stored in DNAnexus, the file will be downloaded locally into 
        the directory given by ``pulsarpy_to_encodedcc.FASTQ_FOLDER`` (the download folder will be
        checked first to see if the file was previously downloaded before attempting to download.
        
        After the file object is created on the ENCODE Portal, it's accession will be stored as the
        upstream identifier in the Pulsar SequencingResult record for the given read.  Thus, if a 
        file object was creatd for a R1 FASTQ file, then the `SequencingResult.read1_upstream_identifier`
        attribute is updated. If instead a file object was created for a R2 FASTQ file, then the
        `SequencingResult.read2_upstream_identifier`` attribute is updated. 

        Some rather complex logic is used to determine the control FASTQ files when submitting an 
        experimental replicate's FASTQ file. If the Biosample associated with the SequencingResult
        is part of a ChipseqExperiment, then the control biosamples consist of the paired input(s) and
        the wild type input, which in Pulsar are given the attribute names 
        `ChipseqExperiment.control_replicates` and `ChipseqExperiment.wild_type_control`. A non-control
        file object on the ENCODE Portal needs to have the ``controlled_by`` property set, which
        points to one or more control FASTQ file accessions on the ENCODE Portal. 
        We normally submit them by matching read numbers, so if the file object we are creating is 
        for a R1 FASTQ file, then all the controlled_by accessions are also R1 FASTQ files. The challenge
        is in knowing which SequencingResult set to use for control FASTQ files. Since a Biosample can
        have multiple Libraries, which can have multiple SequencingRequests, which can have multiple
        SequencingRuns, there can be many sets of SequencingResults. However, since in most cases
        there will only be one of each, the approach taken here is to use the SequencingResults of
        the latest SequencingRun of the latest SequencingRequest. Once this simplicity fails to hold,
        an updated approach will need to be taken. 

        If you have alreay created the file record on the Portal and for some reason the FASTQs didn't
        upload, you can try to reupload the FASTQs by calling this method with patch equal to False. 

        Args:
            pulsar_sres_id: A SequencingResult record in Pulsar. 
            read_num: `int`. being either 1 or 2. Use 1 for the forwrard reads FASTQ file, and 2
                for the reverse reads FASTQ file. A SequencingResult in Pulsar stores the location
                of both files (if paired-end sequening).
            end_replicate_id: `str`. The identifier of the DCC replicate record that the file record
                is to be associated with.

        Returns:
            `dict`. The response from the encode-utils POST or PATCH operation.
        """
        sres = models.SequencingResult(pulsar_sres_id)
        lib = models.Library(sres.library_id)
        bio = models.Biosample(lib.biosample_id)
        srun = models.SequencingRun(sres.sequencing_run_id)
        sreq = models.SequencingRequest(srun.sequencing_request_id)
        platform = models.SequencingPlatform(sreq.sequencing_platform_id)
        payload = {}
        payload["aliases"] = []
        # The Pulsar SequencingPlatform must already have the upstream_identifier attribute set.
        payload["platform"] = platform.upstream_identifier
        payload["file_format"] = "fastq"
        # set flowcell_details
        flowcell_details = {}
        flowcell_details["barcode"] = lib.get_barcode_sequence()
        flowcell_details["lane"] = str(srun.lane)
        payload["flowcell_details"] = [flowcell_details]
        payload["replicate"] = enc_replicate_id
        #if sreq.paired_end:
        #    payload["run_type"] = "paired-ended"
        #else:
        #    payload["run_type"] = "single-ended"

        if read_num == 1:
            payload["run_type"] = "paired-ended"
            payload["paired_end"] = "1"
            payload["output_type"] = "reads"
            file_uri = sres.read1_uri
            upstream_id = sres.read1_upstream_identifier
            read_count = sres.read1_count
            
            if exp_type == "nonSC":
                payload["read_length"] = 100
            elif exp_type in ["snRNA-multiome", "snRNA-nonMultiome"]:
                print("single nucleus RNA-seq, multiome")
                payload["read_length"] = 28
                
                # required by DCC: see github issue 739
                if exp_type == "snRNA-multiome":
                    payload['read_structure'] = [{"sequence_element": "barcode", "start": 1, "end": 16},
                                                    {"sequence_element": "UMI", "start": 17, "end": 28}]

            elif exp_type in ["scAtac-multiome", "scAtac-nonMultiome"]:
                print("single cell Atac-seq read 1 submission.")
                payload["read_length"] = 50

        elif read_num == 2:
            if exp_type == "nonSC":
                payload["run_type"] = "paired-ended"
                payload["paired_end"] = "2"
                payload["read_length"] = 100
                payload["output_type"] = "reads"
                
                file_uri = sres.read2_uri
                upstream_id = sres.read2_upstream_identifier
                read_count = sres.read2_count
                # Need to set paired_with key in the payload. In this case, 
                # it is expected that R1 has been already submitted.
                if not sres.read1_upstream_identifier:
                    raise Exception("Can't set paired_with for SequencingResult {} since R1 doesn't have an upstream set.".format(sres.id))
                payload["paired_with"] = sres.read1_upstream_identifier
            elif exp_type in ["snRNA-multiome", "snRNA-nonMultiome"]:
                print("single nucleus RNA-seq read 2")
                # required by DCC: see github issue 739
                if exp_type == "snRNA-nonMultiome":
                    payload['read_structure'] = [{"sequence_element": "barcode", "start": 1, "end": 10}]
                
                payload["run_type"] = "paired-ended"
                payload["paired_end"] = "2"
                payload["read_length"] = 90
                payload["output_type"] = "reads"
                payload["paired_with"] = sres.read1_upstream_identifier
                file_uri = sres.read2_uri
                upstream_id = sres.read2_upstream_identifier
                read_count = sres.read2_count
            elif exp_type == "scAtac-multiome":
                print("single cell Atac-seq index read (R2) submission.")
                
                # no need for index read to have payload["run_type"] or payload["paired_end"]
                payload["read_length"] = 24
                payload["output_type"] = "index reads"
                
                # read_structure required by DCC
                # For scAtac-nonMultiome, it is attached to index read (read 2)
                payload['read_structure'] = [{"sequence_element": "barcode", "start": 1, "end": 24}]

                if not sres.read1_upstream_identifier:
                    raise Exception("Can't set index_of for SequencingResult {} since read 1 (R1) doesn't have an upstream set.".format(sres.id))
                if not sres.read2_upstream_identifier:
                    raise Exception("Can't set index_of for SequencingResult {} since read 2 (R3) doesn't have an upstream set.".format(sres.id))
                payload["index_of"] = [sres.read1_upstream_identifier, sres.read2_upstream_identifier]
            elif exp_type == "scAtac-nonMultiome":
                print("single cell Atac-seq read 2")
                
                payload["read_length"] = 16
                payload["output_type"] = "index reads"

                # read_structure required by DCC
                payload['read_structure'] = [{"sequence_element": "barcode", "start": 1, "end": 16}]
                
                if not sres.read1_upstream_identifier:
                    raise Exception("Can't set index_of for SequencingResult {} since read 1 (R1) doesn't have an upstream set.".format(sres.id))
                if not sres.read2_upstream_identifier:
                    raise Exception("Can't set index_of for SequencingResult {} since read 2 (R3) doesn't have an upstream set.".format(sres.id))
                payload["index_of"] = [sres.read1_upstream_identifier, sres.read2_upstream_identifier]
        elif read_num == 3:
            if exp_type in ["nonSC", "snRNA-multiome", "snRNA-nonMultiome"]:
                raise Exception("Unimplemented for nonSC or snRNA experiment read 3's submission.")
            elif exp_type in ["scAtac-multiome", "scAtac-nonMultiome"]:
                print("single cell Atac-seq read 2 (R3) submission.")

                payload["run_type"] = "paired-ended"
                payload["paired_end"] = "2"
                payload["read_length"] = 50
                payload["output_type"] = "reads"
                file_uri = sres.read2_uri
                upstream_id = sres.read2_upstream_identifier
                read_count = sres.read2_count
                
                # Need to set paired_with key in the payload. In this case, 
                # it is expected that R1 has been already submitted.
                if not sres.read1_upstream_identifier:
                    raise Exception("Can't set paired_with for SequencingResult {} since R1 doesn't have an upstream set.".format(sres.id))
                payload["paired_with"] = sres.read1_upstream_identifier

        data_storage = models.DataStorage(srun.data_storage_id)
        data_storage_provider = models.DataStorageProvider(data_storage.data_storage_provider_id)

        # Initialize file_path to be empty string.
        file_path = ""

        # extract barcode
        barcode = models.Barcode.find_by({'id': str(lib['barcode_id'])})
        if barcode:
            barcode_sequence = barcode['sequence'].replace(',', '_')
        else:
            pbc = models.PairedBarcode(lib.paired_barcode_id)
            if pbc:
                barcode_sequence = models.Barcode(pbc.index1_id)["sequence"] + "_" + models.Barcode(pbc.index2_id)["sequence"]

        assert(barcode_sequence), "Barcode should never be empty for a library to submit."
            
        # Check file existence
        # Assumption: there is only one fastq file for each read.
        files = []
        sreq_id = lib['sequencing_request_ids']
        if len(sreq_id) == 1:
            #file_dir = '/oak/stanford/scg/prj_ENCODE/Staging/' + str(sreq_id[0])
            file_dir = '/oak/stanford/scg/prj_ENCODE/SREQ/' + str(sreq_id[0])
            if os.path.isdir(file_dir):
                old_dir = os.getcwd()
                os.chdir(file_dir)

                # snRNA-seq and scAtac-seq
                if lane_num:
                    files = glob.glob("**/*" + lane_num + "*" + barcode_sequence + "*R" + str(read_num) + "*.fastq.gz", recursive=True)
                
                # snRNA-seq and scAtac-seq
                if lane_num and not files:
                    files = glob.glob("**/*" + barcode_sequence + "*" + lane_num + "*R" + str(read_num) + "*.fastq.gz", recursive=True)
                
                # snRNA-seq and scAtac-seq
                if lane_num and not files:
                    files = glob.glob("**/*" + barcode_sequence.replace('_', '-') + "*" + lane_num + "*R" + str(read_num) + "*.fastq.gz", recursive=True)
                
                # ChIP and bulkAtac-seq
                if not files:
                    files = glob.glob("**/*" + barcode_sequence + "*R" + str(read_num) + "*.fastq.gz", recursive=True)
        else:
            raise Exception(f'More than one sequencing request exists for library {lib["id"]}')

        # now look up the file's full path
        if len(files) == 1:
            file_path = os.path.abspath(files[0])
            file_ref = "SCG${}".format(read_num)
            #payload["aliases"].append(file_ref)
            payload["aliases"].append(os.path.basename(file_path))
            os.chdir(old_dir)
        elif len(files) > 1:
            raise Exception(f'More than one fastq files exist for sequencing result {pulsar_sres_id} read {read_num}.')
            return None
        else:
            if not file_uri:
                raise NoFastqFile("SequencingResult '{}' for R{} does not have a FASTQ file path set.".format(pulsar_sres_id, read_num))
            elif not upstream_id and patch:
                raise Exception("Can't PATCH file object on the Portal when the SequencingResult {} for read {} doesn't have an upstream ID set.".format(pulsar_sres_id, read_num))

            #raise Exception(f'Fastq file for sequencing result {pulsar_sres_id} read {read_num} does not exist.')
            if data_storage_provider.name == "DNAnexus":
                dx_file = dxpy.DXFile(dxid=file_uri, project=data_storage.project_identifier)
                file_path = os.path.join(FASTQ_FOLDER, dx_file.name)
                # Check if file exists and is non-empty in download directory before attempting to download.
                # file_uri = data_storage.project_identifier + ':' + file_uri
                if not patch:
                    if not os.path.exists(file_path) or not os.path.getsize(file_path):
                        # Download file.
                        dxpy.download_dxfile(dxid=file_uri, project=data_storage.project_identifier, filename=file_path, show_progress=True)
                file_ref = "dnanexus${}".format(file_uri)
                payload["aliases"].append(file_ref)
                payload["aliases"].append(dx_file.name)
            elif data_storage_provider == "AWS S3 Bucket":
                file_path = file_uri
                # md5sum key added to payload by encode-utils.
                payload["aliases"].append(file_uri) # i.e. s3://bucket-name/key
                payload["aliases"].append(os.path.basename(file_uri))
        
        if not file_path:
            raise Exception("Could not locate FASTQ file for SequencingResult {}; read number {}.".format(pulsar_sres_id, read_num))
        payload["submitted_file_name"] = file_path 

        dcc_rep = self.ENC_CONN.get(rec_ids=enc_replicate_id, ignore404=False)
        dcc_exp = dcc_rep["experiment"]
        payload["dataset"] = dcc_exp["accession"]
        dcc_exp_accession = dcc_exp["accession"]
        dcc_exp_type = dcc_exp["assay_term_name"] #i.e. ChIP-seq
        if dcc_exp_type == "ChIP-seq":
            if not bio.control and not bio.wild_type:
                controlled_by = self.get_chipseq_controlled_by(pulsar_biosample=bio, read_num=read_num, dcc_exp_id=dcc_exp_accession)
                if controlled_by:
                    # Will be empty if this is a already a control SequencingResult file.
                    payload["controlled_by"] = controlled_by
                else:
                    pass
                    #raise Exception("No controlled_by could be found for SequencingResult {} read number {}.".format(pulsar_sres_id, read_num))
                    # Instead of raise an Exception, let it slide. Pulsar users aren't setting the 
                    # boolean fields control and wild_type_control as they should be so it's not reliable. 
        else:
            if dcc_exp_type not in ["ATAC-seq", "single-nucleus ATAC-seq", "single-cell RNA sequencing assay"]:
                raise Exception("There isn't yet support to set controlled_by for experiments of type {}.".format(dcc_exp_type))

        # POST to ENCODE Portal. Don't use post() method defined here that is a wrapper over 
        # `encode_utils.connection.Connection.post`, since the wrapper only works if the record we
        # are POSTING has a corresponding record type on the Portal. Since Pulsar doesn't have a 
        # corresponding file model, we can't use the wrapper method. So we'll have to manually 
        # set the upstream identifier in the attribute `SequencingResult.read1_upstream_identifier`
        # or `SequencingResult.read2_upstream_identifier`.
        #patch = True
        if patch:
            upstream_id = self.patch(payload=payload, upstream_id=upstream_id)
        else:
            payload[self.ENC_CONN.PROFILE_KEY] = "file"
            res_json = self.ENC_CONN.post(payload=payload)
            upstream_id = res_json["accession"]
            if read_num == 1:
                sres.patch({"read1_upstream_identifier": upstream_id})
            else:
                sres.patch({"read2_upstream_identifier": upstream_id})
        return upstream_id
        
    def get_chipseq_controlled_by(self, pulsar_biosample, read_num, dcc_exp_id):
        """
        Given a p
        Returns:
            `list`: The upstream identifiers for the control file objects on the ENCODE Portal.
        """
        if pulsar_biosample.control or pulsar_biosample.wild_type:
            return []
        bio_id = pulsar_biosample.id
        controlled_by = []
        chipseq_experiment = pulsarpy.models.ChipseqExperiment(upstream=dcc_exp_id)
        # First add pooled input. Normally one control but could be more. 
        ctl_map = chipseq_experiment.paired_input_control_map() 
        if bio_id in ctl_map:
            for ctl_id in ctl_map[bio_id]:
                ctl = pulsarpy.models.Biosample(ctl_id)
                lib = ctl.get_latest_library()
                controlled_by.extend(self.get_all_seqresult_fastq_file_accessions(lib)[read_num])
        # Next add WT input
        wt_input_id = chipseq_experiment.wild_type_control_id
        if wt_input_id:
            wt_input = pulsarpy.models.Biosample(wt_input_id)
            lib = wt_input.get_latest_library()
            controlled_by.extend(self.get_all_seqresult_fastq_file_accessions(lib)[read_num])
        return controlled_by

    def get_all_seqresult_fastq_file_accessions(self, pulsar_lib):
        res = {1: [], 2: []}
        sequencing_result_ids = pulsar_lib.sequencing_result_ids
        for sres_id in sequencing_result_ids:
            sres = pulsarpy.models.SequencingResult(sres_id)
            r1_accession = sres.get_upstream_identifier(read_num=1)
            if not r1_accession:
                raise Exception("Upstream identifier not set for SequencingResult {}, read number {}.".format(sres_id, 1))
            res[1].append(r1_accession)
            r2_accession = sres.get_upstream_identifier(read_num=2)
            if not r2_accession:
                raise Exception("Upstream identifier not set for SequencingResult {}, read number {}.".format(sres_id, 2))
            res[2].append(r2_accession)
        return res

    def get_barcode_details_for_ssc(self, ssc_id):
        """
        This purpose of this method is to provide a value to the library.barcode_details property
        of the Library profile on the ENCODE Portal. That property taks an array of objects whose
        properties are the 'barcode', 'plate_id', and 'plate_location'. 

        Args:
            ssc_id: The Pulsar ID for a SingleCellSorting record.
        """
        ssc = models.SingleCellSorting(ssc_id)
        lib_prototype_id = ssc.library_prototype_id
        lib = models.Library(lib_prototype_id)
        paired_end = lib.paired_end
        plate_ids = ssc.plate_ids
        plates = [models.Plate(p) for p in plate_ids]
        results = []
        for p in plates:
            for well_id in p.well_ids:
                well = models.Well(well_id)
                details = {}
                details["plate_id"] = p.name
                details["plate_location"] = well.name
                well_biosample = models.Biosample(well.biosample_id)
                lib_id = well_biosample.library_ids[-1]
                # Doesn't make sense to have more than one library for single cell experiments. 
                lib = models.Library(lib_id)
                if not paired_end:
                    barcode_id = lib.barcode_id
                    barcode = models.Barcode(barcode_id).sequence
                else:
                    pbc_id = lib.paired_barcode_id
                    pbc = models.PairedBarcode(pbc_id)
                    barcode = models.Barcode(pbc.index1_id)["sequence"] + "-" + models.Barcode(pbc.index2_id)["sequence"]
                details["barcode"] = barcode
                results.append(details)
        return results

    def post_single_cell_sorting(self, rec_id, patch=False):
        rec = models.SingleCellSorting(rec_id)
        sorting_biosample_id = rec.sorting_biosample_id
        sorting_biosample = models.Biosample(sorting_biosample_id)
        payload = {}
        # Set the explicitly required properties first:
        payload["assay_term_name"] = "single-cell ATAC-seq"
        payload["experiment_classification"] = ["functional genomics"]
        # And now the rest
        payload["biosample_ontology"] = sorting_biosample["biosample_term_name"]["name"]
        desc = rec.description.strip()
        if desc:
            payload["description"] = desc
        payload["documents"] = self.post_documents(rec.document_ids)
        exp_upstream = self.post(payload=payload, dcc_profile="experiment", pulsar_model=models.SingleCellSorting, pulsar_rec_id=rec_id)

        # Submit biosample
        self.post_biosample(rec_id=sorting_biosample.id, patch=patch)
        # Submit library_prototype (which is linked to sorting_biosample when it is created)
        library_prototype_id = rec.library_prototype_id
        library_upstream = self.post_library(rec_id=library_prototype_id, patch=patch)
        # Submit replicate.
        # The experiment will be determined via inspection of associated biosample. 
        replicate_upstream = self.post_replicate(library_upstream=library_upstream, patch=patch)
        # A SingleCellSorting has many SequencingRequests through the Plates association.
        sreq_ids = rec.sequencing_request_ids 
        sreqs = [models.SequencingRequest(s) for s in sreq_ids]
        for sreq in sreqs:
            srun_ids = sreq.sequencing_run_ids
            sruns = [models.SequencingRun(s) for s in srun_ids]
            for run in sruns:
                storage_loc_id = run.storage_location_id
                sres_ids = run.sequencing_result_ids
                # Submit a file record

if __name__ == "__main__":
    s = Submit(dcc_mode="v79x0-test-master.demo.encodedcc.org") 
    s.post_chipseq_exp(rec_id=164, patch=False)
    #s.post_chipseq_control_experiments(164)
   
