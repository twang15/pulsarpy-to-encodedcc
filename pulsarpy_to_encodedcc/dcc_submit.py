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

from pulsarpy import models
import pulsarpy.utils
import encode_utils as eu
import encode_utils.connection as euc
import pdb


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
        return txt.replace("/","-")
                  

    #def patch(self, payload, raise_403=True, extend_array_values=False):
         # Just use self.ENC_CONN.patch() directly.


    def post(self, payload, dcc_profile, pulsar_model, pulsar_rec_id):
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
            upstream identifier if the record already exists. The identifier of either the new record 
            or existing record on the ENCODE Portal is determined to be the 'accession' if that property is present, 
            otherwise it's the first alias if the 'aliases' property is present, otherwise its 
            the value of 'uuid' property.
        """
        payload[self.ENC_CONN.PROFILE_KEY] = dcc_profile
        rec = pulsar_model(pulsar_rec_id)
        upstream = rec.upstream_identifier
        if upstream:
            # No need to POST. 
            return upstream
        aliases = payload.get("aliases", [])
        aliases.append(rec.abbrev_id())
        # Add value of 'name' property as an alias, if this property exists for the given model.
        try:
            name = self.sanitize_prop_val(rec.name)
            if name:
                aliases.append(name)
        except KeyError:
            pass
        aliases = list(set(aliases))
        payload["aliases"] = aliases
    
        # `dict`. The POST response if the record didn't yet exist on the ENCODE Portal, or the
        # record JSON itself if it does already exist. Note that the dict. will be empty if the connection
        # object to the ENCODE Portal has the dry-run feature turned on.
        response_json = self.ENC_CONN.post(payload)
        if "accession" in response_json:
            upstream = response_json["accession"]
        elif "aliases" in response_json:
            upstream = response_json["aliases"][0]
        else:
            upstream = response_json["uuid"]
        # Set value of the Pulsar record's upstream_identifier, but only if we are in prod mode since
        # we don't want to set it to an upstream identifiers from any of the ENOCODE Portal test servers. 
        if self.dcc_mode == eu.DCC_PROD_MODE:
            print("Setting the Pulsar record's upstream_identifier attribute to '{}'.".format(upstream))
            pulsar_rec = pulsar_model(pulsar_rec_id)
            pulsar_rec.patch(payload={"upstream_identifier": upstream})
            print("upstream_identifier attribute set successfully.")
        return upstream
    
    def post_crispr_modification(self, rec_id, patch=False):
        rec = models.CrisprModification(rec_id)
        # CrisprConstruct(s)
        cc_ids = rec.crispr_construct_ids
        ccs = [models.CrisprConstruct(i) for i in cc_ids]
        # DonorConstruct
        dc = models.DonorConstruct(rec.donor_construct_id)

        payload = {}
        payload["category"] = rec.category # Required
        payload["description"] = rec.description
        payload["documents"] = self.post_documents(rec.document_ids)

        guide_seqs = list(c.guide_sequence for c in ccs)
        payload["guide_rna_sequences"] = guide_seqs

        if rec.category in ["insertion", "replacement"]:
            payload["introduced_sequence"] = dc.insert_sequence

        payload["method"] = "CRISPR"       # Required
        payload["modified_site_by_target_id"] = dc.target_id
        payload["purpose"] = rec.purpose   # Required

        # Note that CrisprConstruct can also has_many construct_tags. Those are not part of the donor
        # insert though. 
        construct_tags = [models.ConstructTag(i) for i in dc.construct_tag_ids]
        construct_tag_names = [x.name for x in construct_tags]
        introduced_tags = []
        for tag in construct_tag_names:
            if tag.startswith("eGFP"):
                # Pulsar has two eGFP tags that differ by the linker sequence:
                #    1) eGFP (MH170480)
                #    2) eGFP (MH170481)
                # The Portal, however, only has eGFP and it makes most sense to submit this as 
                # simply eGFP and mention the linker used elsewhere. 
                tag = "eGFP"
            introduced_tags.append({"name": tag, "location": "C-terminal"})
        payload["introduced_tags"] = introduced_tags
        reagents = []
        for i in [*ccs,dc]:
            addgene_id = getattr(i, "addgene_id")
            if addgene_id:
                uri = "http://www.addgene.org/" + addgene_id
                reagents.append({
                    "source": "addgene",
                    "identifier": addgene_id,
                    "url": uri
                })
        payload["reagents"] = reagents
        # ex: ENCGM094ZOS

        if patch: 
            res = self.ENC_CONN.patch(payload=payload, extend_array_values=self.extend_arrays)
        else:
            res = self.post(payload=payload, dcc_profile="genetic_modification", pulsar_model=models.CrisprModification, pulsar_rec_id=rec_id)
            return res
    
    def post_document(self, rec_id, patch=False):
        rec = models.Document(rec_id)
        payload = {}
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
            res = self.ENC_CONN.patch(payload=payload, extend_array_values=self.extend_arrays)
        else:
            res = self.post(payload=payload, dcc_profile="document", pulsar_model=models.Document, pulsar_rec_id=rec_id)
        return res

    def post_documents(self, rec_ids, patch=False):
        upstreams = []
        for i in rec_ids:
            upstreams.append(self.post_document(rec_id=i, patch=patch))
        return upstreams

    def post_treatments(self, rec_ids, patch=False):
        upstreams = []
        for i in rec_ids:
            upstreams.append(self.post_treatment(rec_id=i, patch=patch))
        return upstream

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
            res = self.ENC_CONN.patch(payload=payload, extend_array_values=self.extend_arrays)
        else:
            res = self.post(payload=payload, dcc_profile="treatment", pulsar_model=models.Treatment, pulsar_rec_id=rec_id)
        return res

    
    def post_vendor(self, rec_id, patch=False):
        """
        Vendors must be registered directly by the DCC personel. 
        """
        raise Exception("Vendors must be registered directly by the DCC personel.")

    def post_biosample(self, rec_id, patch=False):
        rec = models.Biosample(rec_id)
        # The alias lab prefixes will be set in the encode_utils package if the DCC_LAB environment
        # variable is set.
        payload = {}
        btn = models.BiosampleTermName(rec.biosample_term_name_id)
        payload["biosample_term_name"] = btn.name.lower() #Portal requires lower-case
        payload["biosample_term_id"] = btn.accession

        bty = models.BiosampleType(rec.biosample_type_id)
        payload["biosample_type"] = bty.name

        date_biosample_taken = rec.date_biosample_taken
        if date_biosample_taken:
            if bty.name == "tissue":
                payload["date_obtained"] = date_biosample_taken
            else:
                payload["culture_harvest_date"] = date_biosample_taken

        desc = rec.description
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
        if nih_cert:
            payload["nih_institutional_certification"] = nih_cert

        payload["organism"] = "human"

        passage_number = rec.passage_number
        if passage_number:
            payload["passage_number"] = passage_number

        starting_amount = rec.starting_amount
        if starting_amount:
            payload["starting_amount"] = starting_amount
            payload["starting_amount_units"] = models.Unit(rec.starting_amount_units_id).name

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
            if not pob_upstream:
                pob_upstream = self.post_biosample(part_of_biosample_id)
            payload["part_of"] = pob_upstream
    
        pooled_from_biosample_ids = rec.pooled_from_biosample_ids
        if pooled_from_biosample_ids:
            pooled_from_biosamples = [models.Biosample(p) for p in pooled_from_biosample_ids]
            payload["pooled_from"] = []
            for p in pooled_from_biosamples:
                p_upstream = p.get_upstream() 
                if not p_upstream:
                    p_upstream = self.post_biosample(p.id)
                payload["pooled_from"].append(p_upstream)
    
        vendor = models.Vendor(rec.vendor_id)
        vendor_upstream = vendor.get_upstream() 
        if not vendor_upstream:
            raise Exception("Biosample '{}' has a vendor without an upstream set: Vendors are requied to be registered by the DCC personel, and Pulsar needs to have the Vendor record's '{}' attribute set.".format(rec_id, models.Model.UPSTREAM_ATTR))
        payload["source"] = vendor_upstream
    
        payload["treatments"] = self.post_treatments(rec.treatment_ids)
   
        if patch:  
            res = self.ENC_CONN.patch(payload=payload, extend_array_values=self.extend_arrays)
        else:
            res = self.post(payload=payload, dcc_profile="biosample", pulsar_model=models.Biosample, pulsar_rec_id=rec_id)
        return res

    def post_library(self, rec_id, patch=False):
        """
        This method will check whether the biosample associated to this library is submitted. If it
        isn't, it will first submit the biosample. 

        Args:
        """
        rec = models.Library(rec_id)
        payload = {}
        biosample = models.Biosample(rec.biosample_id)
        # If this Library record is a SingleCellSorting.library_prototype, then the Biosample it will
        # be linked to is the SingleCellSorting.sorting_biosample.
        payload["biosample"] = self.post_biosample(rec_id=rec.biosample_id, patch=False)
        payload["documents"] = self.post_documents(rec.document_ids)
        fragmentation_method_id = rec.library_fragmentation_method_id
        if fragmentation_method_id:
            fragmentation_method = models.LibraryFragmentationMethod(fragmentation_method_id).name
            payload["fragmentation_method"] = fragmentation_method.name
        payload["lot_id"] = rec.lot_identifier
        payload["nucleic_acid_term_name"] = models.NucleicAcidTerm(rec.nucleic_acid_term_id).name
        payload["product_id"] = rec.vendor_product_identifier
        payload["size_range"] = rec.size_range
        payload["strand_specificity"] = rec.strand_specific
        payload["source"] = self.post_vendor(rec_id=rec.vendor_id, patch=False)
        ssc_id = rec.single_cell_sorting_id
        if ssc_id:
           barcode_details = self.get_barcode_details_for_ssc(ssc_id=ssc_id)
           payload["barcode_details"] = barcode_details

        # Submit payload
        if patch:  
            res = self.ENC_CONN.patch(payload=payload, extend_array_values=self.extend_arrays)
        else:
            res = self.post(payload=payload, dcc_profile="library", pulsar_model=models.Biosample, pulsar_rec_id=rec_id)
        return res

    def post_replicate(self, library_upstream, patch=False):
        """
        Submits a replicated record, linked to the specified library. The associated experiment
        record will be determined via the method :func:`pulsarpy.utils.get_exp_of_biosample`.

        Todo: Check what value to use for technical_replicate_number. 
        Args:
            library_upstream - The identifier of a Library record on the ENCODE Portal, which is
                 stored in Pulsar via the upstream_identifier attribute of a Library record.
        """
        # Required fields to submit to a replicate are:
        #  -biological_replicate_number
        #  -experiment
        #  -technical_replicate_number

        #dcc_lib = self.ENC_CONN.get(ignore404=False, rec_ids=dcc_library_id)       
        lib = models.Library(upstream=library_upstream)
        biosample_id = lib.biosample_id
        biosample = models.Biosample(biosample_id)
        payload = {}
        payload["antibody"] = "ENCAB728YTO" #AB-9 in Pulsar
        #payload["aliases"] = 
        # Set biological_replicate_number and technical_replicate_number. For ssATAC-seq experiments,
        # these two attributes don't really make sense, but they are required to submit, so ...
        brn = biosample.replicate_number
        if not brn:
            # Check if this is a SingleCellSorting.library_prototype - if so, replicate_number isn't
            # really meaningful in Pulsar for the linked biosample, just default it to 1:
            if lib.single_cell_sorting_id:
                brn = 1
            else:
               raise Exception("Can't submit replicate object for library '{}' since the associated biosample doesn't have the replicate_number attribute set.".format(library_upstream))
        payload["biological_replicate_number"] = brn
        payload["technical_replicate_number"] = 1
        # Figure out the experiment that this Biosample is associated to
        exp_rec = pulsarpy.utils.get_exp_of_biosample(biosample)
        if not exp_rec.upstream_identifier:
            raise Exception("Can't submit replicate when the experiment it is linked to hasn't been submitted.")
        payload["experiment"] = experiment_upstream
        payload["library"] = library_upstream
        # Submit payload
        if patch:  
            res = self.ENC_CONN.patch(payload=payload, extend_array_values=self.extend_arrays)
        else:
            res = self.post(payload=payload, dcc_profile="replicate", pulsar_model=models.Biosample, pulsar_rec_id=rec_id)
        return res
        #return payload

    def post_fastq_file(self, pulsar_sres_id, read_num, enc_replicate_id, patch=False):
        """
        Creates a file record on the ENCODE Portal. 

        Args:
            pulsar_sres_id: A SequencingResult record in Pulsar. 
            read_num: `int`. being either 1 or 2. Use 1 for the forwrard reads FASTQ file, and 2
                for the reverse reads FASTQ file. A SequencingResult in Pulsar stores the location
                of both files (if paired-end sequening).
            end_replicate_id: `str`. The identifier of the DCC replicate record that the file record
                is to be associated with.
        """
        sres = models.SequencingResult(pulsar_sres_id)
        sreq = models.SequencingRequest(sres.sequencing_request_id)
        platform = models.Platform(sreq.sequencing_platform_id)
        payload = {}
        # Need to set md5sum and file_size. 
        if read_num == 1:
            payload["paired_end"] = 1
            file_uri = sres.read1_uri
            read_count = sres.read1_count
        elif read_num == 2:
            payload["paired_end"] = 2
            file_uri = sres.read2_uri
            read_count = sres.read2_count
            # Need to set paired_with key in the payload. In this case, 
            # it is expected that R1 has been already submitted.
            payload["paired_with"] = sres.read1_upstream_identifier
        if not file_uri:
            raise NoFastqFile("SequencingResult '{}' for R{read_num} does not have a FASTQ file path set.".format(pulsar_sres_id, read_num))
        data_storage = models.DataStorage(srun.data_storage_id)
        data_storage_provider = models.DataStorageProvider(data_storage.data_storage_provider_id)
        if data_storage_provider.name == "DNAnexus":
            # Download file and set file_size and md5sum payload keys
        dx_file = dxpy.DXFile(dxid=file_uri)
        file_ref = "dnanexus:{file_uri}".format(file_uri)
        aliases = [dx_file.name, file_ref]
        srun = models.Sequencingrun(sres.sequencing_run_id)
        full_path = os.path.join(data_storage.folder_base_path, storage_path, file_uri)
        dcc_rep = self.ENC_CONN.get(rec_ids=enc_replicate_id, ignore404=False)
        payload["aliases"] = aliases
        dcc_rep = self.ENC_CONN.get(rec_ids=enc_replicate_id, ignore404=False)
        payload["file_format"] = "fastq"
        payload["output_type"] = "reads"
        # The Pulsar Platform must already have the upstream_identifier attribute set.
        payload["platform"] = platform.upstream_identifier
        if read_count:
            payload["read_count"] = read_count
        payload["replicate"] = enc_replicate_id
        payload["run_type"] = sreq.paired_end
        payload["submitted_file_name"] = file_ref
        
        
        

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
                    barcode = pbc.index1["sequence"] + "-" + pbc.index2["sequence"]
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
        payload["biosample_type"] = sorting_biosample["biosample_type"]["name"]
        payload["experiment_classification"] = "functional genomics"
        # And now the rest
        payload["biosample_term_name"] = sorting_biosample["biosample_term_name"]["name"]
        payload["biosmple_term_id"] = sorting_biosample["biosample_term_name"]["accession"]
        payload["description"] = rec.description
        payload["documents"] = self.post_documents(rec.document_ids)
        exp_upstream = self.post(payload=payload, dcc_profile="experiment", pulsar_model="SingleCellSorting", pulsar_rec_id=rec_id)

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
