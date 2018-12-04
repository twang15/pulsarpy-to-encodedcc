# -*- coding: utf-8 -*-

###Author
#Nathaniel Watson
#2017-09-18
#nathankw@stanford.edu
###

import json
import pdb
import re

import encode_utils.connection as euc
import encode_utils.utils as euu
import pulsarpy.models as models
import pulsarpy.utils

protocol_regx = re.compile(r'protocol', re.IGNORECASE)

ACCESSION_PROP = "accession"
ALIASES_PROP = "aliases"
UPSTREAM_PROP = "upstream_identifier"
UUID_PROP = "uuid"
ENC_CONN = euc.Connection("prod")
#ADMIN_USER_ID = models.User.find_by({"email": "admin@enc.com"})["id"]
ADMIN_USER_ID = 1

# Biosamples to import for Jessika:
# https://www.encodeproject.org/search/?type=Biosample&lab.title=Michael+Snyder%2C+Stanford&award.rfa=ENCODE4&biosample_type=tissue
# which sum up to the following 30 accessions:
#['ENCBS443KFH', 'ENCBS251BGN', 'ENCBS558OUC', 'ENCBS319AUC', 'ENCBS895UDJ', 'ENCBS303LEQ', 'ENCBS208HYM', 'ENCBS273ENZ', 'ENCBS704NZI', 'ENCBS444XVA', 'ENCBS924ALU', 'ENCBS892KPS', 'ENCBS729ENA', 'ENCBS268DOV', 'ENCBS157OMX', 'ENCBS441BYJ', 'ENCBS858UHJ', 'ENCBS577DQE', 'ENCBS655YSD', 'ENCBS064HXH', 'ENCBS790AKV', 'ENCBS437TFK', 'ENCBS465UNR', 'ENCBS278ABD', 'ENCBS131SFJ', 'ENCBS605RWM', 'ENCBS722FKO', 'ENCBS603AMH', 'ENCBS222RAG', 'ENCBS649EIC']

def set_name(rec):
    """
    Most of the models in Pulsar have a name attribute, and most of the time it is required.
    When backporting a record from the ENCODE Portal, we need some value to use as the record's name,
    and records in the Portal don't have a name prop, so we need to use some other propery value.

    The approach taken here is to use the record's first alias in the 'aliases' property as a name
    if that is present.  The alias prefix (the part that includes the ':' and preceeding text) is
    stripped off. Otherwise, use the 'acccession' property if that is present.  Otherwise,
    use the record's "uuid" property. If for some reason none of these properties are set, an
    Exception is raised.

    Args:
        rec_id: `str`. An identifier for some record on the ENCODE Portal.

    Returns:
        `str` designating the value to use for a Pulsar LIMS record's 'name' attribute.

    Raises:
        `Exception`: A value for name couldn't be set.

    """
    if ALIASES_PROP in rec:
        return rec[ALIASES_PROP][0].split(":")[1]
    elif ACCESSION_PROP in rec:
        return rec[ACCESSION_PROP]
    elif rec[UUID_PROP] in rec:
        return rec[UUID_PROP]
    raise Exception("Can't set name for record {}".format(json.dumps(rec, indent=4)))

def target(rec_id):
    """
    Backports a source record belonging to https://www.encodeproject.org/profiles/source/json into
    the Pulsar model called Target.

    Identifying properties on the Portal are "aliases", "label-organism.name", and "uuid".

    Args:
        rec_id: `str`. An identifier for a source record on the ENCODE Portal.

    Returns:
        `dict`: The JSON representation of the existing Target if it already exists in
        in Pulsar, otherwise the POST response.
    """
    rec = ENC_CONN.get("targets/" + rec_id, ignore404=False)
    aliases = rec[ALIASES_PROP]
    label = rec["label"]
    label_org_name = label + "-" + rec["organism"]["name"]
    #check if upstream exists already in Pulsar:
    pulsar_rec = models.Target.find_by({UPSTREAM_PROP: [*aliases, rec[UUID_PROP], label_org_name, rec["@id"]]})
    if pulsar_rec:
        return pulsar_rec
    print("Importing target '{}'".format(rec_id))
    payload = {}
    payload[UPSTREAM_PROP] = label_org_name
    payload["name"] = label
    xrefs = rec["dbxref"]
    for ref in xrefs:
      prefix, ref = ref.split(":")
      if prefix == "ENSEMBL":
        payload["ensembl"] = ref
      elif prefix == "UniProtKB":
        payload["uniprotkb"] = ref
      elif prefix == "RefSeq":
        payload["refseq"] = ref
    return models.Target.post(payload)

def biosample(rec_id, patch=False):
    """
    Backports a biosample record belonging to
    https://www.encodeproject.org/profiles/biosample.json into the Pulsar model called
    Biosample.

    Identifying properties on the Portal are "accession", "aliases", and "uuid".
    Portal's required props are: award, biosample_term_id, biosample_term_name, biosample_type lab,
    organism, source

    Args:
        rec_id: `str`. An identifier for a document record on the ENCODE Portal.

    Returns:
        `dict`: The JSON representation of the existing Biosample if it already exists
        in Pulsar, otherwise the POST response.
    """
    dcc_rec = ENC_CONN.get(rec_id, ignore404=False)
    aliases = dcc_rec[ALIASES_PROP]
    accession = dcc_rec[ACCESSION_PROP]
    # Check if upstream exists already in Pulsar:
    pulsar_rec = models.Biosample.find_by({UPSTREAM_PROP: [*aliases, accession, dcc_rec[UUID_PROP], dcc_rec["@id"]]})
    # If the record was found, then don't post it, but may still need to backport any additional
    # "has" and "has_many" relations that may be new. 
    payload = {}
    payload[UPSTREAM_PROP] = accession
    payload["name"] = set_name(dcc_rec)
    btn = dcc_rec["biosample_term_name"]
    bti = dcc_rec["biosample_term_id"]
    pulsar_btn_rec = biosample_term_name(biosample_term_name=btn, biosample_term_id=bti)
    payload["biosample_term_name_id"] = pulsar_btn_rec["id"]
    # biosample_type should already be in Pulsar.biosample_type, so won't check to add it first.
    payload["biosample_type_id"] = models.BiosampleType.find_by({"name": dcc_rec["biosample_type"]})["id"]
    date_obtained = dcc_rec.get("date_obtained")
    if not date_obtained:
        date_obtained = dcc_rec.get("culture_harvest_date")
    payload["date_biosample_taken"] = date_obtained
    payload["description"] = dcc_rec.get("description")
    payload["donor_id"] = donor(dcc_rec["donor"]["@id"])["id"]
    payload["lot_identifier"] = dcc_rec.get("lot_id")
    payload["nih_institutional_certification"] = dcc_rec.get("nih_institutional_certification")
    part_of_biosample = dcc_rec.get("part_of")
    if part_of_biosample:
       # Backport the parent.
       pulsar_parent = biosample(part_of_biosample)
       payload["part_of_biosample_id"] = pulsar_parent["id"]
    payload["tissue_preservation_method"] = dcc_rec.get("preservation_method")
    payload["passage_number"] = dcc_rec.get("passage_number")
    payload["starting_amount"] = dcc_rec.get("starting_amount")
    payload["starting_amount_units"] = dcc_rec.get("starting_amount_units")
    payload["submitter_comments"] = dcc_rec.get("submitter_comments")
    payload["vendor_id"] = vendor(dcc_rec["source"]["@id"])["id"]
    payload["vendor_product_identifier"] = dcc_rec.get("product_id")

    post = False
    if not pulsar_rec:
        post = True
        pulsar_rec = models.Biosample.post(payload)

    pulsar_obj = models.Biosample(pulsar_rec["id"])
    del pulsar_rec
    if patch and not post: #don't do both
        pulsar_obj.patch(payload)

    # Patch in any documents (has_many relationship)
    dcc_document_ids = dcc_rec["documents"] # list of DCC identifiers
    linked_document_ids = [x["id"] for x in pulsar_obj.documents]
    payload = {}
    payload["document_ids"] = []
    for doc in dcc_document_ids:
        pulsar_doc_rec_id = document(doc)["id"]
        if not pulsar_doc_rec_id in linked_document_ids:
            payload["document_ids"].append(pulsar_doc_rec_id)
    if payload["document_ids"]:
       pulsar_obj.patch(payload)
    payload = {}
        
    # Patch in any treatments (has_many relationship)
    dcc_treatment_ids = dcc_rec["treatments"] # list of DCC identifiers
    linked_treatment_ids = [x["id"] for x in pulsar_obj.treatments]
    payload = {}
    payload["treatment_ids"] = []
    for treat in dcc_treatment_ids:
        pulsar_treat_rec_id = treatment(treat)["id"]
        if not pulsar_treat_rec_id in linked_treatment_ids:
            payload["treatment_ids"].append(pulsar_treat_rec_id)
    if payload["treatment_ids"]:
       pulsar_obj.patch(payload)
    payload = {}

    # Check if any CRISPR genetic_modifications and if so, associate with biosample.
    # In Pulsar, a biosample can only have one CRISPR genetic modification, so if there are
    # several here specified from the Portal, than that is a problem.
    # CANT backport GMS, see note in the method genetic_modification.
#    dcc_genetic_modification_ids = dcc_rec["genetic_modifications"] # list of DCC identifiers
#    crispr_gm_count = 0
#    for g in dcc_genetic_modifications_ids:
#        gm = ENC_CONN.get(g, ignore404=False)
#        method = gm["method"]
#        if method != "CRISPR":
#            continue
#        crispr_gm_count += 1
#        if crispr_gm_count > 1:
#            raise Exception("Biosample {} has more than 1 genetic modification when only one is expected in Pulsary.".format(pulsar_obj.id))
#        crispr_modification(dcc_rec=gm, pulsar_biosample_id=pulsar_obj.id)
    return pulsar_obj.attrs

def crispr_modification(dcc_rec, pulsar_biosample_id, patch=False):
    """
    Backports a CRISPR genetic_modification record belonging to
    https://www.encodeproject.org/profiles/genetic_modification.json into the Pulsar model called
    CripsrModification. A CRISPR genetic_modification
    has the "method" property set to "CRISPR".

    Identifying properties on the Portal are "accession", "aliases", and "uuid".
    Required properties on the Portal include "category", "method", and "purpose".

    Args:
        dcc_rec: `dict`. The JSON serialization of a genetic_modification record from
            the ENCODE Portal.
        pulsar_biosample_id: `int`. The ID of the Biosample record in Pulsar with which to
            associate the genetic modification.

    Returns:
        `dict`: The JSON representation of the existing Document if it already exists in
        in Pulsar, otherwise the POST response.
    """
    raise Exception("Currently can't backport genetic_modifications from the ENCODE Portal since the original submissions don't specify anything for the 'reagents' list (where the CRISPR construct and donor construct addgene links are stored). Since a GM requires that those objects be created first in Pulsar, the entire GM objects will need to be manually created.")
    aliases = dcc_rec[ALIASES_PROP]
    accession = dcc_rec[ACCESSION_PROP]
    # Check if upstream exists already in Pulsar:
    pulsar_rec = models.CrisprModification.find_by({UPSTREAM_PROP: [*aliases, accession, rec[UUID_PROP], rec["@id"]]})
    if pulsar_rec and not patch:
        return pulsar_rec
    payload = {}
    method = dcc_rec["method"]
    if method != "CRISPR":
        raise Exception("Only CRISPR gentetic_modifications can be backported into Pulsar at this time.")
    characterizations = dcc_rec["characterizations"]
    if characterizations:
         raise Exception("Error while backporting genetic_modification '{}': Backporting characterizations for crispr modifications is not yet supported.".format(dcc_rec["@id"]))
    reagents = dcc_rec["reagents"]
    if reagents:
         raise Exception("Error while backporting genetic_modification '{}': Backporting reagents for crispr modifications is not yet supported.".format(dcc_rec["@id"]))

    payload[UPSTREAM_PROP] = accession
    payload["name"] = set_name(dcc_rec)
    payload["category"] = dcc_rec["category"]
    payload["description"] = dcc_rec["description"]
    payload["purpose"] = dcc_rec["purpose"]

def biosample_term_name(biosample_term_name, biosample_term_id):
    """
    On the ENCODE Portal, a biosample record has a biosample_term_name property and a biosample_term_id
    property. These two properties in Pulsar fall under the BiosampleTermName model with names
    'biosample_term_name' and 'accession', respectivly.
    There isn't a corresponding model for BiosampleTermName on the ENCODE Portal, and to be able to
    determine whether Pulsar already has the provided biosample_term_name, a lookup in Pulsar will
    be done to try and match up the provided biosample_term_id with BiosampleTermName.accession.

    Args:
        biosample_term_id: `str`. The value of a biosample's 'biosample_term_id' property on the Portal.
        biosample_term_name: `str`. The value of a biosample's 'biosample_term_name' property on the Portal.

    Returns:
        `dict`: The JSON representation of the existing BiosampleTermName if it already exists in
        in Pulsar, otherwise the POST response.
    """
    payload = {}
    payload["name"] = biosample_term_name
    payload[ACCESSION_PROP] = biosample_term_id
    ontology_name = biosample_term_id.split(":")[0]
    payload["biosample_ontology_id"] = models.BiosampleOntology.find_by({"name": ontology_name})["id"]
    # Check if upstream exists already in Pulsar:
    pulsar_rec = models.BiosampleTermName.find_by({ACCESSION_PROP: biosample_term_id})
    if pulsar_rec:
        return pulsar_rec
    return models.BiosampleTermName.post(payload)


def document(rec_id):
    """
    Backports a document record belonging to
    https://www.encodeproject.org/profiles/document.json.
    Example document: https://www.encodeproject.org/documents/716003cd-3ce7-41ce-b1e3-6f203b6632a0/?format=json

    Identifying properties on the Portal are "aliases" and "uuid".

    Args:
        rec_id: `str`. An identifier for a document record on the ENCODE Portal.

    Returns:
        `dict`: The JSON representation of the existing Document if it already exists in
        in Pulsar, otherwise the POST response.
    """
    raise Exception("Due to a bug in the ENCODE Portal, can't fetch document '{}' via the REST API.".format(rec))
    rec = ENC_CONN.get(rec_id, ignore404=False)
    aliases = rec[ALIASES_PROP]
    # Check if upstream exists already in Pulsar:
    pulsar_rec = models.Document.find_by({UPSTREAM_PROP: [*aliases, rec[UUID_PROP], rec["@id"]]})
    if pulsar_rec:
        return pulsar_rec
    payload = {}
    payload["description"] = rec["description"]
    if aliases:
        upstream_identifier = aliases[0]
    else:
        upstream_identifier = rec[UUID_PROP]
    payload[UPSTREAM_PROP] = upstream_identifier
    document_type = rec["document_type"]
    # Determine whether this is a protocol document and set Document.is_protocol accordingly.
    protocol = False
    if protocol_regx.search(document_type):
        protocol = True
    payload["is_protocol"] = protocol

def donor(rec_id):
    """
    Backports a huma-donor record belonging to
    https://www.encodeproject.org/profiles/human_donor.json.

    The record will be checked for existence in Pulsar by doing a search on the field
    `donor.upstread_identifer`` using as a query value the record's accession on the ENCODE Portal,
    and also its aliases alias.

    Args:
        rec_id: `str`. An identifier (alias or uuid) for a human-donor record on the ENCODE Portal.

    Returns:
        `dict`: The JSON representation of the existing Donor if it already exists in
        in Pulsar, otherwise the POST response.
    """
    rec = ENC_CONN.get(rec_id, ignore404=False)
    accession = rec[ACCESSION_PROP]
    aliases = rec[ALIASES_PROP]
    # Check if upstream exists already in Pulsar:
    pulsar_rec = models.Donor.find_by({UPSTREAM_PROP: [accession, *aliases]})
    if pulsar_rec:
        return pulsar_rec
    payload = {}
    AGE_PROP = "age"
    if AGE_PROP in rec:
        payload[AGE_PROP] = rec[AGE_PROP]
    GENDER_PROP = "sex"
    if GENDER_PROP in rec:
        payload["gender"] = rec[GENDER_PROP]
    payload[UPSTREAM_PROP] = accession
    payload["name"] = euu.strip_alias_prefix(aliases[0])
    return models.Donor.post(payload)

def treatment(rec_id, patch=False):
    """
    Backports a treatement record belonging to
    https://www.encodeproject.org/profiles/treatment.json.
    The required properties in the ENCODE Portal are:

    1. treatment_term_name,
    2. treatment_type.

    An example on the Portal: https://www.encodeproject.org/treatments/933a1ff2-43a2-4a54-9c87-aad228d0033e/.
    Identifying properties on the Portal are 'aliases' and 'uuid'.

    Args:
        rec_id: `str`. An identifier (alias or uuid) for a treatment record on the ENCODE Portal.

    Returns:
        `dict`: The JSON representation of the existing Treatment if it already exists in
        in Pulsar, otherwise the POST response.
    """
    rec = ENC_CONN.get(rec_id, ignore404=False)
    aliases = rec[ALIASES_PROP]
    #check if upstream exists already in Pulsar:
    pulsar_rec = models.Treatment.find_by({UPSTREAM_PROP: [*aliases, rec[UUID_PROP], rec["@id"]]})
    if pulsar_rec and not patch:
        return pulsar_rec

    payload = {}
    payload["aliases"] = aliases
    pulsar_document_ids = []
    for doc_id in rec["documents"]:
        pulsar_document_ids.append(document(doc_id))
    payload["document_ids"] = pulsar_document_ids
    payload[UPSTREAM_PROP] = rec[UUID_PROP] 
    payload["concentration"] = rec["amount"]
    amount_units = rec["amount_units"]
    pulsar_amount_units_rec = models.ConcentrationUnit.find_by(payload={"name": amount_units})
    if not pulsar_amount_units_rec:
        raise Exception("ConcentrationUnit '{}' couldn't be found.".format(amount_units))
    pulsar_amount_units_rec = pulsar_amount_units_rec["concentration_unit"]
    payload["concentration_unit_id"] = pulsar_amount_units_rec["id"]
    payload["duration"] = rec["duration"]
    payload["duration_units"] = rec["duration_units"]
    temp_prop_name = "temperature"
    if temp_prop_name  in rec:
        temp = rec[temp_prop_name]
        temp_units = rec["temperature_units"]
        if temp_units == "Kelvin":
            temp = pulsarpy.utils.kelvin_to_celsius(temp)
        elif temp_units == "Fahrenheit":
            temp = pulsarpy.utils.fahrenheit_to_celsius(temp)
        payload[temp_prop_name] = rec[temp_prop_name]
    ttn = rec["treatment_term_name"]
    tti = rec["treatment_term_id"]
    #Create new TreatmentTermName if not found in Pulsar.
    pulsar_ttn_rec = treatment_term_name(ttn, tti)["treatment_term_name"]
    payload["treatment_term_name_id"] = pulsar_ttn_rec["id"]
    payload["treatment_type"] = rec["treatment_type"]
    #The Portal's treatment model doesn't include a name prop or a description prop.
    # Thus, 'name' shall be set to be equal to 'upstream_identifier'.
    payload["name"] = euu.strip_alias_prefix(upstream)
    if not pulsar_rec:
        return models.Treatment.post(payload)
    else:
        pulsar_obj = models.Treatment(pulsar_rec["id"])
        return pulsar_obj.patch(payload)

def treatment_term_name(treatment_term_name, treatment_term_id):
    """
    On the ENCODE Portal, a treatment record has a treatment_term_name property and a treatment_term_id
    property. These two properties in Pulsar fall under the TreatmentTermName model with names
    'treatment_term_name' and 'accession', respectivly.

    There isn't a corresponding model for TreatmentTermName on the ENCODE Portal, and to be able to
    determine whether Pulsar already has the provided treatment_term_name, a lookup in Pulsar will
    be done to try and match up the provided treatment_term_id with TreatmentTermName.accession.

    Args:
        treatment_term_id: `str`. The value of a treatment's 'treatment_term_id' property on the Portal.
        treatment_term_name: `str`. The value of a treatment's 'treatment_term_name' property on the Portal.

    Returns:
        `dict`: The JSON representation of the existing TreatmentTermName if it already exists in
        in Pulsar, otherwise the POST response.
    """
    payload = {}
    payload["name"] = treatment_term_name
    payload[ACCESSION_PROP] = treatment_term_id
    pulsar_rec = models.TreatmentTermName.find_by({ACCESSION_PROP: treatment_term_id})
    if pulsar_rec:
        return pulsar_rec
    return models.TreatmentTermName.post(payload)

def vendor(rec_id):
    """
    Backports a source record belonging to
    https://www.encodeproject.org/profiles/source.json.

    Identifying properties on the Portal are "name", and "uuid".
    Portal's required props are: "name", and "title".

    Args:
        rec_id: `str`. An identifier (uuid or value of the '@id' property) for a source record on the ENCODE Portal.
            Note that even though a source's 'name' property is specified as identifying on the ENCODE Portal,
            that alone won't work for a GET request to pull down that record from the Portal.

    Returns:
        `dict`: The JSON representation of the existing source if it already exists in
        in Pulsar, otherwise the POST response.
    """
    rec = ENC_CONN.get(rec_id, ignore404=False)
    name = rec["name"]
    uuid = rec["uuid"]
    pulsar_rec = models.Vendor.find_by({UPSTREAM_PROP: [name, uuid, rec["@id"]]})
    if pulsar_rec:
        return pulsar_rec
    payload = {}
    payload["description"] = rec["description"]
    payload["name"] = name
    payload[UPSTREAM_PROP] = rec["@id"]
    payload["url"] = rec["url"]
    return models.Vendor.post(payload)

