###                                                                                                     
# © 2018 The Board of Trustees of the Leland Stanford Junior University                                 
# Nathaniel Watson                                                                                      
# nathankw@stanford.edu                                                                                 
### 

from pulsarpy import models

class ValidateChipseqExperiments():
    def __init__(self, ids, outfile):
        """
        Called before posting ChipseqExperiments. Works by ensuring that:

          1. target.upstream_identifier is present,
          2. Experiment replicates each have a GM.

        Write any issues to a report file. 

        Args:
            ids: `list`. One or more Pulsar ChipseqExperiment identifiers.
            outfile: `str`. The name of the report file to output. If it exists already it will
                be overwritten.
        """
        self.chipseq_experiment_ids = ids
        #: List of Pulsar Target IDs that don't have an upstream_identifier set. 
        self.unregistered_targets = []
        #: The report file to create.
        self.fout = open(outfile, 'w')
        self.validate_chipseq_experiments()

    def validate_chipseq_experiments(self):
        for rec_id in self.chipseq_experiment_ids:
            exp = models.ChipseqExperiment(rec_id)
            str_exp_id = str(exp.id)
            target = models.Target(exp.target_id)
            target_upstream = target.upstream_identifier.strip()
            if not target_upstream:
                if target.id not in self.unregistered_targets:
                    self.unregistered_targets.append(target.id)
                    msg = "Target {} not registred.".format(target.id)
                    self.fout.write(str_exp_id + "\t" + msg + "\n")
            exp_rep_ids = exp.replicate_ids
            if not exp_rep_ids:
                msg = "ChipseqExperiment missing replicates."
                self.fout.write(str_exp_id + "\t" + msg + "\n")
            else:
                for rep_id in exp_rep_ids:
                    self.validate_gm_for_crispr_biosample(chip_exp_id=exp.id, biosample_id=rep_id)
                      
            # Check WT input present
            if not exp.wild_type_control_id:
                msg = "ChipseqExperiment missing WT control."
                self.fout.write(str_exp_id + "\t" + msg + "\n")
            # Check Paired input present and has GM
            pi_ids = [
                models.Library(lib_id).biosample_id
                for lib_id in exp.control_replicate_ids
            ]
            if not pi_ids:
                msg = "ChipseqExperiment missing paired input."
                self.fout.write(str_exp_id + "\t" + msg + "\n")
            else:
                for pi_id in pi_ids:
                    self.validate_gm_for_crispr_biosample(chip_exp_id=exp.id, biosample_id=pi_id)
        self.fout.close()

    def validate_gm_for_crispr_biosample(self, chip_exp_id, biosample_id):
        biosample = models.Biosample(biosample_id)
        gm_id = biosample.crispr_modification_id
        if not gm_id:
            msg = "Exp. biosample {} missing GM.".format(biosample.id)
            self.fout.write((str(chip_exp_id) + "\t" + msg + "\n"))
        else:
            # Verify that GM target is present.
            gm = models.CrisprModification(gm_id)
            donor_construct = models.DonorConstruct(gm.donor_construct_id)
            dc_target = models.Target(donor_construct.target_id)
            if not dc_target.upstream_identifier:
                if dc_target.id not in self.unregistered_targets:
                    self.unregistered_targets.append(dc_target.id)
                    msg = "Target {} not registred.".format(dc_target.id)
                    self.fout.write(str(chip_exp_id) + "\t" + msg + "\n")
            ccs = [models.CrisprConstruct(x) for x in gm.crispr_construct_ids]
            for c in ccs:
                cc_target = models.Target(c.target_id)
                if not cc_target.upstream_identifier:
                    if cc_target.id not in self.unregistered_targets:
                        self.unregistered_targets.append(cc_target.id)
                        msg = "Target {} not registred.".format(cc_target.id)
                        self.fout.write((str(chip_exp_id) + "\t" + msg + "\n"))
