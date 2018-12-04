# -*- coding: utf-8 -*-                                                                                
                                                                                                       
###                                                                                                    
# © 2018 The Board of Trustees of the Leland Stanford Junior University                                
# Nathaniel Watson                                                                                     
# nathankw@stanford.edu                                                                                
### 

import os

#: The folder where FASTQ files are downloaded. Will be created if it doesn't yet exist.
FASTQ_FOLDER = "PulsarpyToEncodedcc_FASTQS"
if not os.path.exists(FASTQ_FOLDER):
    os.mkdir(FASTQ_FOLDER)



