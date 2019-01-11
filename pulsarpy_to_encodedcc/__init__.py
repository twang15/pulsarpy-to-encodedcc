# -*- coding: utf-8 -*-

###
# © 2018 The Board of Trustees of the Leland Stanford Junior University
# Nathaniel Watson
# nathankw@stanford.edu
###

import logging
import os
import sys

#: The folder where FASTQ files are downloaded. Will be created if it doesn't yet exist.
FASTQ_FOLDER = "PulsarpyToEncodeDcc_FASTQS"
if not os.path.exists(FASTQ_FOLDER):
    os.mkdir(FASTQ_FOLDER)

LOG_DIR = "PulsarpyToEncodeDcc_Logs"
if not os.path.exists(LOG_DIR):
    os.mkdir(LOG_DIR)

#: The name of the error ``logging`` instance.
ERROR_LOGGER_NAME = __package__ + "_error"

#: A ``logging`` instance that accepts messages at the ERROR level.
error_logger = logging.getLogger(ERROR_LOGGER_NAME)
error_logger.setLevel(logging.ERROR)
f_formatter = logging.Formatter('%(asctime)s:%(name)s:\t%(message)s')
file_name = os.path.join(LOG_DIR, "log_" + ERROR_LOGGER_NAME + ".txt")
file_handler = logging.FileHandler(filename=file_name, mode="a")
file_handler.setLevel(logging.ERROR)
file_handler.setFormatter(f_formatter)
error_logger.addHandler(file_handler)

def log_error(msg)
    error_logger.error(msg)
    print(msg)
