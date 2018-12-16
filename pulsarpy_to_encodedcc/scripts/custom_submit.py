#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###Author
#Nathaniel Watson
#2018-12-15
#nathankw@stanford.edu
###

import pulsarpy_to_encodedcc.dcc_submit as d
mode = "v79x0-test-master.demo.encodedcc.org" 
s = d.Submit(dcc_mode=mode)
s.post_chipseq_exp(rec_id=165, patch=False)
