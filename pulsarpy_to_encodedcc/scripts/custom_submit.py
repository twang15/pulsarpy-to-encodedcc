#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###Author
#Nathaniel Watson
#2018-12-15
#nathankw@stanford.edu
###

"""
This script is meant to be called in a for-loop context, i.e. when submitting one experiment after
another.
"""

import pulsarpy_to_encodedcc.dcc_submit as d
mode = "prod"
s = d.Submit(dcc_mode=mode)
ids = []
#ids.append(5) # Done
#ids.append(10) # Done

#ids.append(21) # Missing paired-input rep
#ids.append(201) # Missing paired-input rep

ids.append(200)
ids.append(198)
ids.append(197)
ids.append(196)
ids.append(193)
ids.append(194)
ids.append(192)
ids.append(188)
ids.append(191)
ids.append(190)
ids.append(189)
ids.append(22)
ids.append(185)
ids.append(186)
ids.append(184)
ids.append(183)
ids.append(182)
ids.append(181)
ids.append(180)
ids.append(179)
ids.append(178)
ids.append(177)
ids.append(176)
ids.append(175)
ids.append(174)
ids.append(173)
ids.append(172)
ids.append(171)
ids.append(170)
ids.append(169)
ids.append(168)
ids.append(167)
ids.append(166)
ids.append(165)
ids.append(164)
ids.append(163)
ids.append(162)
ids.append(161)
ids.append(160)
ids.append(159)
ids.append(158)
ids.append(157)
ids.append(156)
ids.append(155)
ids.append(154)
ids.append(153)
ids.append(152)
ids.append(151)
ids.append(150)
ids.append(149)
ids.append(148)
ids.append(147)
ids.append(146)
ids.append(145)
ids.append(144)
ids.append(143)
ids.append(142)
ids.append(141)
ids.append(140)
ids.append(139)
ids.append(138)
ids.append(137)
ids.append(136)
ids.append(50)
ids.append(49)
ids.append(48)
ids.append(47)
ids.append(46)
ids.append(45)
ids.append(44)
ids.append(43)
ids.append(42)
ids.append(41)
ids.append(40)
ids.append(39)
ids.append(38)
ids.append(37)
ids.append(36)
ids.append(35)
ids.append(34)
ids.append(33)
ids.append(32)
ids.append(31)
ids.append(30)
ids.append(29)
ids.append(28)
ids.append(27)
ids.append(20)
ids.append(19)
ids.append(18)
ids.append(17)
ids.append(16)
ids.append(15)
ids.append(14)
ids.append(13)
ids.append(12)
ids.append(11)
ids.append(9)
ids.append(8)
ids.append(7)
ids.append(6)
ids.append(199)
ids.append(195)

for i in ids:
    try:
        s.post_chipseq_exp(rec_id=i)
    except d.ExpMissingReplicates:
        # Already logged to error log file 
        continue
