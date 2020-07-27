import pulsarpy_to_encodedcc.dcc_submit as dcc_submit
import traceback

def log_traceback(ex):
    tb_lines = traceback.format_exception(ex.__class__, ex, ex.__traceback__)
    for tb in tb_lines:
        print(tb)

chipSeq = dcc_submit.Submit(dcc_mode='prod')
#chipSeq = dcc_submit.Submit(dcc_mode='dev')

records = [
#210,
211,

#275,
#276,
#278,
#279,

#285,
#286,
#287,
#291,
#292,
#293,
#294,
#296,
        ]

for record in records:
    chipSeq.post_chipseq_exp(rec_id=record, patch=False)
#for record in records:
#    try:
#        chipSeq.post_chipseq_exp(rec_id=record, patch=False)
#    except:
#        print ("Chip %d fails!!!".format(record))
#        continue
