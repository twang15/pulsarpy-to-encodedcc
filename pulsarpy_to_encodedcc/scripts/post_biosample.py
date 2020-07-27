import pulsarpy_to_encodedcc.dcc_submit as dcc_submit
import traceback

def log_traceback(ex):
    tb_lines = traceback.format_exception(ex.__class__, ex, ex.__traceback__)
    for tb in tb_lines:
        print(tb)

submitter = dcc_submit.Submit(dcc_mode='prod')
#chipSeq = dcc_submit.Submit(dcc_mode='dev')

records = [
    12106,
        ]

for record in records:
    #submitter.post_biosample(rec_id=record, patch=False)
    submitter.post_biosample(rec_id=record, patch=True)
#for record in records:
#    try:
#        chipSeq.post_chipseq_exp(rec_id=record, patch=False)
#    except:
#        print ("Chip %d fails!!!".format(record))
#        continue
