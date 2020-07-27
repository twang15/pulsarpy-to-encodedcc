import pulsarpy_to_encodedcc.dcc_submit as dcc_submit
import traceback

def log_traceback(ex):
    tb_lines = traceback.format_exception(ex.__class__, ex, ex.__traceback__)
    for tb in tb_lines:
        print(tb)

bulkAtacSeq = dcc_submit.Submit(dcc_mode='prod')
#chipSeq = dcc_submit.Submit(dcc_mode='dev')

records = [
#186, ENCSR483RKN
#188, ENCSR042AWH
#192, ENCSR422SUG
#191, ENCSR872WGW
#224, ENCSR499ASS
#225, ENCSR095QNB
#226, ENCSR591PIX
#189, ENCSR200OML
        ]

for record in records:
    bulkAtacSeq.post_bulk_atacseq_exp(rec_id=record, patch=False)
