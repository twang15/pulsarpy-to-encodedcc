# patch genetic modification for biosamples
import pulsarpy_to_encodedcc.dcc_submit as dcc_submit

sub = dcc_submit.Submit(dcc_mode='prod')

# patch biosamples for GM
# in dcc_submit.py, we need dont_extend_arrays = True

biosample_ids = [
    12386
        ]

for bid in biosample_ids:
  sub.post_biosample(rec_id=bid, patch=True)
