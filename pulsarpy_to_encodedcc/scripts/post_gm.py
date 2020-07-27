import pulsarpy_to_encodedcc.dcc_submit as dcc_submit

sub = dcc_submit.Submit(dcc_mode='prod')
sub.post_crispr_modification(391)
sub.post_crispr_modification(393)
