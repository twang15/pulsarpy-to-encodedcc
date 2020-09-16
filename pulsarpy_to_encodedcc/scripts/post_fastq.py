from pulsarpy_to_encodedcc.dcc_submit import Submit as submit

fastq_submit = submit(dcc_mode='prod')
replicates = [
    (698,"/replicates/d8858198-521f-4153-a487-b4875b482082/")
    ]

for (res_id, replicate_id) in replicates:
  try:
    fastq_submit.post_fastq_file(res_id, 1, replicate_id)
  except:
    print (res_id)
    continue
