from pulsarpy_to_encodedcc.dcc_submit import Submit as submit

fastq_submit = submit(dcc_mode='prod')
replicates = [
    (883,"/replicates/2991149b-54d8-4c08-a05e-2998455f1157/")
    ]

for (res_id, replicate_id) in replicates:
  try:
    fastq_submit.post_fastq_file(res_id, 1, replicate_id)
  except:
    print (res_id)
    continue
