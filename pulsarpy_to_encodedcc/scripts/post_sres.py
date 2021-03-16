from pulsarpy_to_encodedcc.dcc_submit import Submit as submit

fastq_submit = submit(dcc_mode='prod')
replicates = [
    #(235, "/replicates/2afe98ad-bae9-45ae-8aaf-3fd8ea69f7ce/")
    #(573, "/replicates/fd657a62-a797-4247-bfd2-7994125cd80b/"),
    #(574, "/replicates/b4fc3a5b-6c75-4046-9d61-bbed83ea99ff/")
    (576, "/replicates/9c7f4f16-c9c5-4e1a-a852-652942913042/"),
    (577, "/replicates/71806725-c48e-40e4-8cfa-769270f5a8d3/"),
    (578, "/replicates/b9709d4e-8cc3-487a-b89a-68589e68aeb7/"),
    (579, "/replicates/1cdeaa53-5529-4e33-8682-5afe94f2d8ea/")
    ]

for (res_id, replicate_id) in replicates:
  try:
    fastq_submit.post_sres(res_id, replicate_id)
  except:
    print (res_id)
    continue
