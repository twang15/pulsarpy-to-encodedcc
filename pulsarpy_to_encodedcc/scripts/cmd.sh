# submit pcrs
pudb post_pcr.py -i pcrs -o t

# submit chipseq

# submit atacseq
pudb atac-seq.py -i atac -m "Atacseq"


# derive projects ready to delete on DX
pudb get_sreqs.py -i lib
awk '{print $4}' enc | sort | uniq | awk '{printf "https://pulsar-encode.herokuapp.com/data_storages/%s\n", $1}'
awk '{printf "%s %s\n", $3, $4}' enc | sort -k2 | uniq | awk '{printf "https://pulsar-encode.herokuapp.com/sequencing_requests/%s\n", $1}'

# PCR submission
pudb post_pcr.py -i pcrs -o h

# IP submission
pudb post_ip_biosample_characterizations_byENC.py -i ips
