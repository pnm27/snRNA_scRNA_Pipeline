#!/usr/bin/env python



import yaml, os, re, glob2, pandas as pd, numpy as np, subprocess, shlex

with open(r'../new_config.yaml', 'r') as fout:
	a=yaml.load(fout, Loader=yaml.FullLoader)


all_logs=glob2.glob(os.path.join(b['split_bams_pipeline']['split_bams_proxy_dir'], "*.txt"))
job_id_l=[]

for log in all_logs:
	try:
		with open('log', 'r') as fin:
			for line in fin:
				if line.startswith("Submitted script for donor"):
					job_id_l.append(line.split(':')[-1].strip())

	except IOError as e:
	    if e.errno == errno.EACCES:
	        print("file exists, but isn't readable")
	    elif e.errno == errno.ENOENT:
	        print("files isn't readable because it isn't there")


for jid in job_id_l:
	a=subprocess.run(shlex.split(f"bhist -l {jid}"), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	samp_n = re.search("Job Name <(.*?)>,", a.stdout.decode('utf-8')).group(1).split('_')[1]
	if a.returncode == 0:
		print(f"Successfull: Job id {jid} for the sample {samp_n}")
	else:
		print(f"Exited: Job id {jid} for the pool {samp_n} faced some issues. Job details:")
		print(f"{a.stdout.decode('utf-8')}")