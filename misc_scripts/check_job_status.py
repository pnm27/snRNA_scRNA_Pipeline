#!/usr/bin/env python

import yaml, os, re, glob2, subprocess, shlex, argparse, datetime
import pandas as pd, numpy as np



#Parse Command-Line arguments
parser = argparse.ArgumentParser(description="Check job_status(es) for all per-donor split bam script(s) (without any parameters) or for a given job (need jobid)")

# Optional parameters
parser.add_argument('-j', '--jobid', help="JobID")
parser.add_argument('-c', '--config_file', help="Config file (assumed to be in YAML format). If the keys inn config.yaml file has been changed please change in the script too!!!!")
# parser.add_argument('-l', '--levels', help="Levels of config file indicating the path to the proxy-dir containing all logs of the per each pool processed by the pipeline (sep by comma)")
parser.add_argument('-d', '--dir', help="Directory containing all logs of the per each pool processed by the pipeline")


args = parser.parse_args()


config_file=args.config_file
# lev=args.levels
log_dir=args.dir

# all([ v==None for k,v in vars(args).items() ])
# Determine script's purpose
if args.jobid == None:
	if config_file != None and os.path.isfile(config_file) and log_dir == None:
		with open(config_file, 'r') as fout:
			a=yaml.load(fout, Loader=yaml.FullLoader)
			# Specify where the log_dir value is at
			all_logs=glob2.glob(os.path.join(a['split_bams_pipeline']['split_bams_proxy_dir'], "*.txt"))
	elif config_file == None and log_dir != None:
			all_logs=glob2.glob(os.path.join(log_dir, "*.txt"))

	else:
		raise ValueError(f"Unsupported combination of the parameters to the script!")

	job_ids_l=[]

	for log in all_logs:
		try:
			with open('log', 'r') as fin:
				for line in fin:
					if line.startswith("Submitted script for donor"):
						job_ids_l.append(line.split(':')[-1].strip())


		# IOEror is  part of this, now
		except OSError as e:
			# Permission denied
		    if e.errno == errno.EACCES:
		        print(f"File {e.filename} exists, but isn't readable")

		    elif e.errno == errno.ENOENT:
		        print(f"File {e.filename} isn't readable because it isn't there")

		    # Handle rest all errors
		    else:
		    	print(e)


	for jid in job_ids_l:
		a=subprocess.run(shlex.split(f"bhist -l {jid}"), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		# This job name is assumed to follow the construct "NPSAD_{samp_num}_<uuid>". PLEASE CHANGE IF OTHERWISE
		samp_n = re.search("Job Name <(.*?)>,", a.stdout.decode('utf-8')).group(1).split('_')[1]
		if a.returncode == 0:
			print(f"Successfull: Job id {jid} for the sample {samp_n}")
		else:
			print(f"Exited: Job id {jid} for the pool {samp_n} faced some issues. Job details:")
			print(f"{a.stdout.decode('utf-8')}")

elif args.jobid != None and config_file == None and log_dir == None:
	a=subprocess.run(shlex.split(f"bhist -l {args.jobid}"), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	ct = datetime.datetime.now()
	print(f"The status for the job with joibid {args.jobid} at {ct} is {a.returncode}")
else:
	raise ValueError(f"Unsupported combination of the parameters to the script!")