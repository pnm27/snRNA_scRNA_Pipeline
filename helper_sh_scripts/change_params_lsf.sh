#!/usr/bin/bash

# Run this script to change certain parameters for all jobs with a pattern (contains) in job name
# sh change_params.sh <option_to_change> <value> <job_name_pattern>
# eg: sh change_params.sh -W 03:00 -n 6,10 STARsolo will execute bmod -W 03:00 -n 6,10 on all jobs which contains "STARsolo" in the JOB_NAME column


params="${@:1:${#@}-1}"
pat=${@:${#@}-1}

bjobs -w | awk -v pattern=${pat} 'NR==1{for(i=1;i<=NF;i++) if ($i == "JOBID") {j=i} else if ($i == "JOB_NAME") {n=i}} (NR>1 && $n ~ /$pattern/){print $j}' | while read l
do
{
	echo "Changing the following parameters: ${params} for jobid: ${l} that contains the pattern ${pat} in the JOB_NAME column"
	bmod ${params} ${l}
	echo -e "Changed\n\n"

}
done
