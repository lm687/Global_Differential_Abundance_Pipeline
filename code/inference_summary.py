## Read the output from slurm to see if there were divergent transitions after warmup, high R hat, low ESS, etc.

## Run: python inference_summary.py  | sort > ../data/diagnostics/inference_diagnostics.csv

import os
import re
import numpy as np
from datetime import datetime, date

months = {'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6, 'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12}

fles = [i for i in os.listdir() if "slurm-" in i]

def to_float_time(i):
	return(float(re.sub(r' seconds.*', '', re.sub(r'^.*?Elapsed Time: ', '', i), flags=re.DOTALL)))
def to_float_time2(i):
	return((re.sub(r'Chain [0-9]+:', '', re.sub(r' seconds.*', '',i, flags=re.DOTALL)).strip()))
	# (re.sub(r' seconds.*', '', re.sub(r'^.*?: ', '', i), flags=re.DOTALL))
	# (re.sub(r' seconds.*', '',i, flags=re.DOTALL))

def give_datetime(data, hms):
	if(len(data) == 5): ## short version
       		return(datetime(year=int(data[4]), month=months[data[1]], day=int(data[2]), hour=int(hms[0]), minute=int(hms[1]), second=int(hms[2])))
	elif(len(data) == 6): ## 6 char; there is a space
                return(datetime(year=int(data[5]), month=months[data[1]], day=int(data[3]), hour=int(hms[0]), minute=int(hms[1]), second=int(hms[2])))
	else:
		exit()

# print("00 CT \t features \t divergent transitions \t Rhat high \t ESS")
print("0timestamp, slurm,CT,features,nits,model,divergent transitions,Rhat high,ESS,Cancelled (time limit), warmup time (s), sampling time (s), total time (s), real time (s)")
for fle in fles:
	# print(fle)
	list_outfile = None; #line_elapsed_time = None;
	time_limit = ''; div_trans = ''; Tail_ESS = ''; Rhat = ''
	a = open(fle, "r")
	# contents = a.read()

	# print('Elapsed' in contents)

	## search lines
	warm_up_time = []; sampling_time = []; total_time = []; dates = []
	timestamp = None
	for line in a:
		if re.search("output\:", line):
			list_outfile = os.path.basename(line).rstrip("ROO.RData\n").split('_')
		if re.search("divergent transitions", line):
			div_trans = True
		if re.search("Tail-ESS", line):
			Tail_ESS = True
		if re.search("Rhat", line):
			Rhat = True
		if re.search("DUE TO TIME LIMIT", line):
			time_limit = True
		if re.search("2020", line) and re.search(r"\[*\]", line):
			timestamp = line.strip('\n')
		if re.search("seconds \(Warm-up\)", line):
			warm_up_time.append(line)
		if re.search("seconds \(Sampling\)", line):
			sampling_time.append(line)
		if re.search("seconds \(Total\)", line):
			total_time.append(line)
		if re.search("2020", line):
			if not re.search("slurmstepd", line):
				dates.append(line)

	warm_up_time = np.array([ to_float_time(i) for i in warm_up_time]).astype(np.float)
	sampling_time = np.array([ to_float_time2(i) for i in sampling_time]).astype(np.float)
	total_time = np.array([ to_float_time2(i) for i in total_time]).astype(np.float)

	if len(total_time > 0):
		times_avg = ','.join([str(np.round(np.mean(warm_up_time))), str(np.round(np.mean(sampling_time))), str(np.round(np.mean(total_time)))])
		# print('\t'.join([str(np.mean(warm_up_time)), str(np.mean(sampling_time))]))		times_avg = ','.join([str(np.round(np.mean(warm_up_time))),\
		 # str(np.round(np.mean(sampling_time))), str(np.round(np.mean(total_time)))])
	else:
		times_avg = ','.join(['', '', ''])

	if(len(dates) == 1):
		real_time = ''		
	else:
		#print((dates))
		a = dates[0] #"[Sun Aug  2 16:46:45 2020]"
		b = dates[1] #"[Sun Aug  2 16:47:34 2020]"
	
		dataa = (a.strip('\n').strip('[|]').split(" "))
		if(len(dataa) == 5):
			hmsa = dataa[3].split(":")
		else:
			hmsa = dataa[4].split(":")
	
		datab = (b.strip('\n').strip('[|]').split(" "))
		if(len(datab) == 5):
			hmsb = datab[3].split(":")
		else:
			hmsb = datab[4].split(":")

		#print(dataa)
		#print(fle)
		#print(hmsa)	
		a_datetime = give_datetime(dataa, hmsa)
		b_datetime = give_datetime(datab, hmsb)
		real_time = str( int ((b_datetime-a_datetime).total_seconds()) )

	print("%s,%s,%s,%s,%s,%s,%s,%s,%s" % (timestamp,fle, ','.join(list_outfile), div_trans, Rhat, Tail_ESS, time_limit, times_avg, real_time))



