import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.patches as mpatches

pcawg = pd.read_csv("~/Documents/PhD/GlobalDA/data/restricted/pcawg/pcawg.wg11.final_sample_list_MARCH2019.txt", sep='\t')
x=pd.read_csv("~/Desktop/inference_diagnostics.csv", dtype={'0timestamp': str,\
 'slurm': str, 'CT': str, 'features': str, 'nits': float, 'model': str,
 'divergent transitions': str, 'Rhat high': str, 'ESS': str,\
 'Cancelled (time limit)': str, ' warmup time (s)': float, ' sampling time (s)': float,\
 ' total time (s)': float})


# ct_types =[a.replace('-', '') for a in pcawg['histology_detailed'].tolist()]
ct_types =[a for a in pcawg['histology_detailed'].tolist()]

# set_ct = [a.replace('-', '') for a in set(ct_types)]
set_ct = [a for a in set(ct_types)]
size_dataset_ct_types = ([ct_types.count(a) for a in set_ct])
# print([ x[' sampling time (s)'].tolist()[ct_types.index(i)] for i in x['CT']])

models = [a for a in set(x['model'].tolist())]

''' loop over the entries in the Stan inference file, and find how many samples the
corresponding cancer type has
'''

# print(x[' sampling time (s)'].tolist())


size_dataset_append = []
sampling_time_append = []
total_time_append = []
warmup_time_append = []
model = []
ict = 0
for i in x['CT'].tolist():
	print(i)
	print(set_ct.index(i))
	size_dataset_append.append((size_dataset_ct_types[set_ct.index(i)]))
	# print(x['CT'].tolist()[ict])
	sampling_time_append.append(x[' sampling time (s)'].tolist()[ict])
	total_time_append.append(x[' total time (s)'].tolist()[ict])
	warmup_time_append.append(x[' warmup time (s)'].tolist()[ict])
	model.append(models.index(x['model'].tolist()[ict]))
	ict += 1

# print(sampling_time_append)

class_colours = ['#993838', '#7836ad', '#5ad3b7']
classes =  [a for a in set(x['model'])]
recs = []
for i in range(0,len(class_colours)):
    recs.append(mpatches.Rectangle((0,0),1,1,fc=class_colours[i]))



# plt.scatter(size_dataset_append, sampling_time_append, c=[class_colours[a] for a in model],\
# 	s = 30, marker=4)
# plt.ylabel('Size of dataset')
# plt.ylabel('Running time')
# plt.legend(recs,classes,loc=4)

# for i, txt in enumerate(x['CT'].tolist()):
# 	if size_dataset_append[i] > 120:
# 	    plt.annotate(txt, (size_dataset_append[i], sampling_time_append[i]), size=8)

# plt.show()
# plt.close()



# plt.scatter(size_dataset_append, total_time_append, c=[class_colours[a] for a in model],\
# 	s = 30, marker=4)
# plt.ylabel('Size of dataset')
# plt.ylabel('Running time')
# plt.legend(recs,classes,loc=4)
# for i, txt in enumerate(x['CT'].tolist()):
# 	if size_dataset_append[i] > 120:
# 	    plt.annotate(txt, (size_dataset_append[i], total_time_append[i]), size=8)

# plt.show()
# plt.close()



plt.scatter(size_dataset_append, sampling_time_append, c=[class_colours[a] for a in model],\
	s = 30, marker="^")
plt.ylabel('Size of dataset')
plt.ylabel('Running time')
plt.legend(recs,classes,loc=4)
for i, txt in enumerate(x['CT'].tolist()):
	if size_dataset_append[i] > 120:
	    plt.annotate(txt, (size_dataset_append[i], sampling_time_append[i]), size=8)

plt.show()
plt.close()
