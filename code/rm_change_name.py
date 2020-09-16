import os

fld = "../data/assessing_models_simulation/datasets/"
x = [i for i in os.listdir("../data/assessing_models_simulation/datasets/") if "2020" in i]
for j in x:
   j2 = j.split("_")
   print('mv ' + fld + j + ' ' + fld + '_'.join(j2[:3]) + '_10_' + '_'.join(j2[3:len(j2)]))


