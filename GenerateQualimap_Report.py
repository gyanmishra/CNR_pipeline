#!/usr/bin/python

import os
import sys
import pandas as pd

# Fetch qualimap directory from console with "-d" as identifier.
qualimap_dirIndex = sys.argv.index("-d") if "-d" in sys.argv else None
qualimap_dir = sys.argv[qualimap_dirIndex+1] if qualimap_dirIndex and ((len(sys.argv)-1) > qualimap_dirIndex) else None

if sys.argv[qualimap_dirIndex+1]:
       print("Qualimap input folder is : " + sys.argv[qualimap_dirIndex+1])
else:
       print("Invalid value for argument \"-d\".")
       sys.exit()

# list to store files
res = []
store = {}
report_dict = {}    
# Iterate directory
for file_path in os.listdir(qualimap_dir):
    print(file_path)
    if os.path.isdir(os.path.join(qualimap_dir, file_path)):

       # uncomment line depending upon the type of run
       if os.path.isfile(os.path.join(qualimap_dir, file_path,'rnaseq_qc_results.txt')):
              res.append({file_path : os.path.join(qualimap_dir, file_path,'rnaseq_qc_results.txt')})
              print(file_path)
       if os.path.isfile(os.path.join(qualimap_dir, file_path,'genome_results.txt')):
              res.append({file_path : os.path.join(qualimap_dir, file_path,'genome_results.txt')})
              print(file_path)

for i in res:
       for j in i.values():
              #print(j)
              file = open(j, 'r')
              lines = file.readlines()
              store = {}
              for line in lines:
                    if '=' in line:
                            temp = line.strip("\n").strip().split('=')
                            store[temp[0].strip()]=temp[1].strip()
                            #print(store)
       for k in i.keys():
              report_dict[k]=store

df = pd.DataFrame.from_dict(report_dict, orient='index')
df.to_csv(qualimap_dir+"Qualimap_report.tsv", sep='\t')

                                   

