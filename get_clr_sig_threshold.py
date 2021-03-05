import numpy as np
import pandas as pd
from sys import argv

data_in = argv[1]

sf_out = pd.read_csv(data_in, sep = '\t', names = ['Pos', 'Clr','Alpha'])

sf_out.sort_values(by = ['Clr'],  inplace = True, ignore_index = True)

cut_95 = round(len(sf_out)*0.95)-1
cut_99 = round(len(sf_out)*0.99)-1
cut_999 = round(len(sf_out)*0.999)-1

out_cut_95 = sf_out.iloc[cut_95]['Clr']
out_cut_99 = sf_out.iloc[cut_99]['Clr']
out_cut_999 = sf_out.iloc[cut_999]['Clr']

print("Significance threshold 0.05: " + str(out_cut_95), end = '\n')
print("Significance threshold 0.01: " + str(out_cut_99), end = '\n')
print("Significance threshold 0.001: " + str(out_cut_999), end = '\n')
