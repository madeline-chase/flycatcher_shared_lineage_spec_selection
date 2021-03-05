import numpy as np
import pandas as pd
from sys import argv

data_in = argv[1]

H_out = pd.read_csv(data_in, sep = '\t', names = ['start','stop','win_sites', 's_sum', 'theta_w', 'theta_pi_sum', 'theta_h_sum', 'theta_l_sum', 'var_pi_l', 'h_stat'])

H_out.sort_values(by = ['h_stat'],  inplace = True, ignore_index = True, ascending = False)

cut_95 = round(len(sf_out)*0.95)-1
cut_99 = round(len(sf_out)*0.99)-1
cut_999 = round(len(sf_out)*0.999)-1

out_cut_95 = H_out.iloc[cut_95]['h_stat']
out_cut_99 = H_out.iloc[cut_99]['h_stat']
out_cut_999 = H_out.iloc[cut_999]['h_stat']

print("Significance threshold 0.05: " + str(out_cut_95), end = '\n')
print("Significance threshold 0.01: " + str(out_cut_99), end = '\n')
print("Significance threshold 0.001: " + str(out_cut_999), end = '\n')
