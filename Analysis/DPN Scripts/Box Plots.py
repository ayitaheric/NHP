# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 12:21:41 2024

@author: dylan

Scipt that creates boxplots for various parameters
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
from scipy import stats

plt.close('all')

folder_dir = r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\Hail_data_MATEUSZ/'
df = pd.read_csv(folder_dir + 'Chosen Profiles.csv', encoding='unicode_escape')
var = pd.read_csv(r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\Data\Variables.csv', encoding='unicode_escape')

medians_sc = []
medians_nsc = []
means_sc = []
means_nsc = []
std_sc = []
std_nsc = []
whiskers_low_sc = []
whiskers_high_sc = []
whiskers_low_nsc = []
whiskers_high_nsc = []
ts = []
ps = []

for j in range(len(var.Variable)):
    fig, ax = plt.subplots(1, 2, sharey='row', figsize=(10,8))
    
    if var.Units[j] == 'm/s':
        unit = 'm s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}'
    elif var.Units[j] == 'J/kg':
        unit = 'J kg\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}'
    elif var.Units[j] == 'm2/s2':
        unit = 'm\N{SUPERSCRIPT TWO} s\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT TWO}'
    elif var.Units[j] == 'g/kg':
        unit = 'g kg\N{SUPERSCRIPT MINUS}\N{SUPERSCRIPT ONE}'
    elif var.Units[j] == 'm AGL':
        unit = 'm AGL'
    elif var.Units[j] == 'mm':
        unit = 'mm'
    elif var.Units[j] == 'm':
        unit = 'm'
    
    df2 = df[df['Storm Type'] == 'Y']
    if var.Variable[j] == 'SHIP':
        pass
    else:
        ax[0].set_ylabel(f'{unit}')
    ax[0].set_ylim(var.Lower[j],var.Upper[j])
    ax[0].set_title('Supercell',size=12,weight='bold')
    bp1=ax[0].boxplot(df2[f'{var.Variable[j]}'], whis=[10, 90], showmeans=True, meanline=True, showfliers=False)
    ax[0].set_xticks([])
    ax[0].grid(axis = "y", alpha=0.4)
    ax[0].legend([bp1['medians'][0], bp1['means'][0],Patch(facecolor = 'None')], [f'Median: {np.format_float_positional(bp1["medians"][0].get_ydata()[1], precision=4, unique=False, fractional=False, trim="k")}', f'Mean: {np.format_float_positional(bp1["means"][0].get_ydata()[1], precision=4, unique=False, fractional=False, trim="k")}',f'StDev: {np.format_float_positional(df2[f"{var.Variable[j]}"].std(), precision=4, unique=False, fractional=False, trim="k")}'], loc='upper right')
    ax[0].plot([1]*len(df2[f'{var.Variable[j]}']), df2[f'{var.Variable[j]}'], 'o', alpha=0.5, markeredgecolor = 'blue', mfc='none')
    
    df3 = df[df['Storm Type'] == 'N']
    ax[1].set_title('Non Supercell',size=12,weight='bold')
    bp2=ax[1].boxplot(df3[f'{var.Variable[j]}'], whis=[10, 90], showmeans=True, meanline=True, showfliers=False)
    ax[1].set_xticks([])
    ax[1].grid(axis = "y", alpha=0.4)
    ax[1].legend([bp2['medians'][0], bp2['means'][0],Patch(facecolor = 'None')], [f'Median: {np.format_float_positional(bp2["medians"][0].get_ydata()[1], precision=4, unique=False, fractional=False, trim="k")}', f'Mean: {np.format_float_positional(bp2["means"][0].get_ydata()[1], precision=4, unique=False, fractional=False, trim="k")}',f'StDev: {np.format_float_positional(df3[f"{var.Variable[j]}"].std(), precision=4, unique=False, fractional=False, trim="k")}'], loc='upper right')
    ax[1].plot([1]*len(df3[f'{var.Variable[j]}']), df3[f'{var.Variable[j]}'], 'o', alpha=0.5, markeredgecolor = 'blue', mfc='none')
    
    
    medians_sc.append(bp1['medians'][0].get_ydata()[1])
    medians_nsc.append(bp2['medians'][0].get_ydata()[1])
    means_sc.append(bp1['means'][0].get_ydata()[1])
    means_nsc.append(bp2['means'][0].get_ydata()[1])
    std_sc.append(np.format_float_positional(df2[f"{var.Variable[j]}"].std()))
    std_nsc.append(np.format_float_positional(df3[f"{var.Variable[j]}"].std()))
    whiskers_low_sc.append(bp1['whiskers'][0].get_ydata()[1])
    whiskers_low_nsc.append(bp2['whiskers'][0].get_ydata()[1])
    whiskers_high_sc.append(bp1['whiskers'][1].get_ydata()[1])
    whiskers_high_nsc.append(bp2['whiskers'][1].get_ydata()[1])
    
    ts.append('{:.9f}'.format(stats.ttest_ind(df2[f'{var.Variable[j]}'],df3[f'{var.Variable[j]}'])[0]))
    ps.append('{:.9f}'.format(stats.ttest_ind(df2[f'{var.Variable[j]}'],df3[f'{var.Variable[j]}'])[1]))
    #'{:f}'.format(a)
    """
    outliers = bp2["fliers"][0].get_ydata()
    for w in outliers:
        print(df3['Event Number'][df3[f'{j}'][df3[f'{j}'] == w].index.tolist()[0]])
    print('######',f'{j}')
    """
    fig.suptitle(f'{var.Name[var[var["Variable"]==f"{var.Variable[j]}"].index.values[0]]}')
    
    #plt.savefig(r'C:\Users\dylan\Downloads\NTP Related Documents\2023-2024 Research Projects\NHP Hail Work\Plots\Box Plots/' + f'{var.Variable[j]}.png',dpi=300,bbox_inches='tight')

data_sc = {'Variable': var.Variable, 'Median':medians_sc, 'Mean':means_sc, 'Standard Deviation': std_sc,'Whiskers Lower':whiskers_low_sc, 'Whisker Higher':whiskers_high_sc}
data_nsc = {'Variable': var.Variable, 'Median':medians_nsc, 'Mean':means_nsc, 'Standard Deviation': std_nsc,'Whiskers Lower':whiskers_low_nsc, 'Whisker Higher':whiskers_high_nsc}
df_sc = pd.DataFrame(data=data_sc)
df_nsc = pd.DataFrame(data=data_nsc)

#df_sc.to_csv(r'C:\Users\dylan\Downloads\NHP Hail Work\Data\Stats_SC.csv', index=False)
#df_nsc.to_csv(r'C:\Users\dylan\Downloads\NHP Hail Work\Data\Stats_NSC.csv', index=False)

ttest_data = {'Variable': var.Variable, 'T-test Statistic': ts, 'p value': ps}
df_ttest = pd.DataFrame(data=ttest_data)
#df_ttest.to_csv(r'C:\Users\dylan\Downloads\NHP Hail Work\Data\ttest.csv', index=False)