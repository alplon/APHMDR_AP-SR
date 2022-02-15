#Analysis of data from Analysis Output

#Best run in Jupyter Notebook as blocks. Much of the script is repeated code with slight variations to generate different figures (distance vs rmsd vs lower bulb...). 
# %%
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
from matplotlib import rc
import numpy as np
import scipy

rc('font',**{'family':'sans-serif','sans-serif':['DejaVu Sans'],'size':12})
rc('mathtext',**{'default':'regular'})

# %%

#plot non-normalized distance comparisons
data = pd.read_csv('/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/distance_comparisons_all.csv')

plt.figure(figsize=(15, 5))
plt.subplots_adjust(left=0.05, bottom=0.25, right=0.99, top=0.96)

#sns.violinplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Exp', 'Af-Af', 'Exp-Af'], palette = {'Exp-Exp':'tab:blue','Af-Af':'white', 'Exp-Af':'tab:orange'}, scale='count' ,order=['T1029', 'T1027', 'T1030', 'T1032', 'T1040', 'T1099', 'T1039', 'T1073', 'T1043', 'T1080', 'T1050', 'T1038', 'T1033', 'T1031', 'T1042', 'T1067', 'T1090', 'T1041', 'T1037', 'T1074', 'T1079', 'T1054', 'T1049', 'T1026', 'T1082', 'T1046s2', 'T1035', 'T1046s1', 'T1056'])
sns.boxplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Exp', 'Af-Af', 'Exp-Af'], palette = {'Exp-Exp':'darkseagreen','Af-Af':'salmon', 'Exp-Af':'cornflowerblue'}, order=['45.20\n(T1029)', '61.11\n(T1027)', '63.55\n(T1030)', '69.12\n(T1032)', '72.11\n(T1040)', '80.34\n(T1099)', '82.30\n(T1039)', '83.90\n(T1073)', '84.12\n(T1043)', '85.90\n(T1080)', '86.28\n(T1050)', '87.37\n(T1038)', '87.50\n(T1033)', '87.63\n(T1031)', '89.49\n(T1042)', '90.39\n(T1067)', '90.48\n(T1090)', '90.70\n(T1041)', '90.72\n(T1037)', '91.10\n(T1074)', '92.18\n(T1079)', '92.66\n(T1054)', '93.28\n(T1049)', '93.84\n(T1026)', '96.00\n(T1082)', '96.45\n(T1046s2)', '96.81\n(T1035)', '97.57\n(T1046s1)', '98.08\n(T1056)'])
#sns.boxplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Af'], palette = {'Exp-Exp':'cornflowerblue','Af-Af':'salmon', 'Exp-Af':'darkseagreen'}, order=['45.20\n(T1029)', '61.11\n(T1027)', '63.55\n(T1030)', '69.12\n(T1032)', '72.11\n(T1040)', '80.34\n(T1099)', '82.30\n(T1039)', '83.90\n(T1073)', '84.12\n(T1043)', '85.90\n(T1080)', '86.28\n(T1050)', '87.37\n(T1038)', '87.50\n(T1033)', '87.63\n(T1031)', '89.49\n(T1042)', '90.39\n(T1067)', '90.48\n(T1090)', '90.70\n(T1041)', '90.72\n(T1037)', '91.10\n(T1074)', '92.18\n(T1079)', '92.66\n(T1054)', '93.28\n(T1049)', '93.84\n(T1026)', '96.00\n(T1082)', '96.45\n(T1046s2)', '96.81\n(T1035)', '97.57\n(T1046s1)', '98.08\n(T1056)'])
#sns.violinplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Af'], palette = {'Exp-Exp':'darkseagreen','Af-Af':'salmon', 'Exp-Af':'cornflowerblue'}, scale='count' ,order=['45.20\n(T1029)', '61.11\n(T1027)', '63.55\n(T1030)', '69.12\n(T1032)', '72.11\n(T1040)', '80.34\n(T1099)', '82.30\n(T1039)', '83.90\n(T1073)', '84.12\n(T1043)', '85.90\n(T1080)', '86.28\n(T1050)', '87.37\n(T1038)', '87.50\n(T1033)', '87.63\n(T1031)', '89.49\n(T1042)', '90.39\n(T1067)', '90.48\n(T1090)', '90.70\n(T1041)', '90.72\n(T1037)', '91.10\n(T1074)', '92.18\n(T1079)', '92.66\n(T1054)', '93.28\n(T1049)', '93.84\n(T1026)', '96.00\n(T1082)', '96.45\n(T1046s2)', '96.81\n(T1035)', '97.57\n(T1046s1)', '98.08\n(T1056)'])

plt.xticks(rotation = 90)
plt.xlabel('GDT_TS (Protein)')
plt.legend().set_visible(False)
#plt.savefig('/Users/alexplonski/Desktop/AlphaFoldProject/temp_figures/violin_non-normalized_cross.png', dpi=300)
plt.show()

# %%

#plot normalized distance comparisons
data = pd.read_csv('/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/distance_comparisons_all_normalized.csv')

plt.figure(figsize=(15, 5))
plt.subplots_adjust(left=0.05, bottom=0.22, right=0.99, top=0.96)

#sns.violinplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Exp', 'Af-Af', 'Exp-Af'], palette = {'Exp-Exp':'tab:blue','Af-Af':'white', 'Exp-Af':'tab:orange'}, scale='count' ,order=['T1029', 'T1027', 'T1030', 'T1032', 'T1040', 'T1099', 'T1039', 'T1073', 'T1043', 'T1080', 'T1050', 'T1038', 'T1033', 'T1031', 'T1042', 'T1067', 'T1090', 'T1041', 'T1037', 'T1074', 'T1079', 'T1054', 'T1049', 'T1026', 'T1082', 'T1046s2', 'T1035', 'T1046s1', 'T1056'])
#sns.boxplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Exp', 'Af-Af', 'Exp-Af'], palette = {'Exp-Exp':'cornflowerblue','Af-Af':'salmon', 'Exp-Af':'darkseagreen'}, order=['45.20\n(T1029)', '61.11\n(T1027)', '63.55\n(T1030)', '69.12\n(T1032)', '72.11\n(T1040)', '80.34\n(T1099)', '82.30\n(T1039)', '83.90\n(T1073)', '84.12\n(T1043)', '85.90\n(T1080)', '86.28\n(T1050)', '87.37\n(T1038)', '87.50\n(T1033)', '87.63\n(T1031)', '89.49\n(T1042)', '90.39\n(T1067)', '90.48\n(T1090)', '90.70\n(T1041)', '90.72\n(T1037)', '91.10\n(T1074)', '92.18\n(T1079)', '92.66\n(T1054)', '93.28\n(T1049)', '93.84\n(T1026)', '96.00\n(T1082)', '96.45\n(T1046s2)', '96.81\n(T1035)', '97.57\n(T1046s1)', '98.08\n(T1056)'])
sns.violinplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Af'], palette = {'Exp-Exp':'cornflowerblue','Af-Af':'salmon', 'Exp-Af':'cornflowerblue'}, scale='count' ,order=['45.20\n(T1029)', '61.11\n(T1027)', '63.55\n(T1030)', '69.12\n(T1032)', '72.11\n(T1040)', '80.34\n(T1099)', '82.30\n(T1039)', '83.90\n(T1073)', '84.12\n(T1043)', '85.90\n(T1080)', '86.28\n(T1050)', '87.37\n(T1038)', '87.50\n(T1033)', '87.63\n(T1031)', '89.49\n(T1042)', '90.39\n(T1067)', '90.48\n(T1090)', '90.70\n(T1041)', '90.72\n(T1037)', '91.10\n(T1074)', '92.18\n(T1079)', '92.66\n(T1054)', '93.28\n(T1049)', '93.84\n(T1026)', '96.00\n(T1082)', '96.45\n(T1046s2)', '96.81\n(T1035)', '97.57\n(T1046s1)', '98.08\n(T1056)'])

plt.xticks(rotation = 90)
plt.xlabel('GDT_TS (Protein)')
plt.ylabel('Normalized Distance Change')
#plt.ylim((0,10))
plt.savefig('/Users/alexplonski/Desktop/AlphaFoldProject/temp_figures/violinplot_normalized.png',bbox_inches = 'tight' ,dpi=300)
plt.show()

# %%
#plot lower bulb, for error bars to work, plot with protein instead of GDT
cross_data = pd.read_csv('/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/mean_based_lower_bulb/bulb_stats/mw_bulbstats_Exp-Af.csv')
exp_data = pd.read_csv('/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/mean_based_lower_bulb/bulb_stats/mw_bulbstats_Exp-Exp.csv')
af_data = pd.read_csv('/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/mean_based_lower_bulb/bulb_stats/mw_bulbstats_Af-Af.csv')

plt.figure(figsize=(5, 5))
plt.subplots_adjust(left=0.05, bottom=0.22, right=0.99, top=0.96)

plt.scatter(cross_data['MW'], cross_data['Mean'], label = 'Exp-Af', color = 'cornflowerblue', marker='s',zorder=3)
plt.errorbar(cross_data['MW'], cross_data['Mean'], yerr=cross_data['StDev'], ls='none', color = 'cornflowerblue')
slope, intercept = np.polyfit(cross_data['MW'], cross_data['Mean'], 1)
abline_values = [slope * i + intercept for i in cross_data['MW']]
plt.plot(cross_data['MW'], abline_values, 'b', color = 'cornflowerblue')

plt.scatter(exp_data['MW'], exp_data['Mean'], label = 'Exp-Exp', color = 'darkseagreen', marker='D',zorder=2)
plt.errorbar(exp_data['MW'], exp_data['Mean'], yerr=exp_data['StDev'], ls='none',color = 'darkseagreen')
slope, intercept = np.polyfit(exp_data['MW'], exp_data['Mean'], 1)
abline_values = [slope * i + intercept for i in exp_data['MW']]
plt.plot(exp_data['MW'], abline_values, 'b', color = 'darkseagreen')
pearson = scipy.stats.pearsonr(af_data['MW'], af_data['Mean'])
print(pearson)

plt.scatter(af_data['MW'], af_data['Mean'], label = 'Af-Af', color ='salmon')
plt.errorbar(af_data['MW'], af_data['Mean'], yerr=af_data['StDev'], ls='none', color ='salmon', marker='o',zorder=1)
slope, intercept = np.polyfit(af_data['MW'], af_data['Mean'], 1)
abline_values = [slope * i + intercept for i in af_data['MW']]
plt.plot(af_data['MW'], abline_values, 'b', color = 'salmon')





plt.xlabel('Moelcular Weight [kDa]')
plt.ylabel(r'$\phi$ ')
#plt.xticks(rotation = 90)
plt.legend()

plt.savefig('/Users/alexplonski/Desktop/AlphaFoldProject/temp_figures/mw_bulb_trialcomparison_square.png', bbox_inches = 'tight', dpi=300)
plt.show()

# %%

#plot non-normalized distance comparisons rot bonds 
data = pd.read_csv('/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/distance_rb5.csv')

plt.figure(figsize=(15, 5))
plt.subplots_adjust(left=0.05, bottom=0.25, right=0.99, top=0.96)

#sns.violinplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Exp', 'Af-Af', 'Exp-Af'], palette = {'Exp-Exp':'tab:blue','Af-Af':'white', 'Exp-Af':'tab:orange'}, scale='count' ,order=['T1029', 'T1027', 'T1030', 'T1032', 'T1040', 'T1099', 'T1039', 'T1073', 'T1043', 'T1080', 'T1050', 'T1038', 'T1033', 'T1031', 'T1042', 'T1067', 'T1090', 'T1041', 'T1037', 'T1074', 'T1079', 'T1054', 'T1049', 'T1026', 'T1082', 'T1046s2', 'T1035', 'T1046s1', 'T1056'])
sns.boxplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Exp', 'Af-Af', 'Exp-Af'], palette = {'Exp-Exp':'darkseagreen','Af-Af':'salmon', 'Exp-Af':'cornflowerblue'}, order=['45.20\n(T1029)', '61.11\n(T1027)', '63.55\n(T1030)', '69.12\n(T1032)', '72.11\n(T1040)', '80.34\n(T1099)', '82.30\n(T1039)', '83.90\n(T1073)', '84.12\n(T1043)', '85.90\n(T1080)', '86.28\n(T1050)', '87.37\n(T1038)', '87.50\n(T1033)', '87.63\n(T1031)', '89.49\n(T1042)', '90.39\n(T1067)', '90.48\n(T1090)', '90.70\n(T1041)', '90.72\n(T1037)', '91.10\n(T1074)', '92.18\n(T1079)', '92.66\n(T1054)', '93.28\n(T1049)', '93.84\n(T1026)', '96.00\n(T1082)', '96.45\n(T1046s2)', '96.81\n(T1035)', '97.57\n(T1046s1)', '98.08\n(T1056)'])
#sns.boxplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Af'], palette = {'Exp-Exp':'cornflowerblue','Af-Af':'salmon', 'Exp-Af':'darkseagreen'}, order=['45.20\n(T1029)', '61.11\n(T1027)', '63.55\n(T1030)', '69.12\n(T1032)', '72.11\n(T1040)', '80.34\n(T1099)', '82.30\n(T1039)', '83.90\n(T1073)', '84.12\n(T1043)', '85.90\n(T1080)', '86.28\n(T1050)', '87.37\n(T1038)', '87.50\n(T1033)', '87.63\n(T1031)', '89.49\n(T1042)', '90.39\n(T1067)', '90.48\n(T1090)', '90.70\n(T1041)', '90.72\n(T1037)', '91.10\n(T1074)', '92.18\n(T1079)', '92.66\n(T1054)', '93.28\n(T1049)', '93.84\n(T1026)', '96.00\n(T1082)', '96.45\n(T1046s2)', '96.81\n(T1035)', '97.57\n(T1046s1)', '98.08\n(T1056)'])
#sns.violinplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Af'], palette = {'Exp-Exp':'cornflowerblue','Af-Af':'salmon', 'Exp-Af':'darkseagreen'}, scale='count' ,order=['45.20\n(T1029)', '61.11\n(T1027)', '63.55\n(T1030)', '69.12\n(T1032)', '72.11\n(T1040)', '80.34\n(T1099)', '82.30\n(T1039)', '83.90\n(T1073)', '84.12\n(T1043)', '85.90\n(T1080)', '86.28\n(T1050)', '87.37\n(T1038)', '87.50\n(T1033)', '87.63\n(T1031)', '89.49\n(T1042)', '90.39\n(T1067)', '90.48\n(T1090)', '90.70\n(T1041)', '90.72\n(T1037)', '91.10\n(T1074)', '92.18\n(T1079)', '92.66\n(T1054)', '93.28\n(T1049)', '93.84\n(T1026)', '96.00\n(T1082)', '96.45\n(T1046s2)', '96.81\n(T1035)', '97.57\n(T1046s1)', '98.08\n(T1056)'])

plt.xticks(rotation = 90)
plt.xlabel('GDT_TS (Protein)')
plt.savefig('/Users/alexplonski/Desktop/AlphaFoldProject/temp_figures/boxplot_non-normalized_rb5.png', dpi=300)
plt.show()

# %%
cross_data = pd.read_csv('/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/mean_rb_stats/rb_5/gdt_meanstats_rb5_Exp-Af.csv')
exp_data = pd.read_csv('/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/mean_rb_stats/rb_5/gdt_meanstats_rb5_Exp-Exp.csv')
af_data = pd.read_csv('/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/mean_rb_stats/rb_5/gdt_meanstats_rb5_Af-Af.csv')

plt.figure(figsize=(5, 5))
plt.subplots_adjust(left=0.05, bottom=0.22, right=0.99, top=0.96)

plt.scatter(cross_data['GDT_TS'], cross_data['Mean'], label = 'Exp-Af', color = 'cornflowerblue', marker='s',zorder=3)
plt.errorbar(cross_data['GDT_TS'], cross_data['Mean'], yerr=cross_data['StDev'], ls='none', color = 'cornflowerblue')
slope, intercept = np.polyfit(cross_data['GDT_TS'], cross_data['Mean'], 1)
abline_values = [slope * i + intercept for i in cross_data['GDT_TS']]
plt.plot(cross_data['GDT_TS'], abline_values, 'b', color = 'cornflowerblue', zorder=3)

plt.scatter(exp_data['GDT_TS'], exp_data['Mean'], label = 'Exp-Exp', color = 'darkseagreen', marker='D',zorder=1, s = 50)
plt.errorbar(exp_data['GDT_TS'], exp_data['Mean'], yerr=exp_data['StDev'], ls='none',color = 'darkseagreen')
slope, intercept = np.polyfit(exp_data['GDT_TS'], exp_data['Mean'], 1)
abline_values = [slope * i + intercept for i in exp_data['GDT_TS']]
plt.plot(exp_data['GDT_TS'], abline_values, 'b', color = 'darkseagreen', zorder=1)

plt.scatter(af_data['GDT_TS'], af_data['Mean'], label = 'Af-Af', color ='salmon', marker='o',zorder=2)
plt.errorbar(af_data['GDT_TS'], af_data['Mean'], yerr=af_data['StDev'], ls='none', color ='salmon')
slope, intercept = np.polyfit(af_data['GDT_TS'], af_data['Mean'], 1)
abline_values = [slope * i + intercept for i in af_data['GDT_TS']]
plt.plot(af_data['GDT_TS'], abline_values, 'b', color = 'salmon', zorder=2)
print(slope, intercept)

pearson = scipy.stats.pearsonr(af_data['GDT_TS'], af_data['Mean'])
print(pearson)



plt.xlabel('GDT_TS')
plt.ylabel('Distance Change [$\AA$]')
#plt.xticks(rotation = 90)
plt.legend()

plt.savefig('/Users/alexplonski/Desktop/AlphaFoldProject/temp_figures/median_comparison_rb5.png',bbox_inches = 'tight', dpi=300)
plt.show()

# %%
#plot median difference (cross - self(s))
exp_data = pd.read_csv('/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/median_rb_stats/rb_5/gdt_median_cross-exp.csv')
af_data = pd.read_csv('/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/median_rb_stats/rb_5/gdt_median_cross-af.csv')

plt.figure(figsize=(15, 5))
plt.subplots_adjust(left=0.05, bottom=0.22, right=0.99, top=0.96)



plt.scatter(exp_data['GDT_TS'], exp_data['Median'], label = 'Exp-Exp', color = 'darkseagreen')
slope, intercept = np.polyfit(exp_data['GDT_TS'], exp_data['Median'], 1)
abline_values = [slope * i + intercept for i in exp_data['GDT_TS']]
plt.plot(exp_data['GDT_TS'], abline_values, 'b', color = 'darkseagreen')

plt.scatter(af_data['GDT_TS'], af_data['Median'], label = 'Af-Af', color = 'salmon')
slope, intercept = np.polyfit(af_data['GDT_TS'], af_data['Median'], 1)
abline_values = [slope * i + intercept for i in af_data['GDT_TS']]
plt.plot(af_data['GDT_TS'], abline_values, 'b', color = 'salmon')

plt.xlabel('GDT_TS')
plt.ylabel(r'$\phi$ ')
plt.xticks(rotation = 90)
plt.legend()

plt.savefig('/Users/alexplonski/Desktop/AlphaFoldProject/PUBLICATION/trial_reproducibility/median_difference_5rb_comparison.png', dpi=300)
plt.show()

# %%


#plot non-normalized distance comparisons
data = pd.read_csv('/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/rmsd_rb5.csv')

plt.figure(figsize=(15, 5))
plt.subplots_adjust(left=0.05, bottom=0.25, right=0.99, top=0.96)

#sns.violinplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Exp', 'Af-Af', 'Exp-Af'], palette = {'Exp-Exp':'tab:blue','Af-Af':'white', 'Exp-Af':'tab:orange'}, scale='count' ,order=['T1029', 'T1027', 'T1030', 'T1032', 'T1040', 'T1099', 'T1039', 'T1073', 'T1043', 'T1080', 'T1050', 'T1038', 'T1033', 'T1031', 'T1042', 'T1067', 'T1090', 'T1041', 'T1037', 'T1074', 'T1079', 'T1054', 'T1049', 'T1026', 'T1082', 'T1046s2', 'T1035', 'T1046s1', 'T1056'])
sns.boxplot(x='Protein', y='RMSD [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Exp', 'Af-Af', 'Exp-Af'], palette = {'Exp-Exp':'darkseagreen','Af-Af':'salmon', 'Exp-Af':'cornflowerblue'}, order=['45.20\n(T1029)', '61.11\n(T1027)', '63.55\n(T1030)', '69.12\n(T1032)', '72.11\n(T1040)', '80.34\n(T1099)', '82.30\n(T1039)', '83.90\n(T1073)', '84.12\n(T1043)', '85.90\n(T1080)', '86.28\n(T1050)', '87.37\n(T1038)', '87.50\n(T1033)', '87.63\n(T1031)', '89.49\n(T1042)', '90.39\n(T1067)', '90.48\n(T1090)', '90.70\n(T1041)', '90.72\n(T1037)', '91.10\n(T1074)', '92.18\n(T1079)', '92.66\n(T1054)', '93.28\n(T1049)', '93.84\n(T1026)', '96.00\n(T1082)', '96.45\n(T1046s2)', '96.81\n(T1035)', '97.57\n(T1046s1)', '98.08\n(T1056)'])
#sns.boxplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Af'], palette = {'Exp-Exp':'cornflowerblue','Af-Af':'salmon', 'Exp-Af':'darkseagreen'}, order=['45.20\n(T1029)', '61.11\n(T1027)', '63.55\n(T1030)', '69.12\n(T1032)', '72.11\n(T1040)', '80.34\n(T1099)', '82.30\n(T1039)', '83.90\n(T1073)', '84.12\n(T1043)', '85.90\n(T1080)', '86.28\n(T1050)', '87.37\n(T1038)', '87.50\n(T1033)', '87.63\n(T1031)', '89.49\n(T1042)', '90.39\n(T1067)', '90.48\n(T1090)', '90.70\n(T1041)', '90.72\n(T1037)', '91.10\n(T1074)', '92.18\n(T1079)', '92.66\n(T1054)', '93.28\n(T1049)', '93.84\n(T1026)', '96.00\n(T1082)', '96.45\n(T1046s2)', '96.81\n(T1035)', '97.57\n(T1046s1)', '98.08\n(T1056)'])
#sns.violinplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Comparison',hue_order= ['Exp-Af'], palette = {'Exp-Exp':'darkseagreen','Af-Af':'salmon', 'Exp-Af':'cornflowerblue'}, scale='count' ,order=['45.20\n(T1029)', '61.11\n(T1027)', '63.55\n(T1030)', '69.12\n(T1032)', '72.11\n(T1040)', '80.34\n(T1099)', '82.30\n(T1039)', '83.90\n(T1073)', '84.12\n(T1043)', '85.90\n(T1080)', '86.28\n(T1050)', '87.37\n(T1038)', '87.50\n(T1033)', '87.63\n(T1031)', '89.49\n(T1042)', '90.39\n(T1067)', '90.48\n(T1090)', '90.70\n(T1041)', '90.72\n(T1037)', '91.10\n(T1074)', '92.18\n(T1079)', '92.66\n(T1054)', '93.28\n(T1049)', '93.84\n(T1026)', '96.00\n(T1082)', '96.45\n(T1046s2)', '96.81\n(T1035)', '97.57\n(T1046s1)', '98.08\n(T1056)'])

plt.xticks(rotation = 90)
plt.xlabel('GDT_TS (Protein)')
plt.legend().set_visible(False)
plt.savefig('/Users/alexplonski/Desktop/AlphaFoldProject/temp_figures/boxplot_non-normalized_rb5.png', dpi=300)
plt.show()

# %%

#plot non-normalized distance comparisons by trial
data = pd.read_csv('/Users/alexplonski/Desktop/AlphaFoldProject/CASP-14/Docking_Results/Distance/trial_reproducibility/output_csv_data/distance_comparisons_all_normalized_trialinfo.csv')

plt.figure(figsize=(15, 5))
plt.subplots_adjust(left=0.05, bottom=0.25, right=0.99, top=0.96)

sns.boxplot(x='Protein', y='Distance Change [$\AA$]', data = data, hue = 'Trial',hue_order= ['Af-Af1-2', 'Af-Af1-3', 'Af-Af2-3'], palette = {'Af-Af1-2':'salmon','Af-Af1-3':'salmon', 'Af-Af2-3':'salmon'}, order=['45.20\n(T1029)', '61.11\n(T1027)', '63.55\n(T1030)', '69.12\n(T1032)', '72.11\n(T1040)', '80.34\n(T1099)', '82.30\n(T1039)', '83.90\n(T1073)', '84.12\n(T1043)', '85.90\n(T1080)', '86.28\n(T1050)', '87.37\n(T1038)', '87.50\n(T1033)', '87.63\n(T1031)', '89.49\n(T1042)', '90.39\n(T1067)', '90.48\n(T1090)', '90.70\n(T1041)', '90.72\n(T1037)', '91.10\n(T1074)', '92.18\n(T1079)', '92.66\n(T1054)', '93.28\n(T1049)', '93.84\n(T1026)', '96.00\n(T1082)', '96.45\n(T1046s2)', '96.81\n(T1035)', '97.57\n(T1046s1)', '98.08\n(T1056)'])

plt.xticks(rotation = 90)
plt.xlabel('GDT_TS (Protein)')
plt.legend().set_visible(False)
#plt.legend(loc = 1)
plt.savefig('/Users/alexplonski/Desktop/AlphaFoldProject/temp_figures/trialseparated_af_boxplot_non-normalized_cross.png', dpi=300)
plt.show()

# %%



