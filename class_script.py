# imports
from pylab import meshgrid, cm, imshow, contour, clabel, colorbar, axis, title, show
import pandas as pd
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('default')

sys.path.insert(0, r'MECH_CLASS\parser')
sys.path.insert(0, r'MECH_CLASS\utils_automech')
sys.path.insert(0, r'ROPA')
sys.path.insert(0, r'MECH_CLASS')


from KineticMechanism import KineticMechanism
from RxnClassGroups import ReadRxnGroups
from rxnclass import rxnclass
from rxnclass import plot_heatmap
from ROPA import RateOfProductionAnalysis

# paths
kin_xml_fld = r'C:\Users\lpratalimaffei\Desktop\OpenSmoke_RXN_CLASS\example\kinetics'
class_groups_fld = r'C:\Users\lpratalimaffei\Desktop\OpenSmoke_RXN_CLASS\example'
simul_fld = r'C:\Users\lpratalimaffei\Desktop\OpenSmoke_RXN_CLASS\example\simul\Output'
species_list = ['C6H6','C6H5','C6H5OH','C6H5O','C5H6','C5H5','FULVENE']
sortlist = ['speciestype'] # classgroup, speciestype, subclass # sum if both apply and sort by this criteria
# parse mech
kinetics = KineticMechanism(os.path.join(kin_xml_fld, 'kinetics.xml'))
kinetics.ReadKinetics(os.path.join(kin_xml_fld, 'reaction_names.xml'))
reactions = kinetics.ProcessingReactions()
# parse classes
_, subcl_grp_dct = ReadRxnGroups(
    class_groups_fld, 'rxn_class_groups.txt')

# create dataframe for classes
rxns_sorted = rxnclass(reactions)
rxns_sorted.assign_class_grp(subcl_grp_dct)

# simul output
for sp in species_list:
    ropa = RateOfProductionAnalysis(kin_xml_fld, simul_fld, sp)
    consumption, production = ropa.ComputeRopa(Type = 'Global')
    tot_rop = -np.array(consumption)+np.array(production)
    tot_rop_df = pd.DataFrame(tot_rop, index=np.arange(1,len(tot_rop)+1), columns=['flux_{}'.format(sp)], dtype=float)
    rxns_sorted.assign_flux(tot_rop_df)

# filter flux
rxns_sorted.filter_flux(threshold=1e-2)
# sum same speciestype-classgroup-subclass together
sortdf = rxns_sorted.sortby(sortlist)
plotpath = os.path.join(class_groups_fld, 'plt_sort.png')
plot_heatmap(sortdf, plotpath)

sys.exit()
nonzero_P = np.nonzero(production)[0] 
nonzero_C = np.nonzero(consumption)[0]
ROPA_results = {'id_P': [], 'P_reactions': [], 'ROPA_Coeff_P': [],'id_C': [],
                'C_reactions': [], 'ROPA_Coeff_C': []}

if len(nonzero_P) != 0:
    for i in nonzero_P:
        ROPA_results['id_P'].append(i)
        ROPA_results['P_reactions'].append(ropa.reactionsNames[i])
        ROPA_results['ROPA_Coeff_P'].append(production[i])

if len(nonzero_C) != 0:
    for j in nonzero_C:

        ROPA_results['id_C'].append(j)
        ROPA_results['C_reactions'].append(ropa.reactionsNames[j])
        ROPA_results['ROPA_Coeff_C'].append(-consumption[j])
# print(ROPA_results)

ROPA_results_plot = {'Reactions': [], 'Coefficients': []}
ROPA_results_plot['Reactions'].append(ROPA_results['P_reactions'] + ROPA_results['C_reactions'])
ROPA_results_plot['Coefficients'].append(ROPA_results['ROPA_Coeff_P'] + ROPA_results['ROPA_Coeff_C'])

ROPA_results_plot['Reactions'] = ROPA_results_plot['Reactions'][0]
ROPA_results_plot['Coefficients'] = ROPA_results_plot['Coefficients'][0]

# print(ROPA_results_plot)


df = pd.DataFrame.from_dict(ROPA_results_plot)
df['Coefficients_abs'] = abs(df['Coefficients'])
df = df.sort_values(by=['Coefficients_abs'])

df_last_15 = df.iloc[-15:]
df_last_15

df = df_last_15
plt.rcParams["figure.figsize"]=20,20
plt.barh(df['Reactions'], df['Coefficients'], color = (df['Coefficients'] >=0.).map({True:'red', False:'blue'}))

plt.yticks(fontsize = 20)
plt.xticks(fontsize = 10)

plt.show()