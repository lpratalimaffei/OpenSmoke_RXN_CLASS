# imports
from pylab import meshgrid, cm, imshow, contour, clabel, colorbar, axis, title, show
import pandas as pd
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import copy
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
mechname = 'POLI_2202_CLASSES_TOTB'
kin_xml_fld = r'C:\Users\lpratalimaffei\Desktop\OpenSmoke_RXN_CLASS\{}'.format(mechname)
class_groups_fld = r'C:\Users\lpratalimaffei\Desktop\OpenSmoke_RXN_CLASS'


simul_flds ={
    'JSR-WANG-PHI1': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\JSR_WANG_2019\PHI_1\000_flux\{}'.format(mechname),
   # 'JSR-WANG-PHI05': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\JSR_WANG_2019\PHI_0.5\000_flux\{}'.format(mechname),
   # 'JSR-WANG-PHI18': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\JSR_WANG_2019\PHI_1.8\000_flux\{}'.format(mechname),
    }
simul_flds ={
   # 'PYR-KIM2010': r'D:\PhD\OPENSMOKE\CPD\PYROLYSIS\PFR_Kim2010\flux\{}'.format(mechname),
    'PYR-KIM2010-1500': r'D:\PhD\OPENSMOKE\CPD\PYROLYSIS\PFR_Kim2010\flux_1500\{}'.format(mechname),
    'PYR-DJOKICLD2014': r'D:\PhD\OPENSMOKE\CPD\PYROLYSIS\PFR_Djokic2014\lowD\1073\{}'.format(mechname),
    'PYR-C2H42019': r'D:\PhD\OPENSMOKE\CPD\PYROLYSIS\PFR_C2H4_2019\PFR5\{}'.format(mechname),
    }
simul_flds ={
#    'ST-IND-LASKIN': r'D:\PhD\OPENSMOKE\INDENE\ST_LASKIN_1998\000_flux\{}'.format(mechname),
    'PFR760-IND-JIN': r'D:\PhD\OPENSMOKE\INDENE\PFR_JIN_2019\760Torr\1255K\{}'.format(mechname),
#    'PFR30-IND-JIN': r'D:\PhD\OPENSMOKE\INDENE\PFR_JIN_2019\30Torr\1394K\{}'.format(mechname),
    }
# species_list = ['C5H6','C5H5CH3','FULVENE','C6H6','C6H5OH','C5H5','C6H5','C6H5O']
species_list = ['C5H6','C5H5','C6H6','FULVENE','C7H8','INDENE','C10H8','C6H5C2H','C6H5C2H3']
species_list = ['C5H6','C5H5','C6H6','C6H5C2H3','INDENE','C10H8','FLUORENE','C14H10']
#species_list = ['C10H8','C6H5C2H3','C14H10','FLUORENE']
species_list = ['INDENE','INDENYL','C5H6','C5H5','C6H6','C7H8','C6H5C2H','C6H5C2H3','C10H8','C12H8','FLUORENE','C14H10','C16H10','C18H14']

sortlists = [['classgroup'],['speciestype','subclass'],['bimoltype']] # classgroup, speciestype, subclass, bimoltype (R+R, RSR+RSR, M+M, ETC) # sum if both apply and sort by this criteria
filter_dcts = [{},{},{'classgroup': ['GROWTH', 'ENLARGE']}] # filter according to selected criteria in name
threshs = [1e-3, 1e-1, 1e-3]
# filter_dct = {'classgroup': ['GROWTH', 'ENLARGE']} # filter according to selected criteria in name
# FAI CRITERIO PER SCRIVERE TUTTE LE REAZIONI APPARTENENTI ALLA CLASSE SELEZIONATA
# parse mech
kinetics = KineticMechanism(os.path.join(kin_xml_fld, 'kinetics.xml'))
kinetics.ReadKinetics(os.path.join(kin_xml_fld, 'reaction_names.xml'))
reactions = kinetics.ProcessingReactions()
# parse classes
_, subcl_grp_dct = ReadRxnGroups(
    class_groups_fld, 'rxn_class_groups.txt')

# create dataframe for classes

for simul_name, simul_fld in simul_flds.items():
    print('processing simul {}'.format(simul_fld))
    rxns_sorted = rxnclass(reactions)
    rxns_sorted.assign_class_grp(subcl_grp_dct)
    # simul output
    for sp in species_list:
        ropa = RateOfProductionAnalysis(kin_xml_fld, simul_fld, sp)
        consumption, production = ropa.ComputeRopa(Type = 'Global')
        tot_rop = -np.array(consumption)+np.array(production)
        tot_rop_df = pd.DataFrame(tot_rop, index=np.arange(1,len(tot_rop)+1), columns=['flux_{}'.format(sp)], dtype=float)
        rxns_sorted.assign_flux(tot_rop_df)

    # filter rxns
    for i, sortlist in enumerate(sortlists):
        filter_dct = filter_dcts[i]
        rxns_sorted_i = copy.deepcopy(rxns_sorted)
        THRESH = threshs[i]
        if filter_dct:
            rxns_sorted_i.filter_class(filter_dct)
        # filter flux
        rxns_sorted_i.filter_flux(threshold=THRESH)
        # sum same speciestype-classgroup-subclass together
        sortdf = rxns_sorted_i.sortby(sortlist)
        # drop unsorted cols
        col_names = sortdf.columns
        for col in col_names:
            if 'UNSORTED' in col:
                sortdf = sortdf.drop(col, axis=1)
        if len(sortlist) > 1:
            criteria_str = '-'.join(sortlist)
        else:
            criteria_str = sortlist[0]
        plotpath = os.path.join(class_groups_fld, '{}_{}.png'.format(simul_name, criteria_str))
        plot_heatmap(sortdf, plotpath)

