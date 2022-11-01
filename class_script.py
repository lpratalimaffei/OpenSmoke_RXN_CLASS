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
mechname = 'POLI_2202_CLASSES_TOTN'
kin_xml_fld = r'C:\Users\lpratalimaffei\Desktop\OpenSmoke_RXN_CLASS\{}'.format(mechname)
class_groups_fld = r'C:\Users\lpratalimaffei\Desktop\OpenSmoke_RXN_CLASS'


simul_flds ={
   'PYR-KIM2010': r'D:\PhD\OPENSMOKE\CPD\PYROLYSIS\PFR_Kim2010\flux\{}'.format(mechname),
   # 'PYR-KIM2010-1500': r'D:\PhD\OPENSMOKE\CPD\PYROLYSIS\PFR_Kim2010\flux_1500\{}'.format(mechname),
    #'PYR-DJOKICLD2014': r'D:\PhD\OPENSMOKE\CPD\PYROLYSIS\PFR_Djokic2014\lowD\1073\{}'.format(mechname),
    #'PYR-C2H42019': r'D:\PhD\OPENSMOKE\CPD\PYROLYSIS\PFR_C2H4_2019\PFR5\{}'.format(mechname),
    }
species_list = ['C5H6','C5H5','C6H6','C6H5C2H3','INDENE','C10H8','FLUORENE','C14H10']
sortlists = [['classgroup'],['speciestype','subclass'],['bimoltype']] # classgroup, speciestype, subclass, bimoltype (R+R, RSR+RSR, M+M, ETC) # sum if both apply and sort by this criteria
filter_dcts = [{},{},{'classgroup': ['GROWTH', 'ENLARGE']}] # filter according to selected criteria in name
threshs = [1e-3, 1e-1, 1e-2]

simul_flds ={
    'ST-IND-LASKIN': r'D:\PhD\OPENSMOKE\INDENE\ST_LASKIN_1998\000_flux\{}'.format(mechname),
#    'PFR760-IND-JIN': r'D:\PhD\OPENSMOKE\INDENE\PFR_JIN_2019\760Torr\1255K\{}'.format(mechname),
#    'PFR30-IND-JIN': r'D:\PhD\OPENSMOKE\INDENE\PFR_JIN_2019\30Torr\1394K\{}'.format(mechname),
    }
sortlists = [['classgroup'],['speciestype','subclass']] # classgroup, speciestype, subclass, bimoltype (R+R, RSR+RSR, M+M, ETC) # sum if both apply and sort by this criteria
filter_dcts = [{},{'classgroup': ['GROWTH', 'ENLARGE']}] # filter according to selected criteria in name
threshs = [3e-2, 5e-2, 1e-2]
species_list = ['INDENE','INDENYL','C6H6','C6H5','C5H6','C5H5','C7H8','C6H5C2H','C10H8','C12H8','FLUORENE','C14H10','C16H10','C18H14']



simul_flds ={
 #   'ST-ORME-PHI0.5': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\ST_Orme\PHI_0.5_FLUX\{}'.format(mechname),
 # 'ST-ORME-PHI1': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\ST_Orme\PHI_1_FLUX\{}'.format(mechname),
    'ST-ORME-PHI2': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\ST_Orme\PHI_2_FLUX\{}'.format(mechname),
 #   'ST-ORME-PHI0.5_1550': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\ST_Orme\PHI_0.5_FLUX_1550K\{}'.format(mechname),
 # 'ST-ORME-PHI1_1550': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\ST_Orme\PHI_1_FLUX_1550K\{}'.format(mechname),
   'ST-ORME-PHI2_1550': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\ST_Orme\PHI_2_FLUX_1550K\{}'.format(mechname),
    }
simul_flds ={
    'LFS-PHI0.9': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\LFS\PHI_0.9\{}'.format(mechname),
    'LFS-PHI1.1': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\LFS\PHI_1.1\{}'.format(mechname),
    'LFS-PHI1.4': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\LFS\PHI_1.4\{}'.format(mechname),
    }
species_list = ['C5H6','C5H5']
sortlists = [['speciestype','subclass']] # classgroup, speciestype, subclass, bimoltype (R+R, RSR+RSR, M+M, ETC) # sum if both apply and sort by this criteria
#filter_dcts = [{'speciestype': ['C5'],'classgroup': ['ADD', 'DECO', 'ABS', 'IPSO','ISSION']}] # filter according to selected criteria in name
filter_dcts = [{},{}] 
threshs = [5e-2]


simul_flds ={
   #'PFR-BUTLER-PHI1.6': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\PFR_Buttler_2009\CASES\5\{}'.format(mechname),
   'PFR-BUTLER-PHI0.6': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\PFR_Buttler_2009\CASES\4\{}'.format(mechname),
      # 'JSR-WANG-PHI1': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\JSR_WANG_2019\PHI_1\000_flux\{}'.format(mechname),
   'JSR-WANG-PHI05': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\JSR_WANG_2019\PHI_0.5\000_flux\{}'.format(mechname),
   # 'JSR-WANG-PHI18': r'D:\PhD\OPENSMOKE\CPD\OXIDATION\JSR_WANG_2019\PHI_1.8\000_flux\{}'.format(mechname),
    }

species_list =['C6H6','C6H5','C5H6','C5H5','C6H5OH','C6H5O','C5H4O','C5H5O','C6H4O2','INDENE','C10H8','BIPHENYL']

simul_flds ={
 #  'JSR-CHAI-PHI1': r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\JSR_Chai_1998\PHI_1.02\000_flux\{}'.format(mechname),
 #   'JSR-CHAI-PHI019': r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\JSR_Chai_1998\PHI_0.19\000_flux\{}'.format(mechname),
 #   'JSR-RISTORI-1': r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\JSR_Ristori\000_flux_1\{}'.format(mechname),
 #  'JSR-RISTORI-0.5':r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\JSR_Ristori\000_flux_05\{}'.format(mechname),
 #    'JSR-RISTORI-1.5': r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\JSR_Ristori\000_flux_15\{}'.format(mechname),
 #  'JSR-RISTORI-0.3':r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\JSR_Ristori\000_flux_03\{}'.format(mechname),
 #  'JSR-MARCHAL-1.5':r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\JSR_Machal\PHI_1.5\000_flux\{}'.format(mechname),
 #  'JSR-MARCHAL-0.5':r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\JSR_Machal\PHI_0.5\000_flux\{}'.format(mechname),
 #  'ST-BURCAT-1':r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\ST_Burcat\A_PHI_1_2.44atm_1450_FLUX\{}'.format(mechname),
 #  'ST-BURCAT-2':r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\ST_Burcat\B_PHI_2_2.44atm_1450_FLUX\{}'.format(mechname),
 #  'ST-BURCAT-0.5':r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\ST_Burcat\C_PHI_0.5_2.44atm_1450_FLUX\{}'.format(mechname),
 # 'ST-DACOSTA-1':r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\ST_DaCosta\C_PHI_1_8.5atm_FLUX\{}'.format(mechname),
#  'PFR-BREZ-0.76':r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\PFR_Brezisky_1988\PHI_0.76\{}'.format(mechname),
 #  'PFR-BREZ-1.0':r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\PFR_Brezisky_1988\PHI_1.0\{}'.format(mechname),
 #  'PFR-BREZ-1.36':r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\PFR_Brezisky_1988\PHI_1.36\{}'.format(mechname),
 #'JSR-DACOSTA-1.9':r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\JSR_DaCosta\PHI_1.9_FLUX\{}'.format(mechname),
 #'JSR-DACOSTA-3.6':r'D:\PhD\OPENSMOKE\BENZENE\OXIDATION\JSR_DaCosta\PHI_3.6_FLUX\{}'.format(mechname),
 #'LFS-DAVIS-0.7':r'D:\PhD\OPENSMOKE\BENZENE\FLAMES\aaa_LFS\sensitivities_fiamme\Davis\FLUX0\{}'.format(mechname),
 #'LFS-DAVIS-1.0':r'D:\PhD\OPENSMOKE\BENZENE\FLAMES\aaa_LFS\sensitivities_fiamme\Davis\FLUX3\{}'.format(mechname),
 #'LFS-DAVIS-1.2':r'D:\PhD\OPENSMOKE\BENZENE\FLAMES\aaa_LFS\sensitivities_fiamme\Davis\FLUX5\{}'.format(mechname),
 # 'LFS-JI-0.7':r'D:\PhD\OPENSMOKE\BENZENE\FLAMES\aaa_LFS\sensitivities_fiamme\Ji\FLUX0\{}'.format(mechname),
 #'LFS-JI-1.1':r'D:\PhD\OPENSMOKE\BENZENE\FLAMES\aaa_LFS\sensitivities_fiamme\Ji\FLUX4\{}'.format(mechname),
 'LFS-JHN-1.1':r'D:\PhD\OPENSMOKE\BENZENE\FLAMES\aaa_LFS\sensitivities_fiamme\Johnston\Output_450_3_3_sens_2202\{}'.format(mechname),
 #'LFS-JI-1.5':r'D:\PhD\OPENSMOKE\BENZENE\FLAMES\aaa_LFS\sensitivities_fiamme\Ji\FLUX8\{}'.format(mechname),
 #'FS-TREG-1.9':r'D:\PhD\OPENSMOKE\BENZENE\FLAMES\Tregrossi1999_077\{}'.format(mechname),
 #'FS-VANDOREEN-2':r'D:\PhD\OPENSMOKE\BENZENE\FLAMES\Vandoreen_2009_Phi2\{}'.format(mechname),
   }
species_list =['C6H6','C6H5','C5H6','C5H5','C6H5OH','C6H5O']

simul_flds ={
 #  'PFR-BREZ-PYR':r'D:\PhD\OPENSMOKE\PHENOL\DB_Phenol\PFR_Brezinsky\Pyrolysis\{}'.format(mechname),
 #'PFR-MANION-HY':r'D:\PhD\OPENSMOKE\PHENOL\DB_Phenol\PFR_Manion\000_flux\{}'.format(mechname),
  #'ST-HORN-20':r'D:\PhD\OPENSMOKE\PHENOL\DB_Phenol\ST_Horn\Simulations_20ppm\{}'.format(mechname),
  # 'PFR-BREZ-PHI1.7':r'D:\PhD\OPENSMOKE\PHENOL\DB_Phenol\PFR_Brezinsky\Ox_1.73\{}'.format(mechname),
 # 'PFR-BREZ-PHI1':r'D:\PhD\OPENSMOKE\PHENOL\DB_Phenol\PFR_Brezinsky\Ox_1.03\{}'.format(mechname),
 #  'PFR-BREZ-PHI0.6':r'D:\PhD\OPENSMOKE\PHENOL\DB_Phenol\PFR_Brezinsky\Ox_0.64\{}'.format(mechname),
 #   'RCM-BUTTG-1-10bar':r'D:\PhD\OPENSMOKE\PHENOL\DB_Phenol\RCM_Aachen\10 bar\Phi_1_flux\{}'.format(mechname),
# 'RCM-BUTTG-05-10bar':r'D:\PhD\OPENSMOKE\PHENOL\DB_Phenol\RCM_Aachen\10 bar\Phi_0.5_flux\{}'.format(mechname),
  #  'PFR-THOMAS-PYR':r'D:\PhD\OPENSMOKE\CATECHOL\PFR_THOMAS\PYR_FLUX\{}'.format(mechname),
   # 'PFR-THOMAS-PHI1':r'D:\PhD\OPENSMOKE\CATECHOL\PFR_THOMAS\PHI_0.92_FLUX\{}'.format(mechname),
  # 'JSR-NOWAK-PYR':r'D:\PhD\OPENSMOKE\ANISOLE\JSR_Nowakoska\pyr_flux\{}'.format(mechname),
 #   'JSR-NOWAK-PYR1000':r'D:\PhD\OPENSMOKE\ANISOLE\JSR_Nowakoska\pyr_flux_1000\{}'.format(mechname),
 #   'JSR-NOWAK-PYR1100':r'D:\PhD\OPENSMOKE\ANISOLE\JSR_Nowakoska\pyr_flux_1100\{}'.format(mechname),
 #  'JSR-NOWAK-OX':r'D:\PhD\OPENSMOKE\ANISOLE\JSR_Nowakoska\ox_flux\{}'.format(mechname),
 # 'JSR-WAGNON-PHI1b':r'D:\PhD\OPENSMOKE\ANISOLE\JSR_Wagnon\flux_phi_1\{}'.format(mechname),
#  'JSR-WAGNON-PHI2':r'D:\PhD\OPENSMOKE\ANISOLE\JSR_Wagnon\flux_phi_2\{}'.format(mechname),
 # 'JSR-WAGNON-PHI0.5':r'D:\PhD\OPENSMOKE\ANISOLE\JSR_Wagnon\flux_phi_0.5\{}'.format(mechname),
 #'ST-AACHEN-PHI1':r'D:\PhD\OPENSMOKE\ANISOLE\ST_Aachen\10 bar flux\{}'.format(mechname),
#'RCM-AACHEN-PHI1':r'D:\PhD\OPENSMOKE\ANISOLE\RCM_Aachen\10 bar flux\{}'.format(mechname),
 #'LFS-PHI1':r'D:\PhD\OPENSMOKE\ANISOLE\LFS_FLUX_phi1\{}'.format(mechname),
# 'GUAIAC-PHI1':r'D:\PhD\OPENSMOKE\GUAIACOL\JSR_NOWAKOWSKA\OXFLUX\{}'.format(mechname),
'GUAIAC-PYR':r'D:\PhD\OPENSMOKE\GUAIACOL\JSR_NOWAKOWSKA\PYRFLUX\{}'.format(mechname),
    }

#PYR
species_list =['C7H8','C6H5C2H3','BIPHENYL','FLUORENE','C14H10']
species_list =['C6H5OH','C6H5O','C5H6','C5H5','C6H6','C6H5','C10H8','INDENE']
species_list =['C6H5OH','C6H5O','C5H6','C5H5']
species_list =['C6H5OH','C6H5O','C5H6','C5H5','C6H6','C10H8']
#OX
species_list =['C6H5OH','CATECHOL','C6H5O','OC6H4OH','C6H4OH','C6H5']
species_list =['C6H5OH','C6H5O','C6H4OH','OC6H4OH','CATECHOL','C6H4O2','C5H6','C5H5','C6H6','C6H5']
species_list =['C6H5OCH3','C6H5O','C6H5OH','OC6H4CH3','CRESOL','C6H6','C6H5CHO','C6H5C2H3','C7H8','C5H6','C5H5CH3','C10H8','BZFUR']
species_list =['C6H5OCH3','C6H5O','C6H5OH','OC6H4CH3','CRESOL','C6H6','C6H5','C6H5CHO']
species_list =['GUAIACOL','OC6H4OH','SALICALD','CRESOL','OC6H4CH3','C6H5OH','C6H5O']
filter_dcts = [{},{},{}] # filter according to selected criteria in name
threshs = [1e-2,1e-2,1e-2]


sortlists = [['classgroup'],['speciestype'],['bimoltype']] # classgroup, speciestype, subclass, bimoltype (R+R, RSR+RSR, M+M, ETC) # sum if both apply and sort by this criteria


###################################################################################################
##################################################################################################
# filter_dct = {'classgroup': ['GROWTH', 'ENLARGE']} # filter according to selected criteria in name

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

    rxns_sorted.sum_fwbw()
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

    del rxns_sorted, ropa, consumption, production, tot_rop, tot_rop_df, sortdf