'''
SCRIPT: rxnclasses
@Authors:
    Luna Pratali Maffei [1]
    [1]: CRECK Modeling Lab, Department of Chemistry, Materials, and Chemical Engineering, Politecnico di Milano
@Contacts:
    luna.pratali@polimi.it
@Additional notes:
    This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
    Please report any bug to: alberto.cuoci@polimi.it
'''

# Import main libraries
from msilib.schema import Class
from operator import index
from pylab import meshgrid, cm, imshow, contour, clabel, colorbar, axis, title, show
import pandas as pd
import sys
import os
import numpy as np
import copy
import matplotlib.pyplot as plt
plt.style.use('default')

R = ['R', 'H', 'OH', 'O2', 'O', 'CH3', 'HO2', 'HCO', 'C2H3']
RSR = ['C5H5', 'C7H7', 'C6H5O']


def check_bimol_type(speciestype, subclass):
    """
    check if rxn is M+M, RSR+M, R+M, R+R, RSR+RSR
    """
    try:
        species_type_1 = speciestype.split('-')[-1]
    except AttributeError:
        return 'UNSORTED'
    if subclass == 'ENLARGE':
        return 'UNIMOL'
    # CHECK SPECIES TYPE 2
    if subclass.split('_')[0] in ['HABS', 'ADD', 'REC', 'ENLARGE'] and '-' not in subclass:
        if len(subclass.split('_')) == 1:
            return 'UNSORTED'
        if subclass.split('_')[1] in R:
            species_type_2 = 'R'
        elif subclass.split('_')[1] in RSR:
            species_type_2 = 'RSR'
        else:
            species_type_2 = 'NA'
    elif subclass.split('_')[0] in ['ADD', 'REC', 'ENLARGE'] and '-' in subclass:
        species_type_2 = subclass.split('-')[-1]
    else:
        return 'UNIMOL'

    type_list = [species_type_1, species_type_2]
    type_list.sort()

    return '+'.join(type_list)


class rxnclass:

    def __init__(self, reactions):
        """ allocate dataframe 
        """

        # turn reactions into a dataframe
        # index = rxn index, columns = [classgroup, class, subclass, flux]
        self.rxn_class_df = pd.DataFrame(index=np.arange(1, len(
            reactions)+1), columns=['name', 'classgroup', 'speciestype', 'subclass', 'bimoltype'], dtype=object)
        self.reactions = reactions

        """ assign class and sublcass
        """
        for rxn in self.reactions:
            idx = rxn['index']
            self.rxn_class_df['name'][idx] = rxn['name']
            self.rxn_class_df['speciestype'][idx] = rxn['class']
            self.rxn_class_df['subclass'][idx] = rxn['subclass']
            self.rxn_class_df['bimoltype'][idx] = check_bimol_type(
                rxn['class'], rxn['subclass'])

        self.rxn_class_df = self.rxn_class_df.replace('None', np.nan)
        self.rxn_class_df['speciestype'] = self.rxn_class_df['speciestype'].replace(
            np.nan, 'UNSORTED')
        self.rxn_class_df['subclass'] = self.rxn_class_df['subclass'].replace(
            np.nan, 'UNSORTED')

    def assign_class_grp(self, subcl_grp_dct):
        """ assign class group if available
        """
        for subcl, subset in self.rxn_class_df.groupby('subclass'):
            # None values are automatically discarded
            rxns = subset.index
            try:
                self.rxn_class_df['classgroup'][rxns] = subcl_grp_dct[subcl]
            except KeyError:
                self.rxn_class_df['classgroup'][rxns] = 'UNSORTED'
                print(
                    '*Warning: subclass {} not found in class groups'.format(subcl))
                continue

    def assign_flux(self, tot_rop_df):
        """ assign flux from the flux analysis
        """

        flux_sp_name = tot_rop_df.columns[0]
        self.rxn_class_df = pd.concat([self.rxn_class_df, tot_rop_df], axis=1)
        # replace missing flux values with 0
        self.rxn_class_df[flux_sp_name] = self.rxn_class_df[flux_sp_name].replace(
            np.nan, 0.0)

        # renormalize
        #renorm_factor = sum(abs(self.rxn_class_df[flux_sp_name]))
        
        # renormalize
        # rxn_class_df = rxn_class_df.sort_values(by='absflux', ascending = False)

    def sum_fwbw(self):
        """
        sum flux for fw and bw rxns
        if '2INDENYL=>C18H14' in self.rxn_class_df['name'].values and 'C18H14=>2INDENYL' in self.rxn_class_df['name'].values:
            idx0 = self.rxn_class_df.index[self.rxn_class_df['name'] == '2INDENYL=>C18H14'][0]
            idx1 = self.rxn_class_df.index[self.rxn_class_df['name'] == 'C18H14=>2INDENYL'][0]
            self.rxn_class_df.loc[idx0, flux_sp_name] += self.rxn_class_df.loc[idx1, flux_sp_name]
            self.rxn_class_df.loc[idx1, flux_sp_name] *= 0
            # self.rxn_class_df = self.rxn_class_df.drop(idx1, axis = 0)
        """     
        self.flux_cols = [
            col for col in self.rxn_class_df.columns if 'flux' in col]
        #list of irrev rxns
        rxns_irrev = [rxn for rxn in self.rxn_class_df['name'] if '=>' in rxn]
        idxs_rxns_irrev = list(set([self.rxn_class_df.index[self.rxn_class_df['name'] == rxn][0] for rxn in rxns_irrev]))
        rxns_irrev_series = pd.DataFrame(index=idxs_rxns_irrev, columns=['name', 'rcts', 'prds'])
        for idx in idxs_rxns_irrev:
            rxn = self.rxn_class_df['name'][idx]
            rxns_irrev_series['name'][idx] = rxn
            rxns_irrev_series['rcts'][idx]  = '+'.join(sorted(rxn.split('=>')[0].split('+')))
            rxns_irrev_series['prds'][idx]  = '+'.join(sorted(rxn.split('=>')[1].split('+')))

        for idx0 in idxs_rxns_irrev:
            # continue if idx was removed
            if idx0 not in rxns_irrev_series.index:
                continue
            rcts, prds = rxns_irrev_series[['rcts', 'prds']].loc[idx0]
            idx1s = list(set(rxns_irrev_series.index[rxns_irrev_series['rcts'] == prds]) & set(rxns_irrev_series.index[rxns_irrev_series['prds'] == rcts]))
            for idx1 in idx1s:
                if self.rxn_class_df['speciestype'][idx1] == 'UNSORTED':
                    idxkeep = copy.deepcopy(idx0)
                    idxrem = copy.deepcopy(idx1)
                elif self.rxn_class_df['speciestype'][idx0] == 'UNSORTED':
                    idxkeep = copy.deepcopy(idx1)
                    idxrem = copy.deepcopy(idx0)
                else:
                    # pick the highest flux
                    fl0 = max(abs(self.rxn_class_df.loc[idx0, self.flux_cols]))
                    fl1 = max(abs(self.rxn_class_df.loc[idx1, self.flux_cols]))
                    maxfl = max([fl0, fl1])
                    idxkeep = [idx0, idx1][[fl0, fl1].index(maxfl)]
                    idxrem = [idx0, idx1][1-[idx0, idx1].index(idxkeep)]
                print('merging flux {} and removing {}'.format(self.rxn_class_df['name'][idxkeep], self.rxn_class_df['name'][idxrem]))
                self.rxn_class_df.loc[idxkeep, self.flux_cols] += self.rxn_class_df.loc[idxrem, self.flux_cols]
                self.rxn_class_df = self.rxn_class_df.drop(idxrem, axis = 0)
                rxns_irrev_series = rxns_irrev_series.drop(idxrem, axis = 0)
                      
        self.rxn_class_df_all = copy.deepcopy(self.rxn_class_df)

    def filter_class(self, filter_dct):
        """ only keep classes according to criteria listed in filter_dct
        """
        dict_indexes = dict.fromkeys(filter_dct.keys())
        for criterion, values in filter_dct.items():
            dict_indexes[criterion] = list(self.rxn_class_df[[any(v in val for v in values) for val in self.rxn_class_df[criterion]]].index)

        indexes_filter = list(set.intersection(*[set(val) for val in dict_indexes.values()]))
        self.rxn_class_df = self.rxn_class_df.loc[list(set(indexes_filter))]

    def filter_flux(self, threshold=1e-3):
        """ delete all reactions with contributions below a threshold 
        """

        indexes_filter = []

        for flux_sp_name in self.flux_cols:
            indexes_filter.extend(
                list(self.rxn_class_df[abs(self.rxn_class_df[flux_sp_name])/max(abs(self.rxn_class_df[flux_sp_name])) > threshold].index))

        self.rxn_class_df = self.rxn_class_df.loc[list(set(indexes_filter))]

    def sortby(self, sortlist):
        """ sum fluxes by criteria in sortlist
        """
        # check that all criteria are columns
        if not all(criterion in self.rxn_class_df.columns for criterion in sortlist):
            print('*Error: criteria not all present in dataframe columns - exiting')
            sys.exit()

        # group and sum
        new_sort_df = pd.DataFrame(index=self.flux_cols)
        for grp_idx, grp_df in self.rxn_class_df.groupby(sortlist):
            if isinstance(grp_idx, str):
                name = grp_idx
            else:
                name = '['+']['.join(grp_idx)+']'
            print(grp_idx, '\n', grp_df, '\n')
            new_sort_df[name] = grp_df[self.flux_cols].sum()

        # renormalize by species
        
        for flux_sp_name in new_sort_df.index:
            renorm_factor = max(abs(new_sort_df.loc[flux_sp_name]))
            weight_factor = sum(abs(self.rxn_class_df[flux_sp_name]))/sum(abs(self.rxn_class_df_all[flux_sp_name]))
            #new_sort_df.loc[flux_sp_name] /= renorm_factor
            new_sort_df.loc[flux_sp_name] *= (weight_factor/renorm_factor)
        
        print(new_sort_df)

        return new_sort_df


def plot_heatmap(sort_df, savepath):
    """ heat maps of x: df.columns, y: df.index """
    # generate the figure
    fig, axes = plt.subplots(
        figsize=[len(sort_df.columns), len(sort_df.index)])
    # valmax = max([np.min(sort_df.values), np.max(sort_df.values)])
    image = axes.imshow(sort_df.values, aspect='auto', cmap='RdBu_r',
    #)
                        vmin=-1, vmax=1)  # coolwarm, RdBu, seismic, bwr
    axes.set_yticks(np.arange(0, len(sort_df.index)))    
    axes.set_xticks(np.arange(0, len(sort_df.columns)))
    size_for_lbl = 'xx-large'*(len(sort_df.index) > 3) + 'medium'*(len(sort_df.index) <= 3)

    axes.set_yticklabels([idx.split('flux_')[1] for idx in sort_df.index], fontsize=size_for_lbl)
    axes.set_xticklabels(sort_df.columns, rotation=90, fontsize=size_for_lbl)
    #image.figure.axes[1].tick_params(axis="y", labelsize=size_for_lbl) # colorbar
    plt.colorbar(image)
    
    

    fig.tight_layout()
    # save the figure

    if os.path.isfile(savepath):
        os.remove(savepath)
    fig.savefig(savepath, dpi=200)

    plt.close()
    """
    colorsrgb = list(matplotlib.cm.gist_rainbow(np.linspace(0, 1, 8)))
    colorsrgb += list(matplotlib.cm.tab10(np.arange(0, 10)))
    colorsrgb += list(matplotlib.cm.Set3(np.arange(0, 12)))
    self.palette_series = pd.Series(
        colorsrgb[:len(SPECIES)], index=SPECIES)
    self.palette_series[REAC] = 'black'
    

    axx.plot(
        t, tW_DF[Si], color=self.palette_series[Si], ls=self.lstyles[kk])

    # labels and formatting
    axes[self.n_rows-1,
            col_plot].set_xlabel(r'$t$ [s]', fontsize=self.fsize_vect[1])
    axes[row_plot, 0].set_ylabel(
        r'$X_{i}$ [-]', fontsize=self.fsize_vect[1])

    axx.set_title((str(T) + ' K'),
                    fontsize=self.fsize_vect[2], fontweight='bold')
    axx.set_xlim([0, min(maxt)])
    axx.set_ylim([0, N_INIT])
    axx.axes.ticklabel_format(axis='both', style='sci', scilimits=(
        0, 0), useLocale=True, useMathText=True)
    axx.yaxis.offsetText.set_fontsize(self.fsize_vect[0])
    axx.xaxis.set_tick_params(labelsize=self.fsize_vect[0])
    axx.xaxis.offsetText.set_fontsize(self.fsize_vect[0])
    axx.yaxis.set_tick_params(labelsize=self.fsize_vect[0])

    axes[0, 0].legend(list(self.SPECIES), loc=4, fontsize=self.legsize)
    # set tight layout
    """
