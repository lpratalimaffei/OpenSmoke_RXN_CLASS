'''
MODULE: KineticMechanism
@Authors:
    Alberto Cuoci [1], Timoteo Dinelli [1]
    [1]: CRECK Modeling Lab, Department of Chemistry, Materials, and Chemical Engineering, Politecnico di Milano
@Contacts:
    alberto.cuoci@polimi.it
    timoteo.dinelli@polimi.it
@Additional notes:
    -This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
     Please report any bug to: alberto.cuoci@polimi.it
    -This is a modified class frome the original one KineticMechanism.py inside PyTools4OpenSMOKE 
     modified by Timoteo Dinelli to handle gas-phase only to post-process data in order to perform 
     Rate Of Production Analysis
'''

import xml.etree.ElementTree as ET
import numpy as np
from scipy import sparse
import sys

class KineticMechanism:

    '''
    Description of the KineticMechanism class
    '''
    def __init__(self, xml_file_name):
    
        tree = ET.parse(xml_file_name)
        root = tree.getroot()
    
        # List of elements
        elements = root.find('NamesOfElements')
        elements = (elements.text).split()
        NumberOfElements = len(elements) # ne

        # List of species
        species = root.find('NamesOfSpecies')
        species = (species.text).split()
        NumberOfSpecies = len(species) # ns

        # Elemental composition
        atomic = root.find('AtomicComposition')
        atomic = (atomic.text).split()
        atomic = np.reshape(atomic, (NumberOfSpecies, NumberOfElements))
        atomic = np.float32(atomic)

        # Indices of relevant elements
        iC = elements.index('C')
        iH = elements.index('H')
        iO = elements.index('O')
        iN = elements.index('N')
        
        # Elements molecular weights
        mwe = []
        for elem in elements:
            if (elem == 'C'): mwe.append(12.010999679565430)
            if (elem == 'H'): mwe.append(1.008000016212463)
            if (elem == 'O'): mwe.append(15.998999595642090)
            if (elem == 'N'): mwe.append(14.0069999694824)
            if (elem == 'HE'): mwe.append(4.002999782562256)
            if (elem == 'AR'): mwe.append(39.948001861572270)
        mwe = np.array(mwe)
        
        # Species molecular weights
        mws = atomic.dot(mwe)
        
        kinetics = root.find('Kinetics')
        
        NumberOfReactions = int( (kinetics.find('NumberOfReactions')).text ) # nr

        reaction_classes = kinetics.find('ReactionClasses')
        self.rxnclass = dict.fromkeys(np.arange(1,NumberOfReactions+1))
        self.rxnsubclass = dict.fromkeys(np.arange(1,NumberOfReactions+1))
        
        if (reaction_classes != None):
            # classes: {mainclass: {subclass: [indices]}}
            classes = {}
            for child in reaction_classes:
                if (child.tag == 'MainClass'):
                    classname = child.attrib['name']
                    if classname not in list(classes.keys()):
                        # initialize
                        classes[classname] = {}
                    for subclass in child:
                        #if (subclass.tag == 'SubClass'): # should only contain subclass type
                        subclassname = subclass.attrib['name']
                        if subclassname not in list(classes[classname].keys()):
                            classes[classname][subclassname] = []
                        for subsec in subclass:
                            if (subsec.tag == 'ReactionIndices'):
                                dummy = (subsec.text).split()
                                
                                for i in dummy: # assign rxn indices
                                    i = int(i) +1
                                    classes[classname][subclassname].append(i)
                                    self.rxnclass[i] = classname
                                    self.rxnsubclass[i] = subclassname

            # print(classes)
            self.classes = classes
            
        # Kinetic parameters
        kinetic_parameters = kinetics.find('KineticParameters')
        direct = kinetic_parameters.find('Direct')
        
        lnA = direct.find('lnA')
        lnA = (lnA.text).split()
        lnA = np.float64(lnA[1:])
        A = np.exp(lnA)
        
        Beta = direct.find('Beta')
        Beta = (Beta.text).split()
        Beta = np.float64(Beta[1:])
        
        E_over_R = direct.find('E_over_R')
        E_over_R = (E_over_R.text).split()
        E_over_R = np.float64(E_over_R[1:])

        
        # Assign internal members
        
        self.elements = elements
        self.atomic = atomic
        self.species = species
        self.NumberOfElements = NumberOfElements 
        self.NumberOfSpecies = NumberOfSpecies
        self.NumberOfReactions = NumberOfReactions
            
        self.iC = iC
        self.iH = iH
        self.iO = iO
        self.iN = iN
        
        self.mwe = mwe
        self.mws = mws
        
        self.groups = []

        self.A = A
        self.Beta = Beta
        self.E_over_R = E_over_R
        
        self.reactions = []
    
   
    def Group(self, name):
  
        for i in range(len(self.groups)):
            if (name == self.groups[i]['name']): return self.groups[i]
    
    def ReadKinetics(self, xml_file_name):
    
        tree = ET.parse(xml_file_name)
        root = tree.getroot()

        sub_root = root.find('reaction-names') 
        reaction_list = sub_root.text.split()
        
        self.reaction_lines = reaction_list
        
    def ProcessingReactions(self):
    
        for i in range(self.NumberOfReactions):

            reaction = {    'index': i+1, 'name': self.reaction_lines[i], \
                            'A': self.A[i], 'Beta': self.Beta[i], 'E_over_R': self.E_over_R[i], \
                            'class': self.rxnclass[i+1], 'subclass': self.rxnsubclass[i+1]
                    }
            self.reactions.append(reaction)
        
        return self.reactions
    
    def ReactionsWithSpecies(self, species_name, flag):
        
        if (species_name in self.species):
        
            index = self.species.index(species_name)
            
            if (flag == 'R'):
                return self.species_in_reaction_as_reactant[index]
            elif (flag == 'P'):
                return self.species_in_reaction_as_product[index]
            elif (flag == 'RP'):
                return self.species_in_reaction_as_reactant_or_product[index]
                
        else:
        
            return []
    
    def SubReactionsWithSpecies(self, species_name, flag, indices):
        
        
        
        if (species_name in self.species):
        
            index = self.species.index(species_name)
            
            reactions = []
            
            if (flag == 'R'):
                n = len(self.species_in_reaction_as_reactant[index])
                for i in range(n):
                    if self.species_in_reaction_as_reactant[index][i] in indices:
                        reactions.append(self.species_in_reaction_as_reactant[index][i])
                        
            elif (flag == 'P'):
                n = len(self.species_in_reaction_as_product[index])
                for i in range(n):
                    if self.species_in_reaction_as_product[index][i] in indices:
                        reactions.append(self.species_in_reaction_as_product[index][i])

            elif (flag == 'RP'):
                n = len(self.species_in_reaction_as_reactant_or_product[index])
                for i in range(n):
                    if self.species_in_reaction_as_reactant_or_product[index][i] in indices:
                        reactions.append(self.species_in_reaction_as_reactant_or_product[index][i])
                        
            return (np.unique(reactions)).tolist()
            
        else:
        
            return []
    
    def ReactionsWithMultipleSpecies(self, species_names, flags, logic_operator):
    
        indices = list(range(0,self.nr))
        return self.SubReactionsWithMultipleSpecies(species_names, flags, indices, logic_operator)
     
    def SubReactionsWithMultipleSpecies(self, species_names, flags, indices, logic_operator):
        
        ns = len(species_names)

        if (logic_operator == 'OR'):
            selection = []
            for i in range(0,ns):
                print(species_names[i])
                local_selection = self.SubReactionsWithSpecies(species_names[i], flags[i], indices)
                selection.extend(local_selection)
            return (np.unique(selection)).tolist()
            
        elif (logic_operator == 'AND'):
            selection = self.SubReactionsWithSpecies(species_names[0], flags[0], indices)
            for i in range(1,ns):
                selection = self.SubReactionsWithSpecies(species_names[i], flags[i], selection)
            return (np.unique(selection)).tolist()
            
        return []
       
    def ReactionsWithoutSpecies(self, species_name, flag):
        
        indices = list(range(0,self.nr))
        return self.SubReactionsWithoutSpecies(species_name, flag, indices)
       
    def SubReactionsWithoutSpecies(self, species_name, flag, indices):  
    
        reactions_with = self.SubReactionsWithSpecies(species_name, flag, indices)
        reactions_without = indices
    
        return np.array(list(set(reactions_without) - set(reactions_with)))
    
    def ReactionsWithoutMultipleSpecies(self, species_names, flags, logical_operator):
        
        indices = list(range(0,self.nr))
        return self.SubReactionsWithoutMultipleSpecies(species_names, flags, indices, logical_operator)  
     
    def SubReactionsWithoutMultipleSpecies(self, species_names, flags, indices, logical_operator):  
    
        reactions_with = self.SubReactionsWithMultipleSpecies(species_names, flags, indices, logical_operator)
        reactions_without = indices
    
        return np.array(list(set(reactions_without) - set(reactions_with)))  
       
    def ReactionsWithListsOfSpecies(self, single_species, flag_single_species, list_of_species, flags_list_of_species, logical_operator):
    
        indices = list(range(0,self.nr))
        return self.SubReactionsWithListsOfSpecies(single_species, flag_single_species, list_of_species, flags_list_of_species, indices, logical_operator)
    
    def SubReactionsWithListsOfSpecies(self, single_species, flag_single_species, list_of_species, flags_list_of_species, indices, logical_operator):
    
        sub_indices = self.SubReactionsWithSpecies(single_species, flag_single_species, indices)
        
        reactions = []
        if (logical_operator == 'AND'):
            for i in range(len(list_of_species)):
                reactions.extend( self.SubReactionsWithSpecies(list_of_species[i], flags_list_of_species[i], sub_indices) )
        
        return (np.unique(reactions)).tolist()
           
    def PrintKineticMechanism(self, file_name, species, reactions):  

        nc = len(self.reaction_class_name)
        nr = len(reactions)

        f = open(file_name, "w")
        
        # ---------------------------------------------------------------------------------------------- #   
        # Write Summary
        # ---------------------------------------------------------------------------------------------- #    
        f.write('! CRECK Modeling Lab @ Politecnico di Milano \n')
        f.write('! Please visit our web-site:  http://creckmodeling.chem.polimi.it/ \n')
        f.write('\n')
        f.write('! Total number of reactions: ' + str(nr) + '\n')
    

        # In case reaction classes are available
        if (nc != 0):
            classes = self.OrganizeReactionsInClasses(reactions)
            f.write('!  * Class [' + 'undefined' + ']: ' + str(len(classes[nc])) + '\n')
            for k in range(nc):
                if (len(classes[k]) != 0):
                    f.write('!  * Class [' + self.reaction_class_name[k] + ']: ' + str(len(classes[k])) + '\n')
                
                
        # ---------------------------------------------------------------------------------------------- #   
        # Write ELEMENTS section
        # ---------------------------------------------------------------------------------------------- #    
        self.PrintElementsSectionInCHEMKIN(f)
        
        # ---------------------------------------------------------------------------------------------- #   
        # Write SPECIES section
        # ---------------------------------------------------------------------------------------------- #
        self.PrintSpeciesSectionInCHEMKIN(f, species)

        # ---------------------------------------------------------------------------------------------- #   
        # Write REACTIONS section
        # ---------------------------------------------------------------------------------------------- #
        f.write('REACTIONS' + '\n')
        
        # No reaction classes
        if (nc == 0):
            
            for i in range(nr):
                index = reactions[i]
                f.write(self.reaction_lines[index])
        
        # Reaction classes
        else:
        
            classes = self.OrganizeReactionsInClasses(reactions)

            # Undefined class
            f.write('! Class: [' + 'undefined' + '] (' + str(len(classes[nc])) + ')\n\n')
            for i in range(len(classes[nc])):
                index = classes[nc][i]
                f.write(self.reaction_lines[index])
            f.write('\n')
                
            # Defined class
            for k in range(nc):
                if (len(classes[k]) != 0):
                    f.write('! Class: [' + self.reaction_class_name[k] + '] (' + str(len(classes[k])) + ')\n\n')
                    for i in range(len(classes[k])):
                        index = classes[k][i]
                        f.write(self.reaction_lines[index])
                    f.write('\n\n')
            
        f.write('\n\n')
        f.write('END' + '\n')
        f.close()
   
    def SpeciesIndicesInSelectedReactions(self, reactions):
        
        flags = [0] * self.ns
        
        nr = len(reactions)
        for j in range(nr):
            i = reactions[j]
            indices = sparse.find(self.nur.getrow(i))[1]
            for k in range(len(indices)):
                flags[indices[k]] = 1
            indices = sparse.find(self.nup.getrow(i))[1]
            for k in range(len(indices)):
                flags[indices[k]] = 1
                
        return flags
    
    def SpeciesInSelectedReactions(self, reactions):
        
        flags = self.SpeciesIndicesInSelectedReactions(reactions)
        
        species = []
        for i in range(self.ns):
            if flags[i] == 1: species.append(self.species[i])
            
        return species
 
    def RemoveCrossingValues(self, v1,v2):
        
        crossing = []
        for i in range(len(v1)):
            for j in range(len(v2)):
                if (v1[i] == v2[j]): crossing.append(v2[j])
                
        crossing = (np.unique(crossing)).tolist()
        
        for i in range(len(crossing)):
            v2.remove(crossing[i])
            
        return v2
            
    def PrintKineticMechanismByClasses(self, file_name, by_groups, species, reaction_indices, reaction_groups, reaction_comments):  

        ngroups = len(reaction_indices)
        nc = len(self.reaction_class_name)
        
        nr_tot = 0
        for i in range(ngroups):
            nr_tot = nr_tot + len(reaction_indices[i])
            
        f = open(file_name, "w")
        
        # ---------------------------------------------------------------------------------------------- #   
        # Write Summary
        # ---------------------------------------------------------------------------------------------- #    
        f.write('! CRECK Modeling Lab @ Politecnico di Milano \n')
        f.write('! Please visit our web-site:  http://creckmodeling.chem.polimi.it/ \n')
        f.write('\n')
        f.write('! Total number of reactions: ' + str(nr_tot) + '\n')
    

        # In case reaction groups are available
        f.write('! --------------------------------------------------------------------------- !\n')
        f.write('! Kinetic mechanism structure by Groups\n')
        f.write('! --------------------------------------------------------------------------- !\n')
        for j in range(ngroups):
            nr = len(reaction_indices[j])
            f.write('!  * Group [' + reaction_groups[j] + ']: ' + str(nr) + '\n')
            if (nc != 0):
                classes = self.OrganizeReactionsInClasses(reaction_indices[j])          
                if (len(classes[nc]) != 0):
                    f.write('!    - Class [' + 'undefined' + ']: ' + str(len(classes[nc])) + '\n')
                for k in range(nc):
                    if (len(classes[k]) != 0):
                        f.write('!    - Class [' + self.reaction_class_name[k] + ']: ' + str(len(classes[k])) + '\n')
        
        
        if (nc != 0):

            f.write('\n')
            f.write('! --------------------------------------------------------------------------- !\n')
            f.write('! Kinetic mechanism structure by Classes\n')
            f.write('! --------------------------------------------------------------------------- !\n')
        
            groups = self.OrganizeReactionsInGroups(reaction_indices)
            
            # Number of reactions per group
            nrc = [0]*(nc+1)
            for k in range(nc+1):
                for j in range(len(groups)):
                    nrc[k] = nrc[k] + len(groups[j][k])        
            
            f.write('!  * Class [' + 'undefined' + ']: ' + str(nrc[nc]) + '\n')
            for j in range(len(groups)):
                if (len(groups[j][nc])!=0): f.write('!    - Group [' + reaction_groups[j] + ']: ' + str(len(groups[j][nc])) + '\n')
            
            for k in range(nc):
                f.write('!  * Class [' + self.reaction_class_name[k] + ']: ' + str(nrc[k]) + '\n')
                for j in range(len(groups)):
                    if (len(groups[j][k])!=0): f.write('!    - Group [' + reaction_groups[j] + ']: ' + str(len(groups[j][k])) + '\n')
                
                
        # ---------------------------------------------------------------------------------------------- #   
        # Write ELEMENTS section
        # ---------------------------------------------------------------------------------------------- #
        self.PrintElementsSectionInCHEMKIN(f)
        
        # ---------------------------------------------------------------------------------------------- #   
        # Write SPECIES section
        # ---------------------------------------------------------------------------------------------- #
        self.PrintSpeciesSectionInCHEMKIN(f, species)

        # ---------------------------------------------------------------------------------------------- #   
        # Write REACTIONS section
        # ---------------------------------------------------------------------------------------------- #
        f.write('REACTIONS' + '\n')
        f.write('! Total number of reactions: ' + str(nr_tot) + '\n\n')
        
        # No reaction classes
        if (nc == 0):
        
            for j in range(ngroups):
            
                nr = len(reaction_indices[j])
                
                f.write('\n')
                f.write('! [' + reaction_groups[j] + ']' + '\n')
                f.write('! ' + reaction_comments[j] + ' (' + str(nr) + ')\n')
                for i in range(nr):
                    index = reaction_indices[j][i]
                    f.write(self.reaction_lines[index])
                f.write('\n')
                f.write('\n')
        
        # Reactions by groups
        elif (by_groups == True and nc != 0):
        
            for j in range(ngroups):
            
                nr = len(reaction_indices[j])
                classes = self.OrganizeReactionsInClasses(reaction_indices[j])
        
                f.write('\n')
                f.write('! [' + reaction_groups[j] + ']' + '\n')
                f.write('! ' + reaction_comments[j] + ' (' + str(nr) + ')\n')
                
                # Undefined class
                if (len(classes[nc]) != 0):
                    f.write('! Class: [' + reaction_groups[j] + '][' + 'undefined' + '] (' + str(len(classes[nc])) + ')\n\n')
                    for i in range(len(classes[nc])):
                        index = classes[nc][i]
                        f.write(self.reaction_lines[index])
                    f.write('\n')
                    
                # Defined class
                for k in range(nc):
                    if (len(classes[k]) != 0):
                        f.write('! Class: [' + reaction_groups[j] + '][' + self.reaction_class_name[k] + '] (' + str(len(classes[k])) + ')\n\n')
                        for i in range(len(classes[k])):
                            index = classes[k][i]
                            f.write(self.reaction_lines[index])
                        f.write('\n\n')
        
                f.write('\n')
                f.write('\n')
                
                
        # Reaction by classes       
        elif (by_groups == False and nc != 0):
        
            groups = self.OrganizeReactionsInGroups(reaction_indices)
            
            # Number of reactions per group
            nrc = [0]*(nc+1)
            for k in range(nc+1):
                for j in range(len(groups)):
                    nrc[k] = nrc[k] + len(groups[j][k])
                    
            if (nrc[nc] != 0):
                f.write('\n')
                f.write('! [CLASS] [' + 'undefined' + '] (' + str(nrc[nc]) + ')\n\n')
                
                for j in range(len(groups)):
                    if (len(groups[j][nc])!=0): 
                        f.write('! Group [' + reaction_groups[j] + ']: ' + str(len(groups[j][nc])) + '\n\n')
                        for i in range(len(groups[j][nc])):
                            index = groups[j][nc][i]
                            f.write(self.reaction_lines[index])
                        f.write('\n\n')        
                    
        
            for k in range(nc):
            
                if (nrc[k] != 0):
                
                    f.write('\n')
                    f.write('! [SOOTCLASS] [' + self.reaction_class_name[k] + '] (' + str(nrc[k]) + ')\n\n')
                    
                    for j in range(len(groups)):
                        if (len(groups[j][k])!=0): 
                            f.write('! Group [' + reaction_groups[j] + ']: ' + str(len(groups[j][k])) + '\n\n')
                            for i in range(len(groups[j][k])):
                                index = groups[j][k][i]
                                f.write(self.reaction_lines[index])
                            f.write('\n\n')
            
        f.write('\n\n')
        f.write('END' + '\n')
        
        f.close()
    
    def OrganizeReactionsInGroups(self, reaction_indices):
    
        ngroups = len(reaction_indices)
        groups = [[] for x in range(ngroups)]
        
        nc = len(self.reaction_class_name)
        for k in range(ngroups):
            classes = [[] for x in range(nc+1)]
            
            nr = len(reaction_indices[k])
            for i in range(nr):
                index = self.reaction_class_belonging[reaction_indices[k][i]]
                if (index == -1): classes[nc].append(reaction_indices[k][i])
                else: classes[ index ].append(reaction_indices[k][i])
                
            groups[k] = classes
            
        return groups
            
    def OrganizeReactionsInClasses(self, reactions):
    
        nc = len(self.reaction_class_name)
        nr = len(reactions)
        classes = [[] for x in range(nc+1)]
        
        for i in range(nr):
            index = self.reaction_class_belonging[reactions[i]]
            if (index == -1): classes[nc].append(reactions[i])
            else: classes[ index ].append(reactions[i])
            
        return classes
    
    def SortAccordingToMolecularWeight(self, species):

        mws = []
        for i in range(len(species)):
            mws.append(self.mws[self.species.index(species[i])])
        return [element for _, element in sorted( zip(mws, species) )]
       
    def SortAndSplitSpecies(self, species, criterion):
    
        species_gas = species
        
        # Inerts
        inerts = []
        if ('HE' in species_gas): inerts.append('HE')
        if ('AR' in species_gas): inerts.append('AR')
        if ('N2' in species_gas): inerts.append('N2')
        species_gas = list(set(species_gas) - set(inerts))
        
        # Carbon
        carbon = []
        if ('CSOLID' in species_gas): carbon.append('CSOLID')
        species_gas = list(set(species_gas) - set(carbon))
        
        # PAHs12
        pahs12 = []
        for i in range(len(self.Group('PAH12')['list'])):
            if (self.Group('PAH12')['list'][i] in species_gas): 
                pahs12.append(self.Group('PAH12')['list'][i])       
        species_gas = list(set(species_gas) - set(pahs12))
            
        # PAHs34
        pahs34 = []
        for i in range(len(self.Group('PAH34')['list'])):
            if (self.Group('PAH34')['list'][i] in species_gas): 
                pahs34.append(self.Group('PAH34')['list'][i])
        species_gas = list(set(species_gas) - set(pahs34))
        
        # PAHsLP
        pahslp = []
        for i in range(len(self.Group('PAHLP')['list'])):
            if (self.Group('PAHLP')['list'][i] in species_gas): 
                pahslp.append(self.Group('PAHLP')['list'][i])
        species_gas = list(set(species_gas) - set(pahslp))   

        # Spherical particles
        sp = []
        for i in range(len(self.Group('SP')['list'])):
            if (self.Group('SP')['list'][i] in species_gas): 
                sp.append(self.Group('SP')['list'][i])
        species_gas = list(set(species_gas) - set(sp))

        # Aggregates
        agg = []
        for i in range(len(self.Group('AGG')['list'])):
            if (self.Group('AGG')['list'][i] in species_gas): 
                agg.append(self.Group('AGG')['list'][i])
        species_gas = list(set(species_gas) - set(agg))
        
        # Fixed gas species
        fixed_species = []
        if ('C2H4' in species_gas): fixed_species.append('C2H4')
        if ('H2' in species_gas):   fixed_species.append('H2')
        if ('C2H2' in species_gas): fixed_species.append('C2H2')
        if ('O2' in species_gas):   fixed_species.append('O2')
        if ('CO2' in species_gas):  fixed_species.append('CO2')
        if ('H2O' in species_gas):  fixed_species.append('H2O')
        if ('CO' in species_gas):   fixed_species.append('CO')
        if ('OH' in species_gas):   fixed_species.append('OH')
        species_gas = list(set(species_gas) - set(fixed_species))
        species_gas = fixed_species + self.SortAccordingToMolecularWeight(species_gas)
                
        return  species_gas, \
                self.SortAccordingToMolecularWeight(pahs12), \
                self.SortAccordingToMolecularWeight(pahs34), \
                sorted(pahslp), \
                sorted(sp), \
                sorted(agg), \
                carbon, inerts