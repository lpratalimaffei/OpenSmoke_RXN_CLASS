'''
MODULE: ROPA.py
@Author:
    Timoteo Dinelli [1]
    [1]: CRECK Modeling Lab, Department of Chemistry, Materials, and Chemical Engineering, Politecnico di Milano
@Contacts:
    timoteo.dinelli@polimi.it
@Additional notes:
    -This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
     Please report any bug to: timoteo.dinelli@polimi.it
    -This is a class to perform Rate Of Production Analysis based on the original OpenSMOKE postprocessor developed by 
     Alberto Cuoci (https://github.com/acuoci/OpenSMOKEppPostProcessor)
'''
# Import External ibraries
import xml.etree.ElementTree as ET
import numpy as np
import os

# Import internal modules
from KineticMechanism import KineticMechanism
from ReadStoichiometry import ReadStoichiometry

class RateOfProductionAnalysis:

    def __init__(self, 
                kineticMechanismFolder: str,
                outputFolder: str,
                Species: str):
        """
        :param kineticsMechanismFolder: path to the folder containing the kinetics.xml and reaction_names.xml files
        :param outputFolder: path to the folder containing output.xml file
        :param Species: species on which it is desired to perform the Rate Of Production Analysis
        """

        self.numberOfSpecies, self.numberOfReactions, self.nur,\
        self.nup, self.reactionsNames, self.speciesNames = self.ProcessKineticMechanism(kineticMechanismFolder = kineticMechanismFolder)
        
        self.indexSpecies = self.GetIndexSpecies(speciesNames = Species)

        self.matrixReactionRates, self.time = self.ReadOutputXML(outputFolder = outputFolder)
        
    def ProcessKineticMechanism(self, kineticMechanismFolder):

        pathkin = os.path.join(kineticMechanismFolder, 'kinetics.xml')
        pathname = os.path.join(kineticMechanismFolder, 'reaction_names.xml')

        kinetics = KineticMechanism(pathkin)
        kinetics.ReadKinetics(pathname)
        Stoichiometry = ReadStoichiometry(pathkin) # leggere la versione e fare check come su opensmoke

        return kinetics.NumberOfSpecies, kinetics.NumberOfReactions, Stoichiometry.nur, Stoichiometry.nup, kinetics.reaction_lines, kinetics.species

    def ReadOutputXML(self, outputFolder):
        
        outputFile = os.path.join(outputFolder, 'Output.xml')

        if not os.path.isfile(outputFile):
            raise ValueError('XML OUTPUT FILE DOES NOT EXIST!')

        # tree = ET.ElementTree(ET.fromstring(string))
        tree_out = ET.parse(outputFile)
        root_out = tree_out.getroot()                                                     
        
        massfractions = root_out.find('mass-fractions')
        massfractions = massfractions.text.split()

        profiles = root_out.find('profiles')
        profiles = profiles.text.split()

        profileSize = root_out.find('profiles-size')
        profileSize = profileSize.text.split()
        
        NumberOfRows = int(profileSize[0])
        NumberOfColumns = int(profileSize[1])
        
        profiles = np.reshape(profiles, (NumberOfRows, int(len(profiles) / NumberOfRows)))
        profiles = profiles.astype(np.float)
    
        time = profiles[:,0]
        
        gas_reaction_rates = root_out.find('reaction-rates') 
        gas_reaction_rates = np.float64(gas_reaction_rates.text.split()) # questo non è un vettore ma bensì una matrice di dimensione
                                                                         # NumberOfReactions x altra dimensione
                                                                         # altra dimensione la trovo in profile size dove primo numero
                                                                         # numero di righe e secondo numero numero di colonne

        matrixReactionRates = np.reshape(gas_reaction_rates,(NumberOfRows, self.numberOfReactions))
    
        gas_formation_rates = root_out.find('formation-rates') 
        gas_formation_rates = np.float64(gas_formation_rates.text.split())

        return matrixReactionRates, time

    def ComputeRopa(self, Type: str, region: list=[], localValue: float = 0):
        
        def ReturnCoefficients():
            Production, Consumption = self.CalculateCoefficientsNotNormalized()

            coefficientProduction = self.ComputeCoefficients(matrix = Production)
            coefficientConsumption = self.ComputeCoefficients(matrix = Consumption)

            threshold = 1e-20
            C = coefficientConsumption
            P = coefficientProduction
            C = [0 if C_ < threshold or C_>= 1 else C_ for C_ in C]
            P = [0 if P_ < threshold or P_>= 1 else P_ for P_ in P]
            return C, P

        if Type == 'Global':
            """
            Description: TODO
            """
            return ReturnCoefficients()
        
        elif Type == 'Region':
            """
            Description: TODO
            """
            tmin = region[0]
            tmax = region[1]

            valuemin, idxmin = self.FindNearest(array = self.time, value = tmin)
            self.time = self.time[idxmin:]
            self.matrixReactionRates = self.matrixReactionRates[idxmin:,:]
            valuemax, idxmax = self.FindNearest(array = self.time, value = tmax)
            self.time = self.time[:idxmax+1]
            self.matrixReactionRates = self.matrixReactionRates[:idxmax+1,:]

            return ReturnCoefficients()
        
        elif Type == 'Local':
            """
            Description: TODO
            """
            timeValue, idx = self.FindNearest(array = self.time, value = localValue)
            self.matrixReactionRates = self.matrixReactionRates[idx]
            Production, Consumption = self.CalculateLocalCoefficientsNotNormalized()

            localProduction = self.ComputeLocalCoefficients(vector= Production)
            localConsumption = self.ComputeLocalCoefficients(vector= Consumption)

            C = localConsumption
            P = localProduction

            return C, P            
        else:
            raise ValueError('The selected type is not supported. Available are: Global | Local | Region')

    def GetIndexSpecies(self, speciesNames):
        for position, species in enumerate(self.speciesNames):
            if species == speciesNames:
                index_species = position
                break
            else:
                continue
        
        return index_species

    def CalculateCoefficientsNotNormalized(self):
        
        Production = np.empty((self.nur.shape[0], len(self.time)))
        Consumption = np.empty((self.nup.shape[0], len(self.time)))

        nup_arr = self.nup.toarray()
        nur_arr = self.nur.toarray()  
        
        for k in range(len(self.time)):
            for j in range(self.numberOfReactions):
                Production[j,k] = nup_arr[j,self.indexSpecies] * self.matrixReactionRates[k,j]
                Consumption[j,k] = nur_arr[j,self.indexSpecies] * self.matrixReactionRates[k,j]
        
        Production, Consumption = self.MatrixSwapping(m1 = Production, m2 = Consumption)
      
        return Production, Consumption

    def CalculateLocalCoefficientsNotNormalized(self):
        
        Production = np.empty((self.nur.shape[0]))
        Consumption = np.empty((self.nup.shape[0]))

        nup_arr = self.nup.toarray()
        nur_arr = self.nur.toarray()  
        
        for j in range(self.numberOfReactions):
            Production[j] = nup_arr[j,self.indexSpecies] * self.matrixReactionRates[j]
            Consumption[j] = nur_arr[j,self.indexSpecies] * self.matrixReactionRates[j]
        
        Production, Consumption = self.VectorSwapping(v1 = Production, v2 = Consumption)
      
        return Production, Consumption

    def MatrixSwapping(self, m1, m2):

        for k in range(len(self.time)):
            for j in range(self.numberOfReactions):
                if m1[j,k]< 0:
                    m2[j,k] = -m1[j,k]
                    m1[j,k] = 0
                elif m2[j,k] < 0:
                    m1[j,k] = -m2[j,k]
                    m2[j,k] = 0
                else:
                    continue

        return m1, m2        

    def VectorSwapping(self, v1, v2):

        for j in range(self.numberOfReactions):
            if v1[j]< 0:
                v2[j] = -v1[j]
                v1[j] = 0
            elif v2[j] < 0:
                v1[j] = -v2[j]
                v2[j] = 0
            else:
                continue
        return v1, v2

    def ComputeCoefficients(self, matrix):

        coefficient = np.empty((self.numberOfReactions, len(self.time)))
        computedCoefficient = np.empty(self.numberOfReactions)
        delta = self.time[-1] - self.time[0]

        for i in range(len(self.time)-1):
            dt = (self.time[i+1] - self.time[i]) / delta
            for j in range(self.numberOfReactions):
                coefficient[j,i] = dt * matrix[j,i]

        for i in range(self.numberOfReactions):
            computedCoefficient[i] = sum(coefficient[i,:])
        
        return computedCoefficient

    def ComputeLocalCoefficients(self, vector):

        coefficient = np.empty((self.numberOfReactions))
        
        for j in range(self.numberOfReactions):
            coefficient[j] = vector[j]
        
        return sum(coefficient)

    def FindNearest(self,array, value):

        idx = (np.abs(array - value)).argmin()

        return array[idx], idx 
