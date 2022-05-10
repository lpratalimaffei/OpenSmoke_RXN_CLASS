'''
MODULE: ReadStoichiometry
@Authors:
    Alberto Cuoci [1], Timoteo Dinelli [1]
    [1]: CRECK Modeling Lab, Department of Chemistry, Materials, and Chemical Engineering, Politecnico di Milano
@Contacts:
    alberto.cuoci@polimi.it
    timoteo.dinelli@polimi.it
@Additional notes:
    This code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
    Please report any bug to: alberto.cuoci@polimi.it
'''

# Import external libraries
import numpy as np
import xml.etree.ElementTree as ET
from scipy import sparse
import sys

# Import internal modules
from KineticMechanism import KineticMechanism

class ReadStoichiometry:

    def __init__(self, kinetics_xml):
        
        self.kinetics = KineticMechanism(kinetics_xml)
        self.ReadMechFile(kinetics_xml)
        self.stoichiometryProcessing()
        self.nur = self.StoichiometryCoeffReactants(self.kinetics)
        self.nup = self.StoichiometryCoeffProducts(self.kinetics)

    def ReadVectorInt(self, vector, pos):
        length = int(vector[pos])
        subvector = vector[pos+1:pos+length+1]
        return np.int32(subvector), pos+length+1

    def ReadVectorFloat64(self, vector, pos):
        length = int(vector[pos])
        subvector = vector[pos+1:pos+length+1]
        return np.float64(subvector), pos+length+1

    def ReadMechFile(self, kinetics_xml):

        tree = ET.parse(kinetics_xml)
        root = tree.getroot()
        kinetics_element = root.find('Kinetics')
        stoichiometry_element = kinetics_element.find('Stoichiometry')
        self.stoichiometry = (stoichiometry_element.text).split()       

    def stoichiometryProcessing(self):
        
        stoichiometry = self.stoichiometry
        
        pos = 0
        numDir1, pos = self.ReadVectorInt(stoichiometry, pos)
        numDir2, pos = self.ReadVectorInt(stoichiometry, pos)
        numDir3, pos = self.ReadVectorInt(stoichiometry, pos)
        numDir4, pos = self.ReadVectorInt(stoichiometry, pos)
        numDir5, pos = self.ReadVectorInt(stoichiometry, pos)

        numRevTot1, pos = self.ReadVectorInt(stoichiometry, pos)
        numRevTot2, pos = self.ReadVectorInt(stoichiometry, pos)
        numRevTot3, pos = self.ReadVectorInt(stoichiometry, pos)
        numRevTot4, pos = self.ReadVectorInt(stoichiometry, pos)
        numRevTot5, pos = self.ReadVectorInt(stoichiometry, pos)

        numRevEq1, pos = self.ReadVectorInt(stoichiometry, pos)
        numRevEq2, pos = self.ReadVectorInt(stoichiometry, pos)
        numRevEq3, pos = self.ReadVectorInt(stoichiometry, pos)
        numRevEq4, pos = self.ReadVectorInt(stoichiometry, pos)
        numRevEq5, pos = self.ReadVectorInt(stoichiometry, pos)

        jDir1, pos = self.ReadVectorInt(stoichiometry, pos)
        jDir2, pos = self.ReadVectorInt(stoichiometry, pos)
        jDir3, pos = self.ReadVectorInt(stoichiometry, pos)
        jDir4, pos = self.ReadVectorInt(stoichiometry, pos)
        jDir5, pos = self.ReadVectorInt(stoichiometry, pos)
        valueDir5, pos = self.ReadVectorFloat64(stoichiometry, pos)

        jDir1 = jDir1-1
        jDir2 = jDir2-1
        jDir3 = jDir3-1
        jDir4 = jDir4-1
        jDir5 = jDir5-1

        jRevTot1, pos = self.ReadVectorInt(stoichiometry, pos)
        jRevTot2, pos = self.ReadVectorInt(stoichiometry, pos)
        jRevTot3, pos = self.ReadVectorInt(stoichiometry, pos)
        jRevTot4, pos = self.ReadVectorInt(stoichiometry, pos)
        jRevTot5, pos = self.ReadVectorInt(stoichiometry, pos)
        valueRevTot5, pos = self.ReadVectorFloat64(stoichiometry, pos)

        jRevTot1 = jRevTot1-1
        jRevTot2 = jRevTot2-1
        jRevTot3 = jRevTot3-1
        jRevTot4 = jRevTot4-1
        jRevTot5 = jRevTot5-1

        jRevEq1, pos = self.ReadVectorInt(stoichiometry, pos)
        jRevEq2, pos = self.ReadVectorInt(stoichiometry, pos)
        jRevEq3, pos = self.ReadVectorInt(stoichiometry, pos)
        jRevEq4, pos = self.ReadVectorInt(stoichiometry, pos)
        jRevEq5, pos = self.ReadVectorInt(stoichiometry, pos)
        valueRevEq5, pos = self.ReadVectorFloat64(stoichiometry, pos)

        jRevEq1 = jRevEq1-1
        jRevEq2 = jRevEq2-1
        jRevEq3 = jRevEq3-1
        jRevEq4 = jRevEq4-1
        jRevEq5 = jRevEq5-1

        changeOfMoles, pos = self.ReadVectorFloat64(stoichiometry, pos)
        explicit_reaction_orders = int(stoichiometry[pos])

        # Elementary reactions only
        if (explicit_reaction_orders == 0):
            
            lambda_numDir1 = numDir1;
            lambda_numDir2 = numDir2;
            lambda_numDir3 = numDir3;
            lambda_numDir4 = numDir4;
            lambda_numDir5 = numDir5;

            lambda_numRevEq1 = numRevEq1;
            lambda_numRevEq2 = numRevEq2;
            lambda_numRevEq3 = numRevEq3;
            lambda_numRevEq4 = numRevEq4;
            lambda_numRevEq5 = numRevEq5;

            lambda_jDir1 = jDir1;
            lambda_jDir2 = jDir2;
            lambda_jDir3 = jDir3;
            lambda_jDir4 = jDir4;
            lambda_jDir5 = jDir5;
            lambda_valueDir5 = valueDir5;

            lambda_jRevEq1 = jRevEq1;
            lambda_jRevEq2 = jRevEq2;
            lambda_jRevEq3 = jRevEq3;
            lambda_jRevEq4 = jRevEq4;
            lambda_jRevEq5 = jRevEq5;
            lambda_valueRevEq5 = valueRevEq5;
            
        else:

            sys.exit("Non-elementary reactions cannot be processed")

        
        self.numDir1 = numDir1
        self.numDir2 = numDir2
        self.numDir3 = numDir3
        self.numDir4 = numDir4
        self.numDir5 = numDir5

        self.numRevTot1 = numRevTot1
        self.numRevTot2 = numRevTot2
        self.numRevTot3 = numRevTot3
        self.numRevTot4 = numRevTot4
        self.numRevTot5 = numRevTot5

        self.numRevEq1 = numRevEq1
        self.numRevEq2 = numRevEq2
        self.numRevEq3 = numRevEq3
        self.numRevEq4 = numRevEq4
        self.numRevEq5 = numRevEq5

        self.jDir1 = jDir1
        self.jDir2 = jDir2
        self.jDir3 = jDir3
        self.jDir4 = jDir4
        self.jDir5 = jDir5

        self.jRevTot1 = jRevTot1
        self.jRevTot2 = jRevTot2
        self.jRevTot3 = jRevTot3
        self.jRevTot4 = jRevTot4
        self.jRevTot5 = jRevTot5
        self.valueRevTot5 = valueRevTot5

        self.jRevEq1 = jRevEq1
        self.jRevEq2 = jRevEq2
        self.jRevEq3 = jRevEq3
        self.jRevEq4 = jRevEq4
        self.jRevEq5 = jRevEq5
        self.valueRevEq5 = valueRevEq5

        self.lambda_numDir1 = lambda_numDir1
        self.lambda_numDir2 = lambda_numDir2
        self.lambda_numDir3 = lambda_numDir3
        self.lambda_numDir4 = lambda_numDir4
        self.lambda_numDir5 = lambda_numDir5

        self.lambda_numRevEq1 = lambda_numRevEq1
        self.lambda_numRevEq2 = lambda_numRevEq2
        self.lambda_numRevEq3 = lambda_numRevEq3 
        self.lambda_numRevEq4 = lambda_numRevEq4 
        self.lambda_numRevEq5 = lambda_numRevEq5

        self.lambda_jDir1 = lambda_jDir1
        self.lambda_jDir2 = lambda_jDir2
        self.lambda_jDir3 = lambda_jDir3
        self.lambda_jDir4 = lambda_jDir4
        self.lambda_jDir5 = lambda_jDir5
        self.lambda_valueDir5 = lambda_valueDir5   

        self.lambda_jRevEq1 = lambda_jRevEq1
        self.lambda_jRevEq2 = lambda_jRevEq2
        self.lambda_jRevEq3 = lambda_jRevEq3
        self.lambda_jRevEq4 = lambda_jRevEq4
        self.lambda_jRevEq5 = lambda_jRevEq5
        self.lambda_valueRevEq5 = lambda_valueRevEq5

    def StoichiometryCoeffReactants(self, kinetics):

        react_species = []
        react_reaction = []
        react_nu = []

        count1=0
        count2=0
        count3=0
        count4=0
        count5=0
        for i in range(kinetics.NumberOfSpecies):

            for k in range(self.numDir1[i]):
                react_species.append(i)
                react_reaction.append(self.jDir1[count1])
                react_nu.append(1.)
                count1 = count1+1
                
            for k in range(self.numDir2[i]):
                react_species.append(i)
                react_reaction.append(self.jDir2[count2])
                react_nu.append(2.)
                count2 = count2+1
                
            for k in range(self.numDir3[i]):
                react_species.append(i)
                react_reaction.append(self.jDir3[count3])
                react_nu.append(3.)
                count3 = count3+1
                
            for k in range(self.numDir4[i]):
                react_species.append(i)
                react_reaction.append(self.jDir4[count4])
                react_nu.append(0.5)
                count4 = count4+1
                
            for k in range(self.numDir5[i]):
                react_species.append(i)
                react_reaction.append(self.jDir5[count5])
                react_nu.append(self.valueDir5[count5])
                count5 = count5+1

        self.nur = sparse.coo_matrix((react_nu,(react_reaction,react_species)),shape=(kinetics.NumberOfReactions,kinetics.NumberOfSpecies))

        return self.nur

    def StoichiometryCoeffProducts(self, kinetics):  

        react_species = []
        react_reaction = []
        react_nu = []

        count1=0
        count2=0
        count3=0
        count4=0
        count5=0
        for i in range(kinetics.NumberOfSpecies):

            for k in range(self.numRevTot1[i]):
                react_species.append(i)
                react_reaction.append(self.jRevTot1[count1])
                react_nu.append(1.)
                count1 = count1+1
                
            for k in range(self.numRevTot2[i]):
                react_species.append(i)
                react_reaction.append(self.jRevTot2[count2])
                react_nu.append(2.)
                count2 = count2+1
                
            for k in range(self.numRevTot3[i]):
                react_species.append(i)
                react_reaction.append(self.jRevTot3[count3])
                react_nu.append(3.)
                count3 = count3+1
                
            for k in range(self.numRevTot4[i]):
                react_species.append(i)
                react_reaction.append(self.jRevTot4[count4])
                react_nu.append(0.5)
                count4 = count4+1
                
            for k in range(self.numRevTot5[i]):
                react_species.append(i)
                react_reaction.append(self.jRevTot5[count5])
                react_nu.append(self.valueRevTot5[count5])
                count5 = count5+1

        self.nup = sparse.coo_matrix((react_nu,(react_reaction,react_species)),shape=(kinetics.NumberOfReactions,kinetics.NumberOfSpecies))

        return self.nup

    def ReactionOrderReactants(self,kinetics):

        react_species = []
        react_reaction = []
        react_lambda = []

        count1=0
        count2=0
        count3=0
        count4=0
        count5=0
        for i in range(kinetics.NumberOfSpecies):

            for k in range(self.lambda_numDir1[i]):
                react_species.append(i)
                react_reaction.append(self.lambda_jDir1[count1])
                react_lambda.append(1.)
                count1 = count1+1
                
            for k in range(self.lambda_numDir2[i]):
                react_species.append(i)
                react_reaction.append(self.lambda_jDir2[count2])
                react_lambda.append(2.)
                count2 = count2+1
                
            for k in range(self.lambda_numDir3[i]):
                react_species.append(i)
                react_reaction.append(self.lambda_jDir3[count3])
                react_lambda.append(3.)
                count3 = count3+1
                
            for k in range(self.lambda_numDir4[i]):
                react_species.append(i)
                react_reaction.append(self.lambda_jDir4[count4])
                react_lambda.append(0.5)
                count4 = count4+1
                
            for k in range(self.lambda_numDir5[i]):
                react_species.append(i)
                react_reaction.append(self.lambda_jDir5[count5])
                react_lambda.append(self.lambda_valueDir5[count5])
                count5 = count5+1

        self.lambdar = sparse.coo_matrix((react_lambda,(react_reaction,react_species)),shape=(kinetics.NumberOfReactions,kinetics.NumberOfSpecies))

        return self.lambdar

    def ReactionOrderProducts(self,kinetics):  

        react_species = []
        react_reaction = []
        react_lambda = []

        count1=0
        count2=0
        count3=0
        count4=0
        count5=0
        for i in range(kinetics.NumberOfSpecies):

            for k in range(self.lambda_numRevEq1[i]):
                react_species.append(i)
                react_reaction.append(self.lambda_jRevEq1[count1])
                react_lambda.append(1.)
                count1 = count1+1
                
            for k in range(self.lambda_numRevEq2[i]):
                react_species.append(i)
                react_reaction.append(self.lambda_jRevEq2[count2])
                react_lambda.append(2.)
                count2 = count2+1
                
            for k in range(self.lambda_numRevEq3[i]):
                react_species.append(i)
                react_reaction.append(self.lambda_jRevEq3[count3])
                react_lambda.append(3.)
                count3 = count3+1
                
            for k in range(self.lambda_numRevEq4[i]):
                react_species.append(i)
                react_reaction.append(self.lambda_jRevEq4[count4])
                react_lambda.append(0.5)
                count4 = count4+1
                
            for k in range(self.lambda_numRevEq5[i]):
                react_species.append(i)
                react_reaction.append(self.lambda_jRevEq5[count5])
                react_lambda.append(self.lambda_valueRevEq5[count5])
                count5 = count5+1

        self.lambdap = sparse.coo_matrix((react_lambda,(react_reaction,react_species)),shape=(kinetics.NumberOfReactions,kinetics.NumberOfSpecies))

        return self.lambdap