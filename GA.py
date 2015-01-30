
#!/usr/bin/env python

__author__      = "Luis E Lara-Gonzalez, David Goode, Andrew Bakshi"
__copyright__   = "Copyright 2015, Tumour Evolution and Drug Resistance"
__credits__     = ["Luis E Lara-Gonzalez", "David Goode"]
__license__     = "GPL"
__version__     = "1.0.0"
__maintainer__  = "Luis E Lara-Gonzlez"
__email__       = "llara@student.unimelb.edu.au"
__status__      = "Prototype"

import pylab as pl
import numpy as np
import math
import main
import sys
import os
import random
import operator
import socket
from itertools import combinations
import itertools
import pygal
from pygal.style import Style
import pandas
from pandas import DataFrame, Series
import scipy, scipy.stats
import datetime
import time
import statsmodels.formula.api as sm
import polygon
#from mpi4py import MPI


'''
	This class is the basic object to kkep track of the inheritance of the
	genetic algorithm.

	The Paret flag will indicate if the element is Parent of children for a given time t.
	The Generation will indicate the generation path.
	The Genome is a dictionary data structure that has the simulation running parametrs.


	If tou want to modify some initial parameters, you can modify those in the the calling of the GA.<object>
	or in the deafault sefinition section.

	A commnad line parser shoul be implemented.

'''	
class Family(object):
	
	'''	#######################################							
		Define Init Parametersa and Constructor
	'''	#######################################

	def __init__(self,
		Parent 		=			True 	, 	#am I Parent or child
		Generation	=		   	   '0' 	, 	#Number of children
		genome 		= 			dict()	   #Simulation Parameters 
 		):

		self.p = Parent
		self.g = Generation
		self.genome = dict()
	

	''' #################
		Function Section
	''' #################

	def writeGenome(self, dict):
		self.genome = dict

	def printGeneration(self):
		return self.g

	def printGenome(self):
		print('\n|-------[ . . . G E N O M E . . . ]-------|')
		print self.genome

''' ##############################################################
	This class will initialise and the fine the structre of the GA
''' ##############################################################
class GA(object):

	'''	#######################################							
		Define Init Parametersa and Constructor
		#######################################

		def __init__(...) will initiliase those relevant parameters if those are not given
	'''
	def __init__(self,
		# This are DEFAULT values if not given			
		ProliferationRate 	   =          0.04 , 	#From 0 to 1 ? 0.01 to 0.1
		DeathRate 		  	   =          0.03 , 	#From 0 to 1 ?
		MutationRate      	   =          0.001, 	#From 0 to 1 ?
		loops 			  	   =      10000    ,   #Let's keep this fixed - 10000000
		MaxSize 		       = 	    10000    ,	#From - 1000000000  
		InitSize 		  	   =         25    ,	#Limits?
		InitialDiversity  	   =          1    ,	#Worthy to randomize?
		SelectivePressure 	   =          0.01 , 	#Maybe if we found good parameters we can vary this
		SelectTime 		 	   =     400000    ,	#This are loops?
		P_BenMutation 		   =          0.01 , 	#We must randomize this parameter 
		P_DelMutation 		   =          0.99 ,   #We must randomize this parameter
		P_IncMutRate 		   =          0.0  ,	#Does this value work?
		P_DecMutRate 		   =          0.0  ,	#Does this value work?
		Permutation 		   =        3    , 	#What is this?
		scale 				   =          0.5  ,	#Scale of what?
		MutScale 			   =          1.0  ,	#??
		PairsofParents 		   =			  1 	, 	#Initial number of parents
		OffspringSize		   =		    10 	, 	#Number of children
		Randomize 			   =  'scratch'	,   #True Randomize from scratch, False randomize from baseline or read from file
		IsolationMode 		   = 		False 	, 	#If true instances will communicate on cloud
		Scheme 				   = 	   'Smooth' , 	#Searching scheme. Agressive or Smooth -> Cumulative
		Exploration			   =	  'basic4'	,	#Define a exploration pattern: basic4 will modify only the P_BEN,P_Del, P_Inc, P_Dec
		NetworkModel		   =	'Individual',	#Defines the network compuation scheme
		MutKernel 			   = 		   	'U',	#option dict
		Kernelparameters 	   = 	  [0.01,0.5,0.1,0.99,0.01,0.5,0.01,0.5], 	#Paraters to search
		indexes 			   = 	   [0,2,4,6],   #indexes of the parameters to search
		P_SNP 				   = 			  0.4 ,   #Probability of introducing a SNP in children
		P_CrossOver			   =			  0.1 ,   #Probability of CrossOver
		P_TransLocation 	   = 			  0.05, 	#Probability to introduce a dramatic mutation
		PopulationSize  	   = 		   100 	, 	#Define a population size
		ConstantPopSize 	   = 		True 	, 	#Force to always have the same population size
		crossoverFunc 		   = 	'uniform'	,	#Function to implement crossover
		MaxIterations 		   = 	   10000 	, 	#Max number of iterations
		HaltFunction		   = 		'smooth',	#Criteria to stop the running
		MatingModel 		   = 	     'HiLow',	#Model to generate more offpring
		SelectionPercentage    = 		  25    ,   #Perecentage of extarction from best members
		Tournament 			   =  		False	,	#Compete within
		NextOffspring 		   = 			2

		):



		""" Constructor """
		self.pr 		= ProliferationRate
		self.dr 		= DeathRate
		self.mr 		= MutationRate
		self.pop 		= PairsofParents
		self.l 			= loops
		self.ms			= MaxSize
		self.insz		= InitSize
		self.indiv 		= InitialDiversity
		self.st 		= SelectTime
		self.sp 		= SelectivePressure
		self.pbm 		= P_BenMutation
		self.pdm 		= P_DelMutation
		self.imr 		= P_IncMutRate
		self.dmr 		= P_DecMutRate
		self.perm 		= Permutation
		self.sc 		= scale
		self.msc 		= MutScale
		self.rand 		= Randomize
		self.expl 		= Exploration
		self.mk 		= MutKernel
		self.kp 		= Kernelparameters
		self.idx 		= indexes
		self.os 		= OffspringSize
		self.crss 		= crossoverFunc
		self.psnp 		= P_SNP
		self.pco 		= P_CrossOver
		self.ptrans 	= P_TransLocation
		self.sch 		= Scheme
		self.topScores 	= dict()
		self.MatMod 	= MatingModel
		self.selper 	= SelectionPercentage
		self.nextGen_id = 0
		self.currentGen_id = 0
		self.HOF = dict()
		self.FamMem = []
		self.Population = dict()
		self.Best = []
		self.Worst = []
		self.Radar_Data = dict()
		self.BestFitRadar = dict()
		self.path_hof = ''
		self.path_Summary = ''
		self.path_Regression = ''
		self.path_Plots = ''
		self.path_df = ''
		self.nextOff = NextOffspring
		self.FitnessPopulation_List = dict()

		#self.createAdan(0,0)
		#self.createEve()

	''' #################
		Function Section
	''' #################

	'''#######################################
		Define Sampling Distribution Functions
	'''#######################################
	def Z(parameter):
		print(parameter)

	def U(value, parameters, idx):
		return round(np.random.uniform(parameters[idx], parameters[idx + 1]), 4) 

	global distribution 
	distribution = {'Z': Z, 'U': U}

	
	'''
		This function print the values of the class
	'''
	def printValues(self):
		print self.pr, self.dr

	
	'''
		This function creates the first Parent, further functions will take
		care of the input considerations (such as randomize, correlation, ...).

		Execution v_1.0: 
		1) Create Parents as a Famility Object
		2) Initilise Adan according with input parameters
		3) Print Adan's genome 
	'''
	def createParents(self):
		print('|-------[ . . . C R E A T I N G  P A R E N T S . . . ]-------|')
		Parents = list()
		for i in range(0, self.pop):
			print('>> Adan_' + str(self.currentGen_id)  )
			Adan = Family(Parent = True, Generation = 'A'+str(self.currentGen_id)+':'+self.child_ID_NextGen()) # 1) 
			self.initFamilyMemeber( Adan )				 	  # 2)
			Adan.printGenome()							      # 3)
			Parents.append(Adan)		
			print('>> Eve_' + str(self.currentGen_id) )
			Eve = Family( Parent = True, Generation = 'E'+str(self.currentGen_id)+':'+self.child_ID_NextGen() )
			self.initFamilyMemeber( Eve )
			Eve.printGenome()
			Parents.append(Eve)

		return Parents

	'''
		This function will initiliase the Familily Object.

		Execution v_1.0:
		1) Initilise genome by calling insertRunningParameters()
		2) Configure and set the mutation scheme

	'''
	def initFamilyMemeber(self, FamilyMemeber):

		self.insertRunningParameters(FamilyMemeber)
		print('\n|-------[ . . . R A N D O M I S A T I O N . . . ]-------|')
		if self.rand == 'baseline':
			print('>> baseline')
			self.mutateBaseline(FamilyMemeber)
		elif self.rand == 'scratch':
			print('>> scratch')
			self.mutateScratch(FamilyMemeber)
		elif self.rand == 'file':
			print('>> from file')

	'''
		This function will create the genome of the Family object and write on it

		Execution v_1.0:
		1) Create dictionary
		2) Write on Object Family
	'''
	def insertRunningParameters(self, FamilyMemeber ):
		genome = {
				'ProliferationRate' : self.pr, 
				'DeathRate' 		: self.dr,
				'MutationRate'		: self.mr,
				'loops' 			: self.l,
				'MaxSize'			: self.ms,
				'InitSize' 			: self.insz,
				'InitialDiversity'	: self.indiv,
				'SelectTime' 		: self.st,
				'SelectivePressure'	: self.sp,
				'P_BenMutation'		: self.pbm,
				'P_DelMutation' 	: self.pdm,
				'P_IncMutRate' 		: self.imr,
				'P_DecMutRate' 		: self.dmr,
				'Permutation' 		: self.perm,
				'scale' 			: self.sc,
				'MutScale'			: self.msc
				}									# 1)

		FamilyMemeber.writeGenome(genome)			# 2)
		
	'''
		This function mutates the Family Member from baseline considering
		the Exploration model

		Execution v_1.0:
		1) 

	'''	
	def mutateBaseline(self, FamilyMemeber):
		if self.expl == 'basic4':
			self.MutateBasic4(FamilyMemeber)
			
	'''
		This function mutates from scratch, this means it will mutate using only as a paramter 
		the distribition and range.

		Execution v_1.0:
		1) Identify exploration pattern.
		2) Call mutation fucntion according to exploration scheme
	'''
	def mutateScratch(self, FamilyMemeber):
		print('\n|-------[ . . . E X P L O R A T I O N . . . ]-------| ')
		if self.expl == 'basic4':
			print('>> Basic 4 Probabilities')
			self.MutateBasic4(FamilyMemeber)
		elif self.expl == 'AllParameters':
			print('>> All Parameters')

	'''
		This function mutates the parameters:
		* P_BenMutation
		* P_DelMutation
		* P_DecMutRate
		* P_IncMutRate

		Execution v_1.0:
		1) Determine randomize pattern
		2) Randomize according to criteria and write into Familymemeber object

	'''
	def MutateBasic4(self, FamilyMemeber):
		print('\n|-------[ . . . R A N D O M I Z I N G . . . ]-------| ')
		print('>> Initial guess: ')
		if self.rand == 'scratch':
			out = True
			P_Ben = distribution[self.mk]( FamilyMemeber.genome['P_BenMutation'], self.kp, self.idx[0]  )
			P_Del = distribution[self.mk]( FamilyMemeber.genome['P_DelMutation'], self.kp, self.idx[1]  )
			P_Dec = distribution[self.mk]( FamilyMemeber.genome['P_DecMutRate'], self.kp, self.idx[2]  )
			P_Inc = distribution[self.mk]( FamilyMemeber.genome['P_IncMutRate'], self.kp, self.idx[3]  )
			if P_Ben + P_Del < 1 and P_Dec + P_Inc < 1:
				print('--> Good to Go', P_Ben, P_Del, P_Dec, P_Inc)
				FamilyMemeber.genome['P_BenMutation'] = P_Ben
				FamilyMemeber.genome['P_DelMutation'] = P_Del
				FamilyMemeber.genome['P_DecMutRate'] = P_Dec
				FamilyMemeber.genome['P_IncMutRate'] = P_Inc
			else:
				print('--> Invalid: Looking for new parameters', P_Ben, P_Del, P_Dec, P_Inc)
				while(out):
					P_Ben = distribution[self.mk]( FamilyMemeber.genome['P_BenMutation'], self.kp, self.idx[0]  )
					P_Del = distribution[self.mk]( FamilyMemeber.genome['P_DelMutation'], self.kp, self.idx[1]  )
					P_Dec = distribution[self.mk]( FamilyMemeber.genome['P_DecMutRate'], self.kp, self.idx[2]  )
					P_Inc = distribution[self.mk]( FamilyMemeber.genome['P_IncMutRate'], self.kp, self.idx[3]  )
					out = P_Ben + P_Del > 1 and P_Dec +  P_Dec > 1
					print ('Values:', P_Ben, P_Del,P_Dec,P_Inc, out)
				FamilyMemeber.genome['P_BenMutation'] = P_Ben
				FamilyMemeber.genome['P_DelMutation'] = P_Del
				FamilyMemeber.genome['P_DecMutRate'] = P_Dec
				FamilyMemeber.genome['P_IncMutRate'] = P_Inc

		elif self.rand == 'baseline':
			pass

	'''
		This function will generate the offspring of size self.os 

	'''
	def createOffspring(self, Parents):
		print('\n|-------[ . . . M A T I N G . . . ]-------| ')
		print('>> Generating ' + str(self.os) + ' children ')
		print('>> ' + Parents[0].g + ' with ' + Parents[1].g)
		c_id = self.newChildrenID(Parents)
		offspring = self.newChildren(c_id)
		mutations = self.mutationLists()
		self.mate(mutations[3], offspring, Parents)
		self.crossOver(mutations[1], offspring, Parents)
		self.Translocation(mutations[2], offspring)
		self.PointMutations(mutations[0], offspring)
		self.inc_CurrentGen()
		return offspring

	'''
		This function prints the basic4 elements in the Offspring
	'''
	def printOffspringB4(self, Children):
		print('\n|-------[ . . . CHILDREN BASIC 4 . . . ]-------| ')
		for child,i in  zip(Children, range(len(Children))):
			print('--------'+str(i)+'--------')
			print('P_BEN: ',child.genome['P_BenMutation'] )	
			print('P_DEL: ',child.genome['P_DelMutation'] )	
			print('P_INC: ',child.genome['P_IncMutRate'] )	
			print('P_DEC: ',child.genome['P_DecMutRate'])
			print('\n|--------------| ')

	'''
		This function is the first layer on the Translocation event,
		this function calls the TranslocateB4(...)
	'''
	def Translocation(self, MutationList, Children):
		if self.expl == 'basic4':
			self.TranslocateB4(MutationList, Children)
		else: # implement the rest of the schmes
			pass

	'''
		This function will call doTranslocationB4(...) to instroduce a
		translocation on selected memebers.
	'''
	def TranslocateB4(self, MutationList, Children):
		for mutation, child in zip(MutationList, Children):
			if mutation:
				print ('Translocating ' + child.g)
				self.doTranslocationB4(child)
		
	'''
		This function will introduce the translocation in the Basic 4 scheme
		by sampling from aa binomial distribution. If the sample form 
		the Bionomial is one  the translocation is in the
		P_Beneficial <-> P_Deleterious.
		If the value is 0 the the translocation is 
		P_Inc_Mutation_Rate <-> P_Dec_Mutation_Rate

	'''
	def doTranslocationB4(self, FamilyMemeber):
		if np.random.binomial(1, 0.5) == 1:
			p_ben = FamilyMemeber.genome['P_BenMutation']
			p_del = FamilyMemeber.genome['P_DelMutation']
			FamilyMemeber.genome['P_BenMutation'] = p_del
			FamilyMemeber.genome['P_DelMutation'] = p_ben
			print ('P_BEN: ' + str(p_ben) + ' -> P_DEL ' + str(p_del))
		else:
			p_inc = FamilyMemeber.genome['P_IncMutRate']
			p_dec = FamilyMemeber.genome['P_DecMutRate']
			FamilyMemeber.genome['P_IncMutRate'] = p_dec
			FamilyMemeber.genome['P_DecMutRate'] = p_inc
			print ('P_INC: ' + str(p_inc) + ' -> P_DEC ' + str(p_dec))
		print 'Done'


	'''
		This Function handles identifies which is the corresponfig case of exploration then,
		call the appropiate function.
	'''
	def PointMutations(self, mutationList, Children):
		if self.expl == 'basic4':
			self.SNP_B4(mutationList, Children)
		elif self.expl == 'other': # Modigy this
			pass

	'''
		This function is specific for the Basic 4 exploration, 
		if a child has a mutation on genome then mutate.
		Execution v_1.0:
		1) loop the SNP and children lists
		2) if the children has a snp
		3) Introduce mutation calling doSNPs(...)
	'''
	def SNP_B4(self, mutationList, Children):
		for snp, child in zip(mutationList, Children):
			if snp:
				self.doSNPs_B4(child)

	'''
		This function forces to SNPs in child's genome.
		Execution v_1.0:
		1) Generate a lsit of two arrays of 0s and 1s of size 4 (due to Basic 4 exploration)
		2) loop over those arrays and if 1 is present call introduce snp
		
		** Note: The calling of introduceSNP_B4(...) can be kernelizable, BUT, I think mutating from
		a Gaussian distribution is the best approach. Consider this in further modifications!

	'''
	def doSNPs_B4(self, child):
		print('Introducing SNPS in ' + child.g)
		snps_list =  np.random.multinomial(1, [1/4.]*4, size=2)
		for loci in snps_list:
			for locus,i in zip(loci, range(len(loci))):
				if locus == 1:
					self.introduceSNP_B4(child, 'Gaussian', i)
		print('Done')

	'''
		This function mutate a loci using different sampling distributions, in the case of the Gaussian
		the mean of the distributin is associated with the parameter value and the variance will depend
		on the GA scheme. Smooth only induces gentle mutations wheares an Aggressive scheme can induce
		more dramatic effects.
		Execution v_1.0:
		1) Identify Kernel
		2) Mutate with respective kernel
			2.a) For Gaussian:
				* round the new value to 4 digits -> obtain a R.V. with mean = child.genome['Parameter_Value']
				* and variance -> according to Scheme [0.1, 0.3]
		3) Validate
	'''
	def introduceSNP_B4(self, FamilyMemeber, Kernel, Parameter):
		L = {'0':'P_BenMutation', '1':'P_DelMutation', '2':'P_IncMutRate', '3':'P_DecMutRate'}
		if Kernel == 'Gaussian':
			if self.sch == 'Smooth':
				mutation =  round( math.fabs( np.random.normal( FamilyMemeber.genome[L[str(Parameter)]], 0.1 ) ), 4)
				FamilyMemeber.genome[L[str(Parameter)]] = mutation
			else:
				print np.random.normal(FamilyMemeber.genome[L[str(Parameter)]], 0.3)# to be implemented
			#print(FamilyMemeber.genome[L[str(Parameter)]])
		elif Kernel == 'Uniform':
			pass
		self.ValidateChildGenomeB4(FamilyMemeber)
		
	'''
		This Function will generate a child template.
	'''
	def newChildrenID(self, Parents):
		return self.getID_FtomString( Parents[0].g ) + 1

	'''
		This Function extracts the ID pattern tag from the Generation string 
		example A0:9 ->> extracts the 9 and cast it into int.
	'''
	def getID_FtomString(self, Generation):
		return int( Generation.split(':',1)[1] )

	'''
		This function creates a default children by calling the function,
		insertRunningParameters used in the mating of parents.
	'''
	def newChildren(self, GenerationID):
		offspring = list()
		for i in xrange(self.os):
			child = Family(Parent = False, Generation = 'C'+str(self.currentGen_id)+':'+self.child_ID_NextGen())
			self.insertRunningParameters(child)
			offspring.append(child)
		return offspring

	'''
		This will create 3 lists indicating teh mutation pattern
		per child.

		*** Consider to change this into a dict to increase readeability ***
		Executioni v_1.0:
		1) Generate a list with mutations sampled from binomial distribution 
		(We only consider succes or fail)
		2) Generate the paretan inheritance calling inheritance()
	'''
	def mutationLists(self):
		children_Mtations = [ np.random.binomial( 1, self.psnp, self.os ).astype(bool),   # SNPS
							  np.random.binomial( 1, self.pco, self.os ).astype(bool),    # Cross Over
							  np.random.binomial( 1, self.ptrans, self.os ).astype(bool), # Translocation
							 ]
		children_Mtations.append( self.inheritance() ) 
							      # Adand or Eve inheritance in parameters
							
		return children_Mtations

	'''
		This function creates the parental chromosomal seggreagation as a binomial distributon.
		It will generate a vector of of 0s and 1s to then be casted to boolean.
		1 Will represent Adan amd 0 Eve.
		
		Execution v_1.0:
		1) Loop from the nomber of children
		2) Generate Bionomial Distribution considereing the parameters

		Notes: len(self.idx) -> indcates how many parameters are we modifying in the selected
				scheme.
	'''
	def inheritance(self):
		inh = list()
		for i in xrange(self.os):
			inh.append( np.random.binomial( 1, 0.5, len(self.idx) ).astype(bool) )
		return inh

	'''
		This function will mate Adan and Eve. Is scheme specific
		
		Execution v_1.0:
		1) Determine scheme
		2) Call the proper function 
	'''
	def mate(self, ParentMutation_List, offspring, Parents):
		if self.expl == 'basic4':
			self.mateBasic4(ParentMutation_List, offspring, Parents)
		else:
			pass# implement rest of the functions

	'''
		This function generates the parental segreggation in the Basic4 schem.
		Execution v_1.0:
		1) Ierate two lists
		2) Iterate indivial elements from two lists
		3) If the mutation is Strand is true it will seggreagate from Adan, Eve otherwise
		4) Validate the criteria of P_Ben + P_Del < 1 & P_Inc + P_Dec < 1
	'''
	def mateBasic4(self, ParentMutation_List, offspring, Parents):
		Basic4 = ['P_BenMutation', 'P_DelMutation','P_IncMutRate','P_DecMutRate']
		for child, ParentStrand in zip(offspring, ParentMutation_List): ## Loop in children and Parental list
			print(child.g) # Print child ID for tracking purposes
			for Strand, parameter in zip(ParentStrand, Basic4):  # Now
				if Strand:
					child.genome[parameter] = Parents[0].genome[parameter]
					print('Adan', parameter, child.genome[parameter])
				else:
					child.genome[parameter] = Parents[1].genome[parameter]
					print('Eve', parameter, child.genome[parameter])
			self.ValidateChildGenomeB4(child)

	'''
		This function will validate if the adjusted parameters in the Basic4 scheme is Valid.
		If the values are  not valid, it will decrease them by a factor of 0.05
		until convergance is achieved.
	'''
	def ValidateChildGenomeB4(self, Child):
		Condtion = True
		P_BEN_DEL = False
		P_INC_DEC = False
		while Condtion:
			if Child.genome['P_BenMutation'] + Child.genome['P_DelMutation'] > 1 :
				Child.genome['P_BenMutation'] = round( math.fabs( Child.genome['P_BenMutation'] - 0.05), 4) 
				Child.genome['P_DelMutation'] = round( math.fabs( Child.genome['P_DelMutation'] - 0.05), 4) 
				print('Invalid: Updating P_BenMutation ', Child.genome['P_BenMutation'], ' & P_DelMutation ->', Child.genome['P_DelMutation'])
			else:
				P_BEN_DEL = True
			if Child.genome['P_IncMutRate'] + Child.genome['P_DecMutRate'] > 1 :
				Child.genome['P_IncMutRate'] = round(math.fabs( Child.genome['P_IncMutRate'] - 0.05), 4)  
				Child.genome['P_DecMutRate'] = round(math.fabs( Child.genome['P_DecMutRate'] - 0.05), 4) 
				print('Invalid: Updating P_IncMutRate ', Child.genome['P_IncMutRate'], ' & P_DecMutRate ->', Child.genome['P_DecMutRate'])
			else:
				P_INC_DEC = True

			if P_BEN_DEL and P_INC_DEC:
				Condtion = False

	'''
		This function generates introuces crossover in selected children. 
		1) Select the proper exploration ptattern
		2) Call the function crossOverB4(...)

	'''
	def crossOver(self, mutationList, children, Parents):
		if self.expl == 'basic4':
			self.crossOverB4(mutationList, children, Parents)
		elif self.expl == 'other': # Modigy this
			pass


	'''
		If the children has a crossover mutation, then introduce it.
	'''
	def crossOverB4(self, mutationList, children, Parents):
		for crossover, child in zip(mutationList, children):
			if crossover:
				self.doCrossOverB4(child, Parents)

	'''
		This function will introduce the cross over mutation in child's genome.
		Execution v_1.0:
		1) Generete breakpoint consideeing 4 locations.
			1.a) if location is the last one 3 then copy idx 0,1,2 from Eve and 3 from Adan
			1.b) if breakpoint location is in the middle then loop after the breakpoint loci 
				 copy from Adan, and befir from Eve
			1.c) if location is in the begining then copy idx 1,2,3 from Eve and 0 from Adan.
		2) Validate genome modifications
	'''
	def doCrossOverB4(self, child, Parents):
		L = {'0':'P_BenMutation', '1':'P_DelMutation', '2':'P_IncMutRate', '3':'P_DecMutRate'}
		print('Generating CrossOver in ' + child.g)
		break_point = random.randint( 0, len(self.idx) - 1 )
		print ('Break Point is '+ str(break_point))
		if break_point == 3:
			print ('CASE  break point is in extreme -> 3')
			for i in range(0,3):
				child.genome[L[str(i)]] = Parents[1].genome[L[str(i)]]
			child.genome[L['3']] = Parents[0].genome[L['3']]
		elif break_point == 0:
			print ('CASE  break point is in extreme -> 0')
			for i in range(1,3):
				child.genome[L[str(i)]] = Parents[1].genome[L[str(i)]]
			child.genome[L['0']] = Parents[0].genome[L['0']]
		else:
			print ('CASE The BP is != 3 & != 0', break_point)
			for i in range(0, break_point):
				child.genome[L[str(i)]] = Parents[0].genome[L[str(i)]]
			for i in range(break_point, len(self.idx)):
				child.genome[L[str(i)]] = Parents[1].genome[L[str(i)]]
		print 'Validating mutations'
		self.ValidateChildGenomeB4(child)
		print 'Done'

	'''
		This function will generate the html files for accesing in the Nectar Allocation
		This function will generate the usmmary file and write it on the 
		../Summrary/Summrary.html path.
		This function is called constatnly.

		Note: Consider create a separate object for this function.
	'''
	def toHtml(self, Hall_of_Fame, Status, FamilyMemeber):
		end_line = "</p>\n"
		start_line = "<p>"
		generation = str(self.currentGen_id - 1 )
		end_file = "</body>\n</html>"
		table = "<table border=1><td>ID</td><td>Fitness</td>\n"
		doc = """<!DOCTYPE html>\n<html>\n<body>\n<h1>Peter MacCalum Bioinformatics Division: GA Algortihm v 1.0</h1>\n<p>Status:"""
		start_row = "<tr><td>"
		mid_element = "<td>"
		end_row = "</td>\n"
		end_table = "</table>"
		if Status:
			doc = doc + 'Running' + end_line
		else:
			doc = doc + 'Done' + end_line

		host = str(socket.gethostname())
		doc = doc + start_line + 'Host: ' + host + end_line
		doc = doc + start_line + 'Exploration: ' + self.expl + end_line
		doc = doc + start_line + 'Scheme: ' + self.sch + end_line
		doc = doc + start_line + 'Generation: ' + generation + end_line
		doc = doc + '<h2> Generation Table of Results </p>\n'
		doc = doc + table 

		
		for element in sorted(Hall_of_Fame, key=Hall_of_Fame.get, reverse=True):
			doc = doc + start_row + element + mid_element +  str(Hall_of_Fame[element]) + end_row # averge 

		doc = doc + end_table + end_file
		
		Html_file= open(self.path_Summary + "/Summary.html","w")
		Html_file.write(doc)
		Html_file.close()
		#self.createHall_of_Fame_html(Hall_of_Fame, Status, FamilyMemeber)

	'''
		This function is a setter.
		It introduces the Fitness into the <dict>.self.HOF[<str>]
	'''
	def insert_Hall_of_Fame(self, Fitness, FamilyMemeber):
		self.HOF[FamilyMemeber.g] = Fitness
		self.FamMem.append(FamilyMemeber)
		#self.createHall_of_Fame_html(self.HOF, True, FamilyMemeber)

	'''
		This function calls the hall of fame function to create
		the html file called Hall_of_Fame.html 

		The Boolean indicator indicates if the code is running.

		Note: The off command is not coded yet. It could be implemented
		by keeping track of the number of iterations and in the las one set
		the value to False.
	'''
	def Hall_of_Fame_toHtml(self, FamilyMemeber):
		self.createHall_of_Fame_html(self.HOF, True, FamilyMemeber)


	'''
		This function is a getter for the Hall of Fame, getHall_of_Fame()
	'''
	def getHall_of_Fame(self):
		return self.HOF

	'''
		This function creates the hall of fame Html file.
		This function saves the file in .../Hall_of_fame/Hall_of_Fame.html
		Note: Consider create a separate object for this function.
	'''
	def createHall_of_Fame_html(self, Hall_of_Fame, Status, FamilyMemeber):
		
		end_line = "</p>\n"
		start_line = "<p>"
		generation = str(self.currentGen_id - 1)
		end_file = "</body>\n</html>"
		table = "<table border=1><td>ID</td><td>P_Ben</td><td>P_Del</td><td>P_Inc</td><td>P_Dec</td><td>Avg_Fitness</td><td>Std</td><td>min</td><td>25%</td><td>50%</td><td>75%</td><td>max</td>\n"
		doc = """<!DOCTYPE html>\n<html>\n<body>\n<h1>Peter MacCalum Bioinformatics Division: GA Algortihm v 1.0</h1>\n<p>Status:"""
		start_row = "<tr><td>"
		col = "</td><td>"
		mid_element = "<td>"
		end_row = "</td>\n"
		end_table = "</table>"
		if Status:
			doc = doc + 'Running' + end_line
		else:
			doc = doc + 'Done' + end_line

		host = str(socket.gethostname())
		doc = doc + start_line + 'Host: ' + host + end_line
		doc = doc + start_line + 'Exploration: ' + self.expl + end_line
		doc = doc + start_line + 'Scheme: ' + self.sch + end_line
		doc = doc + start_line + 'Generation: ' + generation + end_line
		doc = doc + start_line + 'Repetitions: ' + str(self.perm) + end_line
		doc = doc + '<h2> Hall of Fame </p>\n'
		doc = doc + table 
		for element in Hall_of_Fame:
			if not element in self.topScores:
				print 'Inserting '+ FamilyMemeber.g + ' in Hall of Fame'
				self.topScores[element] = FamilyMemeber

		for element in sorted(Hall_of_Fame, key=Hall_of_Fame.get, reverse=True):
			#print (element)
			df = self.PopulationList_to_DataFrame(element)
			doc = doc + start_row + element + col \
				+  str(self.topScores[element].genome['P_BenMutation']) + col \
				+  str(self.topScores[element].genome['P_DelMutation']) + col \
				+  str(self.topScores[element].genome['P_IncMutRate']) + col \
				+  str(self.topScores[element].genome['P_DecMutRate']) + col \
				+  str(Hall_of_Fame[element]) + col \
				+  str(df['Fitness']['std']) + col \
				+  str(df['Fitness']['min']) + col \
				+  str(df['Fitness']['25%']) + col \
				+  str(df['Fitness']['50%']) + col \
				+  str(df['Fitness']['75%']) + col \
				+  str(df['Fitness']['max']) + col \
				+ end_row

		doc = doc + end_table + end_file
		Html_file = open(self.path_hof + "/Hall_of_Fame.html","w")
		Html_file.write(doc)
		Html_file.close()

	'''
		This function transform the Population list to a pandas Data Frame,
		consider implementing this function to only insert new values in 
		the Data  Frame.

		This function returns the basic statistical parameters as a <dict>

	'''
	def PopulationList_to_DataFrame(self, FamilyMemeber):
		df = DataFrame(self.FitnessPopulation_List[FamilyMemeber], columns = ['Fitness'])
		#print(df.describe())
		return df.describe().to_dict()

	'''
		Not implemented!!!
		This function can be deleted
	'''
	def RandomizeAll(self, FamilyMemeber):
		print random.uniform(0,1)	


	'''
		This function sorts the Population Dictionary and sends
		to perform a selection on the desired members.
		Execution v_1.0:
		1) Sort Population dictionary
		2) Call function Select Offspring
	'''
	def UpdateGeneration(self, Parents, Children, Population):
		sorted_Population = sorted(Population, key=Population.get,  reverse=True)
		NewGeneration = self.Select_Offspring(Parents, Children, sorted_Population)
		self.inc_CurrentGen()
		return(NewGeneration)
		

	'''
		This function determines the memeber selection sheme to update the generation
		Execution v_1.0:
		1) Determine Mating Model 
		2) Call the corresponfig function
	'''
	def Select_Offspring(self, Parents, Children, Sorted_Population):
		if self.MatMod == 'HiLow':
			print("Total Population", len(Parents) + len(Children))
			return(self.HiLow_Selection(Parents, Children, Sorted_Population))
		elif self.MatMod == 'Hi':
			pass
		else:
			pass

	'''
		This is the core function of the High and Low Mating scheme
		Execution v_1.0:
		1) Determine the Highest and lowest fitted memeber by taking the
		SelectionPercentage of the Sorted Population. 
		2) Once extracted the ID's of the Best and Worst Fitted memebers, 
		call getNewParents to obtain the genomes and mate the parents of
		the new generation.
	'''
	def HiLow_Selection(self, Parents, Children, Sorted_Population):
		BestFit = []
		WorstFit = []
		Total_Size = len(Parents) + len(Children)
		number_of_selected = int(round( ( Total_Size * self.selper )/ 100.0 ))
		print("Parents: ", len(Parents),"Children: ",len(Children), "Sorted: ", len(Sorted_Population))
		print("Number of SELECTED ", number_of_selected)
		selectedMemebers = dict()
		print("Sorted:",Sorted_Population)

		if number_of_selected == 1:
			BestFit.append(Sorted_Population[0])
			WorstFit.append(Sorted_Population[Total_Size - 1])
	
		else:
			for i in range(int(number_of_selected)):
				BestFit.append(Sorted_Population[i])
				WorstFit.append(Sorted_Population[Total_Size - i - 1])
		print(Parents, Children, BestFit, WorstFit)
		return(self.getNewParents(Parents, Children, BestFit, WorstFit))

	'''
		This function split the data in Best and Worst members,
		and call the function MateHiLow().
		Execution v_1.0:
		1) If the generation is greater than 0, all memebers are children
		 --> First merge Parents and Children data.
		 Note: From generations greater than 0 all memebers have ID as Children,
		 but the execution will focus on the Worst and lowest fited memebers.
		 This will keep the seraching in a ciruclar basis.

		 If the genearion is not greater than 0 then look for children and Parents
		 2) Once the best and worst survivors are found call MateHiLow(...) function. 
	'''
	def getNewParents(self, Parents, Children, BestFit, WorstFit):
		BestSurvivors = []
		WorstSurvivors = []
		
		if self.currentGen_id > 0:
			print("Case gen greater that 1")
			MergedList = Parents + Children
			for best in BestFit:
				BestSurvivors.append( self.getMember_by_ID(best, MergedList) )
			for worst in WorstFit:
				WorstSurvivors.append( self.getMember_by_ID(worst, MergedList) )
		else:
			print("Case gen 0")
			for best in BestFit:
				if self.getMemeberTag(best) == 'C':
					print("Best C found")
					BestSurvivors.append( self.getMember_by_ID(best, Children) )
				else:
					BestSurvivors.append( self.getMember_by_ID(best, Parents) )

		
			for worst in WorstFit:
				print("Worst C found")
				if self.getMemeberTag(worst) == 'C':
					WorstSurvivors.append( self.getMember_by_ID(worst, Children) )
				else:
					WorstSurvivors.append( self.getMember_by_ID(worst, Parents) )

		return (self.MateHiLow(BestSurvivors, WorstSurvivors))

	'''
		This function gets the selected memeber surviving in the  current 
		generation.
		Execution v_1.0:
		1) From all Memebers list look the ID, if found return that value.
		ID and memeber.g are strings
	'''
	def getMember_by_ID(self, ID, Members):
		for member in Members:
			if ID == member.g:
				return member

	'''
		This function is the first step in mating High and Low memebers of
		the current generation. It will select the proper fuction according 
		to the search scheme. This function calls CorrelateBasic4(...)

		Note: This function will require further development if more schemes
		are implemented.
	'''
	def MateHiLow(self, BestSurvivors, WorstSurvivors):
		children = []
		if self.expl == 'basic4':
			children.append(self.CorrelateBasic4(BestSurvivors))
			children.append(self.CorrelateBasic4(WorstSurvivors))
		else:
			pass # implement different schmes
		return(children)
	
	'''
		This function will generate the possible combinations in the list
		of survivors. The output will be a tuple list wchich will be mated.
		The output of this function will generate "self.nextOff" children per
		tuple combinaton.
	'''
	def CorrelateBasic4(self, Survivors):
		couples = []
		children = []
		for combination in combinations(Survivors, 2):
   			couples.append(combination)

   		print("Couples: ", couples)

   		self.os = self.nextOff

   		for couple in couples:
   			children.append(self.CoupleDifference(couple))
   		return(children)

	'''
		This function computes the differentes in the basic 4 scheme parameters,
		then it complements. This approach inteds to keep the GA in a circular fashion.
		Aiming to interpoalte all points in the fitness landscape. The complement of the
		absolute difference of the tuple will be the initial parameter values for the 
		mutaton scheme in generaating the ofsspring.

		This function calls createOffspring_NextGen(...)
	'''	
	def CoupleDifference(self, Couple):
		children = []
		genome = dict()
		#self.CrossOver_Survivors(Couple)
		P_Ben = self.ValidProbablities( 1 - abs( Couple[0].genome['P_BenMutation'] - Couple[1].genome['P_BenMutation'] ) )
		P_Del = self.ValidProbablities( 1 - abs( Couple[0].genome['P_DelMutation'] - Couple[1].genome['P_DelMutation'] ) )
		P_Inc = self.ValidProbablities( 1 - abs( Couple[0].genome['P_IncMutRate'] - Couple[1].genome['P_IncMutRate'] ) )
		P_Dec = self.ValidProbablities( 1 - abs( Couple[0].genome['P_DecMutRate'] - Couple[1].genome['P_DecMutRate'] ) )
		genome['P_BenMutation'] = P_Ben
		genome['P_DelMutation'] = P_Del
		genome['P_IncMutRate'] = P_Inc
		genome['P_DecMutRate'] = P_Dec
		children = self.createOffspring_NextGen(Couple, genome)
		
		return children

	'''
		This function creates new offspring for a couple memeber.
		It follows the same structre as when the parents are mated.
		Execution v_1.0:
		1) Assign a new ID to the children.
		2) Perform crossOver if required
		3) Perform translocation if required
		4) Introduce SNPS if required
	'''
	def createOffspring_NextGen(self, Couple, genome):
		print('\n|-------[ . . . M A T I N G   G E N ' +str(self.currentGen_id-1) + ' . . . ]-------| ')
		print('>> Generating ' + str(self.os) + ' children ')
		print('>> ' + Couple[0].g + ' with ' + Couple[1].g)
		c_id = self.newChildrenID_NextGen(Couple)
		offspring = self.newChildren_NextGen(c_id) # This creates 5 children
		offspring = self.doCrossOver_NextGen(offspring, Couple)
		self.Translocation(np.random.binomial( 1, 0.25, self.os ).astype(bool), offspring)
		self.PointMutations_NextGen(offspring, genome)
		return(offspring)
		
	'''
		This function introduces point mutations in the ofsspring.
		The point mutation induction is sampled from a normal distribution with
		mean assitiated to the value previously defined and variance of +/- 0.1. The
		value is truncated to have only 4 decimal places.
		Execution v_1.0:
		1) For the selected children in offspring, introduce a SNP with mean
		defined in the parameter value and variace of 0.1.
		2) Once mutations are generated, introduce them in child's genome.
	'''
	def PointMutations_NextGen(self, offspring, genome):
		children = []
		L = {'0':'P_BenMutation', '1':'P_DelMutation', '2':'P_IncMutRate', '3':'P_DecMutRate'}
		for child in offspring:
			for i in xrange(0,4):
				mutation = round( math.fabs( np.random.normal(genome[L[str(i)]], 0.1) ) ,4)
				child.genome[L[str(i)]] = mutation
			children.append( self.ValidateChildGenomeB4(child) )
		
		return(children)

	'''
		This function introduces a crossover mutation in child's genome 
		by calling the function Mating_NextGen(...)
	'''
	def doCrossOver_NextGen(self, offspring, Couple):
		new_children = []
		for child in offspring:
			new_children.append (self.Mating_NextGen(child, Couple) )

		return(new_children)
			

	'''
		This function will create the new ofsspring of children using the
		function insertRunningParameters(...). The latter function is implemented
		to have a default init fpor childrne's genome.

	'''
	def newChildren_NextGen(self, GenerationID):
		offspring = list()
		print("Total to create", self.os)
		for i in xrange(self.os):
			child = Family(Parent = False, Generation = (GenerationID + self.child_ID_NextGen()))
			self.insertRunningParameters(child)
			offspring.append(child)


		return offspring

	'''
		This function return the part of the children generation ID,
		The ID is in the following way,
		[C-P][0-9]:[0-000]
	'''
	def newChildrenID_NextGen(self, Couple):
		#member_A = int(Couple[0].g[1])
		#member_B = int(Couple[1].g[1])
		#if member_A >= member_B:
		#	return 'C' + str(member_A + 1) + ':'
		#else:
		#	return 'C' + str(member_B + 1) + ':'
		return ('C'+str(self.currentGen_id)+':')

	'''
		Update the id of via local variable and update, this wil control
		the memeber increment vi local value of of the class.
	'''
	def child_ID_NextGen(self):
		child_id = self.nextGen_id
		self.nextGen_id = self.nextGen_id + 1
		return( str(child_id) )

	'''
		This function increment the current generation plus one for new
		elements in the current generation.
	'''
	def inc_CurrentGen(self):
		self.currentGen_id = self.currentGen_id + 1

	'''
		This function generates the mating scheme, rather if it is
		parent 0 or one. It throws a coin-flip to chose which parameter will
		be inherited from wchi parent.
		If value less than 0.5 it takes the value fome Parent[0], Parent[1] otherwise.
		Finally, the genome is validated by the function ValidateChildGenomeB4(...) 
	'''
	def Mating_NextGen(self, child, Couple):
		L = {'0':'P_BenMutation', '1':'P_DelMutation', '2':'P_IncMutRate', '3':'P_DecMutRate'}
		#print("Original P0 Genome: ", Couple[0].genome[L['0']], Couple[0].genome[L['1']], Couple[0].genome[L['2']], Couple[0].genome[L['3']])
		#print("Original P1 Genome: ", Couple[1].genome[L['0']], Couple[1].genome[L['1']], Couple[1].genome[L['2']], Couple[1].genome[L['3']])
		#print("Original Child's Genome: ", child.genome[L['0']],child.genome[L['1']], child.genome[L['2']], child.genome[L['3']])
		for i in xrange(0,4):
			mating = random.uniform(0,1)
			print(mating)
			if mating < 0.5:
				#print("Parent 0 in ", i)
				child.genome[L[str(i)]] = Couple[0].genome[L[str(i)]]
			else:
				#print("Parent 1 in ", i)
				child.genome[L[str(i)]] = Couple[1].genome[L[str(i)]]
		#print 'Validating mutations'
		self.ValidateChildGenomeB4(child)
		#print("CrossOver Child's Genome: ", child.genome[L['0']],child.genome[L['1']], child.genome[L['2']], child.genome[L['3']])
		return child

	'''
		This function checks for extreme valu cases, if values are extreme
		correct them using a samples from uniform distribution
	'''
	def ValidProbablities(self, Value):
		print(Value)
		if Value == 0.0 or Value == 1.0 :
			Value = random.uniform(0,1)
			print( "Value is 0 or 1 changing to: ",Value )
		else:
			print("Good values")
		return(Value)


	'''
		This function gets the id of the memeber and takes the first character,
		this character is the identifier of the memeber. It can be a Child or 
		a Parent.
		The ID encodign is: [C-P][0-00]:[0-00] for instance C1:0 or E0:0
	'''
	def getMemeberTag(self, ID):
		return( ID[0] )
		
	'''
		This function computes the average fitness, to asses reproducibility it is 
		necesary to simulate 'N' number of times.
	'''
	def AverageFitness(self, Fitness):
		Fitness = self.remove_values_from_list(Fitness, 0)
		validValues = 0
		Total = 0.0
		for value in Fitness:
			if value != 0:
				validValues += 1
				Total += value
		if validValues != 0:
			meanValue = Total / validValues
		else:
			meanValue = 0 # Modify this
		return meanValue

	'''
		This function removes the selected value  form the list.
	'''
	def remove_values_from_list(self, numlist, value):
		clean_list = filter(lambda data: data != value, numlist)
		return clean_list

	'''
		This funtion updates the population list by removing zeros
	'''
	def insert_to_FitnessPopulation_List(self, Fitness, member):
		self.FitnessPopulation_List[member.g] = self.remove_values_from_list(Fitness,0)
 
	'''
		This function is a getter, returns the Fitness Population list
	'''
	def get_FitnessPopulation_list(self):
		return self.FitnessPopulation_List
 
	'''
		This function will return the values forem the Fitness Population,
		usign the memeber as value.
	'''
	def get_FitnessPopulation_list_by_memeber(self, memeber):
		return self.FitnessPopulation_List[memeber.g]

	'''
		This function insert the fitness value of FamilyMemeber class dict()
	'''
	def insertPopulation(self, Fitness, FamilyMemeber):
		self.Population[FamilyMemeber.g] = Fitness
		
	'''
		This function is a getter
	'''
	def getPopulation(self):
		return (self.Population)

	'''
		This function appends in the best list the new value
	'''
	def insertBest(self, Best):
		self.Best = self.Best + Best
		
	'''
		This function is a getter for the best memebers list.
	'''
	def getBest(self):
		return self.Best

	'''
		This function appends to trhe worst list.
	'''
	def insertWorst(self, Worst):
		self.Worst = self.Worst + Worst

	'''
		This function is a getter for the Worst list
	'''
	def getWorst(self):
		return self.Worst

	'''
		This function is scheme specific for the Radar plot,
		it calls the specific function according to a shceme.
		This function calls GenerateRadar_Basic4(...)
	'''
	def insert_to_RadarPlot(self, FamilyMemeber, AverageFitness):
		if self.expl == 'basic4':
			self.GenerateRadar_Basic4(FamilyMemeber, AverageFitness)
		else:
			pass

	'''
		This function will generate a the best fit memebers for the radar plot.
		This function will compare the fitness of the memeber and if it is greater 
		than the mutation rate, then it will insert the value in <dict>.BestFitRadar(...).
		It also insert the value in <self>.Radar_Data(...)

	'''
	def GenerateRadar_Basic4(self, FamilyMemeber, AverageFitness):
		if AverageFitness > self.mr:
			self.BestFitRadar[FamilyMemeber.g] = ( FamilyMemeber.genome['P_BenMutation'],\
												 FamilyMemeber.genome['P_DelMutation'],\
												 FamilyMemeber.genome['P_IncMutRate'], \
												 FamilyMemeber.genome['P_DecMutRate'], \
												 AverageFitness
												 )

		self.Radar_Data[FamilyMemeber.g] = ( FamilyMemeber.genome['P_BenMutation'],\
												 FamilyMemeber.genome['P_DelMutation'],\
												 FamilyMemeber.genome['P_IncMutRate'], \
												 FamilyMemeber.genome['P_DecMutRate'], \
												 AverageFitness
												 )

	'''
		This function generates the multivariate linear regression ussing the Basic 4 parameters.
		Tis function will compute as well the regression of the log transformation.
	'''
	def multivariate_RegressionB4(self):
		Raw_DF = []
		log_DF = []
		if len(self.Radar_Data) > 0:
			for member in self.Radar_Data:
				if self.Radar_Data[member][4] != 0:
					Raw_DF.append([ self.Radar_Data[member][0],\
								    self.Radar_Data[member][1], \
								    self.Radar_Data[member][2], \
								    self.Radar_Data[member][3], \
								    self.Radar_Data[member][4] ])
					log_DF.append([ math.log(self.Radar_Data[member][0]),\
								    math.log(self.Radar_Data[member][1]), \
								    math.log(self.Radar_Data[member][2]), \
								    math.log(self.Radar_Data[member][3]), \
								    math.log(self.Radar_Data[member][4]) ])
			df     = DataFrame(Raw_DF, columns = ['P_Ben','P_Del','P_Inc','P_Dec','AvgFitness'])
			log_df = DataFrame(log_DF, columns = ['P_Ben','P_Del','P_Inc','P_Dec','AvgFitness'])
			result = sm.ols(formula = "AvgFitness ~ P_Ben + P_Del + P_Inc + P_Dec", data = df).fit()
			logRes = sm.ols(formula = "AvgFitness ~ P_Ben + P_Del + P_Inc + P_Dec", data = log_df).fit()
		
			self.regression_to_html(result.summary(), logRes.summary())

	'''
		This function converts the ols ouytup to an html file
	'''
	def regression_to_html(self, results, logResults):
		end_line   = "</p>\n"
		start_line = "<p>"
		generation = str(self.currentGen_id-1)
		end_file   = "</body>\n</html>"
		doc        = """<!DOCTYPE html>\n<html>\n<body>\n<h1>Peter MacCalum Bioinformatics Division: GA Algortihm v 1.0</h1>\n<p>Multivariate Regression for Generation:"""+generation+"</p>"
		host       = str(socket.gethostname())
		doc        = doc + start_line + 'Host: ' + host + end_line
		doc        = doc + start_line +'Exploration: ' + self.expl + end_line
		doc        = doc + start_line + 'Scheme: ' + self.sch + end_line		
		doc        = doc + '<p>Model of Form: Fitness = P_Beneficial + P_Deleterious + P_IncMutRate + P_DecMutRate</p>'
		doc        = doc + '<p> Multivariate Regression </p>\n'
		doc        = doc + '<pre>\n' + str(results) + '\n</pre>\n'
		doc         = doc + '<p> Mulivariate log Regression </p>\n'
		doc = doc + '<pre>\n' + str(logResults) + '\n</pre>\n'+end_file
		print doc
		Html_file = open(self.path_Regression + '/MVRegression_gen_' + generation +'.html',"w")
		Html_file.write(doc)
		Html_file.close()

	'''
		This function generates the radar plots for all generation memebers, 
		it generates the pygal version an the png version.
	'''
	def plot_BestFitRadar(self, Remove):
		bestFit = []
		IDs = []
		radar_chart = pygal.Radar()
		radar_chart.title = 'Best Fit Members in Generation '+ str(self.currentGen_id-1) + ' when end Avg. Mutation Rate is Greater than ' + str(self.mr)
		radar_chart.x_labels = ['P(Beneneficial)', 'P(Deleterious)', 'P(Increment_MR)', 'P(Decrement_MR)', 'Avg. Fitness']
		if len(self.BestFitRadar) > 0:
			for member in self.BestFitRadar:
				IDs.append(member)
				bestFit.append([ self.BestFitRadar[member][0],\
						  self.BestFitRadar[member][1],\
						  self.BestFitRadar[member][2],\
						  self.BestFitRadar[member][3],\
						  self.BestFitRadar[member][4]])
				radar_chart.add(member, [ self.BestFitRadar[member][0],\
										  self.BestFitRadar[member][1],\
										  self.BestFitRadar[member][2],\
										  self.BestFitRadar[member][3],\
										  self.BestFitRadar[member][4]])
			radar_chart.render_to_file(self.path_Plots +'/Radar_Gen_'+str(self.currentGen_id-1) +'.svg') 
		else:
			print("Values do not meet criteria")
		# Delete all to only siplay best in current generation
		if Remove:
			self.BestFitRadar = dict()

		self.save_matplot_img(IDs, bestFit)

	'''  
		This function saves current generation data into a png file, it calls the polygon class.
	''' 
	def save_matplot_img(self, IDs, bestFit):
		save_path = self.path_Plots + '/Radar_Best_Memebers_Gen_'+str(self.currentGen_id-1) +'.png'
		polygon.generate_Single_BestFitB4( bestFit, str(self.currentGen_id-1), save_path )

	'''
		This function generates all radar data using the pygal class and calls the
		plot_Best_and_Worst(...) function,
	'''
	def plot_AllRadar_Data(self):
		radar_chart = pygal.Radar()
		radar_chart.title = 'All Members in Generation '+ str(self.currentGen_id-1) + ' when end Avg. Mutation Rate is Greater than ' + str(self.mr)
		radar_chart.x_labels = ['P(Beneneficial)', 'P(Deleterious)', 'P(Increment_MR)', 'P(Decrement_MR)', 'Avg. Fitness']
		if len(self.Radar_Data) > 0:
			for member in self.Radar_Data:
				radar_chart.add(member, [ self.Radar_Data[member][0],\
										  self.Radar_Data[member][1],\
										  self.Radar_Data[member][2],\
										  self.Radar_Data[member][3],\
										  self.Radar_Data[member][4]])
			radar_chart.render_to_file(self.path_Plots + '/All_Radar_Memebers_Gen_'+str(self.currentGen_id-1) +'.svg') 
		else:
			print("Values do not meet criteria")
		self.plot_Best_and_Worst()

	'''
		This function handles the radar plot for the best and worst elements using two different colors.
		This functions calls the class polygon. The class polygon can be easily modified to make further
		changes in radar plots
	'''
	def plot_Best_and_Worst(self):
		save_path = self.path_Plots + '/Radar_Best_and_Worst_Memebers_Gen_'+str(self.currentGen_id-1) +'.png'
		BestFit = []
		WorstFit = []
		elements = []
		fitness_dict = dict()
		IDs = []
		if len(self.Radar_Data) > 0:
			for element in self.Radar_Data:
				fitness_dict[element] = self.Radar_Data[element][4]
			for element in sorted(fitness_dict, key= fitness_dict.get, reverse = True  ):
				IDs.append(element)
				elements.append( [ self.Radar_Data[element][0],\
								   self.Radar_Data[element][1],\
								   self.Radar_Data[element][2],\
								   self.Radar_Data[element][3],
								   self.Radar_Data[element][4]])

		total = math.ceil( (len(elements) / 2.0 )/2.0 )
		print "elements",elements
		exit_flag = 0
		for element_up, element_down in zip(elements, reversed(elements) ):
			if exit_flag != total:
				exit_flag += 1
				BestFit.append(element_up)
				WorstFit.append(element_down)
			else:
				break
		polygon.generate_Best_and_Worst(BestFit, WorstFit, str(self.currentGen_id-1), save_path)

	'''
		This function saves the generation data as a pandas Drata Frame, 
		separated by ','
	'''
	def save_DataFrame(self):
		Raw_DF = []
		if len(self.Radar_Data) > 0:
			for member in self.Radar_Data:
				if self.Radar_Data[member][4] != 0:
					Raw_DF.append([ self.Radar_Data[member][0],\
								    self.Radar_Data[member][1], \
								    self.Radar_Data[member][2], \
								    self.Radar_Data[member][3], \
								    self.Radar_Data[member][4] ])
			df = DataFrame(Raw_DF, columns = ['P_Ben','P_Del','P_Inc','P_Dec','AvgFitness'])
			df.to_csv(self.path_df + '/DataFrame_gen_'+ str(self.currentGen_id-1)+'.csv')

	'''
		This function creates the different data folders used in the simulation. 
		This function sets the path locations for differnt folders to be used.
		The main path should be modified for the Necatar Allocations
	'''
	def createFolders(self):
		ts = time.time()
		st = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d')
		#modify this to nectar folder /var/www/DarwinGA.com/
		current_Path = os.getcwd()
		print "Creating Folders in: " + current_Path 
		current_Path += "/Simulation_"+ st
		print "Results will be stored in: " + current_Path
		path_hof = current_Path + "/Hall_of_Fame"
		path_Summary = current_Path + "/Summary"
		path_Plots = current_Path + "/Radar"
		path_df = current_Path + "/Data_Frame"
		path_MVR = current_Path + "/MVRegression"
		self.createFolder(path_hof)
		self.createFolder(path_Summary)
		self.createFolder(path_Plots)
		self.createFolder(path_df)
		self.createFolder(path_MVR)
		self.path_hof = path_hof
		self.path_Summary = path_Summary
		self.path_Plots = path_Plots
		self.path_df = path_df
		self.path_Regression = path_MVR

	'''
		This function creates a foder in the input path location
	'''
	def createFolder(self, Path):
		if not os.path.exists(Path):
			os.makedirs(Path)


'''
	Main Testing Section, functional but consider 
	creating this a separate main object.
'''			
# Testing Section
if __name__ == '__main__':
	
	os.system('clear')
	i = GA()
	i.createFolders()
	Population = dict()
	#i.printValues()
	Local_Hall_of_Fame = dict()
	Best = i.createParents()
	print('\nParents created: ' + str(len(Best)))
	Worst = i.createOffspring(Best)
	print('Children created: ' + str(len(Worst)))

	i.insertBest(Best)
	i.insertWorst(Worst)

	for generation in xrange(0, 3):
		print("Running Generation " + str(generation))
		
		for member in Best:
			print('\n-------------------------------------------------\n')
			Fitness = main.ComputeFitness(member)
			i.insert_to_FitnessPopulation_List(Fitness, member)
			avgFitness = i.AverageFitness(Fitness)
			print ('Avg. Fitness ', avgFitness)
			i.insert_Hall_of_Fame(avgFitness, member)
			i.insertPopulation(avgFitness, member)
			i.toHtml(i.getHall_of_Fame(), True, member)
			i.insert_Hall_of_Fame(avgFitness, member)
			i.Hall_of_Fame_toHtml(member)
			i.insert_to_RadarPlot(member, avgFitness)
		
		for child in Worst:
			print('\n-------------------------------------------------\n')
			Fitness = main.ComputeFitness(child)
			i.insert_to_FitnessPopulation_List(Fitness, child)
			avgFitness = i.AverageFitness(Fitness)
			i.insert_Hall_of_Fame(avgFitness, child)
			i.insertPopulation(avgFitness, child)
			print('Avg. Fitness', avgFitness)
			i.toHtml(i.getHall_of_Fame(), True, child)
			i.Hall_of_Fame_toHtml(child)
			i.insert_to_RadarPlot(child, avgFitness)


		i.plot_BestFitRadar(False)
		i.plot_AllRadar_Data()
		i.save_DataFrame()
		i.multivariate_RegressionB4()

			

		print("Creating New Generation")
		NewGeneration = i.UpdateGeneration(i.getBest(), i.getWorst(), i.getPopulation())
		Best = list(itertools.chain(*NewGeneration[0]))
		Worst = list(itertools.chain(*NewGeneration[1]))
		i.insertBest(Best)
		i.insertWorst(Worst)
		print("Reseting")