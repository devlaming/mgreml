import pandas as pd

class MgremlReader:
    
    def __init__(self, logger, args):
        '''
        MgremlReader is provided with a logger and the input arguments
        from the call of the main mgreml script e.g.
        -> python ./mgreml.py --grm mygrm --pheno mypheno --out results
        
        Currently, this initialisation method only reports the attributes
        an object in this class is required to have once it's initialised
        
        An instance of an MgremlReader will be used as input
        when initialising instances of other classes involved in
        MGREML estimation; the attributes of the MgremlReader will be 
        used to determine what needs to be done
        '''
        self.dfA            # DataFrame based on --grm mygrm
        self.sPrefix        # string based on --out results
        self.dfY            # DataFrame based on --pheno mypheno.txt [nolabelpheno]
        self.dfX            # DataFrame based on --covar mycovar.txt; DEFAULT = None
        self.dfBinXY        # DataFrame based on --covar-model mycovmodel.txt [nolabelpheno] [nolabelcovar]; DEFAULT = None
        self.dfGenBinFY     # DataFrame based on --genetic-model mygenmodel.txt [nolabelpheno] [nolabelfactor]; DEFAULT = None
        self.dfEnvBinFY     # DataFrame based on--environment-model myenvmodel.txt [nolabelpheno] [nolabelfactor]; DEFAULT = None
        self.bPerfectRhoG   # Boolean based on --rho-genetic 1; DEFAULT = False
        self.bNoRhoG        # Boolean based on --rho-genetic 0; DEFAULT = False
        self.bNoRhoE        # Boolean based on --rho-environment 0; DEFAULT = False
        self.dfGenBinFY0    # DataFrame based on --restricted-genetic-model mygenmodel0.txt [nolabelpheno] [nolabelfactor]; DEFAULT = None
        self.dfEndBinFY0    # DataFrame based on --restricted-environment-model myenvmodel0.txt [nolabelpheno] [nolabelfactor]; DEFAULT = None
        self.bPerfectRhoG0  # Boolean based on --restricted-rho-genetic 1; DEFAULT = False
        self.bNoRhoG0       # Boolean based on --restricted-rho-genetic 0; DEFAULT = False
        self.bNoRhoE0       # Boolean based on --restricted-rho-environment 0; DEFAULT = False
        self.bSE            # Boolean based on --no-se; DEFAULT = True (so reporting SEs is the default!)
        self.bRelCutoff     # Boolean based on --rel-cutoff 0.025; DEFAULT = False
        self.dRelCutoff     # number based on --rel-cutoff 0.025; DEFAULT = None
        self.iDropLeadPCs   # integer based on --ignore-pcs 40 [1000]; DEFAULT = 20
        self.iDropTrailPCs  # integer based on --ignore-pcs 40 [1000]; DEFAULT = 0
        self.bStoreIter     # Boolean based on --store-iter 50; DEFAULT = False
        self.iStoreIterFreq # integer based on --store-iter 50; DEFAULT = None
        self.sInitValsFile  # string based on --reinitialise results.iter.250.pkl
        self.bAllCoeffs     # Boolean based on --all-coefficients
        self.bBFGS          # Boolean based on --newton; DEFAULT = True
        self.dGradTol       # number based on --grad-tol; DEFAULT = 1E-5
        # note: when an instance of MgremlReader has all the right
        # attributes, after the __init__() method has been invoked,
        # there is no need to return anything! that instance of MgremlReader
        # can then be used by other classes to retrieve the necessary
        # settings/specifications/data
        
        # now, within this class e.g. self.dfA can be set in the __init__()
        # method using some method like ReadGRM. I.e. we can then say:
        sGRMfilename = '' # based on --grm mygrm
        self.ReadGRM(sGRMfilename)
        # here, in the __init__ method, we thus call on the instance method
        # ReadGRM(), where that method in turn is defined within that class
        # e.g. along the following lines (notice it's higher up in indentation):
        
    def ReadGRM(self, sGRMfilename):
        # first step to read in the grm
        # ......
        self.dfA = None # based on result from preceding steps
