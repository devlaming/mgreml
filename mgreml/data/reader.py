import pandas as pd
import numpy as np
import logging

class MgremlReader:
    
    def __init__(self, parser, logger, dStart):
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
        # store logger and parser as attributes of instance
        self.logger = logger
        self.parser = parser
        # initialise arguments
        self.InitialiseArgsAndLogger()
        #initialize out directory
        if self.args.out is not None:
            self.sPrefix=self.args.out
        else:
            raise ValueError('The --out flag is required.')
        # get DataFrame for GRM based on --grm mygrm option
        if self.args.grm is not None:
            self.read_grm()
        # get DataFrame for Pheno based on --pheno mypheno.txt [nolabelpheno]    
        if self.args.pheno is not None:
            self.read_pheno()
        # get Dataframe for Covar based on --covar mycovar.txt; DEFAULT = None     
        if self.args.covar is not None:
            self.read_covar()
        else:
            self.dfX=None
        # get DataFrame for covar model based on --covar-model mycovmodel.txt [nolabelpheno] [nolabelcovar]; DEFAULT = None
        if self.args.covar_model is not None:
            self.read_covar_model
        else:
            self.dfBinXY=None
        # get DataFrame for genetic model based on --genetic-model mygenmodel.txt [nolabelpheno] [nolabelfactor]; DEFAULT = None
        if self.args.genetic_model is not None:
            self.read_model(component='genetic')
        else:
            self.dfGenBinFY=None
        # get DataFrame for environment model based on --environment-model myenvmodel.txt [nolabelpheno] [nolabelfactor]; DEFAULT = None
        if self.args.environment_model is not None:
            self.read_model(component='environment')
        else:
            self.dfEnvBinFY=None
        
        # check if boolean perfect rho g is given based on --rho-genetic 1; DEFAULT = False
        # check if boolean no rho g is given based on  --rho-genetic 0; DEFAULT = False
        self.bPerfectRhoG=False 
        self.bNoRhoG=False 
        if self.args.rho_genetic is not None:
            if self.args.rho_genetic==1:
                self.bPerfectRhoG=True
            elif self.args.rho_genetic==0:
                self.bNoRhoG=True
            else:
                raise ValueError('The argument of the --rho-genetic flag is incorrectly specified, did you mean 0 or 1?')
        # see if rho environment is passed based on  --rho-environment 0; DEFAULT = False
        self.bNoRhoE=False # Boolean based on --rho-environment 0; DEFAULT = False
        if self.args.rho_environment is not None:
            if self.self.args.rho_environment==0:
                self.bNoRhoE=True
            else:
                raise ValueError('The argument of the --rho-environement flag is incorrectly specified, did you mean 0?')
            
       
        # get DataFrame for restricted genetic model based on --restricted-genetic-model mygenmodel0.txt [nolabelpheno] [nolabelfactor]; DEFAULT = None
        if self.args.restricted_genetic_model is not None:
            self.read_model(component='genetic',restricted=True)
        else:
            self.dfGenBinFY0=None
        # get DataFrame for restricted environment model based on --restricted-environment-model myenvmodel0.txt [nolabelpheno] [nolabelfactor]; DEFAULT = None
        
        if self.args.environment_model is not None:
            self.read_model(component='environment',restricted=True)
        else:
            self.dfEnvBinFY0=None
        
        # check if boolean restricted perfect rho g is given based on --restricted-rho-genetic 1; DEFAULT = False
        # check if boolean restricted no rho g is given based on  --restricted-rho-genetic 0; DEFAULT = False
        self.bPerfectRhoG0=False 
        self.bNoRhoG0=False 
        if self.args.restricted_rho_genetic is not None:
            if self.args.restricted_rho_genetic==1:
                self.bPerfectRhoG0=True
            elif self.args.restricted_rho_genetic==0:
                self.bNoRhoG0=True
            else:
                raise ValueError('The argument of the --restricted-rho-genetic flag is incorrectly specified, did you mean 0 or 1?')
        # check if rho environment is passed based on  --rho-environment 0; DEFAULT = False
        self.bNoRhoE0=False # Boolean based on --rho-environment 0; DEFAULT = False
        if self.args.restricted_rho_environment is not None:
            if self.self.args.restricted_rho_environment==0:
                self.bNoRhoE0=True
            else:
                raise ValueError('The argument of the --restricted-rho-environement flag is incorrectly specified, did you mean 0?')
        
        # check if boolean se is passed based on  --no-se 0; DEFAULT = True
        self.bSE=True
        if self.args.no_se is not None:
            self.bSE=False
        
        # check if boolean rel cutoff is passed based on --rel-cutoff 0.025; DEFAULT = None
        self.dRelCutoff =None
        if self.args.rel_cutoff is not None:
            self.dRelCutoff=self.args.rel_cutoff     # number based on --rel-cutoff 0.025; DEFAULT = None
        
        # set the number of lead and trail pcs to be dropped based on --ignore-pcs 40 [1000]; DEFAULT = 20 0
        self.iDropLeadPCs=20
        self.iDropTrailPCs=0
        if self.args.ignore_pcs is not None:
            if len(self.args.ignore_pcs)==1:
                self.iDropLeadPCs=self.args.ignore_pcs[0]
            elif len(self.args.ignore_pcs)==2:
                self.iDropLeadPCs=self.args.ignore_pcs[0]
                self.iDropTrailPCs=self.args.ignore_pcs[1]
            else:
                raise ValueError('The --ignore-pcs flag is incorrectly specified, and requires only 1 or 2 arguments?')
        
        #check if --store-iter flag is given and set iter freq if given 
        self.bStoreIter=False
        self.iStoreIterFreq=None
        if self.args.store_iter is not None:
            self.bStoreIter=True
            self.iStoreIterFreq=self.args.store_iter
        
        #check if needed to reinitialise based on --reinitialise results.iter.250.pkl
        self.sInitValsFile=None
        if self.args.reinitialise is not None:
            self.sInitValsFile=self.args.reinitialise
            
            
        #check if all coefficients are needed based on --all-coefficients
        self.bAllCoeffs=False
        if self.args.all_coefficients is not None:
            self.bAllCoeffs=True
        
        #check if newton approach is required based on --newton; DEFAULT = True
        self.bBFGS=True;
        if self.args.newton is not None:
            self.bBFGS=False
        
        # check if grad tol flag is passed based on  --grad-tol; DEFAULT = 1E-5
        self.dGradTol=10^(-5)
        if self.args.grad_tol is not None:
            self.dGradTol=self.args.dGradTol
            
        # note: when an instance of MgremlReader has all the right
        # attributes, after the __init__() method has been invoked,
        # there is no need to return anything! that instance of MgremlReader
        # can then be used by other classes to retrieve the necessary
        # settings/specifications/data
    
    def InitialiseArgsAndLogger(self):
        #create mutually exclusive groups
        #self.test=self.parser.add_argument_group('foo options', 'various (mutually exclusive) ways to do foo')
        groupGenetic = self.parser.add_mutually_exclusive_group()
        groupRestrictedGenetic = self.parser.add_mutually_exclusive_group()
        groupEnvironment = self.parser.add_mutually_exclusive_group()
        groupRestrictedEnvironment = self.parser.add_mutually_exclusive_group()
        #create arguments
        self.parser.add_argument('--out', default=None, type=str, 
                            help="Output filename prefix.")
        self.parser.add_argument('--grm', default=None, type=str,
                            help='prefix of the binary GRM from GCTA. NB: case insensitive.')
        self.parser.add_argument('--pheno', default=None, type=str,nargs='*',
                            help='Name of the phenotypes file. NB: case insensitive.')
        self.parser.add_argument('--covar', default=None, type=str,nargs='*',
                            help='Name of the covariates file. NB: case insensitive.')
        self.parser.add_argument('--covar-model', default=None, type=str,nargs='*',
                            help='Name of a file that has the specification for which covariate affects which phenotype. Possible to add the flags nolabelpheno and/or nolabelcovar.')
       
        self.parser.add_argument('--no-se', default=None,action='store_true',  
                            help='optional flag to indicate whether to skip calculating SEs (e.g. when only interested in doing a likelihood-ratio test).')
        self.parser.add_argument('--rel-cutoff', default=None, type=float,
                            help='optional flag to specify whether to apply a relatedness cut-off and the corresponding value between 0 and 1 using a greedy algorithm')
        self.parser.add_argument('--ignore-pcs', default=None, type=int, nargs='*',
                            help='Optional flag to specify how many leading PCs and trailing PCs to skip, to control for population stratifcation and remove redundant pcs, e.g. 40 10000 here. Possible to provide only one argument which is for the leading pcs.')
        self.parser.add_argument('--store-iter', default=None, type=int, 
                            help='Optional flag to specify every how many iterations you want to store results, e.g. every 50 itertions.')
        self.parser.add_argument('--reinitialise', default=None, type=str,  
                            help='optional flag to reinitialise MGREML estimation from a .pkl file generated by the --store-iter option; make sure this .pkl file has been generated from earlier iterations for estimating the model under consideration, and not some other model.')                    
        self.parser.add_argument('--all-coefficients', default=None, action='store_true', 
                            help='optional flag to report the all factor coefficients (i.e. from each factor to each trait), including the sampling covariance of these estimates.')
        self.parser.add_argument('--newton', default=None, action='store_true', 
                            help='optional flag to perform Newton instead of BFGS (not recommended, unless the model is well-identified and the number of traits is small).')
        self.parser.add_argument('--grad-tol', default=None, type=float,
                            help='A float to define the threshold on the length of the gradient vector over the number of observations and parameters, before we say MGREML has converged.')
        groupGenetic.add_argument('--genetic-model', default=None, type=str, nargs='*',
                    		help='Name of a file that has the specification for which  genetic factors apply to which phenotypes. Possible to add the flags nolabelpheno and/or nolabelfactor.')
        groupGenetic.add_argument('--rho-genetic', default=None, type=int,
                    		help='Integer 0 or 1 to identify no genetic correlation (0) and perfect genetic correlation with 1 factor (1).')
        groupRestrictedGenetic.add_argument('--restricted-genetic-model', default=None, type=str, nargs='*',
                    		help='Name of a file that has the restricted specification for which  genetic factors apply to which phenotypes. Possible to add the flags nolabelpheno and/or nolabelfactor.')
        groupRestrictedGenetic.add_argument('--restricted-rho-genetic', default=None, type=int,
                            help='Integer 0 or 1 to identify no genetic correlation (0) and perfect genetic correlation with 1 factor (1).')
        groupEnvironment.add_argument('--environment-model', default=None, type=str,nargs='*',
                    		help='Name of a file that has the specification for which  environmental factors apply to which phenotypes. Possible to add the flags nolabelpheno and/or nolabelfactor.')
        groupEnvironment.add_argument('--rho-environment', default=None, type=int,
                    		help='Integer 0 to identify no environmental correlation (0).')
        groupRestrictedEnvironment.add_argument('--restricted-environment-model', default=None, type=str,nargs='*',
                    		help='Name of a file that has the restricted specification for which  environmental factors apply to which phenotypes. Possible to add the flags nolabelpheno and/or nolabelfactor.')
        groupRestrictedEnvironment.add_argument('--restricted-rho-environment', default=None, type=int,
                    		help='Integer 0 to identify no environmental correlation (0).')					
							
        self.args = self.parser.parse_args()
        # customise the logger e.g. usings args.out
        c_handler = logging.StreamHandler()
        f_handler = logging.FileHandler(self.args.out + '.log','w+',encoding="utf-8")
        c_handler.setLevel(logging.DEBUG)
        f_handler.setLevel(logging.DEBUG)
        self.logger.setLevel(logging.DEBUG)
        # Create formatters and add it to handlers
        c_format = logging.Formatter('%(message)s')
        f_format = logging.Formatter('%(message)s - %(asctime)s ')
        c_handler.setFormatter(c_format)
        f_handler.setFormatter(f_format)
        # Add handlers to the logger
        self.logger.addHandler(c_handler)
        self.logger.addHandler(f_handler)

    def read_grm(self):
        # read binary GRM
        # names
        BinFileName = self.args.grm + ".grm.bin"
        NFileName   = self.args.grm + ".grm.N.bin"
        IDFileName  = self.args.grm + ".grm.id"
        # read IDs and sample size
        ids = pd.read_csv(IDFileName, sep = '\t', header = None)
        ids.columns=['FID','IID']
        arrays=[ids['FID'],ids['IID']]  
        tuples = list(zip(*arrays))
        index = pd.MultiIndex.from_tuples(tuples, names=['FID', 'IID'])
        iN   = len(ids.index)
        self.logger.info('{iN} individuals in {f}'.format(iN=iN, f=IDFileName))
        # read GRM from bin file
        dt  = np.dtype('f4') # Relatedness is stored as a float of size 4 in the binary file
        a = np.fromfile(BinFileName, dtype = dt)
        # create GRM as 2D-array
        mA = np.zeros((iN,iN))
        # set counter for elements of GRM read thus far
        k = 0
        # for each row
        for i in range(0,iN):
            # for each column
            for j in range(0,i+1):
                dThisVal = a[k]
                mA[i,j] = dThisVal
                mA[j,i] = dThisVal
                k += 1
        #convert to pd dataframe        
        dfA = pd.DataFrame(mA, columns=index, index=index)
        # store dfA as attribute of instance
        self.dfA = dfA
    def read_pheno(self):
        if len(self.args.pheno)==1:
            bHeader=True
        elif len(self.args.pheno)==2:
            if self.args.pheno[1]!='nolabels':
                raise ValueError('The second argument of the --pheno flag is incorrectly specified, did you mean nolabels?')
            else:
                bHeader=False
        else:
            raise ValueError('The --pheno flag requires either 1 or 2 arguments.')
        if bHeader:
            logger.info('Header specified in Phenofile')
            sPhenofile=self.args.pheno[0]
            logger.info('Reading pheno file {f}'.format(f=sPhenofile))
            dfY=pd.read_csv(sPhenofile, header = 0, sep = "\s+|\t+|,",engine='python')
            dfY=dfY.apply(pd.to_numeric, errors='coerce')
            [iN,iK]=dfY.shape
            logger.info('Found {N} Individuals and {K} Phenotypes'.format(N=iN,K=iK))
            #find number of NA  write number of missings
            iMiss=dfY.isnull().sum().sum()
            logger.info('Encountered {M} missings'.format(M=iMiss))
            
        else:
            logger.info('No header specified in Phenofile')
            sPhenofile=self.args.pheno[0]
            logger.info('Reading pheno file {f}'.format(f=sPhenofile))
            dfY=pd.read_csv(sPhenofile, header = None, sep = "\s+|\t+|,",engine='python')
            dfY=dfY.apply(pd.to_numeric, errors='coerce')
            [iN,iK]=dfY.shape
            logger.info('Found {N} Individuals and {K} Phenotypes (including FID IID)'.format(N=iN,K=iK))
            #find number of NA  write number of missings
            iMiss=dfY.isnull().sum().sum()
            logger.info('Encountered {M} missings'.format(M=iMiss))
            logger.warning('No varNames provided, creating phenotype columns labeled as 0,1,2,...')
            lNames=['FID','IID']+list(map(str,range(0,iK-2)))
            dfY.columns=lNames
                
            
                
            
        arrays=[dfY['FID'],dfY['IID']]  
        tuples = list(zip(*arrays))
        index = pd.MultiIndex.from_tuples(tuples, names=['FID', 'IID'])   
        dfY.index=index
        dfY.drop(['FID','IID'], axis = 1)   
            
        self.dfY=dfY
        
    def read_covar(self):
            #read covariance data
        if len(self.args.covar)==1:
            bHeader=True
        elif len(self.args.covar)==2:
            if self.args.covar[1]!='nolabels':
                raise ValueError('The second argument of the --covar flag is incorrectly specified, did you mean nolabels?')
            bHeader=False
        else:
            raise ValueError('The --covar flag requires either 1 or 2 arguments.')
        if bHeader:
            logger.info('Header specified in covar file')
            sCovarfile=self.args.covar[0]
            logger.info('Reading covar file {f}'.format(f=sCovarfile))
            mX=pd.read_csv(sCovarfile, header = 0, sep = "\s+|\t+|,",engine='python')
            mX=mX.apply(pd.to_numeric, errors='coerce')
            [iN,iK]=mX.shape
            logger.info('Found {N} Individuals and {K} covariates'.format(N=iN,K=iK))
            #find number of NA  write number of missings
            iMiss=mX.isnull().sum().sum()
            logger.info('Encountered {M} missings'.format(M=iMiss))
            
        else:
            logger.info('No header specified in covar file')
            sCovarfile=self.args.covar[0]
            logger.info('Reading covar file {f}'.format(f=sCovarfile))
            mX=pd.read_csv(sCovarfile, header = None, sep = "\s+|\t+|,",engine='python')
            mX=mX.apply(pd.to_numeric, errors='coerce')
            [iN,iK]=mX.shape
            logger.info('Found {N} Individuals and {K} covariates (including FID IID)'.format(N=iN,K=iK))
            #find number of NA  write number of missings
            iMiss=mX.isnull().sum().sum()
            logger.info('Encountered {M} missings'.format(M=iMiss))
            logger.warning('No varNames provided, creating covariates labeled as 0,1,2,...')
            lNames=['FID','IID']+list(map(str,range(0,iK-2)))
            mX.columns=lNames
            
            
        arrays=[mX['FID'],mX['IID']]  
        tuples = list(zip(*arrays))
        index = pd.MultiIndex.from_tuples(tuples, names=['FID', 'IID'])   
        mX.index=index
        mX.drop(['FID','IID'], axis = 1)  
            
        self.dfX=mX
    
    def read_covar_model(self):
        #function to read the covarmodel, required input is T+1 x C+1, where T = #phenotypes and C=#covariates, and save it as a pickled file
    
        #check number of input arguments given
        if len(self.args.covar_model)==1:
            bHeader=True
            bIndex=True
        elif len(self.args.covar_model)==2:
            if self.args.covar[1]=='nolabelpheno':
                bHeader=True
                bIndex=False
            elif self.args.covar[1]=='nolabelcovar':
                bHeader=False
                bIndex=True
            else:
                raise ValueError('The second argument of the --covar-model flag is incorrectly specified, did you mean nolabelpheno or nolabelcovar?')
        elif len(self.args.covar_model)==3:
            if 'nolabelpheno' in self.args.covar_model and 'nolabelcovar' in self.args.covar_model:
                bHeader=False
                bIndex=False
            else:
                raise ValueError('The second or third argument of the --covar-model flag is incorrectly specified, did you mean nolabelpheno or nolabelcovar?')
        else:
            raise ValueError('The --covar-model flag requires either 1, 2, or 3 arguments.')
            
        index_col=False
        if bIndex:
            index_col=1
        if bHeader:  
            sCovarModelfile=self.args.covar_model[0]
            logger.info('Reading covar model file {f}'.format(f=sCovarModelfile))
            dfBinXY=pd.read_csv(sCovarModelfile, header = 0, sep = "\s+|\t+|,",engine='python',index_col=index_col)
        else:
            sCovarModelfile=self.args.covar_model[0]
            logger.info('Reading covar model file {f}'.format(f=sCovarModelfile))
            dfBinXY=pd.read_csv(sCovarModelfile, header = None, sep = "\s+|\t+|,",engine='python',index_col=index_col)
        dfBinXY=dfBinXY.apply(pd.to_numeric, errors='coerce')
        [iP,iK]=dfBinXY.shape
        logger.info('Found {P} Phenotypes and {K} covariates'.format(P=iP,K=iK))
        sCovarModelOut=args.out+".covar-model.pkl"
        logger.info('Writing pickled covar-model file to {f}'.format(f=sCovarModelOut))
        self.dfBinXY=dfBinXY
        
    def read_model(self,component,restricted=False):
        #function to read the factor model, required input is T+1 x F+1, where T = #phenotypes and F=#factors, and save it as a pickled file
        #sort out if the model is an genetic or environmental and if restricted or not for reading correct file
        if component=='genetic':
            if restricted:
                lInput=self.args.restricted_genetic_model
                sFlag='restricted-genetic-model'
            else:
                lInput=self.args.genetic_model
                sFlag='genetic-model'
        elif component=='environment':   
            if restricted:
                lInput=self.args.restricted_environment_model
                sFlag='restricted-environment-model'
            else:
                lInput=self.args.environment_model
                sFlag='environment-model'
        else:
            raise ValueError('The model type is not recognized')
        
        #check number of input arguments given
        if len(lInput)==1:
            bHeader=True
            bIndex=True
        elif len(lInput)==2:
            if lInput[1]=='nolabelpheno':
                bHeader=True
                bIndex=False
            elif lInput[1]=='nolabelfactor':
                bHeader=False
                bIndex=True
            else:
                raise ValueError('The second argument of the --{S} flag is incorrectly specified, did you mean nolabelpheno or nolabelfactor?'.format(S=sFlag))
        elif len(lInput)==3:
            if 'nolabelpheno' in lInput and 'nolabelfactor' in lInput:
                bHeader=False
                bIndex=False
            else:
                raise ValueError('The second or third argument of the --{S} flag is incorrectly specified, did you mean nolabelpheno or nolabelfactor?'.format(S=sFlag))
        else:
            raise ValueError('The --{S} flag requires either 1, 2, or 3 arguments.'.format(S=sFlag))
            
        index_col=False
        if bIndex:
            index_col=0
        if bHeader:  
            sModelfile=lInput[0]
            logger.info('Reading {S} file {f}'.format(S=sFlag,f=sModelfile))
            dfBinFY=pd.read_csv(sModelfile, header = 0, sep = "\s+|\t+|,",engine='python',index_col=index_col)
        else:
            sModelfile=lInput[0]
            logger.info('Reading {S} file {f}'.format(S=sFlag,f=sModelfile))
            dfBinFY=pd.read_csv(sModelfile, header = None, sep = "\s+|\t+|,",engine='python',index_col=index_col)
        dfBinFY=dfBinFY.apply(pd.to_numeric, errors='coerce')
        [iP,iK]=dfBinFY.shape
        logger.info('Found {P} Phenotypes and {K} covariates'.format(P=iP,K=iK))
        
        if component=='genetic':
            if restricted:
                self.dfGenBinFY0=dfBinFY
            else:
                self.dfGenBinFY=dfBinFY
        elif component=='environment':   
            if restricted:
                self.dfEnvBinFY0=dfBinFY
            else:
                self.dfEnvBinFY=dfBinFY
        else:
            raise ValueError('The model type is not recognized')
