import numpy as np
import pandas as pd
import networkx as nx
import logging
import os.path
from numpy.matlib import repmat
from tqdm import tqdm
pd.options.mode.chained_assignment = None

__version__ = '0.03'
MASTHEAD  = "\n"
MASTHEAD += "########################################################################\n"
MASTHEAD += "##                                                                    ##\n"
MASTHEAD += "##                                                                    ##\n"
MASTHEAD += "##    ##      ##                                                      ##\n"
MASTHEAD += "##    ###    ###  ######   ######    ######### ##     ## ##           ##\n"
MASTHEAD += "##    ####  #### ##    ##  ##    ##  ##        ###   ### ##           ##\n"
MASTHEAD += "##    ## #### ## ##        ##    ##  ##        #### #### ##           ##\n"
MASTHEAD += "##    ##  ##  ## ##   #### ######    #######   ## ### ## ##           ##\n"
MASTHEAD += "##    ##      ## ##    ##  ##   ##   ##        ##  #  ## ##           ##\n"
MASTHEAD += "##    ##      ## ##    ##  ##    ##  ##        ##     ## ##           ##\n"
MASTHEAD += "##    ##      ##  ######   ##     ## ######### ##     ## #########    ##\n"
MASTHEAD += "##                                                                    ##\n"
MASTHEAD += "##                                                                    ##\n"
MASTHEAD += "########################################################################\n"
MASTHEAD += "##                                                                    ##\n"
MASTHEAD += "##  BETA VERSION {V}                                                 ##\n".format(V=__version__)
MASTHEAD += "##                                                                    ##\n"
MASTHEAD += "##  (C) 2021 Ronald de Vlaming and Eric Slob                          ##\n"
MASTHEAD += "##                                                                    ##\n"
MASTHEAD += "##  Vrije Universiteit Amsterdam and University of Cambridge          ##\n"
MASTHEAD += "##  GNU General Public License v3                                     ##\n"
MASTHEAD += "##                                                                    ##\n"
MASTHEAD += "########################################################################\n"

class MgremlReader:

    # set list of indices for FID and IID columns
    # in pheno and covar files to [0,1]
    # and give labels
    lIndFID_IID = [0,1]
    lLabelsFID_IID = ['FID', 'IID']
    
    # abbreviated labels for phenos, covars,
    # genetic factors and environment factors
    sPhe = 'phe'
    sCov = 'cov'
    sGen = 'gen'
    sEnv = 'env'
    
    # set boolean for restricted models
    bNull = True
    
    # set default values for dropping leading and trailing PCs
    iDropLeadPCsDefault = 20
    iDropTrailPCsDefault = 0
    
    # set default convergence treshold
    dGradTol = 1E-5
    
    # highest condition number per trait allowed in phenotypic
    # correlation matrix; threshold corresponds all correlations >=0.95
    dCondThreshold = 19 
    
    # define what many dummies are
    iManyDummies = 1000
    
    # define prefix of label for each factor
    sLabelFactor = 'factor '
    
    # define intercept position (column), value (1), and label
    iIntPos = 0
    iIntVal = 1
    sIntLab = 'intercept'
    
    # define seed for random-number generator used for relatedness pruning
    iSeed = 1809234
    
    @staticmethod
    def SetPerfectRho(lLabels):
        lFactorName = [MgremlReader.sLabelFactor + str(0)]
        dfBinFY = pd.DataFrame(data=1, index=pd.Index(lLabels), columns=pd.Index(lFactorName))
        return dfBinFY
        
    @staticmethod
    def SetNoRho(lLabels):
        iT = len(lLabels)
        lFactors = [MgremlReader.sLabelFactor + str(x) for x in range(0,iT)]
        dfBinFY = pd.DataFrame(data=np.eye(iT), index=pd.Index(lLabels), columns=pd.Index(lFactors))
        return dfBinFY
        
    @staticmethod
    def SetNoVar(lLabels):
        lFactorName = [MgremlReader.sLabelFactor + str(0)]
        dfBinFY = pd.DataFrame(data=0, index=pd.Index(lLabels), columns=pd.Index(lFactorName))
        return dfBinFY
    
    def __init__(self, parser, logger, process, bCopy=False, mdData=None):
        if not(bCopy):
            # store logger, parser, and process as attributes of instance
            self.logger = logger
            self.parser = parser
            self.process = process
            # initialise arguments and logger
            self.InitialiseArgumentsAndLogger()
            # print welcome screen and given input args
            self.PrintWelcome()
            # assume we will carry out an analysis
            self.bAnalyse = True
            # if GRM, pheno or output has not been specified: stop analysing
            if self.args.grm is None:
                self.logger.error('Error: no GRM was specified.')
                self.bAnalyse = False
            if self.args.pheno is None:
                self.logger.error('Error: no phenotype file was specified.')
                self.bAnalyse = False            
            if self.args.out is None:
                self.logger.error('Error: no prefix for the output files was specified.')
                self.bAnalyse = False
            # if we can analyse, read in data
            if self.bAnalyse:
                self.logger.info('1. READING IN ALL DATA, MODELS, AND INPUT OPTIONS')
                self.logger.info('Current memory usage is ' + str(int((self.process.memory_info().rss)/(1024**2))) + 'MB')
                self.logger.info('READING OPTIONS')
                # determine whether we will drop missings
                self.DetermineIfDropMissings()
                # determine whether --grm-cutoff has been used
                self.SetRelCutOff()
                # determine how many PCs to drop
                self.SetNumberOfPCs()
                # determine whether SEs are desired
                self.NeedSEs()
                # determine whether intercept needs to be added
                self.NeedIntercept()
                # determine whether all coeffs are desired
                self.NeedAllCoeffs()
                # determine whether all variance components are desired
                self.NeedVarianceComponents()
                # determine whether we do BFGS or Newton
                self.DoBFGS()
                # set convergence threshold
                self.SetGradTol()
                # determine whether we story anything during iterations
                self.SetIterStoreFreq()
                # assert if one or two models have been specified
                self.IsNested()
                # determine if we need to reinitialise anything
                self.NeedToReinitialise()
                self.NeedToReinitialise(MgremlReader.bNull)
                # determine if we want a mediation analysis
                self.IsMediation()
                # determine if we want pairwise bivariate analyses
                self.IsPairwise()
                # determine if we can ignore collinearity
                self.IgnoreCollinearity()
                # print update
                self.logger.info('READING MODELS')
                # determine whether any correlations are fixed to zero or one
                # or if genetic variances are fixed to zero
                self.FindFixedRhoVar()
                # if covar model has been specified: read
                if self.args.covar_model is not None:
                    self.ReadModel(MgremlReader.sCov)
                else:
                    self.dfBinXY = None
                # if genetic model has been specified: read
                if self.args.genetic_model is not None:
                    if self.bPairwise:
                        raise SyntaxError('--genetic-model cannot be combined with --pairwise')
                    elif self.bMediation:
                        raise SyntaxError('--mediation cannot be combined with --genetic-model')
                    self.ReadModel(MgremlReader.sGen)
                else:
                    self.dfGenBinFY = None
                # if genetic model under null has been specified: read
                if self.args.restricted_genetic_model is not None:
                    if self.bPairwise:
                        raise SyntaxError('--restricted-genetic-model cannot be combined with --pairwise')
                    self.ReadModel(MgremlReader.sGen, MgremlReader.bNull)
                else:
                    self.dfGenBinFY0 = None
                # if environment model has been specified: read
                if self.args.environment_model is not None:
                    if self.bPairwise:
                        raise SyntaxError('--environment-model cannot be combined with --pairwise')
                    elif self.bMediation:
                        raise SyntaxError('--mediation cannot be combined with --environment-model')
                    self.ReadModel(MgremlReader.sEnv)
                else:
                    self.dfEnvBinFY = None
                # if genetic model under null has been specified: read
                if self.args.restricted_environment_model is not None:
                    if self.bPairwise:
                        raise SyntaxError('--restricted-environment-model cannot be combined with --pairwise')
                    self.ReadModel(MgremlReader.sEnv, MgremlReader.bNull)
                else:
                    self.dfEnvBinFY0 = None
                # print update
                self.logger.info('READING DATA')
                # read phenotype file
                self.ReadData(MgremlReader.sPhe)
                # if covariate file specified: read
                if self.args.covar is not None:
                    self.ReadData(MgremlReader.sCov)
                else: # else set covariates to NoneType
                    self.dfX = None
                # find out if we have covariates
                self.DetermineIfCovsAreGiven()
                # read GRM
                self.ReadGRM()
                # print update
                self.logger.info('2. CLEANING YOUR DATA')
                self.logger.info('Current memory usage is ' + str(int((self.process.memory_info().rss)/(1024**2))) + 'MB')
                # add intercept if required
                self.AddIntercept()
                # clean up the specification of the part of the model on covariates
                self.CleanSpecificationCovariates()
                # check if there are duplicates and/or rank defficiency
                self.CheckDuplicatesAndRank()
                # find overlapping individuals and sort data
                self.FindOverlapAndSort()
                # drop observations where missingness affects all traits or any trait
                self.DropMissings()
                # apply relatedness pruning
                self.PruneByRelatedness()
                # contruct pheno-specific dummies to address remaining missingness
                self.CreateDummies()
                # store variable names
                self.lPhenos = self.dfY.columns.tolist()
                if self.bCovs:
                    self.lCovs = self.dfX.columns.tolist()
                # set constrained models if needed
                self.SetConstrainedModels()
                # finalise Mgreml data using canonical transformation
                self.FinaliseData()
                # print update
                self.logger.info('Data cleaning completed')
                self.logger.info('Current memory usage is ' + str(int((self.process.memory_info().rss)/(1024**2))) + 'MB\n')
        else:
            self.logger = mdData.logger
            self.parser = mdData.parser
            self.process = mdData.process
            self.sPrefix = mdData.sPrefix
            self.bDropMissings = mdData.bDropMissings
            self.bRelCutoff = mdData.bRelCutoff
            self.dRelCutoff = mdData.dRelCutoff
            self.bPairwise = mdData.bPairwise
            self.bMediation = mdData.bMediation
            self.bIgnoreCollinearity = mdData.bIgnoreCollinearity
            self.bCovs = mdData.bCovs
            self.bSameCovs = mdData.bSameCovs
            self.bNested = mdData.bNested
            self.bAnalyse = mdData.bAnalyse
            self.dLogDetXTX = mdData.dLogDetXTX
            self.iT = mdData.iT
            self.iN = mdData.iN
            self.iDropLeadPCs = mdData.iDropLeadPCs
            self.iDropTrailPCs = mdData.iDropTrailPCs
            self.bSEs = mdData.bSEs
            self.bIntercept = mdData.bIntercept
            self.bAllCoeffs = mdData.bAllCoeffs
            self.bVarComp = mdData.bVarComp
            self.bBFGS = mdData.bBFGS
            self.dGradTol = mdData.dGradTol
            self.bStoreIter = mdData.bStoreIter
            self.iStoreIterFreq = mdData.iStoreIterFreq
            self.bReinitialise = mdData.bReinitialise
            self.bReinitialise0 = mdData.bReinitialise0
            self.sInitValsFile = mdData.sInitValsFile
            self.sInitValsFile0 = mdData.sInitValsFile0
            self.bPerfectRhoG = mdData.bPerfectRhoG
            self.bNoRhoG = mdData.bNoRhoG
            self.bPerfectRhoG0 = mdData.bPerfectRhoG0
            self.bNoRhoG0 = mdData.bNoRhoG0
            self.bNoRhoE = mdData.bNoRhoE
            self.bNoRhoE0 = mdData.bNoRhoE0
            self.bNoVarG = mdData.bNoVarG
            self.bNoVarG0 = mdData.bNoVarG0
            self.mYT = mdData.mYT.copy()
            self.mY = self.mYT.T
            self.vD = mdData.vD.copy()
            self.vDSq = mdData.vDSq.copy()
            if mdData.dfGenBinFY is not None:
                self.dfGenBinFY = mdData.dfGenBinFY.copy()
            else:
                self.dfGenBinFY = None
            if mdData.dfEnvBinFY is not None:
                self.dfEnvBinFY = mdData.dfEnvBinFY.copy()
            else:
                self.dfEnvBinFY = None
            if mdData.dfGenBinFY0 is not None:
                self.dfGenBinFY0 = mdData.dfGenBinFY0.copy()
            else:
                self.dfGenBinFY0 = None
            if mdData.dfEnvBinFY0 is not None:
                self.dfEnvBinFY0 = mdData.dfEnvBinFY0.copy()
            else:
                self.dfEnvBinFY0 = None
            self.mCovY = mdData.mCovY.copy()
            self.lPhenos = mdData.lPhenos.copy()
            if self.bCovs:
                self.iK = mdData.iK
                self.mXT = mdData.mXT.copy()
                self.mX = self.mXT.T
                self.lCovs = mdData.lCovs.copy()
                if not(self.bSameCovs):
                    self.mBinXY = mdData.mBinXY.copy()
                    self.vLogDetXTX = mdData.vLogDetXTX.copy()
                    self.vIndCovs = np.array(np.where(np.array(self.mBinXY).ravel()==1)).ravel()
                    self.iKtotal = self.mBinXY.sum()
    
    def InitialiseArgumentsAndLogger(self):
        #create mutually exclusive groups
        groupGenetic = self.parser.add_mutually_exclusive_group()
        groupRestrictedGenetic = self.parser.add_mutually_exclusive_group()
        groupEnvironment = self.parser.add_mutually_exclusive_group()
        groupRestrictedEnvironment = self.parser.add_mutually_exclusive_group()
        #create arguments
        self.parser.add_argument('--grm', metavar = 'PREFIX', default = None, type = str,
                            help = 'prefix of binary GRM')
        self.parser.add_argument('--grm-cutoff', metavar = 'THRESHOLD', default = None, type = float,
                            help = 'option to drop individuals using greedy algorithm, such that there is no relatedness in GRM in excess of threshold for remaining individuals')
        self.parser.add_argument('--adjust-pcs', metavar = '', default = None, type = int, nargs = '+',
                            help = '\b\b\b\b\b\b\b\b\bINTEGER [INTEGER] option to specify for how many leading principal components (PCs) from genetic data to adjust (to control for population stratification) and for how many trailing PCs to adjust (for computational efficiency); if just one non-negative integer is specified this is taken as the number of leading PCs to adjust for')
        self.parser.add_argument('--pheno', metavar = '', default = None, type = str, nargs = '+',
                            help = '\b\b\b\b\b\b\b\b\b\b\b\b\b\bFILENAME [nolabelpheno] phenotype file: should be comma-, space-, or tab-separated, with one row per individual, with FID and IID as first two fields, followed by a field per phenotype; can be followed by optional flag nolabelpheno, e.g. --pheno mypheno.txt nolabelpheno, but we recommend to label phenotypes')
        self.parser.add_argument('--mediation', action = 'store_true',
                            help = 'option to perform a mediation analysis, in line with the structural equations model proposed by Rietveld et al. (2021) and based on estimates from a saturated bivariate model; the first phenotype in the phenotype file is assumed to act as mediator for the genetic component of the second phenotype in the phenotype file; all further phenotypes are ignored; this flag cannot be combined with --(restricted-)genetic-model, --(restricted-)rho-genetic, --(restricted-)no-var-genetic, --(restricted-)environment-model, --(restricted-)rho-environment, and --(restricted-)reinitialise')
        self.parser.add_argument('--drop-missings', action = 'store_true',
                            help = 'option to drop all observations from data with at least one missing phenotype or at least one missing covariate')
        self.parser.add_argument('--no-intercept', action = 'store_true',
                            help = 'option to indicate an intercept should not be included automatically as covariate')
        self.parser.add_argument('--covar', metavar = '', default = None, type = str, nargs = '+',
                            help = '\b\b\b\b\b\b\b\b\b\b\b\b\b\bFILENAME [nolabelcovar] optional covariate file: should be comma-, space-, or tab-separated, with one row per individual, with FID and IID as first two fields, followed by a field per covariate; can be followed by optional flag nolabelcovar, e.g. --covar mycovar.txt nolabelcovar, but we recommend to label covariates; WARNING: do not include principal components from genetic data as covariates, use --adjust-pcs instead')
        self.parser.add_argument('--covar-model', metavar = '', default = None, type = str, nargs = '+',
                            help = 'FILENAME [nolabelpheno] [nolabelcovar] optional covariate model file: should be comma-, space-, or tab-separated, with one row per phenotype and one column per covariate; can be followed by optional flags nolabelpheno and/or nolabelcovar, but we recommend to label phenotypes and covariates; without --covar-model, all covariates are assumed to apply to all traits')
        groupGenetic.add_argument('--genetic-model', metavar = '', default = None, type = str, nargs = '+',
                            help = 'FILENAME [nolabelpheno] [nolabelfactor] optional genetic model file: should be comma-, space-, or tab-separated, with one row per phenotype and one column per genetic factor; can be followed by optional flags nolabelpheno and/or nolabelfactor, but we recommend to label phenotypes and genetic factors')
        groupGenetic.add_argument('--rho-genetic', metavar = '0 or 1', choices = [0, 1], default = None, type = int,
                            help = 'option followed by 0 or 1, forcing all genetic correlations to take on the specified value; this flag cannot be combined with --genetic-model')
        groupGenetic.add_argument('--no-var-genetic', action = 'store_true',
                            help = 'option to force all genetic variances to equal zero; this flag cannot be combined with --genetic-model and/or --rho-genetic')
        groupRestrictedGenetic.add_argument('--restricted-genetic-model', metavar = '', default = None, type = str, nargs = '+',
                            help = 'FILENAME [nolabelpheno] [nolabelfactor] optional restricted genetic model file: should be comma-, space-, or tab-separated, with one row per phenotype and one column per genetic factor; can be followed by optional flags nolabelpheno and/or nolabelfactor, but we recommend to label phenotypes and genetic factors')
        groupRestrictedGenetic.add_argument('--restricted-rho-genetic', metavar = '0 or 1', choices = [0, 1], default = None, type = int,
                            help = 'option followed by 0 or 1, forcing all genetic correlations in the restricted model to take on the specified value; this flag cannot be combined with --restricted-genetic-model')
        groupRestrictedGenetic.add_argument('--restricted-no-var-genetic', action = 'store_true',
                            help = 'option to force all genetic variances in the restricted model to equal zero; this flag cannot be combined with --restricted-genetic-model and/or --restricted-rho-genetic')
        groupEnvironment.add_argument('--environment-model', metavar = '', default = None, type = str, nargs = '+',
                            help = 'FILENAME [nolabelpheno] [nolabelfactor] optional environment model file: should be comma-, space-, or tab-separated, with one row per phenotype and one column per environment factor; can be followed by optional flags nolabelpheno and/or nolabelfactor, but we recommend to label phenotypes and environment factors')
        groupEnvironment.add_argument('--rho-environment', metavar = '0', choices = [0], default = None, type = int,
                            help = 'option followed by 0, forcing all environment correlations to zero; this flag cannot be combined with --environment-model')
        groupRestrictedEnvironment.add_argument('--restricted-environment-model', metavar = '', default = None, type = str, nargs = '+',
                            help = 'FILENAME [nolabelpheno] [nolabelfactor] optional restricted environment model file: should be comma-, space-, or tab-separated, with one row per phenotype and one column per environment factor; can be followed by optional flags nolabelpheno and/or nolabelfactor, but we recommend to label phenotypes and environment factors')
        groupRestrictedEnvironment.add_argument('--restricted-rho-environment', metavar = '0', choices = [0], default = None, type = int,
                            help = 'option followed by 0, forcing all environment correlations in the restricted model to zero; this flag cannot be combined with --restricted-environment-model')
        self.parser.add_argument('--no-se', action = 'store_true',
                            help = 'option to skip calculation of standard errors and covariance matrix of estimates')
        self.parser.add_argument('--factor-coefficients', action = 'store_true', 
                            help = 'option to report estimated factor coefficients')
        self.parser.add_argument('--variance-components', action = 'store_true', 
                            help = 'option to report estimated variance components')
        self.parser.add_argument('--newton', action = 'store_true',
                            help = 'option to use Newton method instead of BFGS; not recommended, unless the model is well-defined, starting values are of good quality, and the number of traits is small')
        self.parser.add_argument('--grad-tol', metavar = 'THRESHOLD', default = None, type = float,
                            help = 'option to set convergence threshold on the length of the gradient vector per parameter, per observation, different from the default value of 1E-5')
        self.parser.add_argument('--store-iter', metavar = 'INTEGER', default = None, type = int, 
                            help = 'option to specify every how many iterations you want to store results')
        self.parser.add_argument('--reinitialise', metavar = 'FILENAME', default = None, type = str,  
                            help = 'option to reinitialise mgreml for a model and its estimates from a .pkl file stored by --store-iter')
        self.parser.add_argument('--restricted-reinitialise', metavar = 'FILENAME', default = None, type = str,  
                            help = 'option to reinitialise mgreml for a restricted model and its estimates from a .pkl file generated by --store-iter')
        self.parser.add_argument('--pairwise', action = 'store_true',
                            help = 'option to perform pairwise bivariate estimation instead of multivariate estimation; cannot be combined with --mediation, --(restricted-)genetic-model, --(restricted-)no-var-genetic, --(restricted-)environment-model, --factor-coefficients, --variance-components, --store-iter, and/or --(restricted-)reinitialise; WARNING: the number of output files can be very large')
        self.parser.add_argument('--ignore-collinearity', action = 'store_true',
                            help = 'option to ignore multicollinearity between phenotypes; please use this option only as a last resort; model may be poorly identified when your phenotype data is perfectly collinear; preferred route to solve collinearity is to consider e.g. a subset of phenotypes')
        self.parser.add_argument('--out', metavar = 'PREFIX', default = None, type = str,
                            help = 'prefix of output files')
        # try to parse the input arguments
        try:
            # parse the input options
            self.args = self.parser.parse_args()
        except Exception:
            # if doesn't work: print error
            raise SyntaxError('you specified incorrect input options.')
        # customise the logger using the prefix for output-files
        c_handler = logging.StreamHandler()
        # if no --out option has not been specified
        # use generic name for log-file and so on
        if self.args.out is not None:
            # get directory name if present within prefix
            sDir = os.path.dirname(self.args.out)
            # check if output holds directory name at all, and if so whether it doesn't exist
            if not(sDir == '') and not(os.path.isdir(sDir)):
                # if so, raise an error
                raise ValueError('prefix specified using --out may start with a directory name; this directory must exist however. ' + sDir + ' is not a directory.')
            self.sPrefix = self.args.out + '.'
            sFileOut = self.sPrefix + 'log'
        else:
            self.sPrefix = 'output.'
            sFileOut = self.sPrefix + 'log'
        # set the name for the log file
        f_handler = logging.FileHandler(sFileOut,'w+',encoding="utf-8")
        c_handler.setLevel(logging.DEBUG)
        f_handler.setLevel(logging.DEBUG)
        self.logger.setLevel(logging.DEBUG)
        # for output to the console, just print warnings, info, errors, etc.
        c_format = logging.Formatter('%(message)s')
        # for output to the log file also add timestamps
        f_format = logging.Formatter('%(asctime)s: %(message)s')
        c_handler.setFormatter(c_format)
        f_handler.setFormatter(f_format)
        # add handlers to the logger
        self.logger.addHandler(c_handler)
        self.logger.addHandler(f_handler)
        
    def PrintWelcome(self):
        # try to print welcome screen
        try:
            defaults = vars(self.parser.parse_args(''))
            opts = vars(self.args)
            non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
            header = MASTHEAD
            header += "\nYour call: \n"
            header += './mgreml.py \\\n'
            options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults]
            header += '\n'.join(options).replace('True','').replace('False','').replace("', \'", ' ').replace("']", '').replace("['", '').replace('[', '').replace(']', '').replace(', ', ' ').replace('  ', ' ')
            header = header[0:-1]+'\n'
            self.logger.info(header)
        except Exception:
            raise SyntaxError('you specified incorrect input options.')
    
    def DetermineIfDropMissings(self):
        # if --drop-missings option used
        if self.args.drop_missings:
            # do so
            self.bDropMissings = True
            self.logger.info('Observations with any missing data will be dropped from analysis')
        else:
            self.bDropMissings = False
            self.logger.info('MGREML will construct phenotype-specific dummies to control for missing data') 
            
    def SetRelCutOff(self):
        # if --grm-cutoff used
        if self.args.grm_cutoff is not None:
            if (self.args.grm_cutoff < 0):
                raise ValueError('--grm-cutoff should be followed by non-negative value.')
            else:
                self.bRelCutoff = True
                self.dRelCutoff = self.args.grm_cutoff
                self.logger.info('A relatedness cutoff of ' + str(self.dRelCutoff) + ' will be applied to your GRM.')
        else:
            self.bRelCutoff = False
            self.dRelCutoff = None
            self.logger.info('No relatedness cutoff will be applied to your GRM.')
    
    def SetNumberOfPCs(self):
        # if no. of PCs specified
        if self.args.adjust_pcs is not None:
            if len(self.args.adjust_pcs)==1:
                if (self.args.adjust_pcs[0] < 0):
                    raise ValueError('--adjust-pcs should be followed by one or two non-negative integers')
                self.iDropLeadPCs = self.args.adjust_pcs[0]
                self.iDropTrailPCs = MgremlReader.iDropTrailPCsDefault
            elif len(self.args.adjust_pcs)==2:
                if (self.args.adjust_pcs[0] < 0) or (self.args.adjust_pcs[1] < 0):
                    raise ValueError('--adjust-pcs should be followed by one or two non-negative integers')
                self.iDropLeadPCs = self.args.adjust_pcs[0]
                self.iDropTrailPCs = self.args.adjust_pcs[1]
                self.logger.info('Adjusting for ' + str(self.iDropTrailPCs) + ' trailing eigenvectors from GRM to improve computational efficiency')
            else:
                raise SyntaxError('--adjust-pcs should be followed by one or two non-negative integers')
        else:
            self.iDropLeadPCs = MgremlReader.iDropLeadPCsDefault
            self.iDropTrailPCs = MgremlReader.iDropTrailPCsDefault
        self.logger.info('Adjusting for ' + str(self.iDropLeadPCs) + ' leading eigenvectors from GRM to control for population stratification')
    
    def NeedSEs(self):
        # if --no-se option used
        if self.args.no_se:
            # report no SEs
            self.bSEs = False
            self.logger.info('Your results will not include standard errors.')
        else:
            self.bSEs = True
            self.logger.info('Your results will include standard errors.')
            
    def NeedIntercept(self):
        # if --no-intercept option used
        if self.args.no_intercept:
            # do not automatically add intercept
            self.bIntercept = False
        else:
            self.bIntercept = True
            self.logger.info('An intercept will be added automatically to your set of covariates.')
            self.logger.info('This intercept will apply to all phenotypes in your data.')
            
    def NeedAllCoeffs(self):
        # if --factor-coefficients option used
        if self.args.factor_coefficients:
            # report them
            self.bAllCoeffs = True
            self.logger.info('Your results will include estimates of all factor coefficients')
        else:
            self.bAllCoeffs = False
    
    def NeedVarianceComponents(self):
        # if --variance-components option used
        if self.args.variance_components:
            # report them
            self.bVarComp = True
            self.logger.info('Your results will include estimates of all variance components')
        else:
            self.bVarComp = False
    
    def DoBFGS(self):
        # if --newton option used
        if self.args.newton:
            # don't do BFGS but print warning
            self.bBFGS = False
            self.logger.info('MGREML will use a Newton algorithm instead of BFGS for estimation')
            self.logger.warning('Warning: Newton not recommended for a large number of traits')
        else:
            self.bBFGS = True
            self.logger.info('MGREML will use a BFGS algorithm for estimation')
    
    def SetGradTol(self):
        if self.args.grad_tol is not None:
            if (self.args.grad_tol <= 0):
                raise ValueError('--grad-tol should be followed by a positive number e.g. 1E-6, 1e-5, or 0.0001')
            self.dGradTol = self.args.grad_tol
        else:
            self.dGradTol = MgremlReader.dGradTol
        self.logger.info('Setting convergence threshold for length of gradient per observation per parameter to ' + str(self.dGradTol))
    
    def SetIterStoreFreq(self):
        if self.args.store_iter is not None:
            if (self.args.store_iter < 1):
                raise ValueError('--store-iter should be followed by a positive integer')
            self.bStoreIter = True
            self.iStoreIterFreq = self.args.store_iter
            self.logger.info('Storing parameter estimates every ' + str(self.iStoreIterFreq) + ' iterations')
        else:
            self.bStoreIter = False
            self.iStoreIterFreq = None
            self.logger.info('Not storing parameter estimates during iterations')
    
    def IsNested(self):
        self.bNested = (self.args.restricted_genetic_model is not None) \
                    or (self.args.restricted_environment_model is not None) \
                    or (self.args.restricted_rho_genetic is not None) \
                    or (self.args.restricted_rho_environment is not None) \
                    or (self.args.restricted_no_var_genetic) \
                    or (self.args.restricted_reinitialise is not None)
        if self.bNested:
            self.logger.info('You specified two models for comparison. Results will include a likelihood-ratio test comparing the restricted model (null hypothesis) to the main model (alternative hypothesis).')
        else:
            self.logger.info('You specified only the main model. No likelihood-ratio test will be performed.')
    
    def NeedToReinitialise(self, bNull = False):
        if bNull:
            sFile = self.args.restricted_reinitialise
        else:
            sFile = self.args.reinitialise
        if sFile is not None:
            if not(os.path.isfile(sFile)):
                raise TypeError('iteration file ' + sFile + ' does not exist')
            elif sFile[-4:] != '.pkl':
                raise TypeError('file ' + sFile + ' is not a .pkl file')
            if bNull:
                self.bReinitialise0 = True
                self.sInitValsFile0 = sFile
                self.logger.info('MGREML will reinitialise restricted model using estimates in ' + self.sInitValsFile0)
            else:
                self.bReinitialise = True
                self.sInitValsFile = sFile
                self.logger.info('MGREML will reinitialise using estimates in ' + self.sInitValsFile)
        else:
            if bNull:
                self.bReinitialise0 = False
                self.sInitValsFile0 = None
            else:
                self.bReinitialise = False
                self.sInitValsFile = None
    
    def IsMediation(self):
        if self.args.mediation:
            if self.bNested:
                raise SyntaxError('--mediation cannot be combined with any flag of the type --restricted-...')
            elif self.bReinitialise:
                raise SyntaxError('--mediation cannot be combined with --reinitialise')
            self.bMediation = True
            self.logger.info('MGREML will perform a mediation analysis based on results from a bivariaite model for the first two phenotypes in your phenotype file.')
        else:
            self.bMediation = False
    
    def IsPairwise(self):
        if self.args.pairwise:
            if self.bAllCoeffs:
                raise SyntaxError('--pairwise cannot be combined with --factor-coefficients')
            elif self.bVarComp:
                raise SyntaxError('--pairwise cannot be combined with --variance-components')
            elif self.bStoreIter:
                raise SyntaxError('--pairwise cannot be combined with --store-iter')
            elif self.bReinitialise:
                raise SyntaxError('--pairwise cannot be combined with --reinitialise')
            elif self.bReinitialise0:
                raise SyntaxError('--pairwise cannot be combined with --restricted-reinitialise')
            elif self.bMediation:
                raise SyntaxError('--pairwise cannot be combined with --mediation')
            self.bPairwise = True
            self.logger.info('MGREML will perform pairwise bivariate estimation instead of multivariate estimation.')
        else:
            self.bPairwise = False
    
    def IgnoreCollinearity(self):
        if self.args.ignore_collinearity:
            self.bIgnoreCollinearity = True
            self.logger.info('MGREML will ignore multicollinearity between phenotypes.')
        else:
            self.bIgnoreCollinearity = False
    
    def FindFixedRhoVar(self):
        # set all booleans for fixed rhoG and rhoE to False
        self.bPerfectRhoG = False
        self.bNoRhoG = False
        self.bPerfectRhoG0 = False
        self.bNoRhoG0 = False
        self.bNoRhoE = False
        self.bNoRhoE0 = False
        # set booleans for varG fixed to zero to False
        self.bNoVarG = False
        self.bNoVarG0 = False
        # assess whether rhoG is perfect or zero
        if self.args.rho_genetic is not None:
            if self.bReinitialise:
                raise SyntaxError('--rho-genetic cannot be combined with --reinitialise, as the .pkl file is used to set the model')
            elif self.bMediation:
                raise SyntaxError('--mediation cannot be combined with --rho-genetic')
            if self.args.rho_genetic==1:
                self.logger.info('Genetic correlations in the main model all set to one.')
                self.bPerfectRhoG = True
            else:
                self.logger.info('Genetic correlations in the main model all set to zero.')
                self.bNoRhoG = True
        # assess whether rhoG is perfect or zero in the restricted model
        if self.args.restricted_rho_genetic is not None:
            if self.bReinitialise0:
                raise SyntaxError('--restricted-rho-genetic cannot be combined with --restricted-reinitialise, as the .pkl file is used to set the restricted model')
            if self.args.restricted_rho_genetic==1:
                self.logger.info('Genetic correlations in the null model all set to one.')
                self.bPerfectRhoG0 = True
            else:
                self.logger.info('Genetic correlations in the null model all set to zero.')
                self.bNoRhoG0 = True
        # assess whether rhoE is zero
        if self.args.rho_environment is not None:
            if self.bReinitialise:
                raise SyntaxError('--rho-environment cannot be combined with --reinitialise, as the .pkl file is used to set the model')
            elif self.bMediation:
                raise SyntaxError('--mediation cannot be combined with --rho-environment')
            self.logger.info('Environment correlations in the main model all set to zero.')
            self.bNoRhoE = True
        # assess whether rhoE is zero in the restricted model  
        if self.args.restricted_rho_environment is not None:
            if self.bReinitialise0:
                raise SyntaxError('--restricted-rho-environment cannot be combined with --restricted-reinitialise, as the .pkl file is used to set the restricted model')
            self.logger.info('Environment correlations in the null model all set to zero.')
            self.bNoRhoE0 = True
        # assess whether genetic variance is zero
        if self.args.no_var_genetic:
            if self.bReinitialise:
                raise SyntaxError('--no-var-genetic cannot be combined with --reinitialise, as the .pkl file is used to set the model')
            elif self.bPairwise:
                raise SyntaxError('--no-var-genetic cannot be combined with --pairwise')
            elif self.bMediation:
                raise SyntaxError('--mediation cannot be combined with --no-var-genetic')
            self.logger.info('Genetic variance in the main model set to zero.')
            self.bNoVarG = True
        # assess whether genetic variance is zero in the restricted model
        if self.args.restricted_no_var_genetic:
            if self.bReinitialise0:
                raise SyntaxError('--restricted-no-var-genetic cannot be combined with --restricted-reinitialise, as the .pkl file is used to set the restricted model')
            elif self.bPairwise:
                raise SyntaxError('--restricted-no-var-genetic cannot be combined with --pairwise')
            self.logger.info('Genetic variance in the null model set to zero.')
            self.bNoVarG0 = True
        # assess if --rho-genetic 1 and --restricted-rho-genetic 0: non-nested
        if self.bPerfectRhoG and self.bNoRhoG0:
            raise ValueError('--rho-genetic 1 and --restricted-rho-genetic 0 are not nested')
        # assess if --rho-genetic 0 and --restricted-rho-genetic 1: non-nested
        if self.bNoRhoG and self.bPerfectRhoG0:
            raise ValueError('--rho-genetic 0 and --restricted-rho-genetic 1 are not nested')
    
    def ReadModel(self, sType, bNull = False):
        # set no label string for phenotypes (in rows of all models)
        sNoLabelIndex = 'nolabelpheno'
        # figure out what type of input data we have
        # and set strings appropriately
        if sType == MgremlReader.sCov:
            sNoLabelHeader = 'nolabelcovar'
            sData = 'covariate model'
            sDescr = 'covariates'
            sDescrShort = 'covariate'
            sOption = '--covar-model'
            lArgs = self.args.covar_model
        elif sType == MgremlReader.sGen:
            sNoLabelHeader = 'nolabelfactor'
            sDescr = 'genetic factors'
            sDescrShort = 'factor'
            if bNull:
                sData = 'restricted genetic model'
                sOption = '--restricted-genetic-model'
                lArgs = self.args.restricted_genetic_model
            else:
                sData = 'genetic model'
                sOption = '--genetic-model'
                lArgs = self.args.genetic_model
        elif sType == MgremlReader.sEnv:
            sNoLabelHeader = 'nolabelfactor'
            sDescr = 'environment factors'
            sDescrShort = 'factor'
            if bNull:
                sData = 'restricted environment model'
                sOption = '--restricted-environment-model'
                lArgs = self.args.restricted_environment_model
            else:
                sData = 'environment model'
                sOption = '--environment-model'
                lArgs = self.args.environment_model
        else:
            raise TypeError('MgremlReader.ReadModel() can only be used to read a model specification for factors or covariates.')
        if sType != MgremlReader.sCov:
            if bNull:
                if self.bReinitialise0:
                    raise SyntaxError(sOption + ' cannot be combined with --restricted-reinitialise, as the .pkl file is used to set the restricted model')
            else:
                if self.bReinitialise:
                    raise SyntaxError(sOption + ' cannot be combined with --reinitialise, as the .pkl file is used to set the model')
        # if no. of input args. is one
        if len(lArgs) == 1:
            # we have a header and index
            bHeader = True
            bIndex = True
        elif len(lArgs)==2: # if it's 2
            if lArgs[1] == sNoLabelIndex:
                # we have a header but no index
                bHeader = True
                bIndex = False
            elif lArgs[1] == sNoLabelHeader:
                # we have an index but no header
                bHeader = False
                bIndex = True
            else:
                # raise an error
                raise ValueError('the second argument of ' + sOption + ' is incorrect. Did you mean ' + sNoLabelIndex + ' or ' + sNoLabelHeader + '?')
        elif len(lArgs)==3: # if it's 3
            if (sNoLabelIndex in lArgs[1:]) and (sNoLabelHeader in lArgs[1:]):
                # we have neither index nor header
                bHeader = False
                bIndex = False
            else:
                # raise an error
                raise ValueError('the second or third argument of ' + sOption + ' is incorrect. Did you mean ' + sNoLabelIndex + ' ' + sNoLabelHeader + '?')
        else:
            # if we have neither 1, 2, nor 3 input args: raise an errror
            raise SyntaxError(sOption + ' requires a filename and can have the optional arguments ' + sNoLabelIndex + ' and/or ' + sNoLabelHeader + '.')
        # assign False/None as column/row where indices/headers are found
        iIndCol = False
        iHeadRow = None
        # if either available, set to first column resp. row
        if bIndex:
            iIndCol = 0
        if bHeader:
            iHeadRow = 0
        # check if the data file is missing:
        if not(os.path.isfile(lArgs[0])):
            # if so raise an error
            raise TypeError('specified file for the ' + sData + ' does not exist')
        self.logger.info('Reading {S} file {f}'.format(S=sData,f=lArgs[0]))
        # read the file
        dfBin = pd.read_csv(lArgs[0], header = iHeadRow, sep = None, engine = 'python', index_col = iIndCol)
        # set missings to NaN
        dfBin = dfBin.apply(pd.to_numeric, errors='coerce')
        # get the number of traits and factors/covariates involved and report
        (iT,iK) = dfBin.shape
        self.logger.info('Found a ' + sData + ' on ' + str(iT) + ' phenotypes and ' + str(iK) + ' ' + sDescr + '.')
        # if no index provided
        if not(bIndex):
            dfBin.index = ['phenotype ' + str(x) for x in range(0,iT)]
        # if no header provided
        if not(bHeader):
            dfBin.columns = [sDescrShort + ' ' + str(x) for x in range(0,iK)]
        # store dataframe as appropriate attribute of MgremlReader
        if sType == MgremlReader.sCov:
            self.dfBinXY = dfBin
        elif sType == MgremlReader.sGen:
            if bNull:
                self.dfGenBinFY0 = dfBin
            else:
                self.dfGenBinFY = dfBin
        else:
            if bNull:
                self.dfEnvBinFY0 = dfBin
            else:
                self.dfEnvBinFY = dfBin
        # report reading data is done
        self.logger.info('Finished reading in the ' + sData)
        self.logger.info('Current memory usage is ' + str(int((self.process.memory_info().rss)/(1024**2))) + 'MB')
    
    def ReadData(self, sType):
        # figure out what type of input data we have
        # and set strings appropriately
        if sType == MgremlReader.sPhe:
            sNoLabel = 'nolabelpheno'
            sData = 'phenotype'
            sOption = '--pheno'
            lArgs = self.args.pheno
            bTwoOnly = self.bMediation
        elif sType == MgremlReader.sCov:
            sNoLabel = 'nolabelcovar'
            sData = 'covariate'
            sOption = '--covar'
            lArgs = self.args.covar
            bTwoOnly = False
        else:
            raise TypeError('MgremlReader.ReadData() can only be used to read a phenotype file or covariate file.')
        # if no. of input args. is one
        if len(lArgs) == 1:
            # we have a header
            bHeader = True
        elif len(lArgs)==2: # if it's 2
            # but the 2nd arg. != no-label string
            if lArgs[1] != sNoLabel:
                # raise an error
                raise ValueError('the second argument of ' + sOption + ' is incorrect. Did you mean ' + sNoLabel + '?')
            # otherwise: we do not have a header
            bHeader=False
        else:
            # if we have neither 1 nor 2 input args: raise an errror
            raise SyntaxError(sOption + ' requires a filename and can have the optional argument ' + sNoLabel + '.')
        # prepare string to let user know whether header has been set, yes or no
        if bHeader:
            sInfoString1 = 'Headers specified in ' + sData + ' file'
            sInfoString2 = 'Using headers to set labels'
            iHeader = 0
        else:
            sInfoString1 = 'No headers specified in ' + sData + ' file'
            sInfoString2 = 'Your ' + sData + 's will be labeled as ' + sData + ' 0, ' + sData + ' 1, and so on.'
            iHeader = None
        # check if the data file is missing:
        if not(os.path.isfile(lArgs[0])):
            # if so raise an error
            raise TypeError('specified ' + sData + ' file does not exist')
        # report whether header has been set, yes or no
        self.logger.info(sInfoString1)
        self.logger.info(sInfoString2)
        # indicate data is being read
        self.logger.info('Reading ' + sData + ' file {f}'.format(f=lArgs[0]))
        # read data
        dfData = pd.read_csv(lArgs[0], index_col = MgremlReader.lIndFID_IID, header = iHeader, sep = None, engine='python')
        # force missings to NaN
        dfData.apply(pd.to_numeric, errors='coerce')
        # get number of observations and traits, and report
        (iN,iT) = dfData.shape
        self.logger.info('The ' + sData + ' file contains data on {N} individuals and {T} '.format(N=iN,T=iT) + sData + 's')
        # if we only conisder first two phenotypes because of mediation analysis
        if bTwoOnly:
            # do so, and report
            self.logger.info('Selecting only first two '.format(N=iN,T=iT) + sData + 's for mediation analysis')
            dfData=dfData.iloc[:,0:2]
            (iN,iT) = dfData.shape
            self.logger.info('We now have data on {N} individuals and {T} '.format(N=iN,T=iT) + sData + 's')
        # find number of NaNs and report
        iMiss = dfData.isnull().sum().sum()
        self.logger.info('Encountered {M} missings'.format(M=iMiss))
        # if no header has been specified
        if not(bHeader):
            # set labels of the multiindex to FID, IID
            dfData.index.names = MgremlReader.lLabelsFID_IID
            # set labels of trait/covariates to phenotype/covariate 1,
            # phenotype/covariate 2, etc.
            dfData.columns = [sData + ' ' + str(x) for x in range(0,iT)]
        # store data as appropriate attribute of MgremlReader instance
        if sType == MgremlReader.sPhe:
            self.dfY = dfData
        else:
            self.dfX = dfData
        # report reading data is done
        self.logger.info('Finished reading in ' + sData + ' data')
        self.logger.info('Current memory usage is ' + str(int((self.process.memory_info().rss)/(1024**2))) + 'MB')
    
    def DetermineIfCovsAreGiven(self):
        # assert whether we have regular covariates
        self.bCovs = isinstance(self.dfX, pd.DataFrame)
        self.bSameCovs = not(isinstance(self.dfBinXY, pd.DataFrame))
        # if we have no covariates, yet different covariates have been specified
        if not(self.bCovs) and not(self.bSameCovs):
            raise SyntaxError('you have specified different covariates to apply to different traits, while not providing any data on those covariates using --covar')
        # for final assertion if covariates are present, also consider presence of intercept
        self.bCovs = self.bCovs or self.bIntercept
    
    def ReadGRM(self):
        # set filenames for reading binary GRM
        BinFileName = self.args.grm + ".grm.bin"
        IDFileName = self.args.grm + ".grm.id"
        # check if one or more files are missing:
        if not(os.path.isfile(BinFileName)) or not(os.path.isfile(IDFileName)):
            raise TypeError('specified set of GRM files either incomplete or non-existent.')
        # read IDs and sample size
        ids = pd.read_csv(IDFileName, sep = None, engine='python', header = None)
        ids.columns= MgremlReader.lLabelsFID_IID
        arrays=[ids['FID'],ids['IID']]
        tuples = list(zip(*arrays))
        index = pd.MultiIndex.from_tuples(tuples, names=MgremlReader.lLabelsFID_IID)
        iN = len(ids.index)
        # raise error when trying to control for too many PCs
        if iN < (self.iDropLeadPCs+self.iDropTrailPCs):
            raise ValueError('you cannot use --adjust-pcs to adjust for more eigenvectors than there are individuals in your GRM')
        self.logger.info('{iN} individuals in {f}'.format(iN=iN, f=IDFileName))
        self.logger.info('Reading GRM files {f}.*'.format(f=self.args.grm))
        self.logger.info('This may take some time...')
        # read GRM from bin file
        dt = np.dtype('f4') # Relatedness is stored as a float of size 4 in the binary file
        vA = np.fromfile(BinFileName, dtype = dt)
        # create GRM as 2D-array
        mA = np.zeros((iN,iN))
        # set counter for elements of GRM read thus far
        k = 0
        # for each row
        for i in range(0,iN):
            # for each column
            for j in range(0,i+1):
                dThisVal = vA[k]
                mA[i,j] = dThisVal
                mA[j,i] = dThisVal
                k += 1
        # convert to pd dataframe        
        self.dfA = pd.DataFrame(mA, columns=index, index=index)
        # print update
        self.logger.info('Finished reading in GRM')
        self.logger.info('Current memory usage is ' + str(int((self.process.memory_info().rss)/(1024**2))) + 'MB\n')
    
    def AddIntercept(self):
        if self.bIntercept:
            if isinstance(self.dfX, pd.DataFrame):
                self.dfX.insert(MgremlReader.iIntPos, MgremlReader.sIntLab, MgremlReader.iIntVal)
            else:
                self.dfX = pd.DataFrame(MgremlReader.iIntVal, index=self.dfY.index, columns=[MgremlReader.sIntLab])
            if not(self.bSameCovs):
                self.dfBinXY.insert(MgremlReader.iIntPos, MgremlReader.sIntLab, MgremlReader.iIntVal)
    
    def CleanSpecificationCovariates(self):
        # if we have covariates, and different covariates apply to different traits
        if self.bCovs and not(self.bSameCovs):
            # print update
            self.logger.info('INSPECTING YOUR COVARIATE MODEL')
            # get indices phenotypes from dfY and dfBinXY
            indY_Y  = self.dfY.columns
            indXY_Y = self.dfBinXY.index
            # get labels of covariates from dfX and dfBinXY
            indX_X  = self.dfX.columns
            indXY_X = self.dfBinXY.columns
            # print update
            self.logger.info('You specified which covariates apply to ' + str(self.dfBinXY.shape[0]) + ' phenotypes')
            self.logger.info('There are ' + str(self.dfY.shape[1]) + ' phenotypes in your data')
            # all phenotypes in dfY should refer to phenotypes in dfBinXY; abort if problematic
            if not(indY_Y.isin(indXY_Y).all()):
                raise ValueError('there is at least one phenotype for which you have not specified which covariates apply to it') 
            # all covariates in dfBinXY should refer to covariates dfX; abort if problematic
            if not(indXY_X.isin(indX_X).all()):
                raise ValueError('there is at least one phenotype-specific covariate for which you have not supplied the underlying data') 
            if (self.dfBinXY.shape[0] - self.dfY.shape[1]) != 0:
                self.logger.info('The ' + str(self.dfBinXY.shape[0] - self.dfY.shape[1]) + ' redundant phenotypes will be removed from your covariate model')
            # eliminate phenotypes from dfBinXY that are not in dfY
            self.dfBinXY = self.dfBinXY.loc[indY_Y]
            self.logger.info('You specified for ' + str(self.dfBinXY.shape[1]) + ' covariates to which phenotypes they apply')
            self.logger.info('There are ' + str(self.dfX.shape[1]) + ' covariates in your data')
            if (self.dfX.shape[1] - self.dfBinXY.shape[1]) != 0:
                self.logger.info('The ' + str(self.dfX.shape[1] - self.dfBinXY.shape[1]) + ' redundant covariates will be removed from your covariate data')
            # eliminate covariates from dfX that are not used according to dfBinXY
            self.dfX = self.dfX[indXY_X]
            # if dfBinXY is not binary: abort
            if ((self.dfBinXY==0).sum().sum() + (self.dfBinXY==1).sum().sum()) != (self.dfBinXY.shape[0]*self.dfBinXY.shape[1]):
                raise ValueError('your model indicating which covariate affects which phenotype does not only comprise zeros and ones')
            # if dfBinXY is all ones: print warning and drop
            if ((self.dfBinXY==1).sum().sum()) == (self.dfBinXY.shape[0]*self.dfBinXY.shape[1]):
                self.logger.warning('Warning: your model indicating which covariate affects which phenotype now comprises only ones')
                self.logger.warning('Assuming all covariates apply to all traits.')
                self.dfBinXY = None
                self.bSameCovs = True
    
    def CheckDuplicatesAndRank(self):
        # print update
        self.logger.info('CHECKING FOR DUPLICATES')
        # if duplicates in index or columns of dfA
        if (self.dfA.index.duplicated().sum() > 0) or (self.dfA.columns.duplicated().sum() > 0):
            raise ValueError('you have individuals with the same FID-IID combination in your GRM')
        # if duplicates in columns of dfY
        if self.dfY.columns.duplicated().sum() > 0:
            raise ValueError('you have phenotypes with duplicate labels in your data')
        # if duplicates in index of dfY
        if self.dfY.index.duplicated().sum() > 0:
            raise ValueError('you have individuals with the same FID-IID combination in your phenotype data')
        # if we have covariates
        if self.bCovs:
            # if duplicates in columns of dfX
            if self.dfX.columns.duplicated().sum() > 0:
                raise ValueError('you have covariates with duplicate labels in your data')
            # if duplicates in index of dfX
            if self.dfX.index.duplicated().sum() > 0:
                raise ValueError('you have covariates with the same FID-IID combination in your covariate data')
            # if we have different covariates across traits
            if not(self.bSameCovs):
                # if duplicates in columns of dfBinXY
                if self.dfBinXY.columns.duplicated().sum() > 0:
                    raise ValueError('in your specification which covariates apply to which phenotypes, you have specified duplicate covariates')
                # if duplicates in index of dfBinXY
                if self.dfBinXY.index.duplicated().sum() > 0:
                    raise ValueError('in your specification which covariates apply to which phenotypes, you have specified duplicate phenotypes')
            # print update
            self.logger.info('CHECKING FOR MULTICOLLINEARITY IN YOUR FIXED-EFFECT COVARIATES')
            # get matrix of covariates and set missings to zero
            mX = np.array(self.dfX) # note: changes to np.array(df) do not affect df
            mX[np.isnan(mX)] = 0
            # compute the rank of the matrix of covariates
            dRankX = np.linalg.matrix_rank(mX)
            # count the number of covariates
            iK = mX.shape[1]
            # if rank is below the number of covariates
            if dRankX < iK:
                raise ValueError('your matrix of covariates does not have full rank, i.e. there is perfect multicollinearity')
    
    def FindOverlapAndSort(self):
        # print update
        self.logger.info('JOINING YOUR DATA')
        # construct empty list of dataframes, for multi index of each relevant df
        lIDs = list()
        # print counts
        self.logger.info('There are ' + str(self.dfY.shape[0]) + ' individuals in your phenotype data')
        self.logger.info('There are ' + str(self.dfA.shape[0]) + ' individuals in your GRM')
        # get FID-IID combos for pheno and GRM and append to lIDs
        lIDs.append(pd.DataFrame(data=None, columns=None, index=self.dfY.index))
        lIDs.append(pd.DataFrame(data=None, columns=None, index=self.dfA.index))
        # if we have covariates
        if self.bCovs:
            # print count
            self.logger.info('There are ' + str(self.dfX.shape[0]) + ' individuals in your covariate data')
            # get FID-IID combos for covs and append to lIDs
            lIDs.append(pd.DataFrame(data=None, columns=None, index=self.dfX.index))
        # find the FID-IID combos that are present in all relevant dataframes
        miIDs = pd.concat(lIDs, axis=1, join='inner').index
        # print update
        self.logger.info('Keeping only the ' + str(len(miIDs)) + ' overlapping individuals')
        # select and order phenotypes and GRM by those individuals
        self.dfY = self.dfY.loc[miIDs]
        self.dfA = self.dfA.loc[miIDs,miIDs]
        # double-check if everything is now lined up
        if not(all(self.dfY.index == self.dfA.index)) or not(all(self.dfY.index == self.dfA.columns)):
            raise ValueError('the FID-IID combinations in your GRM cannot be lined up properly with those in your phenotype data')
        # if we have covariates
        if self.bCovs:
            # select and order covariates by those individuals
            self.dfX = self.dfX.loc[miIDs]
            # double-check if everything is now lined up
            if not(all(self.dfY.index == self.dfX.index)):
                raise ValueError('the FID-IID combinations in your covariates cannot be lined up properly with those in your phenotype data')
            # if we have different covariates across traits
            if not(self.bSameCovs):
                # get labels of phenotypes and covariates
                # according to order in dfY and dfX resp.
                indY_Y = self.dfY.columns
                indX_X = self.dfX.columns
                # put dfBinXY in the same order
                self.dfBinXY = (self.dfBinXY.loc[indY_Y])[indX_X]
                # double-check if everything is now lined up
                if not(all(self.dfBinXY.index == indY_Y)):
                    raise ValueError('your specification which covariates applies to which phenotype cannot be properly lined up with the phenotypes')
                if not(all(self.dfBinXY.columns == indX_X)):
                    raise ValueError('your specification which covariates applies to which phenotype cannot be properly lined up with the covariates')
    
    def DropMissings(self):
        # print update
        self.logger.info('DROPPING PROBLEMATIC INDIVIDUALS BECAUSE OF MISSING PHENOTYPES AND/OR COVARIATES')
        self.logger.info('Current memory usage is ' + str(int((self.process.memory_info().rss)/(1024**2))) + 'MB')
        # if we have covariates
        if self.bCovs:
            iCount = 0
            # and all covariates apply to all traits
            if self.bSameCovs:
                tDrop = tqdm(total=self.dfY.shape[0])
                # and for a given observations
                for i in self.dfY.index:
                    tDrop.update(1)
                    # if at least one covariate is  missing
                    if (self.dfX.loc[i].isnull().sum() > 0) or (self.dfX.loc[i].isna().sum() > 0):
                        # set all phenotypes for that individual to missing
                        self.dfY.loc[i] = None
                        iCount += 1
                tDrop.close()
                self.logger.info('Found ' + str(iCount) + ' individuals with missing data on one or more covariates')
                self.logger.info('Setting all phenotypes to missing for those individuals')
            else: # if not all covariates apply to all traits
                tDrop = tqdm(total=self.dfY.shape[0])
                # for a given observation
                for i in self.dfY.index:
                    tDrop.update(1)
                    # find indices of missing covariates
                    vIndMissingCovs = np.array(np.where(self.dfX.loc[i].isnull() | self.dfX.loc[i].isna())).ravel()
                    if len(vIndMissingCovs) > 0:
                        iCount += 1
                    # for each missing covariate
                    for j in vIndMissingCovs:
                        # find phenotypes affected by the missing covariate
                        vIndPhenotypes = np.array(np.where(self.dfBinXY.iloc[:,j] == 1)).ravel()
                        # for each of those phenotypes
                        for k in vIndPhenotypes:
                            # set phenotypic value to missing
                            self.dfY.loc[i].iloc[k] = None
                tDrop.close()
                self.logger.info('Found ' + str(iCount) + ' individuals with missing data on one or more covariates')
                self.logger.info('Setting all phenotypes affected by missing covariates to missing for those individiuals')
        # count the number of traits and observations in dfY
        iT = self.dfY.shape[1]
        iN = self.dfY.shape[0]
        # count number of observations to be dropped
        if self.bDropMissings:
            # i.e. with any pheno missing
            iM = ((self.dfY.isnull() | self.dfY.isna()).sum(axis=1) >= 1).sum()
            self.logger.info('Dropping ' + str(iM) + ' out of ' + str(iN) + ' individuals from data for whom at least one phenotype is now missing')
            # keep only observations that have no missing phenotypes
            self.dfY = self.dfY[((self.dfY.isnull() | self.dfY.isna()).sum(axis=1) < 1)]
        else:
            # i.e. with all pheno missing
            iM = ((self.dfY.isnull() | self.dfY.isna()).sum(axis=1) == iT).sum()
            self.logger.info('Dropping ' + str(iM) + ' out of ' + str(iN) + ' individuals from data for whom all phenotypes are now missing')
            # keep only observations that have at least one pheno non-missing
            self.dfY = self.dfY[((self.dfY.isnull() | self.dfY.isna()).sum(axis=1) < iT)]
        # get list of individuals that remain
        miIDs = self.dfY.index
        # keep only those individuals in GRM
        self.dfA = self.dfA.loc[miIDs,miIDs]
        # if we have covariates
        if self.bCovs:
            # keep only those individuals in data on covariates
            self.dfX = self.dfX.loc[miIDs]
    
    def PruneByRelatedness(self):
        # apply relatedness pruning if desired
        if self.bRelCutoff:
            # print update
            self.logger.info('APPLYING RELATEDNESS CUTOFF')
            self.logger.info('Current memory usage is ' + str(int((self.process.memory_info().rss)/(1024**2))) + 'MB')
            self.logger.info('Removing individuals such that there is no relatedness in excess of ' + str(self.dRelCutoff) + ' in GRM')
            # set random-number generator
            rng = np.random.default_rng(MgremlReader.iSeed)
            # create binary matrix == 1 when in excess of cutoff
            mB = ((self.dfA > self.dRelCutoff).values).astype(int)
            # count no. of observations
            iN = mB.shape[0]
            # set vector of indices we start with
            vIDstart = np.arange(iN)
            # set diagonal elements binary marix to zero
            for i in range(0,iN):
                mB[i,i] = 0
            # count missingness per observation
            vMiss = np.array(((self.dfY.isnull() | self.dfY.isna()).sum(axis=1))).ravel()
            # count no. of elements per row in excess of threshold
            vCount = mB.sum(axis=0)
            # if there is any unwanted relatedness at all
            if vCount.sum() > 0:
                # temporarily keep only rows and cols with 1 entry > threshold
                # i.e. keep the singletons
                vSelect = (vCount == 1)
                vID = vIDstart[vSelect]
                mB = mB[vSelect,:][:,vSelect]
                # count no. of elements per row in excess of threshold for singletons
                vCount = mB.sum(axis=0)
                # if there is any unwanted relatedness amongst singletons
                if vCount.sum() > 0:
                    # temporarily keep only rows and cols with 1 entry > threshold
                    # i.e. keep the singletons with singletons, and not
                    # singletons with multiples
                    vSelect = (vCount == 1)
                    vID = vID[vSelect]
                    mB = mB[vSelect,:][:,vSelect]
                    # find upper matrix of subset of rows and columns
                    mBU = np.triu(mB)
                    # find cols and rows where mBU indicates relatedness > threshold
                    (vSelectRow,vSelectCol) = np.where(mBU)
                    # relate them back to original IDs in vID
                    mPairIDs = np.stack((vID[vSelectRow],vID[vSelectCol]))
                    # for those IDs, use vMiss to determine missingness
                    vMissID0 = vMiss[mPairIDs[0,:]]
                    vMissID1 = vMiss[mPairIDs[1,:]]
                    # for each pair: choose whether to drop col or row
                    vDrop1 = (vMissID0<vMissID1).astype(int) + ((vMissID0 == vMissID1)*(rng.uniform(size=vMissID0.shape)>0.5)).astype(int)
                    vDrop0 = 1 - vDrop1
                    vDrop = np.hstack((mPairIDs[0,vDrop0>0],mPairIDs[1,vDrop1>0]))
                    # take full set of IDs and remove the ones to drop
                    vKeepIDs = np.sort(np.array(list(set(vIDstart) - set(vDrop))))
                    # select corresponding subset of rows and cols from DataFrame
                    self.dfA = self.dfA.iloc[vKeepIDs,vKeepIDs]
                    # count no. of dropped rows and columns
                    iDropped = vDrop.shape[0]
                else:
                    # if no relatedness amongst singletons: nothing dropped in 1st pass
                    iDropped = 0
                # print no. of observations dropped in 1st pass
                self.logger.info('First pass: dropped ' + str(iDropped) + ' individuals with only one relatedness value in excess of ' + str(self.dRelCutoff))
                # get data from DataFrame
                mA = self.dfA.values
                iN = mA.shape[0]
                vID = np.arange(iN)
                # create row and column indices
                vR = np.array(repmat(vID,iN,1).ravel())
                vC = np.array(repmat(vID,iN,1).T.ravel())
                # find all entries where row index < column index
                (vIndRltC,) = np.where(vR<vC)
                # keep only those entries
                vR = vR[vIndRltC]
                vC = vC[vIndRltC]
                # get values
                vV = mA[vR,vC]
                # find all entries with value in excess of cutoff
                (vIndVgtTau,) = np.where(vV>self.dRelCutoff)
                # if there is any remaining unwanted relatedness at all
                if len(vIndVgtTau) > 0:
                    # keep only those entries
                    vR = vR[vIndVgtTau]
                    vC = vC[vIndVgtTau]
                    vV = vV[vIndVgtTau]    
                    # get matrix of index pairs with 1 row per pair with relatedness > dRelCutoff
                    mRelated = np.c_[vR,vC]
                    # get all unique IDs appearing in mRelated and turn into list
                    vUnique = np.unique(np.array(mRelated.ravel()))
                    # create a graph
                    myGraph = nx.Graph()
                    # add an edge for each related pair
                    myGraph.add_edges_from(tuple(map(tuple, mRelated)))
                    # use Boppana & Halldorsso 1992 algorithm: keeping large subset with no edges
                    lKeepIDs = nx.maximal_independent_set(myGraph)
                    # find IDs to drop
                    setDropIDs = set(vUnique) - set(lKeepIDs)
                    # take full set of IDs and remove the ones to drop
                    vFinalKeepIDs = np.sort(np.array(list(set(vID) - setDropIDs)))
                    # select corresponding subset of rows and cols from DataFrame
                    self.dfA = self.dfA.iloc[vFinalKeepIDs,vFinalKeepIDs]
                    # print no. of observations dropped in 2nd pass
                    self.logger.info('Second pass: dropped ' + str(len(setDropIDs)) + ' individuals with relatedness values in excess of ' + str(self.dRelCutoff))
                else:
                    self.logger.info('Second pass not required, as there are no individuals with relatedness in excess of ' + str(self.dRelCutoff) + ' left after first pass')
                self.logger.info('Relatedness cutoff has been applied')
            else:
                self.logger.info('No cutoff has been applied, as there are no individuals with relatedness in excess of ' + str(self.dRelCutoff))
            self.logger.info('Remaining sample size is ' + str(self.dfA.shape[0]))
            # get index of observations to keep
            miIDs = self.dfA.index
            # keep appropriate parts of dfY
            self.dfY = self.dfY.loc[miIDs]
            # if covariates specified
            if self.bCovs:
                # keep appropriate parts of dfX
                self.dfX = self.dfX.loc[miIDs]
    
    def CreateDummies(self):
        # if there are any missings at all
        if ((self.dfY.isnull() | self.dfY.isna()).sum()[0] > 0):
            self.logger.info('CREATING PHENOTYPE-SPECIFIC DUMMY VARIABLES TO CONTROL FOR REMAINING MISSINGNESS')
            self.logger.info('Current memory usage is ' + str(int((self.process.memory_info().rss)/(1024**2))) + 'MB')
            # if there are no covariates yet
            if not(self.bCovs):
                # initialise dfX and dfBinXY
                self.dfX = pd.DataFrame(data=None, index=self.dfY.index)
                self.dfBinXY = pd.DataFrame(data=None, index=self.dfY.columns)
                # set covariates being present and not the same across traits
                self.bCovs = True
                self.bSameCovs = False
            # if there are covariates
            else:
                # but same covariates for all traits so far
                if self.bSameCovs:
                    # initialise dfBinXY as bunch of ones
                    self.dfBinXY = pd.DataFrame(data=1, index=self.dfY.columns, columns=self.dfX.columns)
                    # set covariates to being not the same across traits
                    self.bSameCovs = False
            iCount = 0
            # for each trait
            for t in self.dfY.columns:
                # get all observations with missing values
                miIDs = self.dfY.loc[self.dfY[t].isnull() | self.dfY[t].isna(),t].index
                # for each observation with missing value
                for i in miIDs:
                    # set missing phenotype to zero
                    self.dfY.loc[i,t] = 0
                    # construct label for new dummy variable
                    lLabel = ['dummy_trait_' + str(t) + '_obs_' + str(i)]
                    # create dummy variable
                    dfXadd = pd.DataFrame(data=0, index=self.dfX.index, columns=lLabel)
                    dfXadd.loc[i] = 1
                    # append to existing set of covariates
                    self.dfX = pd.concat([self.dfX, dfXadd], axis=1, join='inner')
                    # create corresponding entry for dfBinXY
                    dfBinXYadd = pd.DataFrame(data=0, index=self.dfY.columns, columns=lLabel)
                    dfBinXYadd.loc[t] = 1
                    # append to existing specification which covs apply to which phenos
                    self.dfBinXY = pd.concat([self.dfBinXY, dfBinXYadd], axis=1, join='inner')
                    iCount += 1
            # replace missing in dfX by 0
            self.dfX = self.dfX.fillna(0)
            self.logger.info('Added ' + str(iCount) + ' phenotype-specific dummies to your covariate model')
            if iCount*self.dfY.shape[1] > MgremlData.iManyDummies:
                self.logger.warning('This is a large number of phenotype-specific covariates, given you have ' + str(self.dfY.shape[1]) + ' traits in your data')
                self.logger.warning(str(iCount*self.dfY.shape[1]) + ' additional fixed-effect covariates implied, of which ' + str(iCount*self.dfY.shape[1] - iCount) + ' are set to zero')
                self.logger.warning('CPU time of MGREML may increase dramatically')
                self.logger.warning('Consider running MGREML on a subset of your data with a much lower degree of missingness')
    
    def SetConstrainedModels(self):
        # for standard constrained models: set binary matrices accordingly
        if self.bPerfectRhoG:
            self.dfGenBinFY = MgremlReader.SetPerfectRho(self.lPhenos)
        if self.bNoRhoG:
            self.dfGenBinFY = MgremlReader.SetNoRho(self.lPhenos)
        if self.bNoRhoE:
            self.dfEnvBinFY = MgremlReader.SetNoRho(self.lPhenos)
        if self.bNoVarG:
            self.dfGenBinFY = MgremlReader.SetNoVar(self.lPhenos)
        # if we have a nested model
        if self.bNested:
            # for standard constrained models: set binary matrices accordingly
            if self.bPerfectRhoG0:
                self.dfGenBinFY0 = MgremlReader.SetPerfectRho(self.lPhenos)
            if self.bNoRhoG0:
                self.dfGenBinFY0 = MgremlReader.SetNoRho(self.lPhenos)
            if self.bNoRhoE0:
                self.dfEnvBinFY0 = MgremlReader.SetNoRho(self.lPhenos)
            if self.bNoVarG0:
                self.dfGenBinFY0 = MgremlReader.SetNoVar(self.lPhenos)
    
    def FinaliseData(self):
        # print update
        self.logger.info('FINALISING DATA BEFORE MGREML ANALYSIS')
        self.logger.info('Current memory usage is ' + str(int((self.process.memory_info().rss)/(1024**2))) + 'MB')
        self.logger.info('Computing eigenvalue decomposition of the GRM')
        self.logger.info('WARNING! This may take a long time and double current memory usage...')
        # compute EVD of GRM
        (vD,mP) = np.linalg.eigh(self.dfA.values)
        self.logger.info('Eigenvalue decomposition completed')
        self.logger.info('Current memory usage is ' + str(int((self.process.memory_info().rss)/(1024**2))) + 'MB')
        self.logger.info('Dropping GRM from memory')
        # free up memory
        self.dfA = None
        self.logger.info('Current memory usage is ' + str(int((self.process.memory_info().rss)/(1024**2))) + 'MB')
        # if we have covariates and not same covariates across traits
        if self.bCovs and not(self.bSameCovs):
            # convert appropriate dataframe to numpy array
            self.mBinXY = (self.dfBinXY.values).astype(int)
            # free up memory
            self.dfBinXY = None
        # print update
        self.logger.info('Applying the canonical transformation to your data')
        self.logger.info('Sample size prior to the canonical transformation is ' + str(self.dfY.shape[0]))
        self.logger.info('Adjusting for ' + str(self.iDropLeadPCs) + ' leading eigenvectors to control for population stratification')
        if self.iDropTrailPCs > 0:
            self.logger.info('Adjusting for ' + str(self.iDropTrailPCs) + ' trailing eigenvectors to improve computational efficiency')
        if self.iDropLeadPCs == 0:
            if self.iDropTrailPCs > 0:
                # ignore trailing columns and leading columns from eigenvector matrix
                mP = mP[:,self.iDropTrailPCs:]
                # ignore trailing and leading values from eigenvalue vector
                self.vD = vD[self.iDropTrailPCs:]
        else:
            if self.iDropTrailPCs == 0:
                # ignore trailing columns and leading columns from eigenvector matrix
                mP = mP[:,:-self.iDropLeadPCs]
                # ignore trailing and leading values from eigenvalue vector
                self.vD = vD[:-self.iDropLeadPCs]
            else:
                # ignore trailing columns and leading columns from eigenvector matrix
                mP = mP[:,self.iDropTrailPCs:-self.iDropLeadPCs]
                # ignore trailing and leading values from eigenvalue vector
                self.vD = vD[self.iDropTrailPCs:-self.iDropLeadPCs]
        # store squared eigenvalues
        self.vDSq = (self.vD)**2
        # apply canonical transform to phenotypes and store that and its transpose
        self.mYT = ((self.dfY.values).T@mP)
        self.mY = self.mYT.T
        # free up memory
        self.dfY = None
        # count sample size and number of phenotypes, and store
        self.iN = self.mY.shape[0]
        self.iT = self.mY.shape[1]
        # if less than three traits and --pairwise used
        if (self.iT < 3) and self.bPairwise:
            # print warning and adjust bPairwise
            self.logger.warning('Warning: --pairwise ignored when analysing less than three traits')
            self.bPairwise = False
        # initialise log|X'X|
        self.dLogDetXTX = 0
        # if we have covariates
        if self.bCovs:
            # apply canonical transform to covariates and store
            self.mXT = ((self.dfX.values).T@mP)
            self.mX = self.mXT.T
            # free up memory
            self.dfX = None
            # count number of covariates, and store
            self.iK = self.mX.shape[1]
            # compute XTX
            mXTX = np.matmul(self.mXT,self.mX)
            # if same covariates across traits
            if self.bSameCovs:
                # compute EVD of X'X
                (vDXTX,mPXTX) = np.linalg.eigh(mXTX)
                mInvDPXTXT = repmat(np.reshape(1/vDXTX,(len(vDXTX),1)),1,len(vDXTX))*(mPXTX.T)
                # if any eigenvalue is too close to zero or negative
                if any(vDXTX < abs(np.finfo(float).eps)):
                    # raise an error with a proper explanation of the likely cause
                    raise ValueError('your covariates are rank deficient after the canonical transformation (i.e. perfectly multicollinear). Likely reason: you specified principal components (PCs) from your genetic data as fixed-effect covariates. MGREML already controls for population stratification in the canonical transformation. Please do not control for PCs manually as well. Rather, use --adjust-pcs INTEGER, to indicate for how many PCs you want to control via the canonical transformation.')
                # compute log|X'X| and store
                self.dLogDetXTX = (self.iT)*np.log(vDXTX).sum()
                # compute OLS residual of Y w.r.t. X
                mR = self.mY - self.mX@(mPXTX@(mInvDPXTXT@(self.mXT@self.mY)))
            # if not same across traits
            else:
                # initialise matrix of OLS residuals of Y w.r.t. X
                mR = np.zeros((self.iN,self.iT))
                # get indices of these covariates in terms of Z matrix
                self.vIndCovs = np.array(np.where(np.array(self.mBinXY).ravel()==1)).ravel()
                # count total number of covariates across traits
                self.iKtotal = self.mBinXY.sum()
                # initialise vector of log|X'X| per trait
                self.vLogDetXTX = np.zeros(self.iT)
                # for each trait
                for it in range(0,self.iT):
                    # get binary vector indicating which covariates affect current trait
                    vBinaryX = self.mBinXY[it,:]
                    # find indices
                    vIndBinaryX = np.array(np.where(vBinaryX==1)).ravel()
                    # compute EVD of X'X
                    (vDXTX,mPXTX) = np.linalg.eigh(mXTX[vIndBinaryX,:][:,vIndBinaryX])
                    mInvDPXTXT = repmat(np.reshape(1/vDXTX,(len(vDXTX),1)),1,len(vDXTX))*(mPXTX.T)
                    # compute OLS residual of Y w.r.t. X
                    mR[:,it] = self.mY[:,it] - (self.mX[:,vIndBinaryX])@(mPXTX@(mInvDPXTXT@((self.mXT[vIndBinaryX,:])@(self.mY[:,it]))))
                    # if any eigenvalue is too close to zero or negative
                    if any(vDXTX < abs(np.finfo(float).eps)):
                        # raise an error with a proper explanation of the likely cause
                        raise ValueError('Your covariates are rank deficient after the canonical transformation (i.e. perfectly multicollinear). Likely reason: you specified principal components (PCs) from your genetic data as fixed-effect covariates. MGREML already controls for population stratification in the canonical transformation. Please do not control for PCs manually as well. Rather, use --adjust-pcs INTEGER, to indicate for how many PCs you want to control via the canonical transformation.')
                    # compute log|X'X| and insert in vector
                    self.vLogDetXTX[it] = np.log(vDXTX).sum()
                # compute log|X'X| and store
                self.dLogDetXTX = self.vLogDetXTX.sum()
            # use residuals to initialise empirical covariance of Y
            self.mCovY = (mR.T@mR)/self.iN
        # otherwise get raw covariance of Y
        else:
            # initialise empirical covariance of Y
            self.mCovY = (self.mYT@self.mY)/self.iN
        # compute phenotypic correlation matrix
        vVarY = np.diag(self.mCovY)
        vNegStdY = np.power(vVarY,-0.5)
        mCorrY = np.outer(vNegStdY,vNegStdY)*self.mCovY
        # compute eigenvalues and condition number of correlation matrix
        vEVRY = np.linalg.eigvalsh(mCorrY)
        dCond = abs(max(vEVRY)/min(vEVRY)) # min(vEVRY) can be negative within numerical precision; thus abs() needed
        # if condition number exceeds threshold
        if dCond > (self.iT*MgremlReader.dCondThreshold):
            # if collinearity ignored
            if self.bIgnoreCollinearity:
                # only print warning
                self.logger.warning('Warning: your phenotype data is highly collinear after initial residualisation with respect to fixed-effect covariates (if any) and the canonical transformation.')
                self.logger.warning('This may lead to poorly identified models.')
            else:
                # else raise error
                raise ValueError('your phenotype data is too collinear after initial residualisation with respect to fixed-effect covariates (if any) and the canonical transformation')
        # if there is at least one trait with no variance
        if (vVarY == 0).sum() > 0:
            raise ValueError('you have specified one or more phenotypes without any variance at all')
        self.logger.info('Sample size after the canonical transformation is ' + str(self.iN)) 
        self.logger.info('Final sample size, N = ' + str(self.iN))
        self.logger.info('Final number of traits, T = ' + str(self.iT))
        if self.bCovs:
            self.logger.info('Final number of unique covariates, k = ' + str(self.iK))
            if self.bSameCovs:
                self.logger.info('Final number of fixed effects, K = ' + str(self.iK*self.iT))
            else:
                self.logger.info('Final number of fixed effects, K = ' + str(self.iKtotal))

class PairwiseMgremlReader(MgremlReader):
    
    def __init__(self, mdData, i, j):
        self.bCopy = True
        # create a copy of the the mgregml reader data
        super().__init__(None,None,None,self.bCopy,mdData)
        # set prefix
        self.sPrefix = self.sPrefix + str(i) + '.' + str(j) + '.'
        # get right phenotypes
        self.mYT = self.mYT[[i,j],:]
        self.mY = self.mYT.T
        self.lPhenos = [self.lPhenos[i], self.lPhenos[j]]
        # get right logL constant
        if self.bSameCovs:
            self.dLogDetXTX = (self.dLogDetXTX)*(2/self.iT)
        else:
            self.vLogDetXTX = self.vLogDetXTX[[i,j]]
            self.dLogDetXTX = self.vLogDetXTX.sum()
        # set number of traits now considered
        self.iT = 2
        # if covs, but not same across traits
        if self.bCovs and not(self.bSameCovs):
            # get right active covariates
            self.mBinXY = self.mBinXY[[i,j],:]
            self.vIndCovs = np.array(np.where(np.array(self.mBinXY).ravel()==1)).ravel()
            self.iKtotal = self.mBinXY.sum()
        if self.dfGenBinFY is not None:
            if self.bPerfectRhoG:
                # get genetic model for active traits; selecting all factors;
                # this perfectly reflects rhoG = 1
                self.dfGenBinFY = self.dfGenBinFY.iloc[[i,j],:]
            elif self.bNoRhoG:
                # get genetic model for active traits; select corresponding factors
                # this perfectly reflects rhoG = 0
                self.dfGenBinFY = (self.dfGenBinFY.iloc[[i,j],:]).iloc[:,[i,j]]
            else:
                raise SyntaxError('--pairwise cannot be combined with genetic model other than fully saturated, perfect genetic correlations, or no genetic correlations')
        if self.dfEnvBinFY is not None:
            if self.bNoRhoE:
                # get environment model for active traits; select corresponding factors
                # this perfectly reflects rhoE = 0
                self.dfEnvBinFY = (self.dfEnvBinFY.iloc[[i,j],:]).iloc[:,[i,j]]
            else:
                raise SyntaxError('--pairwise cannot be combined with environment model other than fully saturated or no environment correlations')
        if self.bNested:
            if self.dfGenBinFY0 is not None:
                if self.bPerfectRhoG0:
                    # get genetic model for active traits; selecting all factors;
                    # this perfectly reflects rhoG = 1
                    self.dfGenBinFY0 = self.dfGenBinFY0.iloc[[i,j],:]
                elif self.bNoRhoG0:
                    # get genetic model for active traits; select corresponding factors
                    # this perfectly reflects rhoG = 0
                    self.dfGenBinFY0 = (self.dfGenBinFY0.iloc[[i,j],:]).iloc[:,[i,j]]
                else:
                    raise SyntaxError('--pairwise cannot be combined with restricted genetic model other than fully saturated, perfect genetic correlations, or no genetic correlations')
            if self.dfEnvBinFY0 is not None:
                if self.bNoRhoE0:
                    # get environment model for active traits; select corresponding factors
                    # this perfectly reflects rhoE = 0
                    self.dfEnvBinFY0 = (self.dfEnvBinFY0.iloc[[i,j],:]).iloc[:,[i,j]]
                else:
                    raise SyntaxError('--pairwise cannot be combined with restricted environment model other than fully saturated or no environment correlations')