import pandas as pd
import numpy as np
import logging
import os.path
from data import tools

__version__ = '0.01'
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
MASTHEAD += "##  (C) 2020 Ronald de Vlaming and Eric Slob                          ##\n"
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
    
    def __init__(self, parser, logger):
        # store logger and parser as attributes of instance
        self.logger = logger
        self.parser = parser
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
            self.logger.info('READING OPTIONS')
            # determine whether we will drop missings
            self.DropMissings()
            # determine whether rel-cutoff has been used
            self.SetRelCutOff()
            # determine how many PCs to drop
            self.SetNumberOfPCs()
            # determine whether SEs are desired
            self.NeedSEs()
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
            self.NeedToReinitialise(bNull = False)
            self.NeedToReinitialise(bNull = True)
            # determine whether any correlations are fixed to zero or one
            # or if genetic variances are fixed to zero
            self.FindFixedRhoVar()
            self.logger.info('READING MODELS')
            # if covar model has been specified: read
            if self.args.covar_model is not None:
                self.ReadModel(MgremlReader.sCov)
            else:
                self.dfBinXY = None
            # if genetic model has been specified: read
            if self.args.genetic_model is not None:
                self.ReadModel(MgremlReader.sGen)
            else:
                self.dfGenBinFY = None
            # if genetic model under null has been specified: read
            if self.args.restricted_genetic_model is not None:
                self.ReadModel(MgremlReader.sGen, MgremlReader.bNull)
            else:
                self.dfGenBinFY0 = None
            # if environment model has been specified: read
            if self.args.environment_model is not None:
                self.ReadModel(MgremlReader.sEnv)
            else:
                self.dfEnvBinFY = None
            # if genetic model under null has been specified: read
            if self.args.restricted_environment_model is not None:
                self.ReadModel(MgremlReader.sEnv, MgremlReader.bNull)
            else:
                self.dfEnvBinFY0 = None
            self.logger.info('READING DATA')
            # read phenotype file
            self.ReadData(MgremlReader.sPhe)
            # if covariate file specified: read
            if self.args.covar is not None:
                self.ReadData(MgremlReader.sCov)
            else:
                self.dfX = None
            # read GRM
            self.ReadGRM()
    
    def InitialiseArgumentsAndLogger(self):
        #create mutually exclusive groups
        groupGenetic = self.parser.add_mutually_exclusive_group()
        groupRestrictedGenetic = self.parser.add_mutually_exclusive_group()
        groupEnvironment = self.parser.add_mutually_exclusive_group()
        groupRestrictedEnvironment = self.parser.add_mutually_exclusive_group()
        #create arguments
        self.parser.add_argument('--grm', metavar = 'mygrm', default = None, type = str,
                            help = 'Prefix of the binary GRM e.g. from GCTA or PLINK.')
        self.parser.add_argument('--pheno', metavar = 'myphen.txt [nolabelpheno]', default = None, type = str, nargs = '+',
                            help = 'Name of your phenotype file. Possible to add the flags nolabelpheno e.g. --pheno mypheno.txt nolabelpheno; not recommended, please label your phenotypes! File should be comma-, space-, or tab-separated, with one row per individual, with FID and IID as first two fields, followed by a field per phenotype.')
        self.parser.add_argument('--covar', metavar = 'mycov.txt [nolabelcovar]', default = None, type = str, nargs = '+',
                            help = 'Optional flag naming your covariate file. Possible to add the flags nolabelcovar e.g. --covar mycovar.txt nolabelcovar; not recommended, please label your covariates! File should be comma-, space-, or tab-separated, with one row per individual, with FID and IID as first two fields, followed by a field per covariate. WARNING: do not include principal components from your genetic data as covariates. MGREML automatically controls for 20 principal components using the canonical transformation. If you want to control for more or fewer principal components, use the --adjust-pcs option instead.')
        self.parser.add_argument('--covar-model', metavar = 'mycovmodel.txt [nolabelpheno] [nolabelcovar]', default = None, type = str, nargs = '+',
                            help = 'Optional flag naming the file that specifies which covariates (columns) apply to which phenotypes (rows). Possible to add the flags nolabelpheno and/or nolabelcovar; not recommended, please label your phenotypes and covariates! If no file is specified, all covariates are assumed to apply to all traits.')
        groupGenetic.add_argument('--genetic-model', metavar = 'mygenmodel.txt [nolabelpheno] [nolabelfactor]', default = None, type = str, nargs = '+',
                            help = 'Optional flag naming the file that specifies which genetic factors (columns) affect which phenotypes (rows). Possible to add the flags nolabelpheno and/or nolabelfactor; not recommended, please label your phenotypes and genetic factors!')
        groupGenetic.add_argument('--rho-genetic', choices = [0, 1], default = None, type = int,
                            help = 'Optional flag followed by integer equal to zero or one, forcing all genetic correlations to take on the specified value. This flag cannot be combined with --genetic-model.')
        groupGenetic.add_argument('--no-var-genetic', action = 'store_true',
                            help = 'Optional flag, forcing all genetic variances to be equal to zero. This flag cannot be combined with --genetic-model and/or --rho-genetic.')
        groupRestrictedGenetic.add_argument('--restricted-genetic-model', metavar = 'mygenmodel0.txt [nolabelpheno] [nolabelfactor]', default = None, type = str, nargs = '+',
                            help = 'Optional flag naming the file that specifies for a restricted model which genetic factors (columns) affect which phenotypes (row). Possible to add the flags nolabelpheno and/or nolabelfactor; not recommended, please label your phenotypes and genetic factors!')
        groupRestrictedGenetic.add_argument('--restricted-rho-genetic', choices = [0, 1], default = None, type = int,
                            help = 'Optional flag followed by integer equal to zero or one, forcing all genetic correlations in the restricted model to take on the specified value. This flag cannot be combined with --restricted-genetic-model.')
        groupRestrictedGenetic.add_argument('--restricted-no-var-genetic', action = 'store_true',
                            help = 'Optional flag, forcing all genetic variances in the restricted model to be equal to zero. This flag cannot be combined with --restricted-genetic-model and/or --restricted-rho-genetic.')
        groupEnvironment.add_argument('--environment-model', metavar = 'myenvmodel.txt [nolabelpheno] [nolabelfactor]', default = None, type = str, nargs = '+',
                            help = 'Name of a file that specifies which environment factors (columns) affect which phenotypes (rows). Possible to add the flags nolabelpheno and/or nolabelfactor; not recommended, always name things.')
        groupEnvironment.add_argument('--rho-environment', choices = [0], default = None, type = int,
                            help = 'Optional flag followed by integer equal to zero, forcing all environment correlations to zero. This flag cannot be combined with --environment-model.')
        groupRestrictedEnvironment.add_argument('--restricted-environment-model', metavar = 'myenvmodel0.txt [nolabelpheno] [nolabelfactor]', default = None, type = str, nargs = '+',
                            help = 'Name of a file that specifies for a restricted model which environment factors (columns) affect which phenotypes (row). Possible to add the flags nolabelpheno and/or nolabelfactor; not recommended, always name things.')
        groupRestrictedEnvironment.add_argument('--restricted-rho-environment', choices = [0], default = None, type = int,
                            help = 'Optional flag followed by integer equal to zero, forcing all environment correlations in the restricted model to zero. This flag cannot be combined with --restricted-environment-model.')
        self.parser.add_argument('--adjust-pcs', metavar = '20 [0]', default = None, type = int, nargs = '+',
                            help = 'optional flag to specify for how many leading principal components (PCs) to adjust (as method to control for population stratification) and for how many trailing PCs to adjust (for computational efficiency); if just one non-negative integer is specified this is taken as the number of leading PCs to adjust for')
        self.parser.add_argument('--drop-missings', action = 'store_true',
                            help = 'optional flag to drop all observations from data with at least one missing phenotype or at least one missing covariate')
        self.parser.add_argument('--rel-cutoff', metavar = '0.025', default = None, type = float,
                            help = 'optional flag followed by a value above which overly related individuals are removed from the GRM using a greedy algorithm')
        self.parser.add_argument('--no-se', action = 'store_true',
                            help = 'optional flag to indicate calculation of standard errors and covariance matrix of estimates should be skipped (e.g. when only interested in a likelihood-ratio test for a nested model)')
        self.parser.add_argument('--factor-coefficients', action = 'store_true', 
                            help = 'optional flag to report all estimated factor coefficients')
        self.parser.add_argument('--variance-components', action = 'store_true', 
                            help = 'optional flag to report all estimated variance components')
        self.parser.add_argument('--newton', action = 'store_true',
                            help = 'optional flag to perform Newton instead of BFGS; not recommended, unless the model is well-identified and the number of traits is small')
        self.parser.add_argument('--grad-tol', metavar = '1E-5', default = None, type = float,
                            help = 'optional flag to set convergence threshold on the length of the gradient vector per parameter, per observation, different from the default value of 1E-5')
        self.parser.add_argument('--store-iter', metavar = '50', default = None, type = int, 
                            help = 'optional flag to specify every how many iterations you want to store results (e.g. every 50 itertions)')
        self.parser.add_argument('--reinitialise', metavar = 'myoutput.estimates.iter.100.bfgs.pkl', default = None, type = str,  
                            help = 'optional flag to reinitialise MGREML estimation for a model and its estimates as specified in a .pkl file generated by the --store-iter option')
        self.parser.add_argument('--restricted-reinitialise', metavar = 'myoutput.estimates0.iter.100.bfgs.pkl', default = None, type = str,  
                            help = 'optional flag to reinitialise MGREML estimation for a restricted model and its estimates as specified in a .pkl file generated by the --store-iter option')
        self.parser.add_argument('--out', metavar = 'myoutput', default = None, type = str,
                            help = 'Prefix of the output files.')
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
            header += "\nCall: \n"
            header += 'mgreml \\\n'
            options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults]
            header += '\n'.join(options).replace('True','').replace('False','').replace("', \'", ' ').replace("']", '').replace("['", '').replace('[', '').replace(']', '').replace(', ', ' ').replace('  ', ' ')
            header = header[0:-1]+'\n'
            self.logger.info(header)
        except Exception:
            raise SyntaxError('you specified incorrect input options.')

    def ReadGRM(self):
        # set filenames for reading binary GRM
        BinFileName = self.args.grm + ".grm.bin"
        NFileName = self.args.grm + ".grm.N.bin"
        IDFileName = self.args.grm + ".grm.id"
        # check if one or more files are missing:
        if not(os.path.isfile(BinFileName)) or not(os.path.isfile(NFileName)) or not(os.path.isfile(IDFileName)):
            raise TypeError('specified set of GRM files either incomplete or non-existent.')
        # read IDs and sample size
        ids = pd.read_csv(IDFileName, sep = '\t', header = None)
        ids.columns= MgremlReader.lLabelsFID_IID
        arrays=[ids['FID'],ids['IID']]
        tuples = list(zip(*arrays))
        index = pd.MultiIndex.from_tuples(tuples, names=MgremlReader.lLabelsFID_IID)
        iN = len(ids.index)
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
        self.logger.info('Finished reading in GRM\n')
        
    def ReadData(self, sType):
        # figure out what type of input data we have
        # and set strings appropriately
        if sType == MgremlReader.sPhe:
            sNoLabel = 'nolabelpheno'
            sData = 'phenotype'
            sOption = '--pheno'
            lArgs = self.args.pheno
        elif sType == MgremlReader.sCov:
            sNoLabel = 'nolabelcovar'
            sData = 'covariate'
            sOption = '--covar'
            lArgs = self.args.covar
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
            self.logger.info('Genetic variance in the main model set to zero.')
            self.bNoVarG = True
        # assess whether genetic variance is zero in the restricted model
        if self.args.restricted_no_var_genetic:
            if self.bReinitialise0:
                raise SyntaxError('--restricted-no-var-genetic cannot be combined with --restricted-reinitialise, as the .pkl file is used to set the restricted model')
            self.logger.info('Genetic variance in the null model set to zero.')
            self.bNoVarG0 = True
    
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
    
    def DropMissings(self):
        # if --drop-missings option used
        if self.args.drop_missings:
            # do so
            self.bDropMissings = True
            self.logger.info('Observations with any missing data will be dropped from analysis')
        else:
            self.bDropMissings = False
            self.logger.info('MGREML will construct phenotype-specific dummies to control for missing data')
    
    def SetGradTol(self):
        if self.args.grad_tol is not None:
            if (self.args.grad_tol <= 0):
                raise ValueError('--grad-tol should be followed by a positive number e.g. 1E-6, 1e-5, or 0.0001')
            self.dGradTol = self.args.grad_tol
        else:
            self.dGradTol = MgremlReader.dGradTol
        self.logger.info('Setting convergence threshold for length of gradient per observation per parameter to ' + str(self.dGradTol))
    
    def NeedSEs(self):
        # if --no-se option used
        if self.args.no_se:
            # report no SEs
            self.bSEs = False
            self.logger.info('Your results will not include standard errors.')
        else:
            self.bSEs = True
            self.logger.info('Your results will include standard errors.')
            
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
            
    def SetRelCutOff(self):
        # if --rel-cutoff used
        if self.args.rel_cutoff is not None:
            if (self.args.rel_cutoff < 0):
                raise ValueError('--rel-cutoff should be followed by non-negative value.')
            else:
                self.bRelCutoff = True
                self.dRelCutoff = self.args.rel_cutoff
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
