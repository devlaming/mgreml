import pandas as pd
import numpy as np
import logging
import os.path
from data import tools

__version__ = '0.01'
MASTHEAD  = "████████████████████████████████████████████████████████████████████████\n"
MASTHEAD += "██                                                                    ██\n"
MASTHEAD += "██    ██      ██                                                      ██\n"
MASTHEAD += "██    ███    ███  ██████   ██████    █████████ ██     ██ ██           ██\n"
MASTHEAD += "██    ████  ████ ██    ██  ██    ██  ██        ███   ███ ██           ██\n"
MASTHEAD += "██    ██ ████ ██ ██        ██    ██  ██        ████ ████ ██           ██\n"
MASTHEAD += "██    ██  ██  ██ ██   ████ ██████    ███████   ██ ███ ██ ██           ██\n"
MASTHEAD += "██    ██      ██ ██    ██  ██   ██   ██        ██  █  ██ ██           ██\n"
MASTHEAD += "██    ██      ██ ██    ██  ██    ██  ██        ██     ██ ██           ██\n"
MASTHEAD += "██    ██      ██  ██████   ██     ██ █████████ ██     ██ █████████    ██\n"
MASTHEAD += "██                                                                    ██\n"
MASTHEAD += "██                                                                    ██\n"
MASTHEAD += "████████████████████████████████████████████████████████████████████████\n"
MASTHEAD += "██                                                                    ██\n"
MASTHEAD += "██  VERSION {V}                                                      ██\n".format(V=__version__)
MASTHEAD += "██                                                                    ██\n"
MASTHEAD += "██  (C) 2020 Ronald de Vlaming, Eric Slob, Philip Jansen,             ██\n"
MASTHEAD += "██  Philipp Koellinger, Patrick Groenen, Niels Rietveld               ██\n"
MASTHEAD += "██                                                                    ██\n"
MASTHEAD += "██  Vrij Universiteit Amsterdam and Erasmus University Rotterdam      ██\n"
MASTHEAD += "██  GNU General Public License v3                                     ██\n"
MASTHEAD += "██                                                                    ██\n"
MASTHEAD += "████████████████████████████████████████████████████████████████████████\n"

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
    
    def __init__(self, parser, logger):
        '''
        MgremlReader is provided with a logger and a parser.
        
        The logger is used to print to console and log-file what's happening.
        
        The parser is used to read out e.g. commands along the following lines:
        python mgreml --grm mygrm --pheno mypheno --out results
        
        An instance of an MgremlReader is used as input
        when initialising instances of other classes involved in
        MGREML estimation; the attributes of the MgremlReader will be 
        used to determine what needs to be done
        '''
        # store logger and parser as attributes of instance
        self.logger = logger
        self.parser = parser
        # initialise arguments and logger
        self.InitialiseArgumentsAndLogger()
        # print welcome screen and given input args
        self.PrintStart()
        # assume we will carry out an analysis
        self.bAnalyse = True
        # if either GRM, pheno or output has not been specified
        # do not analyse anything
        if self.args.grm is None:
            self.logger.warning('Error: no GRM has been specified.')
            self.bAnalyse = False
        if self.args.pheno is None:
            self.logger.warning('Error: no phenotype file has been specified.')
            self.bAnalyse = False            
        if self.args.out is None:
            self.logger.warning('Error: no prefix for output files has been specified.')
            self.bAnalyse = False
        # if we can analyse, read in data
        if self.bAnalyse:
            # read GRM
            self.read_grm()
            # read phenotype file
            self.read_data(MgremlReader.sPhe)
            # if covariate file specified: read
            if self.args.covar is not None:
                self.read_data(MgremlReader.sCov)
            else:
                self.dfX = None
            # if covar model has been specified: read
            if self.args.covar_model is not None:
                self.dfBinXY = self.read_model(component = MgremlReader.sCov)
            else:
                self.dfBinXY = None
            # if genetic model has been specified: read
            if self.args.genetic_model is not None:
                self.dfGenBinFY = self.read_model(component = MgremlReader.sGen)
            else:
                self.dfGenBinFY = None
            # if genetic model under null has been specified: read
            if self.args.restricted_genetic_model is not None:
                self.dfGenBinFY0 = self.read_model(component = MgremlReader.sGen, restricted = True)
            else:
                self.dfGenBinFY0 = None
            # if environment model has been specified: read
            if self.args.environment_model is not None:
                self.dfEnvBinFY = self.read_model(component = MgremlReader.sEnv)
            else:
                self.dfEnvBinFY = None
            # if genetic model under null has been specified: read
            if self.args.restricted_environment_model is not None:
                self.dfEnvBinFY0 = self.read_model(component = MgremlReader.sEnv, restricted = True)
            else:
                self.dfEnvBinFY0 = None
        # self.sPrefix        # string based on --out results
        # self.dfY            # DataFrame based on --pheno mypheno.txt [nolabelpheno]
        # self.dfX            # DataFrame based on --covar mycovar.txt; DEFAULT = None
        # self.dfBinXY        # DataFrame based on --covar-model mycovmodel.txt [nolabelpheno] [nolabelcovar]; DEFAULT = None
        # self.dfGenBinFY     # DataFrame based on --genetic-model mygenmodel.txt [nolabelpheno] [nolabelfactor]; DEFAULT = None
        # self.dfEnvBinFY     # DataFrame based on--environment-model myenvmodel.txt [nolabelpheno] [nolabelfactor]; DEFAULT = None
        # self.bPerfectRhoG   # Boolean based on --rho-genetic 1; DEFAULT = False
        # self.bNoRhoG        # Boolean based on --rho-genetic 0; DEFAULT = False
        # self.bNoRhoE        # Boolean based on --rho-environment 0; DEFAULT = False
        # self.dfGenBinFY0    # DataFrame based on --restricted-genetic-model mygenmodel0.txt [nolabelpheno] [nolabelfactor]; DEFAULT = None
        # self.dfEndBinFY0    # DataFrame based on --restricted-environment-model myenvmodel0.txt [nolabelpheno] [nolabelfactor]; DEFAULT = None
        # self.bPerfectRhoG0  # Boolean based on --restricted-rho-genetic 1; DEFAULT = False
        # self.bNoRhoG0       # Boolean based on --restricted-rho-genetic 0; DEFAULT = False
        # self.bNoRhoE0       # Boolean based on --restricted-rho-environment 0; DEFAULT = False
        # self.bSE            # Boolean based on --no-se; DEFAULT = True (so reporting SEs is the default!)
        # self.bRelCutoff     # Boolean based on --rel-cutoff 0.025; DEFAULT = False
        # self.dRelCutoff     # number based on --rel-cutoff 0.025; DEFAULT = None
        # self.iDropLeadPCs   # integer based on --ignore-pcs 40 [1000]; DEFAULT = 20
        # self.iDropTrailPCs  # integer based on --ignore-pcs 40 [1000]; DEFAULT = 0
        # self.bStoreIter     # Boolean based on --store-iter 50; DEFAULT = False
        # self.iStoreIterFreq # integer based on --store-iter 50; DEFAULT = None
        # self.sInitValsFile  # string based on --reinitialise results.iter.250.pkl
        # self.bAllCoeffs     # Boolean based on --all-coefficients
        # self.bBFGS          # Boolean based on --newton; DEFAULT = True
        # self.dGradTol       # number based on --grad-tol; DEFAULT = 1E-5
    
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
                            help = 'Optional flag naming your covariate file. Possible to add the flags nolabelcovar e.g. --covar mycovar.txt nolabelcovar; not recommended, please label your covariates! File should be comma-, space-, or tab-separated, with one row per individual, with FID and IID as first two fields, followed by a field per covariate. WARNING: do not include principal components from your genetic data as covariates. MGREML automatically controls for 20 principal components using the canonical transformation. If you want to control for more or fewer principal components, use the --ignore-pcs option instead.')
        self.parser.add_argument('--covar-model', metavar = 'mycovmodel.txt [nolabelpheno] [nolabelcovar]', default = None, type = str, nargs = '+',
                            help = 'Optional flag naming the file that specifies which covariates (columns) apply to which phenotypes (rows). Possible to add the flags nolabelpheno and/or nolabelcovar; not recommended, please label your phenotypes and covariates! If no file is specified, all covariates are assumed to apply to all traits.')
        groupGenetic.add_argument('--genetic-model', metavar = 'mygenmodel.txt [nolabelpheno] [nolabelfactor]', default = None, type = str, nargs = '+',
                            help = 'Optional flag naming the file that specifies which genetic factors (columns) affect which phenotypes (rows). Possible to add the flags nolabelpheno and/or nolabelfactor; not recommended, please label your phenotypes and genetic factors!')
        groupGenetic.add_argument('--rho-genetic', choices = [0, 1], default = None, type = int,
                            help = 'Optional flag followed by integer equal to zero or one, forcing all genetic correlations to take on the specified value. This flag cannot be combined with --genetic-model.')
        groupRestrictedGenetic.add_argument('--restricted-genetic-model', metavar = 'mygenmodel0.txt [nolabelpheno] [nolabelfactor]', default = None, type = str, nargs = '+',
                            help = 'Optional flag naming the file that specifies for a restricted model which genetic factors (columns) affect which phenotypes (row). Possible to add the flags nolabelpheno and/or nolabelfactor; not recommended, please label your phenotypes and genetic factors!')
        groupRestrictedGenetic.add_argument('--restricted-rho-genetic', choices = [0, 1], default = None, type = int,
                            help = 'Optional flag followed by integer equal to zero or one, forcing all genetic correlations in the restricted model to take on the specified value. This flag cannot be combined with --restricted-genetic-model.')
        groupEnvironment.add_argument('--environment-model', metavar = 'myenvmodel.txt [nolabelpheno] [nolabelfactor]', default = None, type = str, nargs = '+',
                            help = 'Name of a file that specifies which environment factors (columns) affect which phenotypes (rows). Possible to add the flags nolabelpheno and/or nolabelfactor; not recommended, always name things.')
        groupEnvironment.add_argument('--rho-environment', choices = [0], default = None, type = int,
                            help = 'Optional flag followed by integer equal to zero, forcing all environment correlations to zero. This flag cannot be combined with --environment-model.')
        groupRestrictedEnvironment.add_argument('--restricted-environment-model', metavar = 'myenvmodel0.txt [nolabelpheno] [nolabelfactor]', default = None, type = str, nargs = '+',
                            help = 'Name of a file that specifies for a restricted model which environment factors (columns) affect which phenotypes (row). Possible to add the flags nolabelpheno and/or nolabelfactor; not recommended, always name things.')
        groupRestrictedEnvironment.add_argument('--restricted-rho-environment', choices = [0], default = None, type = int,
                            help = 'Optional flag followed by integer equal to zero, forcing all environment correlations in the restricted model to zero. This flag cannot be combined with --restricted-environment-model.')
        self.parser.add_argument('--ignore-pcs', metavar = '20 [0]', default = None, type = int,
                            help = 'optional flag to specify how many leading principal components (PCs) to ignore (as method to control for population stratification) and how many trailing PCs to ignore (for computational efficiency); if just one non-negative integer is supplied this take as the number of leading PCs to ignore')
        self.parser.add_argument('--rel-cutoff', metavar = '0.025', default = None, type = float,
                            help = 'optional flag followed by a value above which overly related individuals are removed from the GRM using a greedy algorithm')
        self.parser.add_argument('--no-se', action = 'store_true',
                            help = 'optional flag to indicate calculation of standard errors should be skipped (e.g. when only interested in a likelihood-ratio test for a nested model)')
        self.parser.add_argument('--all-coefficients', action = 'store_true', 
                            help = 'optional flag to report the all factor coefficients (i.e. from each factor to each trait) including the sampling covariance of these estimates')
        self.parser.add_argument('--newton', action = 'store_true',
                            help = 'optional flag to perform Newton instead of BFGS; not recommended, unless the model is well-identified and the number of traits is small')
        self.parser.add_argument('--grad-tol', metavar = '1E-5', default = None, type = float,
                            help = 'optional flag to set convergence threshold on the length of the gradient vector per parameter, per observation, different from the default value of 1E-5')
        self.parser.add_argument('--store-iter', metavar = '50', default = None, type = int, 
                            help = 'optional flag to specify every how many iterations you want to store results (e.g. every 50 itertions)')
        self.parser.add_argument('--reinitialise', metavar = 'myestimates.iter.100.pkl', default = None, type = str,  
                            help = 'optional flag to reinitialise MGREML estimation from a .pkl file generated by the --store-iter option; make sure this .pkl file has obtained under the same data and model you now want to continue estimation for')
        self.parser.add_argument('--restricted-reinitialise', metavar = 'myestimates0.iter.100.pkl', default = None, type = str,  
                            help = 'optional flag to reinitialise MGREML estimation from a .pkl file generated by the --store-iter option for a restricted model; make sure this .pkl file has obtained under the same data and model you now want to continue estimation for')
        self.parser.add_argument('--out', metavar = 'myoutput', default = None, type = str,
                            help = 'Prefix of the output files.')
        # try to parse the input arguments
        try:
            # parse the input options
            self.args = self.parser.parse_args()
        except Exception:
            self.logger.error('Wrong input options used: MGREML cannot not parse input options')
            raise ValueError
        # customise the logger e.g. usings args.out
        c_handler = logging.StreamHandler()
        # if no --out option has not been specified
        # use generic name for log-file and so on
        if self.args.out is not None:
            self.sPrefix = self.args.out + '.'
            sFileOut = self.sPrefix + 'log'
        else:
            self.sPrefix = 'output.'
            sFileOut = self.sPrefix + 'log'
        f_handler = logging.FileHandler(sFileOut,'w+',encoding="utf-8")
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
        
    def PrintStart(self):
        # try to print welcome screen
        try:
            defaults = vars(self.parser.parse_args(''))
            opts = vars(self.args)
            non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
            header = MASTHEAD
            header += "\nCall: \n"
            header += 'mgreml \\\n'
            options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults]
            header += '\n'.join(options).replace('True','').replace('False','').replace("', \'", ' ').replace("']", '').replace("['", '')
            header = header[0:-1]+'\n'
            self.logger.info(header)
        except Exception:
            self.logger.error('Wrong input values used: MGREML does not understand input options')
            raise ValueError

    def read_grm(self):
        # set filenames for reading binary GRM
        BinFileName = self.args.grm + ".grm.bin"
        NFileName = self.args.grm + ".grm.N.bin"
        IDFileName = self.args.grm + ".grm.id"
        # check if one or more files are missing:
        if not(os.path.isfile(BinFileName)) or not(os.path.isfile(NFileName)) or not(os.path.isfile(IDFileName)):
            self.logger.error('Specified set of GRM files either incomplete or does not exist')
            raise ValueError
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
        print('Printing first and last few rows of GRM')
        print(self.dfA)
        self.logger.info('Finished reading in GRM')
        
    def read_data(self, sType):
        ''' read data '''
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
            raise ValueError('MgremlReader.read_data() can only be used to read a phenotype file or covariate file')
        if len(lArgs) == 1: # if no. of input args. is one
            bHeader = True # we have a header
        elif len(lArgs)==2: # if it's 2
            if lArgs[1] != sNoLabel: # but the 2nd arg. != no-label string
                # raise an error
                raise ValueError('The second argument of the ' + sOption + ' flag is incorrectly specified, did you mean ' + sNoLabel + '?')
            # otherwise: we do not have a header
            bHeader=False
        else:
            # if we have neither 1 nor 2 input args: raise an errror
            raise ValueError('The ' + sOption + ' flag requires a filename and the optional modifier ' + sNoLabel + '.')
        # prepare string to let user know whether header has been set, yes or no
        if bHeader:
            sInfoString = 'Headers specified in ' + sData + ' file\nUsing headers to set labels'
            iHeader = 0
        else:
            sInfoString = 'No headers specified in ' + sData + ' file\nYour ' + sData + 's will be labeled as ' + sData + ' 0, ' + sData + ' 1, and so on.'
            iHeader = None
        # check if the data file is missing:
        if not(os.path.isfile(lArgs[0])):
            # if so raise an error
            self.logger.error('Specified ' + sData + ' file does not exist')
            raise ValueError
        # report whether header has been set, yes or no
        self.logger.info(sInfoString)
        self.logger.info('Reading ' + sData + ' file {f}'.format(f=lArgs[0]))
        dfData = pd.read_csv(lArgs[0], index_col = MgremlReader.lIndFID_IID, header = iHeader, sep = "\s+|\t+|,",engine='python')
        dfData.apply(pd.to_numeric, errors='coerce')
        (iN,iT) = dfData.shape
        self.logger.info(sData + ' file contains data on {N} individuals and {T} '.format(N=iN,T=iT) + sData + 's')
        #find number of NA  write number of missings
        iMiss = dfData.isnull().sum().sum()
        self.logger.info('Encountered {M} missings'.format(M=iMiss))
        if not(bHeader):
            dfData.index.names = MgremlReader.lLabelsFID_IID
            dfData.columns = [sData + str(x) for x in range(0,iT)]
        print('Printing first and last few rows of ' + sData + ' data')
        print(dfData)
        if sType == MgremlReader.sPhe:
            self.dfY = dfData
        else:
            self.dfX = dfData
        self.logger.info('Finished reading in ' + sData + ' data')
        