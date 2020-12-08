import pandas as pd
import numpy as np
import logging
import time
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
        # store logger, parser and starting time as attributes of instance
        self.logger = logger
        self.parser = parser
        self.dStart = dStart
        # initialise arguments
        self.InitialiseArgsAndLogger()
        # print welcome screen and given input args
        self.PrintStart()
        # get DataFrame for GRM based on --grm mygrm option
        self.read_grm()
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
    
    def InitialiseArgsAndLogger(self):
        #create mutually exclusive groups
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
                            help='Name of the phenotypes file. NB: case insensitive. Possible to add the flags nolabelpheno; not recommended, always name things.')
        self.parser.add_argument('--covar', default=None, type=str,nargs='*',
                            help='Name of the covariates file. NB: case insensitive. Possible to add the flags nolabelcovar; not recommended, always name things.')
        self.parser.add_argument('--covar-model', default=None, type=str,nargs='*',
                            help='Name of a file that specifies which covariates (columns) affect which phenotypes (rows). Possible to add the flags nolabelpheno and/or nolabelcovar; not recommended, always name things.')
        self.parser.add_argument('--rel-cutoff', default=None, type=int,
                            help='Relatedness cut-off to apply to GRM. No cut-off is applied by default.')
        self.parser.add_argument('--no-se',
                            help='Option to indicate no standard errors should be reported (e.g. when only interested in a likelihood-ratio test for a nested model)')
        groupGenetic.add_argument('--genetic-model', default=None, type=str,nargs='*',
                            help='Name of a file that specifies which genetic factors (columns) affect which phenotypes (rows). Possible to add the flags nolabelpheno and/or nolabelfactor; not recommended, always name things.')
        groupGenetic.add_argument('--rho-genetic', default=None, type=int,nargs='1',
                            help='Integer equal to zero or one, forcing all genetic correlations to take on the specified value.')
        groupRestrictedGenetic.add_argument('--restricted-genetic-model', default=None, type=str,nargs='*',
                            help='Name of a file that specifies for a restricted model which genetic factors (columns) affect which phenotypes (row). Possible to add the flags nolabelpheno and/or nolabelfactor; not recommended, always name things.')
        groupRestrictedGenetic.add_argument('--restricted-rho-genetic', default=None, type=int,nargs='1',
                            help='Integer equal to zero or one, forcing all genetic correlations in a restricted model to take on the specified value.')
        groupEnvironment.add_argument('--environment-model', default=None, type=str,nargs='*',
                            help='Name of a file that specifies which environment factors (columns) affect which phenotypes (rows). Possible to add the flags nolabelpheno and/or nolabelfactor; not recommended, always name things.')
        groupEnvironment.add_argument('--rho-environment', default=None, type=int,nargs='1',
                            help='Integer that must equal zero if option is used, forcing all environment correlations to zero.')
        groupRestrictedEnvironment.add_argument('--restricted-rho-environment', default=None, type=int,nargs='1',
                            help='Integer that must equal zero if option is used, forcing all environment correlations in a restricted model to zero.')
        groupRestrictedEnvironment.add_argument('--restricted-environment-model', default=None, type=str,nargs='*',
                            help='Name of a file that specifies for a restricted model which environment factors (columns) affect which phenotypes (row). Possible to add the flags nolabelpheno and/or nolabelfactor; not recommended, always name things.')
        # parse the input options
        self.args = self.parser.parse_args()
        # customise the logger e.g. usings args.out
        c_handler = logging.StreamHandler()
        if self.args.out is not None:
            f_handler = logging.FileHandler(self.args.out + '.log','w+',encoding="utf-8")
        else:
            f_handler = logging.FileHandler('output.log','w+',encoding="utf-8")
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
        try:
            defaults = vars(self.parser.parse_args(''))
            opts = vars(self.args)
            non_defaults = [x for x in opts.keys() if opts[x] != defaults[x]]
            header = MASTHEAD
            header += "\n\nCall: \n"
            header += 'mgreml \\\n'
            options = ['--'+x.replace('_','-')+' '+str(opts[x])+' \\' for x in non_defaults]
            header += '\n'.join(options).replace('True','').replace('False','')
            header = header[0:-1]+'\n'
            self.logger.info(header)
            if self.args.grm is not None:
                dfA=read_grm(args,logger)
            if self.args.pheno is not None:
                dfY=read_pheno(args,logger)
            if self.args.covar is not None:
                dfX=read_covar(args,logger)    
            if self.args.covar_model is not None:
                dfBinXY=read_covar_model(args,logger)
            if self.args.genetic_model is not None: 
                dfGenBinFY=read_model(args,logger,component='genetic')
            if self.args.restricted_genetic_model is not None: 
                dfResGenBinFY=read_model(args,logger,component='genetic',restricted=True)    
            if self.args.environment_model is not None: 
                dfGenEnvFY=read_model(args,logger,component='environment')
            if self.args.restricted_environment_model is not None: 
                dfResEnvBinFY=read_model(args,logger,component='environment',restricted=True)   
            if self.args.rho_genetic is not None:
                dfGenBinFY=create_model(args,logger,component='genetic')
            if len(vars(self.args))==1:
                self.logger.warning('Error: no analysis selected.')
                self.logger.info('mgreml_prepare.py -h describes options')
        except Exception:
            ex_type, ex, tb = sys.exc_info()
            self.logger.error( traceback.format_exc(ex) )
            raise   
        finally:
            self.logger.info('Total time elapsed: {T}'.format(
                T=tools.sec_to_str(round(time.time() - self.dStart, 2))))

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
