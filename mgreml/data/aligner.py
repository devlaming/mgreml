import numpy as np
import pandas as pd
from tqdm import tqdm
pd.options.mode.chained_assignment = None

class MgremlData:

    iManyDummies = 1000
    sF = 'factor '
    
    def __init__(self, mReader):
        # read out MgremlReader
        self.logger = mReader.logger
        dfA = mReader.dfA
        dfY = mReader.dfY
        dfX = mReader.dfX
        dfBinXY = mReader.dfBinXY
        iDropLeadPCs = mReader.iDropLeadPCs
        iDropTrailPCs = mReader.iDropTrailPCs
        self.logger.info('2. CLEANING YOUR DATA')
        # find out if we have covariates
        self.DetermineIfCovsAreGiven(dfX, dfBinXY)
        # if we have covariates, and different covariates apply to different traits
        if self.bCovs and not(self.bSameCovs):
            # clean up the specification of the part of the model on covariates
            (dfX, dfBinXY) = self.CleanSpecificationCovariates(dfY, dfX, dfBinXY)
        # check if there are duplicates and/or rank defficiency
        self.CheckDuplicatesAndRank(dfY, dfA, dfX, dfBinXY)
        # find overlapping individuals and sort data
        (dfY, dfA, dfX, dfBinXY) = self.FindOverlapAndSort(dfY, dfA, dfX, dfBinXY)
        # drop observations where missingness affects all traits
        (dfY, dfA, dfX) = self.DropMissings(dfY, dfA, dfX, dfBinXY)
        # apply relatedness pruning if desired
        if mReader.bRelCutoff:
            (dfY, dfA, dfX) = self.PruneByRelatedness(dfY, dfA, dfX, mReader.dRelCutoff)
        # contruct pheno-specific dummies to address remaining missingness
        (dfY, dfX, dfBinXY) = self.CreateDummies(dfY, dfX, dfBinXY)
        # store variable names
        self.lCovs = dfX.columns.tolist()
        self.lPhenos = dfY.columns.tolist()
        # finalise Mgreml data using canonical transformation
        self.FinaliseData(dfY, dfA, dfX, dfBinXY, iDropLeadPCs, iDropTrailPCs)
        
        self.logger.info('3. STORING ALL MGREML SETTINGS')
        # store whether we have a nested model
        self.bNested = mReader.bNested
        # store booleans indicating if any models
        # are of the type perfect or no correlations
        self.bPerfectRhoG = mReader.bPerfectRhoG
        self.bNoRhoG = mReader.bNoRhoG
        self.bNoRhoE = mReader.bNoRhoE
        self.bPerfectRhoG0 = mReader.bPerfectRhoG0
        self.bNoRhoG0 = mReader.bNoRhoG0
        self.bNoRhoE0 = mReader.bNoRhoE0
        self.dfGenBinFY = mReader.dfGenBinFY
        self.dfGenBinFY0 = mReader.dfGenBinFY0
        self.dfEnvBinFY = mReader.dfEnvBinFY
        self.dfEnvBinFY0 = mReader.dfEnvBinFY0
        # if such rho=0,1 models present: set binary matrices accordingly
        if self.bPerfectRhoG:
            self.dfGenBinFY = MgremlData.SetPerfectRho(self.lPhenos, MgremlData.sF)
        if self.bNoRhoG:
            self.dfGenBinFY = MgremlData.SetNoRho(self.lPhenos, MgremlData.sF)
        if self.bNoRhoE:
            self.dfEnvBinFY = MgremlData.SetNoRho(self.lPhenos, MgremlData.sF)
        if self.bNested:
            if self.bPerfectRhoG0:
                self.dfGenBinFY0 = MgremlData.SetPerfectRho(self.lPhenos, MgremlData.sF)
            if self.bNoRhoG0:
                self.dfGenBinFY0 = MgremlData.SetNoRho(self.lPhenos, MgremlData.sF)
            if self.bNoRhoE0:
                self.dfEnvBinFY0 = MgremlData.SetNoRho(self.lPhenos, MgremlData.sF)
        # store all the other stuff
        self.sPrefix = mReader.sPrefix
        self.bBFGS = mReader.bBFGS
        self.dGradTol = mReader.dGradTol
        self.bSEs = mReader.bSEs
        self.bAllCoeffs = mReader.bAllCoeffs
        self.bStoreIter = mReader.bStoreIter
        self.iStoreIterFreq = mReader.iStoreIterFreq
        self.bReinitialise = mReader.bReinitialise
        self.sInitValsFile = mReader.sInitValsFile
        self.bReinitialise0 = mReader.bReinitialise0
        self.sInitValsFile0 = mReader.sInitValsFile0
        self.logger.info('Settings stored\n')
    
    @staticmethod
    def SetPerfectRho(lLabels, sPrefix):
        lFactorName = [sPrefix + str(0)]
        dfBinFY = pd.DataFrame(data=1, index=pd.Index(lLabels), columns=pd.Index(lFactorName))
        return dfBinFY
        
    @staticmethod
    def SetNoRho(lLabels, sPrefix):
        iT = len(lLabels)
        lFactors = [sPrefix + str(x) for x in range(0,iT)]
        dfBinFY = pd.DataFrame(data=np.eye(iT), index=pd.Index(lLabels), columns=pd.Index(lFactors))
        return dfBinFY
        
    def DetermineIfCovsAreGiven(self, dfX, dfBinXY):
        # assert whether we have covariates and whether we have same covs across traits
        self.bCovs = isinstance(dfX, pd.DataFrame)
        self.bSameCovs = not(isinstance(dfBinXY, pd.DataFrame))
        # if we have no covariates, yet different covariates have been specified
        if not(self.bCovs) and not(self.bSameCovs):
            raise SyntaxError('you have specified different covariates to apply to different traits, without supplying data on those covariates using --covar')
    
    def CleanSpecificationCovariates(self, dfY, dfX, dfBinXY):
        self.logger.info('INSPECTING YOUR COVARIATE MODEL')
        # get indices phenotypes from dfY and dfBinXY
        indY_Y  = dfY.columns
        indXY_Y = dfBinXY.index
        # get labels of covariates from dfX and dfBinXY
        indX_X  = dfX.columns
        indXY_X = dfBinXY.columns
        self.logger.info('You specified which covariates apply to ' + str(dfBinXY.shape[0]) + ' phenotypes')
        self.logger.info('There are ' + str(dfY.shape[1]) + ' phenotypes in your data')
        # all phenotypes in dfY should refer to phenotypes in dfBinXY; abort if problematic
        if not(indY_Y.isin(indXY_Y).all()):
            raise ValueError('there is at least one phenotype for which you have not specified which covariates apply to it') 
        # all covariates in dfBinXY should refer to covariates dfX; abort if problematic
        if not(indXY_X.isin(indX_X).all()):
            raise ValueError('there is at least one phenotype-specific covariate for which you have not supplied the underlying data') 
        if (dfBinXY.shape[0] - dfY.shape[1]) != 0:
            self.logger.info('The ' + str(dfBinXY.shape[0] - dfY.shape[1]) + ' redundant phenotypes will be removed from your coviarate model')
        # eliminate phenotypes from dfBinXY that are not in dfY
        dfBinXY = dfBinXY.loc[indY_Y]
        self.logger.info('You specified for ' + str(dfBinXY.shape[1]) + ' covariates to which phenotypes they apply')
        self.logger.info('There are ' + str(dfX.shape[1]) + ' covariates in your data')
        if (dfX.shape[1] - dfBinXY.shape[1]) != 0:
            self.logger.info('The ' + str(dfX.shape[1] - dfBinXY.shape[1]) + ' redundant covariates will be removed from your coviarate data')
        # eliminate covariates from dfX that are not used according to dfBinXY
        dfX = dfX[indXY_X]
        # if dfBinXY is not binary: abort
        if ((dfBinXY==0).sum().sum() + (dfBinXY==1).sum().sum()) != (dfBinXY.shape[0]*dfBinXY.shape[1]):
            raise ValueError('your model indicating which covariate affects which phenotype does not only comprise zeros and ones')
        # if dfBinXY is all ones: print warning and drop
        if ((dfBinXY==1).sum().sum()) == (dfBinXY.shape[0]*dfBinXY.shape[1]):
            self.logger.warning('Warning: your model indicating which covariate affects which phenotype now comprises only ones')
            self.logger.warning('Assuming all covariates apply to all traits.')
            dfBinXY = None
            self.bSameCovs = True
        return dfX, dfBinXY
    
    def CheckDuplicatesAndRank(self, dfY, dfA, dfX, dfBinXY):
        self.logger.info('CHECKING FOR DUPLICATES AND MULTICOLLINEARITY')
        # if duplicates in index or columns of dfA
        if (dfA.index.duplicated().sum() > 0) or (dfA.columns.duplicated().sum() > 0):
            raise ValueError('you have individuals with the same FID-IID combination in your GRM')
        # if duplicates in columns of dfY
        if dfY.columns.duplicated().sum() > 0:
            raise ValueError('you have phenotypes with duplicate labels in your data')
        # if duplicates in index of dfY
        if dfY.index.duplicated().sum() > 0:
            raise ValueError('you have individuals with the same FID-IID combination in your phenotype data')
        # if we have covariates
        if self.bCovs:
            # if duplicates in columns of dfX
            if dfX.columns.duplicated().sum() > 0:
                raise ValueError('you have covariates with duplicate labels in your data')
            # if duplicates in index of dfX
            if dfX.index.duplicated().sum() > 0:
                raise ValueError('you have covariates with the same FID-IID combination in your covariate data')
            # if we have different covariates across traits
            if not(self.bSameCovs):
                # if duplicates in columns of dfBinXY
                if dfBinXY.columns.duplicated().sum() > 0:
                    raise ValueError('in your specification which covariates apply to which phenotypes, you have specified duplicate covariates')
                # if duplicates in index of dfBinXY
                if dfBinXY.index.duplicated().sum() > 0:
                    raise ValueError('in your specification which covariates apply to which phenotypes, you have specified duplicate phenotypes')
            # get matrix of covariates and set missings to zero
            mX = np.array(dfX)
            mX[np.where(np.isnan(mX))] = 0
            # compute the rank of the matrix of covariates
            dRankX = np.linalg.matrix_rank(mX)
            # count the number of covariates
            iK = dfX.shape[1]
            # if rank is below the number of covariates
            if dRankX < iK:
                raise ValueError('your matrix of covariates does not have full rank, i.e. there is perfect multicollinearity')
        # get matrix of phenotypes and set missings to zero
        mY = np.array(dfY.copy())
        mY[np.where(np.isnan(mY))] = 0
        # compute the rank of the matrix of phenotypes
        dRankY = np.linalg.matrix_rank(mY)
        # count the number of phenotypes
        iT = dfY.shape[1]
        # if rank is below the number of covariates
        if dRankY < iT:
            self.logger.warning('Warning: your phenotype data is rank deficient, i.e. there is perfect multicollinearity.')
            self.logger.warning('This may lead to poorly identified models.')
        # compute the phenotype variance of each trait
        vVarY = np.var(mY,axis=0)
        # if there is at least one trait with no variance
        if (vVarY == 0).sum() > 0:
            raise ValueError('you have specified one or more phenotypes without any variance at all')
        self.logger.info('No irregularities found')
    
    def FindOverlapAndSort(self, dfY, dfA, dfX, dfBinXY):
        self.logger.info('JOINING YOUR DATA')
        self.logger.info('There are ' + str(dfY.shape[0]) + ' individuals in your phenotype data')
        self.logger.info('There are ' + str(dfA.shape[0]) + ' individuals in your GRM')
        # if we have covariates
        if self.bCovs:
            self.logger.info('There are ' + str(dfX.shape[0]) + ' individuals in your covariate data')
            # define list of dataframes to include dfX
            ldfs = [dfY,dfA,dfX]
        else:
            # define list of dataframes without dfX
            ldfs = [dfY,dfA]
        # construct empty list of dataframes, for multi index of each df in ldfs
        lIDs = list()
        # for each dataframe in ldfs
        for df in ldfs:
            # get only the multi index for FID-IID combos, ignoring data and columns
            # and append to lIDs
            lIDs.append(pd.DataFrame(data=None, columns=None, index=df.index))
        # find the FID-IID combos that are present in all relevant dataframes
        miIDs = pd.concat(lIDs, axis=1, join='inner').index
        self.logger.info('Keeping only the ' + str(len(miIDs)) + ' overlapping individuals')
        # select and order phenotypes and GRM by those individuals
        dfY = dfY.loc[miIDs]
        dfA = dfA.loc[miIDs,miIDs]
        # double-check if everything is now lined up
        if not(all(dfY.index == dfA.index)) or not(all(dfY.index == dfA.columns)):
            raise ValueError('the FID-IID combinations in your GRM cannot be lined up properly with those in your phenotype data')
        # if we have covariates
        if self.bCovs:
            # select and order covariates by those individuals
            dfX = dfX.loc[miIDs]
            # double-check if everything is now lined up
            if not(all(dfY.index == dfX.index)):
                raise ValueError('the FID-IID combinations in your covariates cannot be lined up properly with those in your phenotype data')
            # if we have different covariates across traits
            if not(self.bSameCovs):
                # get labels of phenotypes and covariates
                # according to order in dfY and dfX resp.
                indY_Y = dfY.columns
                indX_X = dfX.columns
                # put dfBinXY in the same order
                dfBinXY = (dfBinXY.loc[indY_Y])[indX_X]
                # double-check if everything is now lined up
                if not(all(dfBinXY.index == indY_Y)):
                    raise ValueError('your specification which covariates applies to which phenotype cannot be properly lined up with the phenotypes')
                if not(all(dfBinXY.columns == indX_X)):
                    raise ValueError('your specification which covariates applies to which phenotype cannot be properly lined up with the covariates')
        return dfY, dfA, dfX, dfBinXY
    
    def DropMissings(self, dfY, dfA, dfX, dfBinXY):
        self.logger.info('DROPPING PROBLEMATIC INDIVIDUALS BECAUSE OF MISSING PHENOTYPES AND/OR COVARIATES')
        # if we have covariates
        if self.bCovs:
            iCount = 0
            # and all covariates apply to all traits
            if self.bSameCovs:
                tDrop = tqdm(total=dfY.shape[0])
                # and for a given observations
                for i in dfY.index:
                    tDrop.update(1)
                    # if at least one covariate is  missing
                    if (dfX.loc[i].isnull().sum() > 0) or (dfX.loc[i].isna().sum() > 0):
                        # set all phenotypes for that individual to missing
                        dfY.loc[i] = None
                        iCount += 1
                tDrop.close()
                self.logger.info('Found ' + str(iCount) + ' individuals with missing data on one or more covariates')
                self.logger.info('Setting all phenotypes to missing for those individuals')
            else: # if not all covariates apply to all traits
                tDrop = tqdm(total=dfY.shape[0])
                # for a given observation
                for i in dfY.index:
                    tDrop.update(1)
                    # find indices of missing covariates
                    vIndMissingCovs = np.array(np.where(dfX.loc[i].isnull() | dfX.loc[i].isna())).ravel()
                    if len(vIndMissingCovs) > 0:
                        iCount += 1
                    # for each missing covariate
                    for j in vIndMissingCovs:
                        # find phenotypes affected by the missing covariate
                        vIndPhenotypes = np.array(np.where(dfBinXY.iloc[:,j] == 1)).ravel()
                        # for each of those phenotypes
                        for k in vIndPhenotypes:
                            # set phenotypic value to missing
                            dfY.loc[i].iloc[k] = None
                tDrop.close()
                self.logger.info('Found ' + str(iCount) + ' individuals with missing data on one or more covariates')
                self.logger.info('Setting all phenotypes affected by missing covariates to missing for those individiuals')
        # count the number of traits and observations in dfY
        iT = dfY.shape[1]
        iN = dfY.shape[0]
        # count number of observations to be dropped i.e. with all phenos missing
        iM = ((dfY.isnull() | dfY.isna()).sum(axis=1) == iT).sum()
        self.logger.info('Dropping ' + str(iM) + ' out of ' + str(iN) + ' individuals from data for whom all phenotypes are now missing')
        # keep only observations that have at least one pheno non-missing
        dfY = dfY[((dfY.isnull() | dfY.isna()).sum(axis=1) < iT)]
        # get list of individuals that remain
        miIDs = dfY.index
        # keep only those individuals in GRM
        dfA = dfA.loc[miIDs,miIDs]
        # if we have covariates
        if self.bCovs:
            # keep only those individuals in data on covariates
            dfX = dfX.loc[miIDs]
        return dfY, dfA, dfX
    
    def PruneByRelatedness(self, dfY, dfA, dfX, dRelCutoff):
        self.logger.info('APPLYING RELATEDNESS CUTOFF')
        self.logger.warning('Warning: this method is currently still only a placeholder. Cutoff has not been applied.')
        self.logger.info('Removing individuals such that there is no relatedness in excess of ' + str(dRelCutoff) + ' in GRM')
        self.logger.info('Using greedy algorithm as implemented in PLINK v...')
        iN = dfY.shape[0]
        iDropped = 0
        # INSERT CODE ERIC HERE
        self.logger.info('Dropping ' + str(iDropped) + ' out of ' + str(iN) + ' individiuals')
        return dfY, dfA, dfX
        
    def CreateDummies(self, dfY, dfX, dfBinXY):
        # if there are any missings at all
        if any(dfY.isnull() | dfY.isna()):
            self.logger.info('CREATING PHENOTYPE-SPECIFIC DUMMY VARIABLES TO CONTROL FOR REMAINING MISSINGNESS')
            # if there are no covariates yet
            if not(self.bCovs):
                # initialise dfX and dfBinXY
                dfX = pd.DataFrame(data=None, index=dfY.index)
                dfBinXY = pd.DataFrame(data=None, index=dfY.columns)
                # set covariates being present and not the same across traits
                self.bCovs = True
                self.bSameCovs = False
            # if there are covariates
            else:
                # but same covariates for all traits so far
                if self.bSameCovs:
                    # initialise dfBinXY as bunch of ones
                    dfBinXY = pd.DataFrame(data=1, index=dfY.columns, columns=dfX.columns)
                    # set covariates to being not the same across traits
                    self.bSameCovs = False
            iCount = 0
            # for each trait
            for t in dfY.columns:
                # get all observations with missing values
                miIDs = dfY.loc[dfY[t].isnull() | dfY[t].isna(),t].index
                # for each observation with missing value
                for i in miIDs:
                    # set missing phenotype to zero
                    dfY.loc[i,t] = 0
                    # construct label for new dummy variable
                    lLabel = ['dummy_trait_' + str(t) + '_obs_' + str(i)]
                    # create dummy variable
                    dfXadd = pd.DataFrame(data=0, index=dfX.index, columns=lLabel)
                    dfXadd.loc[i] = 1
                    # append to existing set of covariates
                    dfX = pd.concat([dfX, dfXadd], axis=1, join='inner')
                    # create corresponding entry for dfBinXY
                    dfBinXYadd = pd.DataFrame(data=0, index=dfY.columns, columns=lLabel)
                    dfBinXYadd.loc[t] = 1
                    # append to existing specification which covs apply to which phenos
                    dfBinXY = pd.concat([dfBinXY, dfBinXYadd], axis=1, join='inner')
                    iCount += 1
            # replace missing in dfX by 0
            dfX = dfX.fillna(0)
            self.logger.info('Added ' + str(iCount) + ' phenotype-specific dummies to your covariate model')
            if iCount*dfY.shape[1] > MgremlData.iManyDummies:
                self.logger.warning('This is a large number of phenotype-specific covariates, given you have ' + str(dfY.shape[1]) + ' traits in your data')
                self.logger.warning(str(iCount*dfY.shape[1]) + ' additional fixed-effect covariates implied, of which ' + str(iCount*dfY.shape[1] - iCount) + ' are set to zero')
                self.logger.warning('CPU time of MGREML may increase dramatically')
                self.logger.warning('Consider running MGREML on a subset of your data with a much lower degree of missingness')
        return dfY, dfX, dfBinXY
    
    def FinaliseData(self, dfY, dfA, dfX, dfBinXY, iDropLeadPCs, iDropTrailPCs):
        self.logger.info('FINALISING DATA BEFORE MGREML ANALYSIS')
        # convert dataframes to numpy arrays
        mY = np.array(dfY)
        mA = np.array(dfA)
        # if we have covariates
        if self.bCovs:
            # also convert the set of covariates to numpy array
            mX = np.array(dfX)
            # if we do not have same covariates across traits
            if not(self.bSameCovs):
                # also convert appropriate dataframe to numpy array
                # and store
                self.mBinXY = np.array(dfBinXY).astype(int)
        # stabilise GRM
        mA = (mA+mA.T)/2
        self.logger.info('Computing eigenvalue decomposition of the GRM')
        self.logger.info('This may take a long time...')
        # compute its EVD
        (vD,mP) = np.linalg.eigh(mA)
        self.logger.info('Applying the canonical transformation to your data')
        self.logger.info('Sample size prior to the canonical transformation is ' + str(mY.shape[0]))
        self.logger.info('Ignoring ' + str(iDropLeadPCs) + ' leading eigenvectors to control for population stratification')
        if iDropTrailPCs > 0:
            self.logger.info('Ignoring ' + str(iDropTrailPCs) + ' trailing eigenvectors to improve computational efficiency')
        # ignore trailing columns and leading columns from eigenvector matrix
        mPT = (mP[:,iDropTrailPCs:-iDropLeadPCs]).T
        # ignore trailing and leading values from eigenvalue vector
        vD  = vD[iDropTrailPCs:-iDropLeadPCs]
        # store eigenvalues and squared eigenvalues
        self.vD = vD
        self.vDSq = vD**2
        # apply canonical transform to phenotypes and store
        self.mY = mPT@mY
        # count sample size and number of phenotypes, and store
        self.iN = self.mY.shape[0]
        self.iT = self.mY.shape[1]
        # if we have covariates
        if self.bCovs:
            # apply canonical transform to covariates and store
            self.mX = mPT@mX
            # for computational reasons, also store its transpose
            self.mXT = self.mX.T
            # count number of covariates, and store
            self.iK = self.mX.shape[1]
            # compute XTX
            mXTX = np.matmul(self.mXT,self.mX)
            # if same covariates across traits
            if self.bSameCovs:
                # compute log|X'X| across all traits in one go using EVD of X'X
                (vThetaXTX,_) = np.linalg.eigh(mXTX)
                # if any eigenvalue is too close to zero or negative
                if any(vThetaXTX < abs(np.finfo(float).eps)):
                    # raise an error with a proper explanation of the likely cause
                    raise ValueError('your covariates are rank deficient after the canonical transformation (i.e. perfectly multicollinear). Likely reason: you specified principal components (PCs) from your genetic data as fixed-effect covariates. MGREML already controls for population stratification in the canonical transformation. Please do not control for PCs manually as well. Rather, use --ignore-pcs INTEGER, to indicate for how many PCs you want to control via the canonical transformation.')
                # compute log|X'X| and store
                self.dLogDetXTX = (self.iT)*np.log(vThetaXTX).sum()
            # if not same across traits
            else:
                # get indices of these covariates in terms of Z matrix
                self.vIndCovs = np.array(np.where(np.array(self.mBinXY).ravel()==1)).ravel()
                # count total number of covariates across traits
                self.iKtotal = self.mBinXY.sum()
                # initialise log|X'X|
                dLogDetXTX = 0
                # for each trait
                for it in range(0,self.iT):
                    # get binary vector indicating which covariates affect current trait
                    vBinaryX = self.mBinXY[it,:]
                    # find indices
                    vIndBinaryX = np.array(np.where(vBinaryX==1)).ravel()
                    # compute log|X'X| for given trait using EVD
                    (vThetaXTX,_) = np.linalg.eigh(mXTX[vIndBinaryX,:][:,vIndBinaryX])
                    # if any eigenvalue is too close to zero or negative
                    if any(vThetaXTX < abs(np.finfo(float).eps)):
                        # raise an error with a proper explanation of the likely cause
                        raise ValueError('Your covariates are rank deficient after the canonical transformation (i.e. perfectly multicollinear). Likely reason: you specified principal components (PCs) from your genetic data as fixed-effect covariates. MGREML already controls for population stratification in the canonical transformation. Please do not control for PCs manually as well. Rather, use --ignore-pcs INTEGER, to indicate for how many PCs you want to control via the canonical transformation.')
                    # compute log|X'X| and add to grand total
                    dLogDetXTX = dLogDetXTX + np.log(vThetaXTX).sum()
                # store log|X'X|
                self.dLogDetXTX = dLogDetXTX
        # otherwise set log|X'X| to zero
        else:
            self.dLogDetXTX = 0
        self.logger.info('Sample size after the canonical transformation is ' + str(self.iN)) 
        self.logger.info('Final sample size, N = ' + str(self.iN))
        self.logger.info('Final number of traits, T = ' + str(self.iT))
        if self.bCovs:
            self.logger.info('Final number of unique covariates, k = ' + str(self.iK))
            if self.bSameCovs:
                self.logger.info('Final number of fixed effects, K = ' + str(self.iK*self.iT))
            else:
                self.logger.info('Final number of fixed effects, K = ' + str(self.iKtotal))
        self.logger.info('Data cleaning complete\n')