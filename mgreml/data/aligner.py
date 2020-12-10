import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None

class MgremlData:
    
    def __init__(self, mReader):
        # read out MgremlReader
        self.logger = mReader.logger
        dfA = mReader.dfA
        dfY = mReader.dfY
        dfX = mReader.dfX
        dfBinXY = mReader.dfBinXY
        iDropLeadPCs = mReader.iDropLeadPCs
        iDropTrailPCs = mReader.iDropTrailPCs
        # find out if we have covariates
        (bCovs, bSameCovs) = self.DetermineIfCovsAreGiven(dfX, dfBinXY)
        # store resulting booleans as attributes
        self.bCovs     = bCovs
        self.bSameCovs = bSameCovs
        # if we have covariates, and different covariates apply to different traits
        if bCovs and not(bSameCovs):
            # clean up the specification of the part of the model on covariates
            (dfX, dfBinXY) = self.CleanSpecificationCovariates(dfY, dfX, dfBinXY)
        # check if there are duplicates and/or rank defficiency
        self.CheckDuplicatesAndRank(dfY, dfA, dfX, dfBinXY)
        # find overlapping individuals and sort data
        (dfY, dfA, dfX, dfBinXY) = self.FindOverlapAndSort(dfY, dfA, dfX, dfBinXY)
        # drop observations where missingness affects all traits
        (dfY, dfA, dfX) = self.DropMissings(dfY, dfA, dfX, dfBinXY)
        # contruct pheno-specific dummies to address remaining missingness
        (dfY, dfX, dfBinXY) = self.CreateDummies(dfY, dfX, dfBinXY)
        # store variable names
        self.lCovs = dfX.columns.tolist()
        self.lPhenos = dfY.columns.tolist()
        # finalise Mgreml data using canonical transformation
        self.FinaliseData(dfY, dfA, dfX, dfBinXY, iDropLeadPCs, iDropTrailPCs)
        
    def DetermineIfCovsAreGiven(self, dfX, dfBinXY):
        # assert whether we have covariates and whether we have same covs across traits
        bCovs     = isinstance(dfX, pd.DataFrame)
        bSameCovs = not(isinstance(dfBinXY, pd.DataFrame))
        # if we have no covariates, yet different covariates have been specified
        if not(bCovs) and not(bSameCovs):
            self.logger.error('Error: you have specified different covariates to apply to different traits, without supplying data on those covariates using --covar.')
            raise TypeError
        return bCovs, bSameCovs
    
    def CleanSpecificationCovariates(self, dfY, dfX, dfBinXY):
        # get indices phenotypes from dfY and dfBinXY
        indY_Y  = dfY.columns
        indXY_Y = dfBinXY.index
        # get labels of covariates from dfX and dfBinXY
        indX_X  = dfX.columns
        indXY_X = dfBinXY.columns
        # all phenotypes in dfY should refer to phenotypes in dfBinXY; abort if problematic
        if not(indY_Y.isin(indXY_Y).all()):
            raise ValueError('There is at least one phenotype for which you have not specified which covariates apply to it.') 
        # all covariates in dfBinXY should refer to covariates dfX; abort if problematic
        if not(indXY_X.isin(indX_X).all()):
            raise ValueError('There is at least one phenotype-specific covariate for which you have not supplied the underlying data.') 
        # eliminate phenotypes from dfBinXY that are not in dfY
        dfBinXY = dfBinXY.loc[indY_Y]
        # eliminate covariates from dfX that are not used according to dfBinXY
        dfX = dfX[indXY_X]
        # if dfBinXY is not binary: abort
        if ((dfBinXY==0).sum().sum() + (dfBinXY==1).sum().sum()) != (dfBinXY.shape[0]*dfBinXY.shape[1]):
            raise ValueError('Your data indicating which covariate affects which phenotype does not only comprise zeros and ones.')
        # if dfBinXY is all ones: print warning and drop
        if ((dfBinXY==1).sum().sum()) == (dfBinXY.shape[0]*dfBinXY.shape[1]):
            print('Warning! Your data indicating which covariate affects which phenotype comprises only ones.')
            print('Assuming all covariates apply to all traits.')
            dfBinXY        = None
            self.bSameCovs = True
        return dfX, dfBinXY
    
    def CheckDuplicatesAndRank(self, dfY, dfA, dfX, dfBinXY):
        # if duplicates in index or columns of dfA
        if (dfA.index.duplicated().sum() > 0) or (dfA.columns.duplicated().sum() > 0):
            raise ValueError('You have individuals with the same FID-IID combination in your GRM.')
        # if duplicates in columns of dfY
        if dfY.columns.duplicated().sum() > 0:
            raise ValueError('You have phenotypes with duplicate labels in your data.')
        # if duplicates in index of dfY
        if dfY.index.duplicated().sum() > 0:
            raise ValueError('You have individuals with the same FID-IID combination in your phenotype data.')
        # if we have covariates
        if self.bCovs:
            # if duplicates in columns of dfX
            if dfX.columns.duplicated().sum() > 0:
                raise ValueError('You have covariates with duplicate labels in your data.')
            # if duplicates in index of dfX
            if dfX.index.duplicated().sum() > 0:
                raise ValueError('You have covariates with the same FID-IID combination in your covariate data.')
            # if we have different covariates across traits
            if not(self.bSameCovs):
                # if duplicates in columns of dfBinXY
                if dfBinXY.columns.duplicated().sum() > 0:
                    raise ValueError('In your specification which covariates apply to which phenotypes, you have specified duplicate covariates.')
                # if duplicates in index of dfBinXY
                if dfBinXY.index.duplicated().sum() > 0:
                    raise ValueError('In your specification which covariates apply to which phenotypes, you have specified duplicate phenotypes.')
            # get matrix of covariates and set missings to zero
            mX = np.array(dfX)
            mX[np.where(np.isnan(mX))] = 0
            # compute the rank of the matrix of covariates
            dRankX = np.linalg.matrix_rank(mX)
            # count the number of covariates
            iK = dfX.shape[1]
            # if rank is below the number of covariates
            if dRankX < iK:
                raise ValueError('Your matrix of covariates does not have full rank. I.e. there is perfect multicollinearity.')
        # get matrix of phenotypes and set missings to zero
        mY = np.array(dfY)
        mY[np.where(np.isnan(mY))] = 0
        # compute the rank of the matrix of phenotypes
        dRankY = np.linalg.matrix_rank(mY)
        # count the number of phenotypes
        iT = dfY.shape[1]
        # if rank is below the number of covariates
        if dRankY < iT:
            print('Warning! Your data of phenotypes does not have full rank.')
            print('I.e. there is perfect multicollinearity between your phenotypes.')
            print('This may lead to poorly identified models.')
        # compute the phenotype variance of each trait
        vVarY = np.var(mY,axis=0)
        # if there is at least one trait with no variance
        if (vVarY == 0).sum() > 0:
            raise ValueError('You have specified one or more phenotypes without any variance at all.')
    
    def FindOverlapAndSort(self, dfY, dfA, dfX, dfBinXY):
        # if we have covariates
        if self.bCovs:
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
        # select and order phenotypes and GRM by those individuals
        dfY = dfY.loc[miIDs]
        dfA = dfA.loc[miIDs,miIDs]
        # double-check if everything is now lined up
        if not(all(dfY.index == dfA.index)) or not(all(dfY.index == dfA.columns)):
            raise ValueError('The FID-IID combinations in your GRM cannot be lined up properly with those in your phenotype data.')
        # if we have covariates
        if self.bCovs:
            # select and order covariates by those individuals
            dfX = dfX.loc[miIDs]
            # double-check if everything is now lined up
            if not(all(dfY.index == dfX.index)):
                raise ValueError('The FID-IID combinations in your covariates cannot be lined up properly with those in your phenotype data.')
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
                    raise ValueError('Your specification which covariates applies to which phenotype cannot be properly lined up with the phenotypes.')
                if not(all(dfBinXY.columns == indX_X)):
                    raise ValueError('Your specification which covariates applies to which phenotype cannot be properly lined up with the covariates.')   
        return dfY, dfA, dfX, dfBinXY
    
    def DropMissings(self, dfY, dfA, dfX, dfBinXY):
        # if we have covariates
        if self.bCovs:
            # and all covariates apply to all traits
            if self.bSameCovs:
                # and for a given observations
                for i in dfY.index:
                    # if at least one covariate is  missing
                    if (dfX.loc[i].isnull().sum() > 0) or (dfX.loc[i].isna().sum() > 0):
                        # set all phenotypes for that individual to missing
                        dfY.loc[i] = None
            else: # if not all covariates apply to all traits
                # for a given observation
                for i in dfY.index:
                    # find indices of missing covariates
                    vIndMissingCovs = np.array(np.where(dfX.loc[i].isnull() | dfX.loc[i].isna())).ravel()
                    # for each missing covariate
                    for j in vIndMissingCovs:
                        # find phenotypes affected by the missing covariate
                        vIndPhenotypes = np.array(np.where(dfBinXY.iloc[:,j] == 1)).ravel()
                        # for each of those phenotypes
                        for k in vIndPhenotypes:
                            # set phenotypic value to missing
                            dfY.loc[i].iloc[k] = None
        # count the number of traits in dfY
        iT = dfY.shape[1]
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
    
    def CreateDummies(self, dfY, dfX, dfBinXY):
        # if there are any missings at all
        if any(dfY.isnull() | dfY.isna()):
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
            # for each observation
            for i in dfX.index:
                # find columns where dfX has missing values, and set value
                # of those columns to zero
                dfX.loc[i,(dfX.loc[i].isna() | dfX.loc[i].isnull())] = 0
        return dfY, dfX, dfBinXY
    
    def FinaliseData(self, dfY, dfA, dfX, dfBinXY, iDropLeadPCs, iDropTrailPCs):
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
        # compute its EVD
        (vD,mP) = np.linalg.eigh(mA)
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
                    raise ValueError('Your covariates are rank deficient after the canonical transformation (i.e. perfectly multicollinear). Likely reason: you specified principal components (PCs) from your genetic data as fixed-effect covariates. MGREML already controls for population stratification in the canonical transformation. Please do not control for PCs manually as well. Rather, use --control-lead-pcs NUM, to indicate for how many PCs you want to control via the canonical transformation.')
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
                        raise ValueError('Your covariates are rank deficient after the canonical transformation (i.e. perfectly multicollinear). Likely reason: you specified principal components (PCs) from your genetic data as fixed-effect covariates. MGREML already controls for population stratification in the canonical transformation. Please do not control for PCs manually as well. Rather, use --control-lead-pcs NUM, to indicate for how many PCs you want to control via the canonical transformation.')
                    # compute log|X'X| and add to grand total
                    dLogDetXTX = dLogDetXTX + np.log(vThetaXTX).sum()
                # store log|X'X|
                self.dLogDetXTX = dLogDetXTX
        # otherwise set log|X'X| to zero
        else:
            self.dLogDetXTX = 0

