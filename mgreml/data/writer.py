import pandas as pd
from mgreml.analysis import comparison
from mgreml.analysis import estimator

class DataWriter:
    
    sPath = './output/'
    sExtension = 'out'
    sH2 = 'HSq.'
    sRG = 'RhoG.'
    sRE = 'RhoE.'
    sH0 = 'null.'
    sHA = 'alt.'
    sSE = 'SE.'
    
    def __init__(self, estimates, sPrefix = 'results.'):
        self.sPrefix = sPrefix
        if isinstance(estimates, comparison.NestedEstimators):
            self.bNested = True
            self.estimates = estimates
        elif isinstance(estimates, estimator.MgremlEstimator):
            self.bNested = False
            self.estimates = estimates
        else:
            raise TypeError('No estimates have been provided to write result files from')
        if not(self.estimates.IsConverged()):
            raise ValueError('Trying to write output results while estimates have not converged')
        if not(self.estimates.IsDone()):
            raise ValueError('Trying to write output results while final statistic (e.g. genetic correlations) have not been calculated')
       
    def WriteHSq(self):
        if self.bNested:
            # get heritabilities 
            vHSq0 = self.estimates.estimator_res.vHSq
            vHSqA = self.estimates.estimator_unres.vHSq
            # get trait labels
            lPhenos = self.estimates.estimator_unres.mgreml_model.data.lPhenos
            # set dataframes
            dfHSq0 = pd.DataFrame(vHSq0,index=lPhenos)
            dfHSqA = pd.DataFrame(vHSqA,index=lPhenos)
            # set filenames
            sHSq0 = DataWriter.sPath + DataWriter.sH2 + DataWriter.sH0 + self.sPrefix + DataWriter.sExtension
            sHSqA = DataWriter.sPath + DataWriter.sH2 + DataWriter.sHA + self.sPrefix + DataWriter.sExtension
            # write dataframes
            dfHSq0.to_csv(sHSq0)
            dfHSqA.to_csv(sHSqA)
            # if SEs are desired, store them
            if (self.estimates.estimator_unres.bSEs):
                # get heritability standard errors 
                vHSq0SE = self.estimates.estimator_res.vHSqSE
                vHSqASE = self.estimates.estimator_unres.vHSqSE
                # set dataframes
                dfHSq0SE = pd.DataFrame(vHSq0SE,index=lPhenos)
                dfHSqASE = pd.DataFrame(vHSqASE,index=lPhenos)
                # set filenames
                sHSq0SE = DataWriter.sPath + DataWriter.sH2 + DataWriter.sH0 + DataWriter.sSE + self.sPrefix + DataWriter.sExtension
                sHSqASE = DataWriter.sPath + DataWriter.sH2 + DataWriter.sHA + DataWriter.sSE + self.sPrefix + DataWriter.sExtension
                # write dataframes
                dfHSq0SE.to_csv(sHSq0SE)
                dfHSqASE.to_csv(sHSqASE)
        else:
            # get heritability
            vHSq = self.estimates.vHSq
            # get trait labels
            lPhenos = self.estimates.mgreml_model.data.lPhenos
            # set dataframe
            dfHSq = pd.DataFrame(vHSq,index=lPhenos)
            # set filename
            sHSq = DataWriter.sPath + DataWriter.sH2 + self.sPrefix + DataWriter.sExtension
            # write dataframe
            dfHSq.to_csv(sHSq)
            # if SEs are desired, store them
            if (self.estimates.bSEs):
                # get heritability standard errors 
                vHSqSE = self.estimates.vHSqSE
                # set dataframe
                dfHSqSE = pd.DataFrame(vHSqSE,index=lPhenos)
                # set filenames
                sHSqSE = DataWriter.sPath + DataWriter.sH2 + DataWriter.sSE + self.sPrefix + DataWriter.sExtension
                # write dataframes
                dfHSqSE.to_csv(sHSqSE)

    def WriteRho(self):
        if self.bNested:
            # get correlation matrices 
            mRhoG0 = self.estimates.estimator_res.mRhoG
            mRhoE0 = self.estimates.estimator_res.mRhoE
            mRhoGA = self.estimates.estimator_unres.mRhoG
            mRhoEA = self.estimates.estimator_unres.mRhoE
            # get trait labels
            lPhenos = self.estimates.estimator_unres.mgreml_model.data.lPhenos
            # set dataframes
            dfRhoG0 = pd.DataFrame(mRhoG0,index=lPhenos,columns=lPhenos)
            dfRhoE0 = pd.DataFrame(mRhoE0,index=lPhenos,columns=lPhenos)
            dfRhoGA = pd.DataFrame(mRhoGA,index=lPhenos,columns=lPhenos)
            dfRhoEA = pd.DataFrame(mRhoEA,index=lPhenos,columns=lPhenos)
            # set filenames
            sRhoG0 = DataWriter.sPath + DataWriter.sRG + DataWriter.sH0 + self.sPrefix + DataWriter.sExtension
            sRhoE0 = DataWriter.sPath + DataWriter.sRE + DataWriter.sH0 + self.sPrefix + DataWriter.sExtension
            sRhoGA = DataWriter.sPath + DataWriter.sRG + DataWriter.sHA + self.sPrefix + DataWriter.sExtension
            sRhoEA = DataWriter.sPath + DataWriter.sRE + DataWriter.sHA + self.sPrefix + DataWriter.sExtension
            # write dataframes
            dfRhoG0.to_csv(sRhoG0)
            dfRhoE0.to_csv(sRhoE0)
            dfRhoGA.to_csv(sRhoGA)
            dfRhoEA.to_csv(sRhoEA)
            # if SEs are desired, store them
            if (self.estimates.estimator_unres.bSEs):
                # get SE matrices 
                mRhoG0SE = self.estimates.estimator_res.mRhoGSE
                mRhoE0SE = self.estimates.estimator_res.mRhoESE
                mRhoGASE = self.estimates.estimator_unres.mRhoGSE
                mRhoEASE = self.estimates.estimator_unres.mRhoESE
                # set dataframes
                dfRhoG0SE = pd.DataFrame(mRhoG0SE,index=lPhenos,columns=lPhenos)
                dfRhoE0SE = pd.DataFrame(mRhoE0SE,index=lPhenos,columns=lPhenos)
                dfRhoGASE = pd.DataFrame(mRhoGASE,index=lPhenos,columns=lPhenos)
                dfRhoEASE = pd.DataFrame(mRhoEASE,index=lPhenos,columns=lPhenos)
                # set filenames
                sRhoG0SE = DataWriter.sPath + DataWriter.sRG + DataWriter.sH0 + DataWriter.sSE + self.sPrefix + DataWriter.sExtension
                sRhoE0SE = DataWriter.sPath + DataWriter.sRE + DataWriter.sH0 + DataWriter.sSE + self.sPrefix + DataWriter.sExtension
                sRhoGASE = DataWriter.sPath + DataWriter.sRG + DataWriter.sHA + DataWriter.sSE + self.sPrefix + DataWriter.sExtension
                sRhoEASE = DataWriter.sPath + DataWriter.sRE + DataWriter.sHA + DataWriter.sSE + self.sPrefix + DataWriter.sExtension
                # write dataframes
                dfRhoG0SE.to_csv(sRhoG0SE)
                dfRhoE0SE.to_csv(sRhoE0SE)
                dfRhoGASE.to_csv(sRhoGASE)
                dfRhoEASE.to_csv(sRhoEASE)
        else:
            # get correlation matrices 
            mRhoG = self.estimates.mRhoG
            mRhoE = self.estimates.mRhoE
            # get trait labels
            lPhenos = self.estimates.mgreml_model.data.lPhenos
            # set dataframes
            dfRhoG = pd.DataFrame(mRhoG,index=lPhenos,columns=lPhenos)
            dfRhoE = pd.DataFrame(mRhoE,index=lPhenos,columns=lPhenos)
            # set filenames
            sRhoG = DataWriter.sPath + DataWriter.sRG + self.sPrefix + DataWriter.sExtension
            sRhoE = DataWriter.sPath + DataWriter.sRE + self.sPrefix + DataWriter.sExtension
            # write dataframes
            dfRhoG.to_csv(sRhoG)
            dfRhoE.to_csv(sRhoE)
            # if SEs are desired, store them
            if (self.estimates.bSEs):
                # get SE matrices 
                mRhoGSE = self.estimates.mRhoGSE
                mRhoESE = self.estimates.mRhoESE
                # set dataframes
                dfRhoGSE = pd.DataFrame(mRhoGSE,index=lPhenos,columns=lPhenos)
                dfRhoESE = pd.DataFrame(mRhoESE,index=lPhenos,columns=lPhenos)
                # set filenames
                sRhoGSE = DataWriter.sPath + DataWriter.sRG + DataWriter.sSE + self.sPrefix + DataWriter.sExtension
                sRhoESE = DataWriter.sPath + DataWriter.sRE + DataWriter.sSE + self.sPrefix + DataWriter.sExtension
                # write dataframes
                dfRhoGSE.to_csv(sRhoGSE)
                dfRhoESE.to_csv(sRhoESE)

