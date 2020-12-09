import pandas as pd
import numpy as np
from analysis import comparison
from analysis import estimator

class DataWriter:
    
    sExtension = 'out'
    sH2 = 'HSq.'
    sRG = 'RhoG.'
    sRE = 'RhoE.'
    sH0 = 'null.'
    sHA = 'alt.'
    sSE = 'SE.'
    sLL = 'loglik.'
    sLRT = 'LRT.'
    sGLSest = 'GLS.est.'
    sGLSvar = 'GLS.var.'
    sCoeff = 'coeff.'
    sCoeffvar = 'coeff.var.'
    
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
            vHSq0 = self.estimates.estimator0.vHSq
            vHSqA = self.estimates.estimatorA.vHSq
            # get trait labels
            lPhenos = self.estimates.estimatorA.mgreml_model.data.lPhenos
            # set dataframes
            dfHSq0 = pd.DataFrame(vHSq0,index=lPhenos)
            dfHSqA = pd.DataFrame(vHSqA,index=lPhenos)
            # set filenames
            sHSq0 = self.sPrefix + DataWriter.sH2 + DataWriter.sH0 + DataWriter.sExtension
            sHSqA = self.sPrefix + DataWriter.sH2 + DataWriter.sHA + DataWriter.sExtension
            # write dataframes
            dfHSq0.to_csv(sHSq0)
            dfHSqA.to_csv(sHSqA)
            # if SEs are desired, store them
            if (self.estimates.estimatorA.bSEs):
                # get heritability standard errors 
                vHSq0SE = self.estimates.estimator0.vHSqSE
                vHSqASE = self.estimates.estimatorA.vHSqSE
                # set dataframes
                dfHSq0SE = pd.DataFrame(vHSq0SE,index=lPhenos)
                dfHSqASE = pd.DataFrame(vHSqASE,index=lPhenos)
                # set filenames
                sHSq0SE = self.sPrefix + DataWriter.sH2 + DataWriter.sH0 + DataWriter.sSE + DataWriter.sExtension
                sHSqASE = self.sPrefix + DataWriter.sH2 + DataWriter.sHA + DataWriter.sSE + DataWriter.sExtension
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
            sHSq = self.sPrefix + DataWriter.sH2 + DataWriter.sExtension
            # write dataframe
            dfHSq.to_csv(sHSq)
            # if SEs are desired, store them
            if (self.estimates.bSEs):
                # get heritability standard errors 
                vHSqSE = self.estimates.vHSqSE
                # set dataframe
                dfHSqSE = pd.DataFrame(vHSqSE,index=lPhenos)
                # set filenames
                sHSqSE = self.sPrefix + DataWriter.sH2 + DataWriter.sSE + DataWriter.sExtension
                # write dataframes
                dfHSqSE.to_csv(sHSqSE)

    def WriteLRT(self):
        if not(self.bNested):
            raise TypeError('Trying to write results for likelihood-ratio test, while no nested model has been estimated')
        # set filename
        sLRTfile = self.sPrefix + DataWriter.sLRT + DataWriter.sExtension
        with open(sLRTfile, 'w') as oLRTfile:
            oLRTfile.write('Results of likelihood-ratio test with ' + str(self.estimates.iDF) + ' degrees of freedom:\n')
            oLRTfile.write('Chi-square test statistic is ' + str(self.estimates.dTestStat) + '\n')
            oLRTfile.write('with P-value = ' + str(self.estimates.dPval) + '\n')
    
    def WriteModelCoefficients(self):
        if self.bNested:
            # get all parameter estimates, with indices of traits and factors
            (vIndTGA, vIndFGA, vParamGA, vIndTEA, vIndFEA, vParamEA) = self.estimates.estimatorA.mgreml_model.model.GetSplitParamsAndIndices()
            (vIndTG0, vIndFG0, vParamG0, vIndTE0, vIndFE0, vParamE0) = self.estimates.estimator0.mgreml_model.model.GetSplitParamsAndIndices()
            # get labels of the phenotypes and factors
            indPhenosA = pd.Index(self.estimates.estimatorA.mgreml_model.data.lPhenos)
            indPhenos0 = pd.Index(self.estimates.estimator0.mgreml_model.data.lPhenos)
            indFGA = pd.Index(self.estimates.estimatorA.mgreml_model.model.genmod.lFactors)
            indFG0 = pd.Index(self.estimates.estimator0.mgreml_model.model.genmod.lFactors)
            indFEA = pd.Index(self.estimates.estimatorA.mgreml_model.model.envmod.lFactors)
            indFE0 = pd.Index(self.estimates.estimator0.mgreml_model.model.envmod.lFactors)
            # for each active coefficient, find labels of phenos and factors
            indAllPhenosGA = indPhenosA[vIndTGA]
            indAllPhenosG0 = indPhenos0[vIndTG0]
            indAllPhenosEA = indPhenosA[vIndTEA]
            indAllPhenosE0 = indPhenos0[vIndTE0]
            indAllFGA = indFGA[vIndFGA]
            indAllFG0 = indFG0[vIndFG0]
            indAllFEA = indFEA[vIndFEA]
            indAllFE0 = indFE0[vIndFE0]
            # convert to lists
            lAllPhenosA = indAllPhenosGA.to_list()
            lAllPhenos0 = indAllPhenosG0.to_list()
            lAllPhenosA.extend(indAllPhenosEA.to_list())
            lAllPhenos0.extend(indAllPhenosE0.to_list())
            lAllFA = indAllFGA.to_list()
            lAllF0 = indAllFG0.to_list()
            lAllFA.extend(indAllFEA.to_list())
            lAllF0.extend(indAllFE0.to_list())
            # concatenate estimates
            vParamA = np.hstack((vParamGA, vParamEA))
            vParam0 = np.hstack((vParamG0, vParamE0))
            # construct dataframes
            dfParamsA = pd.DataFrame(vParamA, index=[lAllPhenosA,lAllFA])
            dfParams0 = pd.DataFrame(vParam0, index=[lAllPhenos0,lAllF0])
            dfSamplingVA = pd.DataFrame(self.estimates.estimatorA.mSamplingV, index=[lAllPhenosA,lAllFA], columns=[lAllPhenosA,lAllFA])
            dfSamplingV0 = pd.DataFrame(self.estimates.estimator0.mSamplingV, index=[lAllPhenos0,lAllF0], columns=[lAllPhenos0,lAllF0])
            # set output names
            sParams0 = self.sPrefix + DataWriter.sCoeff + DataWriter.sH0 + DataWriter.sExtension
            sParamsA = self.sPrefix + DataWriter.sCoeff + DataWriter.sHA + DataWriter.sExtension
            sSamplingV0 = self.sPrefix + DataWriter.sCoeffvar + DataWriter.sH0 + DataWriter.sExtension
            sSamplingVA = self.sPrefix + DataWriter.sCoeffvar + DataWriter.sHA + DataWriter.sExtension
            # write dataframes
            dfParamsA.to_csv(sParamsA)
            dfParams0.to_csv(sParams0)
            dfSamplingVA.to_csv(sSamplingVA)
            dfSamplingV0.to_csv(sSamplingV0)
        else:
            # get all parameter estimates, with indices of traits and factors
            (vIndTG, vIndFG, vParamG, vIndTE, vIndFE, vParamE) = self.estimates.mgreml_model.model.GetSplitParamsAndIndices()
            # get labels of the phenotypes and factors
            indPhenos = pd.Index(self.estimates.mgreml_model.data.lPhenos)
            indFG = pd.Index(self.estimates.mgreml_model.model.genmod.lFactors)
            indFE = pd.Index(self.estimates.mgreml_model.model.envmod.lFactors)
            # for each active coefficient, find labels of phenos and factors
            indAllPhenosG = indPhenos[vIndTG]
            indAllPhenosE = indPhenos[vIndTE]
            indAllFG = indFG[vIndFG]
            indAllFE = indFE[vIndFE]
            # convert to lists
            lAllPhenos = indAllPhenosG.to_list()
            lAllPhenos.extend(indAllPhenosE.to_list())
            lAllF = indAllFG.to_list()
            lAllF.extend(indAllFE.to_list())
            # concatenate estimates
            vParam = np.hstack((vParamG, vParamE))
            # construct dataframes
            dfParams = pd.DataFrame(vParam, index=[lAllPhenos,lAllF])
            dfSamplingV = pd.DataFrame(self.estimates.mSamplingV, index=[lAllPhenos,lAllF], columns=[lAllPhenos,lAllF])
            # set output names
            sParams = self.sPrefix + DataWriter.sCoeff + DataWriter.sExtension
            sSamplingV = self.sPrefix + DataWriter.sCoeffvar + DataWriter.sExtension
            # write dataframes
            dfParams.to_csv(sParams)
            dfSamplingV.to_csv(sSamplingV)
    
    def WriteLogLik(self):
        if self.bNested:
            # set filenames
            sLL0 = self.sPrefix + DataWriter.sLL + DataWriter.sH0 + DataWriter.sExtension
            sLLA = self.sPrefix + DataWriter.sLL + DataWriter.sHA + DataWriter.sExtension
            with open(sLL0, 'w') as oLLfile:
                oLLfile.write('Log-likelihood of nested model (null hypothesis) = ' + str(self.estimates.estimator0.dLogL) + ',\n')
                oLLfile.write('based on data on ' + str(self.estimates.estimator0.mgreml_model.data.iT) + ' traits and ' + str(self.estimates.estimator0.mgreml_model.data.iN) + ' observations,\n')
                oLLfile.write('with a model consisting of ' + str(self.estimates.estimator0.mgreml_model.model.genmod.iF) + ' genetic factors and ' + str(self.estimates.estimator0.mgreml_model.model.envmod.iF) + ' environment factors,\n')
                oLLfile.write('comprising ' + str(self.estimates.estimator0.mgreml_model.model.iParamsG) + ' free genetic factor coefficients and ' + str(self.estimates.estimator0.mgreml_model.model.iParamsE) + ' free environment factor coefficients in turn.\n')
                if self.estimates.estimator0.bBFGS:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.estimator0.iIter) + ' BFGS iterations \n')
                else:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.estimator0.iIter) + ' Newton iterations \n')
            with open(sLLA, 'w') as oLLfile:
                oLLfile.write('Log-likelihood of parent model (alternative hypothesis) = ' + str(self.estimates.estimatorA.dLogL) + ',\n')
                oLLfile.write('based on data on ' + str(self.estimates.estimatorA.mgreml_model.data.iT) + ' traits and ' + str(self.estimates.estimatorA.mgreml_model.data.iN) + ' observations,\n')
                oLLfile.write('with a model consisting of ' + str(self.estimates.estimatorA.mgreml_model.model.genmod.iF) + ' genetic factors and ' + str(self.estimates.estimatorA.mgreml_model.model.envmod.iF) + ' environment factors,\n')
                oLLfile.write('comprising ' + str(self.estimates.estimatorA.mgreml_model.model.iParamsG) + ' free genetic factor coefficients and ' + str(self.estimates.estimatorA.mgreml_model.model.iParamsE) + ' free environment factor coefficients in turn.\n')
                if self.estimates.estimatorA.bBFGS:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.estimatorA.iIter) + ' BFGS iterations \n')
                else:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.estimatorA.iIter) + ' Newton iterations \n')
        else:
            # set filenames
            sLL = self.sPrefix + DataWriter.sLL + DataWriter.sExtension
            with open(sLL, 'w') as oLLfile:
                oLLfile.write('Log-likelihood of model = ' + str(self.estimates.dLogL) + ',\n')
                oLLfile.write('based on data on ' + str(self.estimates.mgreml_model.data.iT) + ' traits and ' + str(self.estimates.mgreml_model.data.iN) + ' observations,\n')
                oLLfile.write('with a model consisting of ' + str(self.estimates.mgreml_model.model.genmod.iF) + ' genetic factors and ' + str(self.estimates.mgreml_model.model.envmod.iF) + ' environment factors,\n')
                oLLfile.write('comprising ' + str(self.estimates.mgreml_model.model.iParamsG) + ' free genetic factor coefficients and ' + str(self.estimates.mgreml_model.model.iParamsE) + ' free environment factor coefficients in turn.\n')
                if self.estimates.bBFGS:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.iIter) + ' BFGS iterations \n')
                else:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.iIter) + ' Newton iterations \n')
    
    def WriteEstimatesGLS(self):
        if self.bNested:
            # throw error if no covariates in model
            if not(self.estimates.estimatorA.mgreml_model.data.bCovs):
                raise TypeError('Trying to write GLS estimates, while no covariates were given')
            # get labels of the phenotypes and covariates
            indPhenos = pd.Index(self.estimates.estimatorA.mgreml_model.data.lPhenos)
            indCovs = pd.Index(self.estimates.estimatorA.mgreml_model.data.lCovs)
            # if not the same covariates
            if not(self.estimates.estimatorA.mgreml_model.data.bSameCovs):
                # get the binary matrix inidicating which covariate
                # applies to which trait
                mBinXY = self.estimates.estimatorA.mgreml_model.data.mBinXY
            else:
                mBinXY = np.ones((len(indPhenos),len(indCovs))).astype(int)
            # find trait and covariate indices of active covariates
            (vIndT,vIndC) = np.where(mBinXY==1)
            # for each active covariate, find labels of phenos and covs
            indAllPhenos = indPhenos[vIndT]
            indAllCovs = indCovs[vIndC]
            # get GLS estimates
            vBetaGLS0 = self.estimates.estimator0.mgreml_model.vBetaGLS
            mVarGLS0 = self.estimates.estimator0.mgreml_model.mVarGLS
            vBetaGLSA = self.estimates.estimatorA.mgreml_model.vBetaGLS
            mVarGLSA = self.estimates.estimatorA.mgreml_model.mVarGLS
            # construct dataframes
            dfBetaGLS0 = pd.DataFrame(vBetaGLS0, index=[indAllPhenos,indAllCovs])
            dfBetaGLSA = pd.DataFrame(vBetaGLSA, index=[indAllPhenos,indAllCovs])
            dfVarGLS0 = pd.DataFrame(mVarGLS0, index=[indAllPhenos,indAllCovs], columns=[indAllPhenos,indAllCovs])
            dfVarGLSA = pd.DataFrame(mVarGLSA, index=[indAllPhenos,indAllCovs], columns=[indAllPhenos,indAllCovs])
            # set output names
            sBetaGLS0 = self.sPrefix + DataWriter.sGLSest + DataWriter.sH0 + DataWriter.sExtension
            sBetaGLSA = self.sPrefix + DataWriter.sGLSest + DataWriter.sHA + DataWriter.sExtension
            sVarGLS0 = self.sPrefix + DataWriter.sGLSvar + DataWriter.sH0 + DataWriter.sExtension
            sVarGLSA = self.sPrefix + DataWriter.sGLSvar + DataWriter.sHA + DataWriter.sExtension
            # write dataframes
            dfBetaGLS0.to_csv(sBetaGLS0)
            dfBetaGLSA.to_csv(sBetaGLSA)
            dfVarGLS0.to_csv(sVarGLS0)
            dfVarGLSA.to_csv(sVarGLSA)
        else:
            # throw error if no covariates in model
            if not(self.estimates.mgreml_model.data.bCovs):
                raise TypeError('Trying to write GLS estimates, while no covariates were given')
            # get labels of the phenotypes and covariates
            indPhenos = pd.Index(self.estimates.mgreml_model.data.lPhenos)
            indCovs = pd.Index(self.estimates.mgreml_model.data.lCovs)
            # if not the same covariates
            if not(self.estimates.mgreml_model.data.bSameCovs):
                # get the binary matrix inidicating which covariate
                # applies to which trait
                mBinXY = self.estimates.mgreml_model.data.mBinXY
            else:
                mBinXY = np.ones((len(indPhenos),len(indCovs))).astype(int)
            # find trait and covariate indices of active covariates
            (vIndT,vIndC) = np.where(mBinXY==1)
            # for each active covariate, find labels of phenos and covs
            indAllPhenos = indPhenos[vIndT]
            indAllCovs = indCovs[vIndC]
            # get GLS estimates
            vBetaGLS = self.estimates.mgreml_model.vBetaGLS
            mVarGLS = self.estimates.mgreml_model.mVarGLS
            # construct dataframes
            dfBetaGLS = pd.DataFrame(vBetaGLS, index=[indAllPhenos,indAllCovs])
            dfVarGLS = pd.DataFrame(mVarGLS, index=[indAllPhenos,indAllCovs], columns=[indAllPhenos,indAllCovs])
            # set output names
            sBetaGLS = self.sPrefix + DataWriter.sGLSest + DataWriter.sExtension
            sVarGLS = self.sPrefix + DataWriter.sGLSvar + DataWriter.sExtension
            # write dataframes
            dfBetaGLS.to_csv(sBetaGLS)
            dfVarGLS.to_csv(sVarGLS)
    
    def WriteRho(self):
        if self.bNested:
            # get correlation matrices 
            mRhoG0 = self.estimates.estimator0.mRhoG
            mRhoE0 = self.estimates.estimator0.mRhoE
            mRhoGA = self.estimates.estimatorA.mRhoG
            mRhoEA = self.estimates.estimatorA.mRhoE
            # get trait labels
            lPhenos = self.estimates.estimatorA.mgreml_model.data.lPhenos
            # set dataframes
            dfRhoG0 = pd.DataFrame(mRhoG0,index=lPhenos,columns=lPhenos)
            dfRhoE0 = pd.DataFrame(mRhoE0,index=lPhenos,columns=lPhenos)
            dfRhoGA = pd.DataFrame(mRhoGA,index=lPhenos,columns=lPhenos)
            dfRhoEA = pd.DataFrame(mRhoEA,index=lPhenos,columns=lPhenos)
            # set filenames
            sRhoG0 = self.sPrefix + DataWriter.sRG + DataWriter.sH0 + DataWriter.sExtension
            sRhoE0 = self.sPrefix + DataWriter.sRE + DataWriter.sH0 + DataWriter.sExtension
            sRhoGA = self.sPrefix + DataWriter.sRG + DataWriter.sHA + DataWriter.sExtension
            sRhoEA = self.sPrefix + DataWriter.sRE + DataWriter.sHA + DataWriter.sExtension
            # write dataframes
            dfRhoG0.to_csv(sRhoG0)
            dfRhoE0.to_csv(sRhoE0)
            dfRhoGA.to_csv(sRhoGA)
            dfRhoEA.to_csv(sRhoEA)
            # if SEs are desired, store them
            if (self.estimates.estimatorA.bSEs):
                # get SE matrices 
                mRhoG0SE = self.estimates.estimator0.mRhoGSE
                mRhoE0SE = self.estimates.estimator0.mRhoESE
                mRhoGASE = self.estimates.estimatorA.mRhoGSE
                mRhoEASE = self.estimates.estimatorA.mRhoESE
                # set dataframes
                dfRhoG0SE = pd.DataFrame(mRhoG0SE,index=lPhenos,columns=lPhenos)
                dfRhoE0SE = pd.DataFrame(mRhoE0SE,index=lPhenos,columns=lPhenos)
                dfRhoGASE = pd.DataFrame(mRhoGASE,index=lPhenos,columns=lPhenos)
                dfRhoEASE = pd.DataFrame(mRhoEASE,index=lPhenos,columns=lPhenos)
                # set filenames
                sRhoG0SE = self.sPrefix + DataWriter.sRG + DataWriter.sH0 + DataWriter.sSE + DataWriter.sExtension
                sRhoE0SE = self.sPrefix + DataWriter.sRE + DataWriter.sH0 + DataWriter.sSE + DataWriter.sExtension
                sRhoGASE = self.sPrefix + DataWriter.sRG + DataWriter.sHA + DataWriter.sSE + DataWriter.sExtension
                sRhoEASE = self.sPrefix + DataWriter.sRE + DataWriter.sHA + DataWriter.sSE + DataWriter.sExtension
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
            sRhoG = self.sPrefix + DataWriter.sRG + DataWriter.sExtension
            sRhoE = self.sPrefix + DataWriter.sRE + DataWriter.sExtension
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
                sRhoGSE = self.sPrefix + DataWriter.sRG + DataWriter.sSE + DataWriter.sExtension
                sRhoESE = self.sPrefix + DataWriter.sRE + DataWriter.sSE + DataWriter.sExtension
                # write dataframes
                dfRhoGSE.to_csv(sRhoGSE)
                dfRhoESE.to_csv(sRhoESE)

