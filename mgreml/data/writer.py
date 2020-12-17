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
    sGLSest = 'GLS.est.'
    sGLSvar = 'GLS.var.'
    sCoeff = 'coeff.'
    sCoeffvar = 'coeff.var.'
    lHsqSE = ['heritability', 'standard error']
    lHsq = ['heritability']
    lBetaSE = ['beta hat', 'standard error']
    
    def __init__(self, estimates, data):
        self.logger = data.logger
        self.sPrefix = data.sPrefix
        self.bNested = data.bNested
        self.bSEs = data.bSEs
        self.bAllCoeffs = data.bAllCoeffs
        self.bCovs = data.bCovs
        self.estimates = estimates
        self.vScaleY = data.vScaleY
        if self.bCovs:
            self.vScaleX = data.vScaleX
    
    def WriteResults(self):
        self.logger.info('6. WRITING MGREML RESULTS')
        if not(self.estimates.IsConverged()):
            raise RuntimeError('Trying to write output results while estimates have not converged')
        if not(self.estimates.IsDone()):
            raise RuntimeError('Trying to write output results while final statistic (e.g. genetic correlations) have not been calculated')
        self.logger.info('WRITING HERITABILITIES')
        self.WriteHSq()
        self.logger.info('WRITING GENETIC AND ENVIRONMENT CORRELATIONS')
        self.WriteRho()
        self.logger.info('WRITING LOG-LIKELIHOOD AND OTHER DESCRIPTIVES')
        self.WriteLogLik()
        if self.bCovs:
            self.logger.info('WRITING GLS ESTIMATES OF FIXED EFFECTS')
            self.WriteEstimatesGLS()
        if self.bAllCoeffs:
            self.logger.info('WRITING ALL FACTOR COEFFICIENTS')
            self.WriteModelCoefficients()
        if self.bNested:
            self.logger.info('WRITING RESULTS LIKELIHOOD-RATIO TEST')
            self.WriteLRT()
        self.logger.info('Done writing results\n')
       
    def WriteHSq(self):
        if self.bNested:
            # get heritabilities 
            vHSq0 = self.estimates.estimator0.vHSq
            vHSqA = self.estimates.estimatorA.vHSq
            # get trait labels
            lPhenos = self.estimates.estimatorA.mgreml_model.data.lPhenos
            # set filenames
            sHSq0 = self.sPrefix + DataWriter.sH2 + DataWriter.sH0 + DataWriter.sExtension
            sHSqA = self.sPrefix + DataWriter.sH2 + DataWriter.sHA + DataWriter.sExtension
            # if SEs are desired, store them together with heritabilities
            if (self.bSEs):
                # get heritability standard errors 
                vHSq0SE = self.estimates.estimator0.vHSqSE
                vHSqASE = self.estimates.estimatorA.vHSqSE
                # set dataframes
                dfHSq0 = pd.DataFrame(np.stack((vHSq0,vHSq0SE)).T,index=lPhenos,columns=DataWriter.lHsqSE)
                dfHSqA = pd.DataFrame(np.stack((vHSqA,vHSqASE)).T,index=lPhenos,columns=DataWriter.lHsqSE)
            else:
                # set dataframes
                dfHSq0 = pd.DataFrame(vHSq0,index=lPhenos,columns=DataWriter.lHsq)
                dfHSqA = pd.DataFrame(vHSqA,index=lPhenos,columns=DataWriter.lHsq)
            # write dataframes
            dfHSq0.to_csv(sHSq0, sep='\t')
            dfHSqA.to_csv(sHSqA, sep='\t')
        else:
            # get heritability
            vHSq = self.estimates.vHSq
            # get trait labels
            lPhenos = self.estimates.mgreml_model.data.lPhenos
            # set filename
            sHSq = self.sPrefix + DataWriter.sH2 + DataWriter.sExtension
            # if SEs are desired, store them together with heritabilities
            if (self.bSEs):
                # get heritability standard errors 
                vHSqSE = self.estimates.vHSqSE
                # set dataframe
                dfHSq = pd.DataFrame(np.stack((vHSq,vHSqSE)).T,index=lPhenos,columns=DataWriter.lHsqSE)
            else:
                # set dataframe
                dfHSq = pd.DataFrame(vHSq,index=lPhenos,columns=DataWriter.lHsq)
            # write dataframe
            dfHSq.to_csv(sHSq, sep='\t')

    def WriteLRT(self):
        if not(self.bNested):
            raise TypeError('Trying to write results for likelihood-ratio test, while no nested model has been estimated')
        # set filename
        sLL = self.sPrefix + DataWriter.sLL + DataWriter.sExtension
        with open(sLL, 'a') as oLRTfile:
            oLRTfile.write('Results of likelihood-ratio test with ' + str(self.estimates.iDF) + ' degrees of freedom:\n')
            oLRTfile.write('Chi-square test statistic is ' + str(self.estimates.dTestStat) + '\n')
            oLRTfile.write('with P-value = ' + str(self.estimates.dPval) + '\n\n')
    
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
            dfParamsA.to_csv(sParamsA, sep='\t')
            dfParams0.to_csv(sParams0, sep='\t')
            dfSamplingVA.to_csv(sSamplingVA, sep='\t')
            dfSamplingV0.to_csv(sSamplingV0, sep='\t')
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
            dfParams.to_csv(sParams, sep='\t')
            dfSamplingV.to_csv(sSamplingV, sep='\t')
    
    def WriteLogLik(self):
        if self.bNested:
            # set filename
            sLL = self.sPrefix + DataWriter.sLL + DataWriter.sExtension
            with open(sLL, 'w') as oLLfile:
                oLLfile.write('Log-likelihood of nested model (null hypothesis) = ' + str(self.estimates.estimator0.dLogL) + ',\n')
                oLLfile.write('based on data on ' + str(self.estimates.estimator0.mgreml_model.data.iT) + ' traits and ' + str(self.estimates.estimator0.mgreml_model.data.iN) + ' observations,\n')
                oLLfile.write('with a model consisting of ' + str(self.estimates.estimator0.mgreml_model.model.genmod.iF) + ' genetic factors and ' + str(self.estimates.estimator0.mgreml_model.model.envmod.iF) + ' environment factors,\n')
                oLLfile.write('comprising ' + str(self.estimates.estimator0.mgreml_model.model.iParamsG) + ' free genetic factor coefficients and ' + str(self.estimates.estimator0.mgreml_model.model.iParamsE) + ' free environment factor coefficients in turn.\n')
                if self.estimates.estimator0.bBFGS:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.estimator0.iIter) + ' BFGS iterations \n\n')
                else:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.estimator0.iIter) + ' Newton iterations \n\n')
            with open(sLL, 'a') as oLLfile:
                oLLfile.write('Log-likelihood of parent model (alternative hypothesis) = ' + str(self.estimates.estimatorA.dLogL) + ',\n')
                oLLfile.write('based on data on ' + str(self.estimates.estimatorA.mgreml_model.data.iT) + ' traits and ' + str(self.estimates.estimatorA.mgreml_model.data.iN) + ' observations,\n')
                oLLfile.write('with a model consisting of ' + str(self.estimates.estimatorA.mgreml_model.model.genmod.iF) + ' genetic factors and ' + str(self.estimates.estimatorA.mgreml_model.model.envmod.iF) + ' environment factors,\n')
                oLLfile.write('comprising ' + str(self.estimates.estimatorA.mgreml_model.model.iParamsG) + ' free genetic factor coefficients and ' + str(self.estimates.estimatorA.mgreml_model.model.iParamsE) + ' free environment factor coefficients in turn.\n')
                if self.estimates.estimatorA.bBFGS:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.estimatorA.iIter) + ' BFGS iterations \n\n')
                else:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.estimatorA.iIter) + ' Newton iterations \n\n')
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
            vBetaGLS0SE = np.power(np.diag(mVarGLS0),0.5)
            vBetaGLSA = self.estimates.estimatorA.mgreml_model.vBetaGLS
            mVarGLSA = self.estimates.estimatorA.mgreml_model.mVarGLS
            vBetaGLSASE = np.power(np.diag(mVarGLSA),0.5)
            # construct dataframes
            dfBetaGLS0 = pd.DataFrame(np.stack((vBetaGLS0,vBetaGLS0SE)).T, index=[indAllPhenos,indAllCovs], columns=DataWriter.lBetaSE)
            dfBetaGLSA = pd.DataFrame(np.stack((vBetaGLSA,vBetaGLSASE)).T, index=[indAllPhenos,indAllCovs], columns=DataWriter.lBetaSE)
            dfVarGLS0 = pd.DataFrame(mVarGLS0, index=[indAllPhenos,indAllCovs], columns=[indAllPhenos,indAllCovs])
            dfVarGLSA = pd.DataFrame(mVarGLSA, index=[indAllPhenos,indAllCovs], columns=[indAllPhenos,indAllCovs])
            # set output names
            sBetaGLS0 = self.sPrefix + DataWriter.sGLSest + DataWriter.sH0 + DataWriter.sExtension
            sBetaGLSA = self.sPrefix + DataWriter.sGLSest + DataWriter.sHA + DataWriter.sExtension
            sVarGLS0 = self.sPrefix + DataWriter.sGLSvar + DataWriter.sH0 + DataWriter.sExtension
            sVarGLSA = self.sPrefix + DataWriter.sGLSvar + DataWriter.sHA + DataWriter.sExtension
            # write dataframes
            dfBetaGLS0.to_csv(sBetaGLS0, sep='\t')
            dfBetaGLSA.to_csv(sBetaGLSA, sep='\t')
            dfVarGLS0.to_csv(sVarGLS0, sep='\t')
            dfVarGLSA.to_csv(sVarGLSA, sep='\t')
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
            vBetaGLSSE = np.power(np.diag(mVarGLS),0.5)
            # construct dataframes
            dfBetaGLS = pd.DataFrame(np.stack((vBetaGLS,vBetaGLSSE)).T, index=[indAllPhenos,indAllCovs], columns=DataWriter.lBetaSE)
            dfVarGLS = pd.DataFrame(mVarGLS, index=[indAllPhenos,indAllCovs], columns=[indAllPhenos,indAllCovs])
            # set output names
            sBetaGLS = self.sPrefix + DataWriter.sGLSest + DataWriter.sExtension
            sVarGLS = self.sPrefix + DataWriter.sGLSvar + DataWriter.sExtension
            # write dataframes
            dfBetaGLS.to_csv(sBetaGLS, sep='\t')
            dfVarGLS.to_csv(sVarGLS, sep='\t')
    
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
            dfRhoG0.to_csv(sRhoG0, sep='\t')
            dfRhoE0.to_csv(sRhoE0, sep='\t')
            dfRhoGA.to_csv(sRhoGA, sep='\t')
            dfRhoEA.to_csv(sRhoEA, sep='\t')
            # if SEs are desired, store them
            if (self.bSEs):
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
                dfRhoG0SE.to_csv(sRhoG0SE, sep='\t')
                dfRhoE0SE.to_csv(sRhoE0SE, sep='\t')
                dfRhoGASE.to_csv(sRhoGASE, sep='\t')
                dfRhoEASE.to_csv(sRhoEASE, sep='\t')
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
            dfRhoG.to_csv(sRhoG, sep='\t')
            dfRhoE.to_csv(sRhoE, sep='\t')
            # if SEs are desired, store them
            if (self.bSEs):
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
                dfRhoGSE.to_csv(sRhoGSE, sep='\t')
                dfRhoESE.to_csv(sRhoESE, sep='\t')

