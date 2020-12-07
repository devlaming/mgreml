import pandas as pd
import numpy as np
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

    def WriteLRT(self):
        if not(self.bNested):
            raise TypeError('Trying to write results for likelihood-ratio test, while no nested model has been estimated')
        # set filename
        sLRTfile = DataWriter.sPath + DataWriter.sLRT + self.sPrefix + DataWriter.sExtension
        with open(sLRTfile, 'w') as oLRTfile:
            oLRTfile.write('Results of likelihood-ratio test with ' + str(self.estimates.iDF) + ' degrees of freedom:\n')
            oLRTfile.write('Chi-square test statistic is ' + str(self.estimates.dTestStat) + '\n')
            oLRTfile.write('with P-value = ' + str(self.estimates.dPval) + '\n')
    
    def WriteModelCoefficients(self):
        if self.bNested:
            # get all parameter estimates, with indices of traits and factors
            (vIndTG_unres, vIndFG_unres, vParamG_unres, vIndTE_unres, vIndFE_unres, vParamE_unres) = self.estimates.estimator_unres.mgreml_model.model.GetSplitParamsAndIndices()
            (vIndTG_res, vIndFG_res, vParamG_res, vIndTE_res, vIndFE_res, vParamE_res) = self.estimates.estimator_res.mgreml_model.model.GetSplitParamsAndIndices()
            # get labels of the phenotypes and factors
            indPhenos_unres = pd.Index(self.estimates.estimator_unres.mgreml_model.data.lPhenos)
            indPhenos_res = pd.Index(self.estimates.estimator_res.mgreml_model.data.lPhenos)
            indFG_unres = pd.Index(self.estimates.estimator_unres.mgreml_model.model.genmod.lFactors)
            indFG_res = pd.Index(self.estimates.estimator_res.mgreml_model.model.genmod.lFactors)
            indFE_unres = pd.Index(self.estimates.estimator_unres.mgreml_model.model.envmod.lFactors)
            indFE_res = pd.Index(self.estimates.estimator_res.mgreml_model.model.envmod.lFactors)
            # for each active coefficient, find labels of phenos and factors
            indAllPhenosG_unres = indPhenos_unres[vIndTG_unres]
            indAllPhenosG_res = indPhenos_res[vIndTG_res]
            indAllPhenosE_unres = indPhenos_unres[vIndTE_unres]
            indAllPhenosE_res = indPhenos_res[vIndTE_res]
            indAllFG_unres = indFG_unres[vIndFG_unres]
            indAllFG_res = indFG_res[vIndFG_res]
            indAllFE_unres = indFE_unres[vIndFE_unres]
            indAllFE_res = indFE_res[vIndFE_res]
            # convert to lists
            lAllPhenos_unres = indAllPhenosG_unres.to_list()
            lAllPhenos_res = indAllPhenosG_res.to_list()
            lAllPhenos_unres.extend(indAllPhenosE_unres.to_list())
            lAllPhenos_res.extend(indAllPhenosE_res.to_list())
            lAllF_unres = indAllFG_unres.to_list()
            lAllF_res = indAllFG_res.to_list()
            lAllF_unres.extend(indAllFE_unres.to_list())
            lAllF_res.extend(indAllFE_res.to_list())
            # concatenate estimates
            vParam_unres = np.hstack((vParamG_unres, vParamE_unres))
            vParam_res = np.hstack((vParamG_res, vParamE_res))
            # construct dataframes
            dfParams_unres = pd.DataFrame(vParam_unres, index=[lAllPhenos_unres,lAllF_unres])
            dfParams_res = pd.DataFrame(vParam_res, index=[lAllPhenos_res,lAllF_res])
            dfSamplingV_unres = pd.DataFrame(self.estimates.estimator_unres.mSamplingV, index=[lAllPhenos_unres,lAllF_unres], columns=[lAllPhenos_unres,lAllF_unres])
            dfSamplingV_res = pd.DataFrame(self.estimates.estimator_res.mSamplingV, index=[lAllPhenos_res,lAllF_res], columns=[lAllPhenos_res,lAllF_res])
            # set output names
            sParams0 = DataWriter.sPath + DataWriter.sCoeff + DataWriter.sH0 + self.sPrefix + DataWriter.sExtension
            sParamsA = DataWriter.sPath + DataWriter.sCoeff + DataWriter.sHA + self.sPrefix + DataWriter.sExtension
            sSamplingV0 = DataWriter.sPath + DataWriter.sCoeffvar + DataWriter.sH0 + self.sPrefix + DataWriter.sExtension
            sSamplingVA = DataWriter.sPath + DataWriter.sCoeffvar + DataWriter.sHA + self.sPrefix + DataWriter.sExtension
            # write dataframes
            dfParams_unres.to_csv(sParamsA)
            dfParams_res.to_csv(sParams0)
            dfSamplingV_unres.to_csv(sSamplingVA)
            dfSamplingV_res.to_csv(sSamplingV0)
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
            sParams = DataWriter.sPath + DataWriter.sCoeff + self.sPrefix + DataWriter.sExtension
            sSamplingV = DataWriter.sPath + DataWriter.sCoeffvar + self.sPrefix + DataWriter.sExtension
            # write dataframes
            dfParams.to_csv(sParams)
            dfSamplingV.to_csv(sSamplingV)
    
    def WriteLogLik(self):
        if self.bNested:
            # set filenames
            sLL0 = DataWriter.sPath + DataWriter.sLL + DataWriter.sH0 + self.sPrefix + DataWriter.sExtension
            sLLA = DataWriter.sPath + DataWriter.sLL + DataWriter.sHA + self.sPrefix + DataWriter.sExtension
            with open(sLL0, 'w') as oLLfile:
                oLLfile.write('Log-likelihood of nested model (null hypothesis) = ' + str(self.estimates.estimator_res.dLogL) + ',\n')
                oLLfile.write('based on data on ' + str(self.estimates.estimator_res.mgreml_model.data.iT) + ' traits and ' + str(self.estimates.estimator_res.mgreml_model.data.iN) + ' observations,\n')
                oLLfile.write('with a model consisting of ' + str(self.estimates.estimator_res.mgreml_model.model.genmod.iF) + ' genetic factors and ' + str(self.estimates.estimator_res.mgreml_model.model.envmod.iF) + ' environment factors,\n')
                oLLfile.write('comprising ' + str(self.estimates.estimator_res.mgreml_model.model.iParamsG) + ' free genetic factor coefficients and ' + str(self.estimates.estimator_res.mgreml_model.model.iParamsE) + ' free environment factor coefficients in turn.\n')
                if self.estimates.estimator_res.bBFGS:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.estimator_res.iIter) + ' BFGS iterations \n')
                else:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.estimator_res.iIter) + ' Newton iterations \n')
            with open(sLLA, 'w') as oLLfile:
                oLLfile.write('Log-likelihood of parent model (alternative hypothesis) = ' + str(self.estimates.estimator_unres.dLogL) + ',\n')
                oLLfile.write('based on data on ' + str(self.estimates.estimator_unres.mgreml_model.data.iT) + ' traits and ' + str(self.estimates.estimator_unres.mgreml_model.data.iN) + ' observations,\n')
                oLLfile.write('with a model consisting of ' + str(self.estimates.estimator_unres.mgreml_model.model.genmod.iF) + ' genetic factors and ' + str(self.estimates.estimator_unres.mgreml_model.model.envmod.iF) + ' environment factors,\n')
                oLLfile.write('comprising ' + str(self.estimates.estimator_unres.mgreml_model.model.iParamsG) + ' free genetic factor coefficients and ' + str(self.estimates.estimator_unres.mgreml_model.model.iParamsE) + ' free environment factor coefficients in turn.\n')
                if self.estimates.estimator_unres.bBFGS:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.estimator_unres.iIter) + ' BFGS iterations \n')
                else:
                    oLLfile.write('Estimates converged after ' + str(self.estimates.estimator_unres.iIter) + ' Newton iterations \n')
        else:
            # set filenames
            sLL = DataWriter.sPath + DataWriter.sLL + self.sPrefix + DataWriter.sExtension
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
            if not(self.estimates.estimator_unres.mgreml_model.data.bCovs):
                raise TypeError('Trying to write GLS estimates, while no covariates were given')
            # get labels of the phenotypes and covariates
            indPhenos = pd.Index(self.estimates.estimator_unres.mgreml_model.data.lPhenos)
            indCovs = pd.Index(self.estimates.estimator_unres.mgreml_model.data.lCovs)
            # if not the same covariates
            if not(self.estimates.estimator_unres.mgreml_model.data.bSameCovs):
                # get the binary matrix inidicating which covariate
                # applies to which trait
                mBinXY = self.estimates.estimator_unres.mgreml_model.data.mBinXY
            else:
                mBinXY = np.ones((len(indPhenos),len(indCovs))).astype(int)
            # find trait and covariate indices of active covariates
            (vIndT,vIndC) = np.where(mBinXY==1)
            # for each active covariate, find labels of phenos and covs
            indAllPhenos = indPhenos[vIndT]
            indAllCovs = indCovs[vIndC]
            # get GLS estimates
            vBetaGLS0 = self.estimates.estimator_res.mgreml_model.vBetaGLS
            mVarGLS0 = self.estimates.estimator_res.mgreml_model.mVarGLS
            vBetaGLSA = self.estimates.estimator_unres.mgreml_model.vBetaGLS
            mVarGLSA = self.estimates.estimator_unres.mgreml_model.mVarGLS
            # construct dataframes
            dfBetaGLS0 = pd.DataFrame(vBetaGLS0, index=[indAllPhenos,indAllCovs])
            dfBetaGLSA = pd.DataFrame(vBetaGLSA, index=[indAllPhenos,indAllCovs])
            dfVarGLS0 = pd.DataFrame(mVarGLS0, index=[indAllPhenos,indAllCovs], columns=[indAllPhenos,indAllCovs])
            dfVarGLSA = pd.DataFrame(mVarGLSA, index=[indAllPhenos,indAllCovs], columns=[indAllPhenos,indAllCovs])
            # set output names
            sBetaGLS0 = DataWriter.sPath + DataWriter.sGLSest + DataWriter.sH0 + self.sPrefix + DataWriter.sExtension
            sBetaGLSA = DataWriter.sPath + DataWriter.sGLSest + DataWriter.sHA + self.sPrefix + DataWriter.sExtension
            sVarGLS0 = DataWriter.sPath + DataWriter.sGLSvar + DataWriter.sH0 + self.sPrefix + DataWriter.sExtension
            sVarGLSA = DataWriter.sPath + DataWriter.sGLSvar + DataWriter.sHA + self.sPrefix + DataWriter.sExtension
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
            sBetaGLS = DataWriter.sPath + DataWriter.sGLSest + self.sPrefix + DataWriter.sExtension
            sVarGLS = DataWriter.sPath + DataWriter.sGLSvar + self.sPrefix + DataWriter.sExtension
            # write dataframes
            dfBetaGLS.to_csv(sBetaGLS)
            dfVarGLS.to_csv(sVarGLS)
    
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

