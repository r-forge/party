
/**
    S4 classes 
    *\file $RCSfile$
    *\author $Author$
    *\date $Date$
*/

#include "PL2_common.h"

SEXP 
    PL2_expectationSym,
    PL2_covarianceSym,
    PL2_linearstatisticSym,
    PL2_expcovinfSym,
    PL2_expcovinfssSym,
    PL2_sumweightsSym,
    PL2_dimensionSym,
    PL2_MPinvSym,  
    PL2_rankSym,
    PL2_svdSym,
    PL2_svdmemSym,
    PL2_methodSym,
    PL2_jobuSym, 
    PL2_jobvSym, 
    PL2_uSym,
    PL2_vSym,
    PL2_sSym,
    PL2_pSym,
    PL2_teststattypeSym,
    PL2_pvalueSym,
    PL2_tolSym,
    PL2_maxptsSym,
    PL2_absepsSym,
    PL2_relepsSym,
    PL2_minprobSym,
    PL2_minsplitSym,
    PL2_minbucketSym,
    PL2_variablesSym, 
    PL2_transformationsSym, 
    PL2_is_nominalSym, 
    PL2_is_ordinalSym, 
    PL2_is_censoredSym, 
    PL2_orderingSym, 
    PL2_levelsSym, 
    PL2_scoresSym, 
    PL2_has_missingsSym, 
    PL2_whichNASym, 
    PL2_jointtransfSym, 
    PL2_nobsSym, 
    PL2_ninputsSym,
    PL2_linexpcov2sampleSym, 
    PL2_weightsSym, 
    PL2_varmemorySym,
    PL2_varMmemorySym, 
    PL2_MscorematricesSym,
    PL2_responsesSym, 
    PL2_inputsSym,
    PL2_testtypeSym, 
    PL2_nresampleSym,
    PL2_varctrlSym, 
    PL2_splitctrlSym, 
    PL2_gtctrlSym,
    PL2_mincriterionSym,
    PL2_maxsurrogateSym,
    PL2_randomsplitsSym,
    PL2_mtrySym,
    PL2_dontuseSym,
    PL2_dontusetmpSym,
    PL2_stumpSym,
    PL2_tgctrlSym;

SEXP party_init(void) {
    PL2_expectationSym = install("expectation");
    PL2_covarianceSym = install("covariance");
    PL2_linearstatisticSym = install("linearstatistic");
    PL2_expcovinfSym = install("expcovinf");
    PL2_expcovinfssSym = install("expcovinfss");
    PL2_sumweightsSym = install("sumweights");
    PL2_dimensionSym = install("dimension");
    PL2_MPinvSym = install("MPinv");
    PL2_rankSym = install("rank");
    PL2_svdSym = install("svd");
    PL2_svdmemSym = install("svdmem");
    PL2_methodSym = install("method");
    PL2_jobuSym = install("jobu");
    PL2_jobvSym = install("jobv");
    PL2_uSym = install("u");
    PL2_vSym = install("v");
    PL2_sSym = install("s");
    PL2_pSym = install("p"); 
    PL2_teststattypeSym = install("teststattype");
    PL2_pvalueSym = install("pvalue");
    PL2_tolSym = install("tol");
    PL2_maxptsSym = install("maxpts");
    PL2_absepsSym = install("abseps");
    PL2_relepsSym = install("releps");
    PL2_minsplitSym = install("minsplit");
    PL2_minprobSym = install("minprob");
    PL2_minbucketSym = install("minbucket");
    PL2_variablesSym = install("variables"); 
    PL2_transformationsSym = install("transformations"); 
    PL2_is_nominalSym = install("is_nominal"); 
    PL2_is_ordinalSym = install("is_ordinal"); 
    PL2_is_censoredSym = install("is_censored"); 
    PL2_orderingSym = install("ordering"); 
    PL2_levelsSym = install("levels"); 
    PL2_scoresSym = install("scores"); 
    PL2_has_missingsSym = install("has_missings"); 
    PL2_whichNASym = install("whichNA"); 
    PL2_jointtransfSym = install("jointtransf"); 
    PL2_nobsSym = install("nobs"); 
    PL2_ninputsSym = install("ninputs"); 
    PL2_linexpcov2sampleSym = install("linexpcov2sample"); 
    PL2_weightsSym = install("weights"); 
    PL2_varmemorySym = install("varmemory"); 
    PL2_varMmemorySym = install("varMmemory"); 
    PL2_MscorematricesSym = install("Mscorematrices"); 
    PL2_responsesSym = install("responses"); 
    PL2_inputsSym = install("inputs"); 
    PL2_testtypeSym = install("testtype"); 
    PL2_nresampleSym = install("nresample"); 
    PL2_varctrlSym = install("varctrl"); 
    PL2_splitctrlSym = install("splitctrl"); 
    PL2_gtctrlSym = install("gtctrl"); 
    PL2_mincriterionSym = install("mincriterion"); 
    PL2_maxsurrogateSym = install("maxsurrogate"); 
    PL2_randomsplitsSym = install("randomsplits"); 
    PL2_mtrySym = install("mtry"); 
    PL2_dontuseSym = install("dontuse"); 
    PL2_dontusetmpSym = install("dontusetmp"); 
    PL2_stumpSym = install("stump"); 
    PL2_tgctrlSym = install("tgctrl"); 
    return(R_NilValue);
}

int get_dimension(SEXP object) {
    return(INTEGER(GET_SLOT(object, PL2_dimensionSym))[0]);
}

int get_teststattype(SEXP object) {
    return(INTEGER(GET_SLOT(object, PL2_teststattypeSym))[0]);
}

int get_pvalue(SEXP object) {
    return(INTEGER(GET_SLOT(object, PL2_pvalueSym))[0]);
}

double get_tol(SEXP object) {
    return(REAL(GET_SLOT(object, PL2_tolSym))[0]);
}

int get_maxpts(SEXP object) {
    return(INTEGER(GET_SLOT(object, PL2_maxptsSym))[0]);
}

double get_abseps(SEXP object) {
    return(REAL(GET_SLOT(object, PL2_absepsSym))[0]);
}

double get_releps(SEXP object) {
    return(REAL(GET_SLOT(object, PL2_relepsSym))[0]);
}

double get_minsplit(SEXP object) {
    return(REAL(GET_SLOT(object, PL2_minsplitSym))[0]);
}

double get_minprob(SEXP object) {
    return(REAL(GET_SLOT(object, PL2_minprobSym))[0]);
}

double get_minbucket(SEXP object) {
    return(REAL(GET_SLOT(object, PL2_minbucketSym))[0]);
}

SEXP get_transformation(SEXP object, int variable) {
    return(VECTOR_ELT(
               GET_SLOT(object, PL2_transformationsSym), 
               variable - 1));
}

SEXP get_jointtransf(SEXP object) {
    return(GET_SLOT(object, PL2_jointtransfSym));
}

SEXP get_variable(SEXP object, int variable) {
    return(VECTOR_ELT(
               GET_SLOT(object, PL2_variablesSym), 
               variable - 1));
}

int is_nominal(SEXP object, int variable) {
    return(INTEGER(GET_SLOT(object, PL2_is_nominalSym))[variable - 1]);
}

int is_ordinal(SEXP object, int variable) {
    return(INTEGER(GET_SLOT(object, PL2_is_ordinalSym))[variable - 1]);
}

int is_censored(SEXP object, int variable) {
    return(INTEGER(GET_SLOT(object, PL2_is_censoredSym))[variable - 1]);
}

int has_missings(SEXP object, int variable) {
    return(INTEGER(GET_SLOT(object, PL2_has_missingsSym))[variable - 1]);
}

SEXP get_ordering(SEXP object, int variable) {
    if (!is_nominal(object, variable)) {
        return(VECTOR_ELT(
               GET_SLOT(object, PL2_orderingSym), 
               variable - 1));
    } else {
        error("Variable %d is not ordered", variable);
        return(R_NilValue);
    }
}

SEXP get_levels(SEXP object, int variable) {
    if (is_nominal(object, variable) || 
        is_ordinal(object, variable)) {
        return(VECTOR_ELT(
               GET_SLOT(object, PL2_levelsSym), 
               variable - 1));
    } else {
        error("Variable %d is not an (ordered) factor", variable);
        return(R_NilValue);
    }
}

SEXP get_scores(SEXP object, int variable) {
    if (is_ordinal(object, variable)) {
        return(VECTOR_ELT(
               GET_SLOT(object, PL2_scoresSym), 
               variable - 1));
    } else {
        error("Variable %d is not an ordered factor", variable);
        return(R_NilValue);
    }
}

SEXP get_missings(SEXP object, int variable) {
    if (has_missings(object, variable)) {
        return(VECTOR_ELT(
               GET_SLOT(object, PL2_whichNASym), 
               variable - 1));
    } else {
        error("Variable %d has no missing values", variable);
        return(R_NilValue);
    }
}

SEXP get_varmemory(SEXP object, int variable) {
    return(VECTOR_ELT(GET_SLOT(object, PL2_varmemorySym), 
                      variable - 1));
}

SEXP get_varMmemory(SEXP object, int variable) {
    return(VECTOR_ELT(GET_SLOT(object, PL2_varMmemorySym), 
                      variable - 1));
}

SEXP get_Mscorematrix(SEXP object, int variable) {
    return(VECTOR_ELT(GET_SLOT(object, PL2_MscorematricesSym), 
                      variable - 1));
}

int get_nobs(SEXP object) {
    return(INTEGER(GET_SLOT(object, PL2_nobsSym))[0]);
}

int get_ninputs(SEXP object) {
    return(INTEGER(GET_SLOT(object, PL2_ninputsSym))[0]);
}

SEXP get_weights(SEXP object, int variable) {
    return(VECTOR_ELT(GET_SLOT(object, PL2_weightsSym), variable - 1));
}

int get_testtype(SEXP object) {
    return(INTEGER(GET_SLOT(object, PL2_testtypeSym))[0]);
}

int get_nresample(SEXP object) {
    return(INTEGER(GET_SLOT(object, PL2_nresampleSym))[0]);
}

SEXP get_varctrl(SEXP object) {
    return(GET_SLOT(object, PL2_varctrlSym));
}

SEXP get_splitctrl(SEXP object) {
    return(GET_SLOT(object, PL2_splitctrlSym));
}

SEXP get_gtctrl(SEXP object) {
    return(GET_SLOT(object, PL2_gtctrlSym));
}

double get_mincriterion(SEXP object) {
    return(REAL(GET_SLOT(object, PL2_mincriterionSym))[0]);
}

int get_maxsurrogate(SEXP object) {
    return(INTEGER(GET_SLOT(object, PL2_maxsurrogateSym))[0]);
}

int get_randomsplits(SEXP object) {
    return(INTEGER(GET_SLOT(object, PL2_randomsplitsSym))[0]);
}

int get_mtry(SEXP object) {
    return(INTEGER(GET_SLOT(object, PL2_mtrySym))[0]);
}

SEXP get_dontuse(SEXP object) {
    return(GET_SLOT(object, PL2_dontuseSym));
}

SEXP get_dontusetmp(SEXP object) {
    return(GET_SLOT(object, PL2_dontusetmpSym));
}

int get_stump(SEXP object) {
    return(INTEGER(GET_SLOT(object, PL2_stumpSym))[0]);
}
