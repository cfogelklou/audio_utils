/**
  This file was autogenerated by ResampleCodeGen.m
  Author: Chris Fogelklou
  Command: ResampleCodeGen(4, 8, 60, 0, 68, 'ownsampleX4Coeffs')
*/
// Include files
#include "AudResampler.h"


// The FIR coefficients for this filter 
static const ar_coef firCoeffs[] = {
	(ar_coef)0.00005270271781855073200,
	(ar_coef)0.00017906215105847859000,
	(ar_coef)0.00033546314070629220000,
	(ar_coef)0.00035820843995879454000,
	(ar_coef)0.00010304651862362536000,
	(ar_coef)-0.00045113693295956468000,
	(ar_coef)-0.00109890529568568850000,
	(ar_coef)-0.00142304607469103320000,
	(ar_coef)-0.00098206914647798373000,
	(ar_coef)0.00037993993662728255000,
	(ar_coef)0.00226716608215259780000,
	(ar_coef)0.00372122468598593080000,
	(ar_coef)0.00359112947498611530000,
	(ar_coef)0.00120063173198812460000,
	(ar_coef)-0.00300478465159405120000,
	(ar_coef)-0.00726002552347918170000,
	(ar_coef)-0.00903226557963008290000,
	(ar_coef)-0.00625901885399085360000,
	(ar_coef)0.00123650272888142090000,
	(ar_coef)0.01093750289054846500000,
	(ar_coef)0.01810881030077349200000,
	(ar_coef)0.01769667320894969400000,
	(ar_coef)0.00703624175888454760000,
	(ar_coef)-0.01185641381818893600000,
	(ar_coef)-0.03171450596927302000000,
	(ar_coef)-0.04188763961561887300000,
	(ar_coef)-0.03228085874429514200000,
	(ar_coef)0.00208605093976066850000,
	(ar_coef)0.05788890466399749400000,
	(ar_coef)0.12331424224628441000000,
	(ar_coef)0.18133795124189722000000,
	(ar_coef)0.21541921534600122000000,
	(ar_coef)0.21541921534600122000000,
	(ar_coef)0.18133795124189722000000,
	(ar_coef)0.12331424224628441000000,
	(ar_coef)0.05788890466399749400000,
	(ar_coef)0.00208605093976066850000,
	(ar_coef)-0.03228085874429514200000,
	(ar_coef)-0.04188763961561887300000,
	(ar_coef)-0.03171450596927302000000,
	(ar_coef)-0.01185641381818893600000,
	(ar_coef)0.00703624175888454760000,
	(ar_coef)0.01769667320894969400000,
	(ar_coef)0.01810881030077349200000,
	(ar_coef)0.01093750289054846500000,
	(ar_coef)0.00123650272888142090000,
	(ar_coef)-0.00625901885399085360000,
	(ar_coef)-0.00903226557963008290000,
	(ar_coef)-0.00726002552347918170000,
	(ar_coef)-0.00300478465159405120000,
	(ar_coef)0.00120063173198812460000,
	(ar_coef)0.00359112947498611530000,
	(ar_coef)0.00372122468598593080000,
	(ar_coef)0.00226716608215259780000,
	(ar_coef)0.00037993993662728255000,
	(ar_coef)-0.00098206914647798373000,
	(ar_coef)-0.00142304607469103320000,
	(ar_coef)-0.00109890529568568850000,
	(ar_coef)-0.00045113693295956468000,
	(ar_coef)0.00010304651862362536000,
	(ar_coef)0.00035820843995879454000,
	(ar_coef)0.00033546314070629220000,
	(ar_coef)0.00017906215105847859000,
	(ar_coef)0.00005270271781855073200,
	0
};

// The resampler configuration struct
const AudResampleConfig DownsampleX4CoeffsCfg = {
  (sizeof(firCoeffs)/sizeof(ar_coef)-1),
  4,
  -1,
  firCoeffs
};