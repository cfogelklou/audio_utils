/**
  This file was autogenerated by ResampleCodeGen.m
  Author: Chris Fogelklou
  Command: ResampleCodeGen(10, 20, 60, 0, 68, 'ownsampleX10Coeffs')
*/
// Include files
#include "aud_resampler.h"


// The FIR coefficients for this filter 
static const ar_coef firCoeffs[] = {
	(ar_coef)0.00011884022274186547000,
	(ar_coef)0.00014177620987582376000,
	(ar_coef)0.00024850129218651672000,
	(ar_coef)0.00038431473142678716000,
	(ar_coef)0.00053799024989214996000,
	(ar_coef)0.00068863819850353001000,
	(ar_coef)0.00080519473724929207000,
	(ar_coef)0.00084745923215067841000,
	(ar_coef)0.00076900780323513871000,
	(ar_coef)0.00052210591669823178000,
	(ar_coef)0.00006447169166240871600,
	(ar_coef)-0.00063256715714818368000,
	(ar_coef)-0.00157528660675491580000,
	(ar_coef)-0.00273917269803713920000,
	(ar_coef)-0.00406210895352451890000,
	(ar_coef)-0.00544036464115712180000,
	(ar_coef)-0.00672847637631980230000,
	(ar_coef)-0.00774382230853019790000,
	(ar_coef)-0.00827624463227718750000,
	(ar_coef)-0.00810252096270313790000,
	(ar_coef)-0.00700487680048583750000,
	(ar_coef)-0.00479214237535823420000,
	(ar_coef)-0.00132166536495993110000,
	(ar_coef)0.00348023083652750430000,
	(ar_coef)0.00960154794730574270000,
	(ar_coef)0.01693381101612880000000,
	(ar_coef)0.02526887281709741900000,
	(ar_coef)0.03430437135864797300000,
	(ar_coef)0.04365805551438314800000,
	(ar_coef)0.05289027931318632900000,
	(ar_coef)0.06153307130028077400000,
	(ar_coef)0.06912341735634533500000,
	(ar_coef)0.07523784866183171000000,
	(ar_coef)0.07952517327298877500000,
	(ar_coef)0.08173426919691025000000,
	(ar_coef)0.08173426919691025000000,
	(ar_coef)0.07952517327298877500000,
	(ar_coef)0.07523784866183171000000,
	(ar_coef)0.06912341735634533500000,
	(ar_coef)0.06153307130028077400000,
	(ar_coef)0.05289027931318632900000,
	(ar_coef)0.04365805551438314800000,
	(ar_coef)0.03430437135864797300000,
	(ar_coef)0.02526887281709741900000,
	(ar_coef)0.01693381101612880000000,
	(ar_coef)0.00960154794730574270000,
	(ar_coef)0.00348023083652750430000,
	(ar_coef)-0.00132166536495993110000,
	(ar_coef)-0.00479214237535823420000,
	(ar_coef)-0.00700487680048583750000,
	(ar_coef)-0.00810252096270313790000,
	(ar_coef)-0.00827624463227718750000,
	(ar_coef)-0.00774382230853019790000,
	(ar_coef)-0.00672847637631980230000,
	(ar_coef)-0.00544036464115712180000,
	(ar_coef)-0.00406210895352451890000,
	(ar_coef)-0.00273917269803713920000,
	(ar_coef)-0.00157528660675491580000,
	(ar_coef)-0.00063256715714818368000,
	(ar_coef)0.00006447169166240871600,
	(ar_coef)0.00052210591669823178000,
	(ar_coef)0.00076900780323513871000,
	(ar_coef)0.00084745923215067841000,
	(ar_coef)0.00080519473724929207000,
	(ar_coef)0.00068863819850353001000,
	(ar_coef)0.00053799024989214996000,
	(ar_coef)0.00038431473142678716000,
	(ar_coef)0.00024850129218651672000,
	(ar_coef)0.00014177620987582376000,
	(ar_coef)0.00011884022274186547000,
	0
};

// The resampler configuration struct
const AudResampleConfig DownsampleX10CoeffsCfg = {
  (sizeof(firCoeffs)/sizeof(ar_coef)-1),
  10,
  -1,
  firCoeffs
};
