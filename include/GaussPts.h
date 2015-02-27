/*
 * GaussPts.h
 *
 *  Created on: Feb 25, 2015
 *      Author: morpheus
 */

#ifndef INCLUDE_GAUSSPTS_H_
#define INCLUDE_GAUSSPTS_H_


#define QUAD28

#ifdef QUAD12

#	define QUADORD 12
	static double xpt11[QUADORD]={-0.981560634,-0.904117256,-0.769902674,
								  -0.587317954,-0.367831499,-0.125233409,
								   0.125233409,0.367831499,0.587317954,
								   0.769902674,0.904117256,0.981560634},
				  xw11[QUADORD]  ={0.047175336387,0.10693932600,0.16007832854,
								   0.20316742672,0.23349253654,0.24914704581,
								   0.24914704581,0.23349253654,0.20316742672,
								   0.16007832854,0.10693932600,0.047175336387};

	static double xpt01[QUADORD] ={0.00921968288,0.0479413718,0.1150486629,
								   0.2063410229,0.316084251,0.437383296,
								   0.562616704,0.683915749,0.793658977,
								   0.884951337,0.952058628,0.990780317},
				  xw01[QUADORD]  ={0.023587668193,0.053469662998,0.08003916427,
								   0.10158371336,0.11674626827,0.12457352291,
								   0.12457352291,0.11674626827,0.10158371336,
								   0.08003916427,0.053469662998,0.023587668193};
#endif

#ifdef QUAD20

#	define QUADORD 20
	static double xpt11[QUADORD]={-0.9931286,-0.9639719,-0.9122344,
								  -0.9122344,-0.8391170,-0.7463319,
								  -0.6360537,-0.5108670,-0.3737061,
								  -0.2277859,-0.076527,0.076527,
								   0.2277859,0.3737061,0.5108670,
								   0.6360537,0.7463319,0.8391170,
								   0.9639719,0.9931286},
				  xw11[QUADORD]  ={0.0176140071,0.0406014298,0.0626720483,
								   0.083276742,0.101930120,0.118194532,
								   0.131688638,0.142096109,0.149172986,
								   0.152753387,0.152753387,0.149172986,
								   0.142096109,0.131688638,0.118194532,
								   0.101930120,0.083276742,0.0626720483,
								   0.0406014298,0.0176140071};

	static double xpt01[QUADORD] ={0.003435700,0.01801404,0.04388279,
								   0.08044151,0.1268340,0.1819732,
								   0.2445665,0.3131470,0.3861071,
								   0.4617367,0.5382633,0.6138929,
								   0.6868530,0.7554335,0.8180268,
								   0.8731660,0.9195585,0.9561172,
								   0.9819860,0.9965643},
				  xw01[QUADORD]  ={0.0088070036,0.0203007149,0.0313360242,
								   0.0416383708,0.0509650599,0.0590972660,
								   0.065844319,0.071048055,0.074586493,
								   0.076376694,0.076376694,0.074586493,
								   0.071048055,0.065844319,0.0590972660,
								   0.0509650599,0.0416383708,0.0313360242,
								   0.0203007149,0.0088070036};

#endif


#ifdef QUAD22

#	define QUADORD 22
	static double xpt11[QUADORD]={-0.99429459, -0.9700605 , -0.92695677, -0.86581258, -0.78781681,
		       -0.69448726, -0.5876404 , -0.46935584, -0.34193582, -0.20786043,
		       -0.06973927,  0.06973927,  0.20786043,  0.34193582,  0.46935584,
		        0.5876404 ,  0.69448726,  0.78781681,  0.86581258,  0.92695677,
		        0.9700605 ,  0.99429459},
				  xw11[QUADORD]  ={0.014628  ,  0.0337749 ,  0.05229334,  0.06979647,  0.08594161,
					        0.10041414,  0.1129323 ,  0.12325238,  0.1311735 ,  0.1365415 ,
					        0.13925187,  0.13925187,  0.1365415 ,  0.1311735 ,  0.12325238,
					        0.1129323 ,  0.10041414,  0.08594161,  0.06979647,  0.05229334,
					        0.0337749 ,  0.014628 };

	static double xpt01[QUADORD] ={ 0.00285271,  0.01496975,  0.03652161,  0.06709371,  0.1060916 ,
	        0.15275637,  0.2061798 ,  0.26532208,  0.32903209,  0.39606979,
	        0.46513036,  0.53486964,  0.60393021,  0.67096791,  0.73467792,
	        0.7938202 ,  0.84724363,  0.8939084 ,  0.93290629,  0.96347839,
	        0.98503025,  0.99714729},
	              xw01[QUADORD]  ={0.007314  ,  0.01688745,  0.02614667,  0.03489823,  0.0429708 ,
	            	        0.05020707,  0.05646615,  0.06162619,  0.06558675,  0.06827075,
	            	        0.06962594,  0.06962594,  0.06827075,  0.06558675,  0.06162619,
	            	        0.05646615,  0.05020707,  0.0429708 ,  0.03489823,  0.02614667,
	            	        0.01688745,  0.007314 };

#endif


#ifdef QUAD28

#	define QUADORD 28
	static double xpt11[QUADORD]={-0.9964425, -0.98130317,  -0.95425928,
			                      -0.91563303, -0.86589252, -0.80564137,
								  -0.73561088, -0.65665109, -0.56972047,
								  -0.47587422, -0.37625152, -0.27206163,
								  -0.16456928, -0.05507929,  0.05507929,
								  0.16456928,  0.27206163,  0.37625152,
								  0.47587422,  0.56972047,  0.65665109,
								  0.73561088,  0.80564137,  0.86589252,
								  0.91563303,  0.95425928,  0.98130317,
								  0.9964425 },
				  xw11[QUADORD]  ={0.00912428,  0.02113211,  0.03290143,
						           0.04427293,  0.05510735,  0.06527292,
								   0.07464621,  0.08311342,  0.09057174,
								   0.09693066,  0.10211297,  0.10605577,
								   0.10871119,  0.11004701,  0.11004701,
								   0.10871119,  0.10605577,  0.10211297,
								   0.09693066,  0.09057174,  0.08311342,
								   0.07464621,  0.06527292,  0.05510735,
								   0.04427293,  0.03290143,  0.02113211,
								   0.00912428};

	static double xpt01[QUADORD] ={0.00177875,  0.00934842,  0.02287036,
			                       0.04218349,  0.06705374,  0.09717931,
								   0.13219456,  0.17167445,  0.21513976,
								   0.26206289,  0.31187424,  0.36396919,
								   0.41771536,  0.47246036,  0.52753964,
								   0.58228464,  0.63603081,  0.68812576,
								   0.73793711,  0.78486024,  0.82832555,
								   0.86780544,  0.90282069,  0.93294626,
								   0.95781651,  0.97712964,  0.99065158,
								   0.99822125},
	              xw01[QUADORD]  ={0.00456214,  0.01056606,  0.01645071,
	            		           0.02213647,  0.02755367,  0.03263646,
								   0.03732311,  0.04155671,  0.04528587,
								   0.04846533,  0.05105648,  0.05302788,
								   0.0543556 ,  0.05502351,  0.05502351,
								   0.0543556 ,  0.05302788,  0.05105648,
								   0.04846533,  0.04528587,  0.04155671,
								   0.03732311,  0.03263646,  0.02755367,
								   0.02213647,  0.01645071,  0.01056606,
								   0.00456214};

#endif





#endif /* INCLUDE_GAUSSPTS_H_ */