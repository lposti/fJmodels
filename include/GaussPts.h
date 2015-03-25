/*
 * GaussPts.h
 *
 *  Created on: Feb 25, 2015
 *      Author: morpheus
 */

#ifndef INCLUDE_GAUSSPTS_H_
#define INCLUDE_GAUSSPTS_H_


#define QUAD22

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


/*
 *  This is rectangle method
 */
#ifdef QUAD22

#	define QUADORD 22
	static double xpt11[QUADORD]={-1.        , -0.50118723, -0.25118864, -0.12589254, -0.06309573,
		       -0.03162278, -0.01584893, -0.00794328, -0.00398107, -0.00199526,
		       -0.001     ,  0.001     ,  0.00199526,  0.00398107,  0.00794328,
		        0.01584893,  0.03162278,  0.06309573,  0.12589254,  0.25118864,
		        0.50118723,  1. },
				  xw11[QUADORD]  ={0.49881277,  0.37440568,  0.18764735,  0.09404645,  0.04713488,
					        0.0236234 ,  0.01183975,  0.00593393,  0.00297401,  0.00149054,
					        0.00149763,  0.00149763,  0.00149054,  0.00297401,  0.00593393,
					        0.01183975,  0.0236234 ,  0.04713488,  0.09404645,  0.18764735,
					        0.37440568,  0.49881277 };

	// tanh-sinh
	static double xpt01[QUADORD] ={ 1.12614037692e-05 ,
			8.09053908976e-05 ,
			0.000425204351858 ,
			0.00172603419131 ,
			0.0056596260047 ,
			0.015541463911 ,
			0.0367576002085 ,
			0.0764489202523 ,
			0.141870454058 ,
			0.237263531493 ,
			0.360233468529 ,
			0.5 ,
			0.639766531471 ,
			0.762736468507 ,
			0.858129545942 ,
			0.923551079748 ,
			0.963242399792 ,
			0.984458536089 ,
			0.994340373995 ,
			0.998273965809 ,
			0.999574795648 ,
			0.999919094609},
	              xw01[QUADORD]  ={2.42000467048e-05 ,
	            		  0.000146089994545 ,
	            		  0.000647128630152 ,
	            		  0.00222235986647 ,
	            		  0.00618889854184 ,
	            		  0.0144761668119 ,
	            		  0.0291728616842 ,
	            		  0.0514732311361 ,
	            		  0.0801433368009 ,
	            		  0.110279557896 ,
	            		  0.133823374104 ,
	            		  0.142799666072 ,
	            		  0.133823374104 ,
	            		  0.110279557896 ,
	            		  0.0801433368009 ,
	            		  0.0514732311361 ,
	            		  0.0291728616842 ,
	            		  0.0144761668119 ,
	            		  0.00618889854184 ,
	            		  0.00222235986647 ,
	            		  0.000647128630152 ,
	            		  0.000146089994545};

#endif


/*
 *  This is rectangle method
 */
#ifdef QUAD26

#	define QUADORD 26
	static double xpt11[QUADORD]={-1.00000000e+00,  -3.98107171e-01,  -1.58489319e-01,
	        -6.30957344e-02,  -2.51188643e-02,  -1.00000000e-02,
	        -3.98107171e-03,  -1.58489319e-03,  -6.30957344e-04,
	        -2.51188643e-04,  -1.00000000e-04,   1.00000000e-04,
	         2.51188643e-04,   6.30957344e-04,   1.58489319e-03,
	         3.98107171e-03,   1.00000000e-02,   2.51188643e-02,
	         6.30957344e-02,   1.58489319e-01,   3.98107171e-01,
	         1.00000000e+00},
				  xw11[QUADORD]  ={6.01892829e-01,   4.20755340e-01,   1.67505718e-01,
					         6.66852275e-02,   2.65478672e-02,   1.05688963e-02,
					         4.20755340e-03,   1.67505718e-03,   6.66852275e-04,
					         2.65478672e-04,   1.75594322e-04,   1.75594322e-04,
					         2.65478672e-04,   6.66852275e-04,   1.67505718e-03,
					         4.20755340e-03,   1.05688963e-02,   2.65478672e-02,
					         6.66852275e-02,   1.67505718e-01,   4.20755340e-01,
					         6.01892829e-01};

	static double xpt01[QUADORD] ={1.00000000e-04,   1.44543977e-04,   2.08929613e-04,
	         3.01995172e-04,   4.36515832e-04,   6.30957344e-04,
	         9.12010839e-04,   1.31825674e-03,   1.90546072e-03,
	         2.75422870e-03,   3.98107171e-03,   5.75439937e-03,
	         8.31763771e-03,   1.20226443e-02,   1.73780083e-02,
	         2.51188643e-02,   3.63078055e-02,   5.24807460e-02,
	         7.58577575e-02,   1.09647820e-01,   1.58489319e-01,
	         2.29086765e-01,   3.31131121e-01,   4.78630092e-01,
	         6.91830971e-01,   1.00000000e+00},
	              xw01[QUADORD]  ={4.45439771e-05,   5.44648065e-05,   7.87255975e-05,
	            	         1.13793110e-04,   1.64481086e-04,   2.37747504e-04,
	            	         3.43649697e-04,   4.96724939e-04,   7.17985982e-04,
	            	         1.03780549e-03,   1.50008534e-03,   2.16828300e-03,
	            	         3.13412249e-03,   4.53018529e-03,   6.54810998e-03,
	            	         9.46489859e-03,   1.36809409e-02,   1.97749760e-02,
	            	         2.85835368e-02,   4.13157809e-02,   5.97194728e-02,
	            	         8.63209011e-02,   1.24771664e-01,   1.80349925e-01,
	            	         2.60684954e-01,   3.08169029e-01};

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
