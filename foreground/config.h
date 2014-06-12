// Config header.
// Paris Kaimakis 29 Jan 2013


#ifndef _CONFIG_H_
#define _CONFIG_H_


#define NUM_CONFIG_VARS		7

#define E_INT				0
#define E_DOUBLE			1



struct s_settings
{
	int U, V;
	double minVariance, minClusterSize;
	double thSigma, DeltaStDev, sigmaMin;
	double propTopLeftTemplateSky, propTopLeftTemplateSea;

};
typedef struct s_settings t_settings;






#endif //_CONFIG_H_







