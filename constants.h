#include "axi.h"                       // axisymmetric geometry
#include "navier-stokes/centered.h"    // solve NS equations
#define FILTERED                       // Smear density and viscosity jumps
#include "two-phase.h"
#include "tension.h"                   // include surface tension between phases
#include "tag.h"                       // helps track droplet properties
#include "curvature.h"
#include "reduced.h"


#define DIM_NONDIM_EXP			'd' // d: dimension; n: nondimension; e: experimentalization

#if DIM_NONDIM_EXP == 'd' || DIM_NONDIM_EXP == 'D'

/*The transient force profile of low-speed droplet impact: measurements and model doi:10.1017/jfm.2019.141 
TABLE 1. Material properties of the liquids employed  Page 305 
Impacting drop of diameter D = 2.050e-03
Liquid                 Density ρ  (kg m−3)               Viscosity µ              Surface tension σ                           
Glycerin                       1258                 1490.0    [1.49 Pa-s]             [0.063 N/m]  
Glycerin (22.0 °C)             1256.6               1136.5    [1.1365 Pa-s]           [0.0632 N/m] 
1000 cst Silicone oil          970                  970.0     [0.97 Pa-s]             [0.0217 N/m] 
Glycerin (26.0 °C)             1254.8               872.3     [0.8723 Pa-s]           [0.06301 N/m]
95 % Glycerin                  1248                 523       [0.523 Pa-s]            [0.063 N/m]  
500 cst Silicone oil           970                  465.0     [0.465 Pa-s]            [0.01989 N/m]   
90 % Glycerin                  1235                 219       [0.219 Pa-s]            [0.064 N/m]
100 cst Silicone oil           960                  96.0      [0.096 Pa-s]            [0.0211 N/m]
80 % Glycerin                  1208                 60        [0.06 Pa-s]             [0.065 N/m]  
70 % Glycerin                  1181.3               22.5      [0.0225 Pa-s]           [0.0664 N/m] 
60 % Glycerin                  1153.8               10.8      [0.0108 Pa-s]           [0.0676 N/m]  
20 cst Silicone oil            950.0                20        [0.02 Pa-s ]            [0.0189 N/m] 
50 % Glycerin                  1126                  6        [0.006 Pa-s]            [0.068 N/m]   
Water                          998                  1.0       [0.001 Pa-s]            [0.073 N/m]
Ethanol  

Ji, B., Yang, Z. & Feng, J. Compound jetting from bubble bursting at an air-oil-water interface. Nat Commun 12, 6305 (2021). 
https://doi.org/10.1038/s41467-021-26382-w

liquids                     ρ (kg/m3)                  μ (mPa .s)                γwa (mN/m)          γoa (mN/m)           γow (mN/m) 
  
DI water                    998                          0.89                     71.6±1.0           N/A                 N/A

5cSt silicone oil          913                          4.6                       N/A                18.7±0.3           38.1±0.4

20cSt silicone oil         950                          19                        N/A                 19.4±0.7          40.9±0.5

     
100cSt silicone oil        960                          96                        N/A                 20.1±0.2          43.7±0.4

            
500cSt silicone oil        970                          465                       N/A                 19.9±0.1          38.7±0.2

1000cSt silicone oil       970                          970                        N/A                21.7±0.2          40.8±1.0 
*/




#define VELOCITY			0.10                    // Velocity of 500 cst Silicone oil    m/s  Si unit  for         
#define DROP_DIAMETER		2.050e-03               // Diameter of 500 cst Silicone oil   drop  meter Si unit 
#define RHO_L				970.0                  // Density of 500cst Silicone oil   25 degree kg/m^3  Si unit 
#define RHO_G				1.21                    // Density of air at 25 degree degree kg/m^3  Si unit 
#define MU_L				0.465                  // Dynamisc Viscosity of  500 cst Silicone oil   at 25 degree Pa s in Si unit 
#define MU_G				1.81e-5                 // Dynamic Viscosity of air at 25 degree
#define SIGMA				0.01989                 // Surface tension of 500 cst Silicone oil  drop  at 25 degree   N/m  Si unit 
#define GRAVITY				9.81  
//
#define RHO_GL				0.0
#define MU_GL				0.0
#define REYNOLDS			0.0
#define WEBER				0.0
#define FROUDE              0.0

#elif DIM_NONDIM_EXP == 'n' || DIM_NONDIM_EXP == 'N'

#define WEBER				100.0
#define REYNOLDS			100.0
#define FROUDE              70.0
#define RHO_GL				(0.0012) // air-water at 25C: 0.001187503
#define MU_GL				(0.0210) // air-water at 25C: 0.020898876
//
#define VELOCITY			0.0
#define DROP_DIAMETER		0.0
#define RHO_L				0.0
#define MU_L				0.0
#define SIGMA				0.0
#define RHO_G				0.0
#define MU_G				0.0
#define GRAVITY				0.0

#elif DIM_NONDIM_EXP == 'e' || DIM_NONDIM_EXP == 'E'

#define WEBER				300.0
#define REYNOLDS			1000.0
#define FROUDE              70.0
#define DROP_DIAMETER		2.0e-3
#define SIGMA				17.6e-3
#define RHO_L				816.0
#define RHO_G				1.2041
#define MU_G				1.94e-5
//
#define VELOCITY			0.0
#define RHO_GL				0.0
#define MU_GL				0.0
#define MU_L				0.0
#define GRAVITY				0.0

#endif

#define INITAL_GRID_LEVEL		9
#define MAX_GRID_LEVEL			11
#define DOMAIN_WIDTH			4.00
#define POOL_DEPTH				0.00
#define INITIAL_DISTANCE		0.04
#define BUBBLE_DIAMETER		    0.00      
#define DBDELTA       		    0.00    
#define REFINE_GAP				0.02
#define MAX_TIME				1.00
#define SAVE_FILE_EVERY			0.01 

#define REFINE_VAR				{f, u.x, u.y} 
#define REFINE_VAR_TEXT			"f, u.x, u.y" 
#define REFINE_VALUE_0			-6
#define REFINE_VALUE_1			-3
#define REFINE_VALUE_2			-3

#define REMOVE_DROP_YESNO		'n'
#define REMOVE_DROP_SIZE		4.0 // equivalent diameter base on the maximum refinement
#define REMOVE_DROP_PERIOD		4
#define REMOVE_BUBBLE_YESNO		'n'
#define REMOVE_BUBBLE_SIZE		4.0 // equivalent diameter base on the maximum refinement
#define REMOVE_BUBBLE_PERIOD	4


#define FILENAME_DATA			"data"
#define FILENAME_DURATION		"duration"
#define FILENAME_PARAMETERS		"parameters.txt"
#define FILENAME_ENDOFRUN		"endofrun"
#define FILENAME_LASTFILE		"lastfile"

#define R_VOFLIMIT				1.0e-9
#define R_PI					3.1415926535897932384626433832795

int LEVELmin = INITAL_GRID_LEVEL, LEVELmax = MAX_GRID_LEVEL ;
double maxruntime = HUGE;
scalar fdrop[], pressure[];
scalar fb[]; 

struct CFDValues {
	double rhoL, rhoG, muL, muG, Sigma;
	double vel, Reynolds, Weber, Froude,Bond, Oh,GXnormlised;
	double diameter, domainsize, refinegap, pooldepth, initialdis;
	double timecontact, timeend, timestep;
	double bubblediameter;    
	double dbdelta;        
};

void readfromarg(char **argv, int argc, struct CFDValues *bvalues);

int numericalmainvalues(char **argv, int argc, struct CFDValues *bvalues)
{
	double velocity = VELOCITY, mu_l = MU_L;
	bvalues->rhoL = 1.0;
	bvalues->vel = 1.0;
	bvalues->diameter = 1.0;
	;
	bvalues->Reynolds = -1.0;
	bvalues->Weber = -1.0;
	bvalues->pooldepth = -1.0;
	bvalues->timeend = -1.0;
	bvalues->timestep = -1.0;
	;
	readfromarg(argv, argc, bvalues);
	switch(DIM_NONDIM_EXP)
	{
	case 'd':
	case 'D':
	{
		bvalues->Reynolds = (RHO_L * VELOCITY * DROP_DIAMETER / mu_l);
		bvalues->Weber = (RHO_L * VELOCITY * VELOCITY * DROP_DIAMETER / SIGMA);
		bvalues->Froude = (VELOCITY / sqrt (GRAVITY * DROP_DIAMETER));
		;
		bvalues->Bond = (RHO_L * GRAVITY * DROP_DIAMETER  * DROP_DIAMETER / SIGMA);
		bvalues->Oh = (mu_l/ sqrt(RHO_L * SIGMA * DROP_DIAMETER));
		;
		bvalues->muL = (bvalues->rhoL * bvalues->vel * bvalues->diameter / bvalues->Reynolds);
		bvalues->Sigma = (bvalues->rhoL * bvalues->vel * bvalues->vel * bvalues->diameter / bvalues->Weber);
		bvalues->rhoG = (RHO_G / RHO_L) * bvalues->rhoL; 
		bvalues->muG = (MU_G / mu_l)  * bvalues->muL; 
		bvalues->GXnormlised = (GRAVITY * DROP_DIAMETER / (VELOCITY*VELOCITY));
		;
	    bvalues->bubblediameter = BUBBLE_DIAMETER * bvalues->diameter;    
	    bvalues->dbdelta = DBDELTA * bvalues->diameter;                 
		break;
	}
	case 'n':
	case 'N':
	{
		if (bvalues->Reynolds < 0.0)
			bvalues->Reynolds = REYNOLDS;
		if (bvalues->Weber < 0.0)
			bvalues->Weber = WEBER;
		if (bvalues->Froude < 0.0)
			bvalues->Froude = FROUDE;
		bvalues->muL = (bvalues->rhoL * bvalues->vel * bvalues->diameter / bvalues->Reynolds);
		bvalues->Sigma = (bvalues->rhoL * bvalues->vel * bvalues->vel * bvalues->diameter / bvalues->Weber);
		bvalues->rhoG = RHO_GL * bvalues->rhoL;
		bvalues->muG = MU_GL * bvalues->muL;
		break;
	}
	case 'e':
	case 'E':
	{
		if (bvalues->Reynolds < 0.0)
			bvalues->Reynolds = REYNOLDS;
		if (bvalues->Weber < 0.0)
			bvalues->Weber = WEBER;
		bvalues->muL = (bvalues->rhoL * bvalues->vel * bvalues->diameter / bvalues->Reynolds);
		bvalues->Sigma = (bvalues->rhoL * bvalues->vel * bvalues->vel * bvalues->diameter / bvalues->Weber);
		velocity = sqrt (bvalues->Weber * SIGMA / (DROP_DIAMETER * RHO_L));
		mu_l = (RHO_L * velocity * DROP_DIAMETER / bvalues->Reynolds);
		bvalues->muG = (MU_G / mu_l) * bvalues->muL;
		bvalues->rhoG = (RHO_G / RHO_L) * bvalues->rhoL;
		break;
	}
	}
	bvalues->domainsize = DOMAIN_WIDTH * bvalues->diameter;
	if (bvalues->pooldepth < 0.0)
		bvalues->pooldepth = POOL_DEPTH * bvalues->diameter;
	bvalues->initialdis = INITIAL_DISTANCE * bvalues->diameter;
	bvalues->refinegap = REFINE_GAP * bvalues->diameter;
	;
	bvalues->timecontact = bvalues->initialdis / bvalues->vel;
	if (bvalues->timeend < 0.0)
		bvalues->timeend = MAX_TIME;
	if (bvalues->timestep < 0.0)
		bvalues->timestep = SAVE_FILE_EVERY;
	;
	switch (pid())
	{
	case 0:
	{
		printf("R: %f --- W: %f --- H: %f\r\n", bvalues->Reynolds, bvalues->Weber, bvalues->pooldepth);
		FILE *fp;
		fp = fopen (FILENAME_PARAMETERS, "w");
	    fprintf (fp, "Name of Liquid : 500 cst Silicone oil   drop  \r\n");
	    fprintf (fp, "Experimental parameters  / Numerical simulation parameters in Basilisk unit\r\n");
		fprintf (fp, "Diameter_Experimental: %.3e / Normalized diameter  of drop (D): %.3e\r\n", DROP_DIAMETER, bvalues->diameter);
		fprintf (fp, "Velocity_Experimental: %.3e / Normalized Velocity  of liquid (V) : %.3e\r\n", velocity, bvalues->vel);
		fprintf (fp, "Rho(L)_Experimental: %.3e /  Normalized density of liquid (rho1): %.3e\r\n", RHO_L, bvalues->rhoL); 
		fprintf (fp, "Rho(G)_Experimental: %.3e / Normalized density  of Air (rho2) : %.3e\r\n", RHO_G, bvalues->rhoG);
		fprintf (fp, "Mu(L)_Experimental: %.3e / Normalized viscosity of liquid (mu1): %.3e\r\n", mu_l, bvalues->muL);
		fprintf (fp, "Mu(G)_Experimental: %.3e / Normalized viscosity of Air (mu2): %.3e\r\n", MU_G, bvalues->muG);
		fprintf (fp, "Sigma(L-G)_Experimental: %.3e / Normalized surface tension (f.sigma): %.3e\r\n", SIGMA, bvalues->Sigma);
		fprintf (fp, "Acceleration due to gravity: %.3e / Normalized gravity (G.x): %.3e\r\n", GRAVITY, bvalues->GXnormlised);
		fprintf (fp, "\r\n");
		fprintf (fp, "Reynolds: %.10f\r\n", bvalues->Reynolds);
		fprintf (fp, "Weber: %.10f\r\n", bvalues->Weber);
		fprintf (fp, "Froude: %.10f\r\n", bvalues->Froude);
		fprintf (fp, "Bond: %.10f\r\n", bvalues->Bond);
		fprintf (fp, "Ohsorge: %.10f\r\n", bvalues->Oh);
		fprintf (fp, "\r\n");
		fprintf (fp, "Level Max: %d\r\n", LEVELmax);
		fprintf (fp, "Level Min: %d\r\n", LEVELmin);
		fprintf (fp, "Domain Size: %.2f\r\n", bvalues->domainsize);
		fprintf (fp, "Pool Depth: %.2f\r\n", bvalues->pooldepth);
		fprintf (fp, "Initial Distance: %.2f\r\n", bvalues->initialdis);
		fprintf (fp, "Refine Gap: %.2f\r\n", bvalues->refinegap);
		fprintf (fp, "Contact Time: %.2f\r\n", bvalues->timecontact);
		fprintf (fp, "Domain Size: %.2f\r\n", bvalues->domainsize);
		fprintf (fp, "\r\n");
		fprintf (fp, "Bubble Diameter: %.2f\r\n", bvalues->bubblediameter);     
		fprintf (fp, "Dbdelta (disatance btw drop and bubble) : %.6f\r\n", bvalues->dbdelta); 
		fprintf (fp, "\r\n");  
		fprintf (fp, "Refine Variables: %s\r\n", REFINE_VAR_TEXT);
		fprintf (fp, "Refine Variables powers: %d, %d, %d\r\n", REFINE_VALUE_0, REFINE_VALUE_1, REFINE_VALUE_2);
		fprintf (fp, "\r\n");
		fprintf (fp, "Remove Drop YesNo: %c\r\n", REMOVE_DROP_YESNO);
		fprintf (fp, "Remove Drop Size: %f\r\n", REMOVE_DROP_SIZE);
		fprintf (fp, "Remove Drop Period: %d\r\n", REMOVE_DROP_PERIOD);
		fprintf (fp, "\r\n");
		fprintf (fp, "Remove Bubble YesNo: %c\r\n", REMOVE_BUBBLE_YESNO);
		fprintf (fp, "Remove Bubble Size: %f\r\n", REMOVE_BUBBLE_SIZE);
		fprintf (fp, "Remove Bubble Period: %d\r\n", REMOVE_BUBBLE_PERIOD);
		fclose (fp);
		break;
	}
	}
	return 1;
}

void readfromarg(char **argv, int argc, struct CFDValues *bvalues)
{
	int i, j;
	char tmp[100];
	if (argc < 2)
		return;
	for (i = 1; i < argc; i++)
	{
		switch(argv[i][0])
		{
		case 'r':
		case 'R':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			bvalues->Reynolds = atof(tmp);
			break;
		}
		case 'w':
		case 'W':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			bvalues->Weber = atof(tmp);
			break;
		}
		case 'F':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			bvalues->Froude = atof(tmp);
			break;
		}
		case 'h':
		case 'H':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			bvalues->pooldepth = atof(tmp);
			break;
		}
		case 'x':
		case 'X':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			LEVELmax = atoi(tmp);
			break;
		}
		case 'n':
		case 'N':
		{
			for(j = 1; j < (int)strlen(argv[i]); j++)
				tmp[j - 1] = argv[i][j];
			tmp[j - 1] = '\0';
			LEVELmin = atoi(tmp);
			break;
		}
		case 't':
		case 'T':
		{
			switch(argv[i][1])
			{
			case 'e':
			case 'E':
			{
				for(j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				bvalues->timeend = atof(tmp);
				break;
			}
			case 's':
			case 'S':
			{
				for(j = 2; j < (int)strlen(argv[i]); j++)
					tmp[j - 2] = argv[i][j];
				tmp[j - 2] = '\0';
				bvalues->timestep = atof(tmp);
				break;
			}
			}
			break;
		}
		}
	}
}

int timecalculation(double t, char *chartime)
{
	int d, h, m, s;
	if(t < 60.0)
	{
		d = 0;
		h = 0;
		m = 0;
		s = (int) t;
	}
	else if(t < 3600.0)
	{
		d = 0;
		h = 0;
		m = (int) (t / 60.0);
		s = (int) (t - m*60.0);
	}
	else if(t < 3600.0*24.0)
	{
		d = 0;
		h = (int) (t / 3600.0);
		m = (int) ((t - h*3600.0) / 60.0);
		s = (int) (t - h*3600.0 - m*60.0);
	}
	else
	{
		d = (int) (t / 3600.0 / 24.0);
		h = (int) ((t - d*3600.0*24.0) / 3600.0);
		m = (int) ((t - d*3600.0*24.0 - h*3600.0) / 60.0);
		s = (int) (t - d*3600.0*24.0 - h*3600.0 - m*60.0);
	}
	sprintf(chartime, "%d:%02d:%02d:%02d", d, h, m, s);
	return 1;
}

