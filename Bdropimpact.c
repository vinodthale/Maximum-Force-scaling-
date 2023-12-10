#include "constants.h"
double drop_time_1file, bubble_time_1file, writefile_time_1file, simulation_time_1file;
double drop_time_total, bubble_time_total, writefile_time_total, simulation_time_total;
clock_t simulation_str_time, simulation_end_time;

struct CFDValues cfdbv;
// Boundary conditions
 u.t[left] = dirichlet(0);  // No slip at surface
 f[left] = 0.;    // non wetting 
 u.n[right] = neumann(0);   // Free flow condition
 p[right] = dirichlet(0);   // 0 pressure far from surface
 u.n[top] = neumann(0);     // Allows outflow through boundary
 p[top] = dirichlet(0);     // 0 pressure far from surface
// Default for bottom is symmetry
int main(int argc, char **argv)
{
	simulation_str_time = clock();
	simulation_time_1file = 0.0;
	writefile_time_1file = 0.0;
	drop_time_1file = 0.0;
	bubble_time_1file = 0.0;
	numericalmainvalues(argv, argc, &cfdbv);
	;
	size(cfdbv.domainsize);
#if AXI
	;
#else
	origin(0, -cfdbv.domainsize / 2., -cfdbv.domainsize / 2.);
#endif
	int initialgrid = pow(2, LEVELmin);
	init_grid(initialgrid);
	;
	char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);
	rho1 = cfdbv.rhoL;
	rho2 = cfdbv.rhoG;
	mu1 = cfdbv.muL;
	mu2 = cfdbv.muG;
	f.sigma = cfdbv.Sigma;
	G.x -= 1.0/sq(cfdbv.Froude);  
  Z.x = 0.0;
  /* Poisson solver constants */
  //DT = 1.0e-4;          // Minimum timestep
  //NITERMIN = 1;         // Min number of iterations (default 1)
  //NITERMAX = 100;       // Max number of iterations (default 100)
  TOLERANCE = 1e-6;       // Possion solver tolerance (default 1e-3)
  run();
	return 1;               
}

event defaults(i = 0)
{
	interfaces = list_add(NULL, f);
	interfaces = list_add(interfaces,fb); 
}

 event initfraction (t = 0)
{
  double x0 = cfdbv.pooldepth + cfdbv.initialdis + cfdbv.diameter*0.50;// This is center drop centre
  double Bubtx0 = cfdbv.pooldepth + cfdbv.initialdis + cfdbv.diameter - ((cfdbv.dbdelta*cfdbv.diameter)) - 0.50*(cfdbv.bubblediameter*cfdbv.diameter);  
  fraction(f,(min(sq(0.50*cfdbv.diameter)-(sq(x - x0) + sq(y) + sq(z)),-(sq(0.50*(cfdbv.bubblediameter*cfdbv.diameter))-(sq(x - Bubtx0) + sq(y) + sq(z))))));  
}

event init(i = 0)
{
	if (restore(file = FILENAME_LASTFILE))
	{
#if AXI
		boundary((scalar *){fm});
		//boundary({p});
#endif
	}
	else
	{
		double x0 = cfdbv.pooldepth + cfdbv.initialdis + cfdbv.diameter*0.50; // This is drop center 
		//double Bubtx0 = cfdbv.pooldepth + cfdbv.initialdis + cfdbv.diameter - ((cfdbv.dbdelta * cfdbv.diameter)) - 0.50 * (cfdbv.bubblediameter * cfdbv.diameter);  //  This is center of Bubble 
		//double Bubtx0 = cfdbv.pooldepth + cfdbv.initialdis + cfdbv.diameter * 0.50;  //  This is center of Bubble
		;
		refine(sq(x - x0) + sq(y) + sq(z) < sq(0.50*cfdbv.diameter + cfdbv.refinegap) && sq(x - x0) + sq(y) + sq(z) > sq(0.50*cfdbv.diameter - cfdbv.refinegap) && level < LEVELmax); // refinement along Dorp 
		//refine(sq(x - x0) + sq(y) + sq(z) < sq(0.50 * cfdbv.diameter + cfdbv.refinegap) &&  level < LEVELmax); // refinement along Dorp 
		foreach ()
		{
      f[] = 0.0;
		  if(sq(x - x0) + sq(y) + sq(z) < sq(0.50*cfdbv.diameter))  // this is for Drop 
			{
        f[] = 1.0;
        u.x[] = -cfdbv.vel;
        u.y[] = 0.0;
			}
      /*if(sq(x - Bubtx0) + sq(y) + sq(z) < sq(0.50*(cfdbv.bubblediameter * cfdbv.diameter))) // this for the Inside the bubble 
			{
				f[] = 0.0;  
        u.x[] = -cfdbv.vel;
        u.y[] = 0.0;
			}*/
		};
		clock_t timestr, timeend;
		timestr = clock();
		FILE *fp;
		char name[100], tmp[50];
		sprintf(name, FILENAME_DURATION);
		sprintf(tmp, "-CPU%02d.plt", pid());
		strcat(name, tmp);
		fp = fopen(name, "w");
		fprintf(fp, "Variables = Iteration DeltaTime CriticalTime PhysicalTime LastDuration DropDuration BubbleDuration FileDuration CellNumber TotalLastDuration TotalDropDuration TotalBubbleDuration TotalFileDuration\r\nzone\r\n");
		fclose(fp);
		timeend = clock();
		writefile_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC;
	}
}



event adapt(i++)
{
	double refine[3];
	refine[0] = pow(10.0, REFINE_VALUE_0);
	refine[1] = pow(10.0, REFINE_VALUE_1);
	refine[2] = pow(10.0, REFINE_VALUE_2);
	adapt_wavelet(REFINE_VAR, (double[]){refine[0], refine[1], refine[2]}, maxlevel = LEVELmax, minlevel = LEVELmin);
}

event showiteration(i++)
{
	switch (pid())
	{
	case 0:
	{
		char name[500], tmp[100];
		if (t - cfdbv.timecontact < 0.0)
			sprintf(name, "i%05d_dt%.2e_tb%.3f_P%02d", i, dt, t - cfdbv.timecontact, (int)(100.0 * t / MAX_TIME));
		else
		sprintf(name, "i%05d_dt%.2e_ta%.3f_P%02d", i, dt, t - cfdbv.timecontact, (int)(100.0 * t / MAX_TIME));
		sprintf(tmp, "_Re%.5f_We%.5f", (double)cfdbv.Reynolds, (double)cfdbv.Weber);
		//sprintf(tmp, "_Re%d_We%d", (int)cfdbv.Reynolds, (int)cfdbv.Weber);
		strcat(name, tmp);
#if AXI
		sprintf(tmp, "_AXI");
		strcat(name, tmp);
#else
#if dimension == 3
		sprintf(tmp, "_3D");
		strcat(name, tmp);
#else
		sprintf(tmp, "_2D");
		strcat(name, tmp);
#endif
#endif
		sprintf(tmp, "_L%02d%02d", LEVELmin, LEVELmax);
		strcat(name, tmp);
		printf("%s\r\n", name);
	}
	}
}

event end(t = cfdbv.timecontact + cfdbv.timeend)
{
	FILE *fp;
	char name[500], tmp[100];
	sprintf(name, FILENAME_ENDOFRUN);
	sprintf(tmp, "-CPU%02d.txt", pid());
	strcat(name, tmp);
	fp = fopen(name, "w");
	fprintf(fp, "SimulationTime %e\r\nRemoveDropTime %e\r\nRemoveBubbleTime %e\r\nWriteFileTime %e\r\n", simulation_time_total, drop_time_total, bubble_time_total, writefile_time_total);
	fclose(fp);
}

event outputfiles (t += SAVE_FILE_EVERY)//remaining the beginning time
{
	clock_t timestr, timeend;
	timestr = clock();
	static FILE *fp;
	char name[500], tmp[100];
	;
	foreach ()
	{
		if (fb[] < 0.0)
			fb[] = 0.0;
		else if (fb[] > 1.0)
			fb[] = 1.0;
	};
	foreach(){
		pressure[] = p[]; 
	}
	p.nodump = false;
	sprintf (name, "intermediate/snapshot-%5.4f", t);
	dump(file = name);
	;
	dump(file = FILENAME_LASTFILE);
	timeend = clock();
	writefile_time_1file += (double)(timeend - timestr) / CLOCKS_PER_SEC;
	int cellnumber = 0;
	foreach ()
		cellnumber++;
	simulation_end_time = clock();
	double estimatetimeleft;
	char LDc[100], TDc[100], ETLc[100];
	simulation_time_1file = (double)(simulation_end_time - simulation_str_time) / CLOCKS_PER_SEC;
	;
	simulation_time_total += simulation_time_1file;
	bubble_time_total += bubble_time_1file;
	drop_time_total += drop_time_1file;
	writefile_time_total += writefile_time_1file;
	if (t == 0.0)
		estimatetimeleft = 0.0;
	else
		estimatetimeleft = simulation_time_total * (cfdbv.timecontact + cfdbv.timeend) / t - simulation_time_total;
	timecalculation(simulation_time_1file, LDc);
	timecalculation(simulation_time_total, TDc);
	timecalculation(estimatetimeleft, ETLc);
	sprintf(name, FILENAME_DURATION);
	sprintf(tmp, "-CPU%02d.plt", pid());
	strcat(name, tmp);
	fp = fopen(name, "a");
	//fprintf(fp, "Variables = Iteration DeltaTime CriticalTime PhysicalTime LastDuration DropDuration BubbleDuration FileDuration CellNumber TotalLastDuration TotalDropDuration TotalBubbleDuration TotalFileDuration\r\nzone\r\n");
	fprintf(fp, "%d %e %e %e %e %e %e %e %d %e %e %e %e\r\n", i, dt, t - cfdbv.timecontact, t, simulation_time_1file, drop_time_1file, bubble_time_1file, writefile_time_1file, cellnumber, simulation_time_total, drop_time_total, bubble_time_total, writefile_time_total);
	fclose(fp);
	simulation_str_time = clock();
	simulation_time_1file = 0.0;
	writefile_time_1file = 0.0;
	drop_time_1file = 0.0;
	bubble_time_1file = 0.0;
	;
	switch (pid())
	{
	case 0:
	{
		printf("\r\nData Files are Written!\r\nDuration Last: %s\r\nTotal: %s, Time Left: %s\r\n\r\n", LDc, TDc, ETLc);
		break;
	}
	}
}



