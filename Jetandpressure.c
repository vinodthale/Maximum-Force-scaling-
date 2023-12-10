#include "axi.h"
#include "navier-stokes/centered.h"
#define FILTERED            // Smear density and viscosity jumps
#include "two-phase.h"
#include "tension.h"
#include "tag.h"
#include "curvature.h"


#define VOFFOLDER "FVFVDb0.00delta0.00V0.10"
scalar pressure[];
//1] qcc -O2 -Wall -D_FORTIFY_SOURCE=0 Jetandpressure.c -o Jetandpressure -lm
//2] ./Jetandpressure
double cfdbvbubblediameter = 0.00;
double PreFactor = 2*pi; // 2*pi; for theta integration okay
int main (int argc, char **argv)
{
  run();
}

event init (t = 0)
{
    double PressureDropMaxima;
    double timebgn = 0.00;
    double timestp = 0.01;
    double timeend = 1.04;
    ;
    char namefile[500];
    char name[100];
    double timeload;
    FILE *ForceonLfet;
    FILE *SeVAt;
    FILE *StresPre;
    char folder[500];
    strcpy(folder, "mkdir ");
    strcat(folder, VOFFOLDER);
    system(folder);
    for (timeload = timebgn; timeload <= timeend; timeload += timestp)
   {
    sprintf (namefile, "intermediate/snapshot-%5.4f", timeload);
    printf ("load the file %s!\r\n", namefile);
    restore (file = namefile);
    sprintf (name, "%s/VOFFVDb0.00delta0.00V0.10-%.2f.gnu", VOFFOLDER, timeload);
    FILE *ip = fopen (name, "w");
    output_facets (f, ip);        // ################################################ tracer 
    fclose (ip);
   //calculate the force on the substrate
    double pleft = 0.;
    double pForce  = 0.;
    static int nff = 0;
      sprintf (name, "%s/ForceFVDb0.00delta0.00V0.10.txt", VOFFOLDER);
    if (!nff)
      ForceonLfet = fopen(name, "w");
    else
    ForceonLfet = fopen(name, "a");
    double pdatum = 0, wt = 0;
    foreach_boundary(top){
    pdatum += 2*pi*y*pressure[]*(Delta);
    wt += 2*pi*y*(Delta);
    } 
    if (wt >0){
    pdatum /= wt;
    }
    foreach_boundary(left)
    {
    pForce += 2*pi*y*(Delta)*(pressure[]-pdatum);
    pleft += pressure[];
    }
    boundary((scalar *){f, u.x, u.y, pressure});
    ;
    fprintf (ForceonLfet, "%f  %.10f %.10f\r\n", timeload, pForce, pleft);
    fclose (ForceonLfet);
    nff++;
 }
}

event end(t = 0.0)
{
    printf("\r\n-------\r\nEND!\r\n");
}



