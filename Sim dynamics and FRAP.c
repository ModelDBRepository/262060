
#define ISAAC

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>
#include <malloc.h>
#include <time.h>
#include <math.h>

#ifdef ISAAC

#ifndef STANDARD
#include "standard.h"
#endif
#ifndef RAND
#include "rand.h"
#endif

#endif

#define	MATRIXDIM		50
#define N_SYNAPSES		4000
#define	STEPS			800
#define P_ON			1.25E-6
#define P_OFF			0.5
#define	ALPHA			(0.0007 / 640000)
#define TOTAL_MOLECULES		640000
#define FRAP_CYCLE		600
#define FRAP_SYNAPSES		200

#define OUTPUTFILENAME		"C:\\Users\\noamz\\Documents\\Data\\testkest1.txt"
#define FRAPOUTPUTFILENAME	"C:\\Users\\noamz\\Documents\\Data\\simfrap.txt"

#define RANDOMSEED

#ifndef ISAAC
#define _CRT_RAND_S
#endif

typedef struct
{
    int			S[MATRIXDIM][MATRIXDIM];
    int		        n;
} SYNAPSE;


int main(int argc, char* argv[])
{
    double		R;
    double		matrix_elements;
    double		k_bind[9], k_unbind[9];
    int			N[MATRIXDIM][MATRIXDIM];
    int			Free_molecules, fr_mol; 
    SYNAPSE		*Synapses;
    int			*SimResult;
    int			*FrapResult;
    int			i, j, k, x, y, n, nbrs;
    int			dt, bleached;
    FILE		*fout, *ffrap;
    time_t		timeval, endval, timeinsec;
    char		*filename = OUTPUTFILENAME;
    char		*frapfilename = FRAPOUTPUTFILENAME;
    unsigned int	randnumber;

    /* initialize ISAAC random number generator */

#ifndef ISAAC
    errno_t		err;
#endif
#ifdef ISAAC
    randctx ctx;

    ctx.randa = ctx.randb = ctx.randc = (ub4) 0;

    /* default seed of pseudo random number generator */

    for (i = 0; i < 256; ++i) 
	ctx.randrsl[i] = (ub4) 0;

    /* unless the seed is changed, the random number streams are identical in all runs.
       The code below generates a different seed on each run, resulting in a different simulation 
       each time it is run. 
       If the seed is not changed, all simulation runs lead to identical results */
#ifdef RANDOMSEED
    for (i = 0; i < 256; ++i) 
    {
	time (&timeinsec);
	randnumber = (timeinsec * i) & 0xFFFFFFF;
	ctx.randrsl[i] = (ub4) randnumber;
    }
#endif
    
    /* if 'flag' is set to 'TRUE', ctx.randrsl is used as a seed */

    randinit(&ctx, TRUE);

#endif

    /* allocate memory for matrices and simulation results */

    matrix_elements = (double) (MATRIXDIM * MATRIXDIM);

    Synapses = (SYNAPSE *) malloc(N_SYNAPSES * sizeof (SYNAPSE));
    if (Synapses == NULL)
    {
	printf ("Memory allocation error\n");
	return 0;
    }
    
    SimResult = (int *) malloc(N_SYNAPSES * (STEPS + 1) * sizeof (int));
    if (SimResult == NULL)
    {
	free ((void *) Synapses);
	printf ("Memory allocation error\n");
	return 0;
    }

    FrapResult = (int *) malloc((STEPS - FRAP_CYCLE + 1) * FRAP_SYNAPSES * sizeof (int));
    if (FrapResult == NULL)
    {
	free ((void *) Synapses);
	free ((void *) SimResult);
	printf ("Memory allocation error\n");
	return 0;
    }
 
    /* clear the memory blocks */

    memset (Synapses, 0, N_SYNAPSES * sizeof (SYNAPSE));
    memset (SimResult, 0, N_SYNAPSES * (STEPS + 1) * sizeof (int));
    memset (FrapResult, 0, (STEPS - FRAP_CYCLE + 1) * FRAP_SYNAPSES * sizeof (int));
  
    /* open the output file */

    if (0 != fopen_s (&fout, filename, "w"))
    {
	printf ("Could not open file %s\n", filename);
	return 0;
    }

    /* open the FRAP output file */

    if (0 != fopen_s (&ffrap, frapfilename, "w"))
    {
	printf ("Could not open FRAP file %s\n", frapfilename);
	return 0;
    }

    /* calculate the cooperativity functions */

    for (i = 0 ; i < 9 ; ++i) 
    {
	k_bind[i] = (double) i * ((double) P_ON / 8.0) + (double) ALPHA;
	k_unbind[8 - i] = (double) i * ((double) P_OFF / 8.0);
    }

    /* initialize the total molecules counter */

    Free_molecules = TOTAL_MOLECULES;

    /* start output to console */

    printf ("Time step:\n");
    time (&timeval);

    /* loop on simulation steps */

    for (dt = 1 ; dt < STEPS ; ++dt)
    {
	/* update the temporary free molecule counter */

	fr_mol = Free_molecules;

	/* loop on synapses */

	for (i = 0 ; i < N_SYNAPSES ; ++i)
	{
	    /* get the number of bound molecules for this synapse from the previous iteration */	

	    n = Synapses[i].n;

	    /* get the number of neighbours for each position in the current matrix */

	    for (j = 0 ; j < MATRIXDIM ; ++j)
	    {
	    	for (k = 0 ; k < MATRIXDIM; ++k)
		{
		    nbrs = 0;

		    /* loop on matrix dimensions, accounting for edges */

		    for (x = max(0, j - 1) ; x <= min(MATRIXDIM - 1, j + 1) ; ++x)
			for (y = max(0, k - 1) ; y <= min(MATRIXDIM - 1, k + 1) ; ++y)
			{
			    if (x == j && y == k)
				continue;
			    nbrs += (Synapses[i].S[x][y] > 0);
			}
		    N[j][k] = nbrs;
		}		
	    }
	    

	    /* photobleach */
	    
	    if (dt == FRAP_CYCLE && i < FRAP_SYNAPSES)
		for (j = 0 ; j < MATRIXDIM ; ++j)
	    	    for (k = 0 ; k < MATRIXDIM; ++k)
			if (Synapses[i].S[j][k])
			    Synapses[i].S[j][k] = 2;

	    /* store # of bleached molecules */

	    if (dt >= FRAP_CYCLE && i < FRAP_SYNAPSES)	
	    {
		bleached = 0;
		for (j = 0 ; j < MATRIXDIM ; ++j)
	    	    for (k = 0 ; k < MATRIXDIM; ++k)
			if (Synapses[i].S[j][k] == 1) //!! reporting non-bleached molecules
			    ++bleached;
		*(FrapResult + i * (STEPS - FRAP_CYCLE) + (dt - FRAP_CYCLE)) = bleached;
	    }
			    
	    /* loop on the matrix of the current synapse */

	    for (j = 0 ; j < MATRIXDIM ; ++j)
	    {
		for (k = 0 ; k < MATRIXDIM; ++k)
		{
		    /* get a random number between 0 and 1 */
#ifndef ISAAC  
		    err = rand_s (&randnumber);
		    if (err != 0)
		    {
			printf ("The rand_s function failed!\n");
			fclose (fout);
    			free ((void *) SimResult);
			return 0;   
		    } 
#endif	
#ifdef ISAAC
		    randnumber = rand_ISAAC (&ctx);
#endif			    
		    R = (double) randnumber / (double) UINT_MAX;
	    		    
		    /* add / remove molecules from the current slot */

		    if (Synapses[i].S[j][k] == 0)	/* empty slot? */
		    {
			/* bind a molecule with a probability based on the cooperativity function 
			   and the number of free molecules and the random number */

			if (R < (k_bind[N[j][k]] * (double) Free_molecules))
			{
			    Synapses[i].S[j][k] = 1;
			    n++;			/* increase the number of bound molecules for the current synapse */
			    fr_mol--;			/* decrease the number of free molecules */
			}
		    }
		    else				/* occupied slot */
		    {
			/* remove a molecule with a probability based on the cooperativity function */

			if (R < k_unbind[N[j][k]])
			{
			    Synapses[i].S[j][k] = 0;
			    n--;			/* decrease the number of bound molecules for the current synapse */
			    fr_mol++;			/* increase the number of free molecules */
			}
		    }
		}
	    }
	    
	    /* store the number of bound molecules for the current synapse */

	    Synapses[i].n = n;
	    
	    /* store the simulation results for the current synapse and time step */

	    *(SimResult + i * STEPS + dt) = n;

	}	/* end of loop on synapses */

	/* update the free molecule counter after the completion of a time step */

	Free_molecules = fr_mol;

	/* print out a result at the end of the current time step iteration */

	printf ("Step %d, effective P_ON %f\n", dt, P_ON * Free_molecules);

    }	    /* end of loop on time steps */

    /* store the data to file  */

    for (i = 0 ; i < N_SYNAPSES ; ++i)
    {
	for (j = 0 ; j < STEPS ; ++j)
	{
	    fprintf (fout, "%d\t", *(SimResult + i * STEPS + j));
	}
	fprintf (fout, "\n");
    }

    /* store the frap data to file  */

    for (i = 0 ; i < FRAP_SYNAPSES ; ++i)
    {
	/* print the number of molecules one step before the FRAP procedure */
	
	fprintf (ffrap, "%d\t", *(SimResult + i * STEPS + FRAP_CYCLE - 1));

	for (j = 0 ; j < (STEPS - FRAP_CYCLE) ; ++j)
	{
	    fprintf (ffrap, "%d\t", *(FrapResult + i * (STEPS - FRAP_CYCLE) + j));
	}
	fprintf (ffrap, "\n");
    }

    /* print the average effective P_ON over the last 10 steps */

    k = n = 0;
    for (i = STEPS - 10 ; i < STEPS; ++i)
    {
	for (j = 0 ; j < N_SYNAPSES ; ++j)
	    k += *(SimResult + j * STEPS + i);
	n++;
    }
    k = k / n;
    printf ("\nAverage P_ON in last 10 steps was %f\n", P_ON * (double)(TOTAL_MOLECULES - k));

    /* close the output file and free the memory */

    fclose (fout);
    fclose (ffrap);
    
    free ((void *) FrapResult);
    free ((void *) SimResult);
    free ((void *) Synapses);

    /* print the completion message and wait for a keystroke before closing */

    time (&endval);
    printf ("\n Done in %d seconds\nPress any key to continue...\n", (int)(endval - timeval));

    _getch();

    return 0;
}

