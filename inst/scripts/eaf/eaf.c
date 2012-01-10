/*************************************************************************

 eaf: Computes the empirical attainment function (EAF) from a number
 of approximation sets.

 ---------------------------------------------------------------------

                    Copyright (c) 2006, 2007, 2008
                  Carlos Fonseca <cmfonsec@ualg.pt>
          Manuel Lopez-Ibanez <manuel.lopez-ibanez@ulb.ac.be>

 This program is free software (software libre); you can redistribute
 it and/or modify it under the terms of the GNU General Public License
 as published by the Free Software Foundation; either version 2 of the
 License, or (at your option) any later version.

 This program is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, you can obtain a copy of the GNU
 General Public License at:
                 http://www.gnu.org/copyleft/gpl.html
 or by writing to:
           Free Software Foundation, Inc., 59 Temple Place,
                 Suite 330, Boston, MA 02111-1307 USA

 ----------------------------------------------------------------------


*************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <limits.h>
#include "eaf.h"

#ifndef DEBUG
#define DEBUG 0
#endif

#ifdef R_PACKAGE
#include <R.h>
#define EAF_MALLOC(WHAT, NMEMB, SIZE)                           \
    WHAT = malloc (NMEMB * SIZE);                               \
    if (!WHAT) { error("malloc failed: %s", #WHAT); }
#define fatalprintf(X) error(X)
#else
#define EAF_MALLOC(WHAT, NMEMB, SIZE)                           \
    WHAT = malloc (NMEMB * SIZE);                               \
    if (!WHAT) { perror (#WHAT ":"); exit (EXIT_FAILURE); }
#define fatalprintf(X) \
    do {                                                                       \
        errprintf (X); exit (EXIT_FAILURE);                                    \
    } while(0)
#endif

static int compare_x_asc (const void *p1, const void *p2)
{
    objective_t x1 = **(objective_t **)p1;
    objective_t x2 = **(objective_t **)p2;
	
    return (x1 != x2)
        ? ((x1 < x2) ? -1 : 1)
        : 0;
}

static int compare_y_desc (const void *p1, const void *p2)
{
	objective_t y1 = *(*(objective_t **)p1+1);
	objective_t y2 = *(*(objective_t **)p2+1);
	
	return (y1 != y2)
            ? ((y1 > y2) ? -1 : 1)
            : 0;
}

eaf_t * eaf_create (int nobj, int nruns, int npoints)
{
    eaf_t *eaf;
    EAF_MALLOC (eaf, 1, sizeof(eaf_t));
    eaf->nobj = nobj;
    eaf->nruns = nruns;
    eaf->size = 0;
    eaf->maxsize = 1 + npoints / 4; /* Maximum is npoints, but normally it
                                       will be smaller, so at most two
                                       realloc will occur.  */
    EAF_MALLOC (eaf->data, nobj * eaf->maxsize, sizeof(objective_t));
    EAF_MALLOC (eaf->attained, nruns * eaf->maxsize, sizeof(int));
    return eaf;
}

void eaf_delete (eaf_t * eaf)
{
    free (eaf->data);
    free (eaf->attained);
    free (eaf);
}

static void
eaf_store_point (eaf_t * eaf, objective_t x, objective_t y, 
                 const int *save_attained)
{
    objective_t * pos;
    const int nruns = eaf->nruns;
    const int nobj = eaf->nobj;

    if (eaf->size == eaf->maxsize) {
        assert (eaf->size < INT_MAX / 2);
        eaf->maxsize = eaf->maxsize * 2;
        eaf->attained = realloc (eaf->attained, 
                                 sizeof(int) * nruns * eaf->maxsize);
        assert(eaf->attained);
        eaf->data = realloc (eaf->data,
                             sizeof(objective_t) * nobj * eaf->maxsize);
        assert(eaf->data);
    }
    pos = eaf->data + nobj * eaf->size;
    pos[0] = x;
    pos[1] = y;
    memcpy (eaf->attained + nruns * eaf->size, save_attained,
            sizeof(int) * nruns);
    eaf->size++;
}

static void
eaf_print_line (FILE *coord_file, FILE *indic_file, FILE *diff_file, 
                objective_t x, objective_t y, 
                const int *save_attained, int nruns)
{
    int count1 = 0;
    int count2 = 0;
    int k;

    if (coord_file) {
        fprintf (coord_file, point_printf_format "\t" point_printf_format,
                 (double)x, (double)y);
        fprintf (coord_file, 
                 (coord_file == indic_file) || (coord_file == diff_file)
                 ? "\t" : "\n");
    }

    if (indic_file) {
        fprintf (indic_file, "%d", 
                 save_attained[0] ? (count1++,1) : 0);
        for (k = 1; k < nruns/2; k++) 
            fprintf (indic_file, "\t%d", 
                    save_attained[k] ? (count1++,1) : 0);
        for (k = nruns/2; k < nruns; k++)
            fprintf (indic_file, "\t%d", 
                     save_attained[k] ? (count2++,1) : 0);

        fprintf (indic_file, (indic_file == diff_file) ? "\t" : "\n");
    } else if (diff_file) {
        for (k = 0; k < nruns/2; k++) 
            if (save_attained[k]) count1++;
        for (k = nruns/2; k < nruns; k++)
            if (save_attained[k]) count2++;
    }

    if (diff_file)
        fprintf (diff_file,"%d\t%d\n", count1, count2);
}

void
eaf_print (eaf_t * eaf, FILE *coord_file,  FILE *indic_file, FILE *diff_file)
{
    int i;
    for (i = 0; i < eaf->size; i++) {
        objective_t *p = eaf->data + i * eaf->nobj;
        eaf_print_line (coord_file, indic_file, diff_file,
                        p[0], p[1], eaf->attained + i * eaf->nruns, eaf->nruns);
    }
}

/* 
   attsurf: compute attainment surfaces from points in objective space,
            using dimension sweeping.

   Input arguments:
        data : a pointer to the data matrix, stored as a linear array of 
               objective_t in row major order.
        nobj : the number of objectives (must be 2 in this implementation).
        cumsize : an array containing the cumulative number of rows in each 
                  non-dominated front (must be non-decreasing).
        nruns :	the number of independent non-dominated fronts.
        attlevel : an array containing the attainment levels to compute.
        nlevel : number of attainment levels to compute.
        coord_file  : stream to write the resulting attainment surfaces.
        indic_file  : stream to write the resulting attainment indices.
        diff_file   : stream to write the difference between the the 
                      first half and the second half of nruns.
*/

eaf_t **
attsurf (const objective_t *data, int nobj, const int *cumsize, int nruns,
         const int *attlevel, const int nlevels)
{
    eaf_t **eaf;
    const objective_t **datax, **datay; /* used to access the data sorted
                                           according to x or y */
    
    int ntotal = cumsize[nruns - 1]; /* total number of points in data */
    int *runtab;	
    int *attained, nattained, *save_attained;
    int k, j, l;

    if (nobj != 2) {
        fatalprintf ("this implementation only supports two dimensions.\n");
    }

    /* Access to the data is made via two arrays of pointers: ix, iy
       These are sorted, to allow for dimension sweeping */

    datax = malloc (ntotal * sizeof(objective_t *));
    datay = malloc (ntotal * sizeof(objective_t *));

    for (k = 0; k < ntotal ; k++)
        datax[k] = datay[k] = data + nobj * k;

#if DEBUG > 1
    fprintf (stderr, "Original data:\n");
    fprint_set (stderr, datax, ntotal);
#endif

    qsort (datax, ntotal, sizeof(*datax), &compare_x_asc);
    qsort (datay, ntotal, sizeof(*datay), &compare_y_desc);

#if DEBUG > 1
    fprintf (stderr, "Sorted data (x):\n");
    fprint_set (stderr, datax, ntotal);
    fprintf (stderr, "Sorted data (y):\n");
    fprint_set (stderr, datay, ntotal);
#endif

    /* Setup a lookup table to go from a point to the approximation
       set (run) to which it belongs.  */

    runtab = malloc (ntotal * sizeof(int));
    for (k = 0, j = 0; k < ntotal; k++) {
        if (k == cumsize[j])
            j++;
        runtab[k] = j;
    }

#if DEBUG > 1
    fprintf (stderr, "Runtab:\n");
    for (k = 0; k < ntotal; k++) 
        fprintf (stderr, "%6d: %6d\n", k, runtab[k]);
#endif

    /* Setup tables to keep attainment statistics. In particular,
       save_attained is needed to cope with repeated values on the same
       axis. */

    attained = malloc (nruns * sizeof(int));
    save_attained = malloc (nruns * sizeof(int));
    eaf = malloc(nlevels * sizeof(eaf_t*));

    for (l = 0; l < nlevels; l++) {
        eaf[l] = eaf_create (nobj, nruns, ntotal);
        int level = attlevel[l];
        int x = 0;
        int y = 0;
        int run;

        nattained = 0;
        for (k = 0; k < nruns; attained[k++] = 0);

        /* Start at upper-left corner */
        run = runtab[(datax[x] - data) / nobj];
        attained[run]++;
        nattained++;

        do {
            /* Move right until desired attainment level is reached */
            while (x < ntotal - 1 && 
                   (nattained < level || datax[x][0] == datax[x+1][0])) {
                x++;
                if (datax[x][1] <= datay[y][1]) {
                    run = runtab[(datax[x] - data)/nobj];
                    if (!attained[run])
                        nattained++;
                    attained[run]++;
                }
            }
#if DEBUG > 1
            for (k = 0; k < nruns; k++)
                fprintf (stderr, "%d ", attained[k]);
            fprintf (stderr, "\n");
#endif

            if (nattained < level)
                continue;

            /* Now move down until desired attainment level is no
               longer reached.  */
            do {
                /* If there are repeated values along the y axis,
                   we need to remember where we are.  */
                /*save_nattained = nattained;*/
                memcpy (save_attained, attained, nruns * sizeof(*attained));

                do {
                    if (datay[y][0] <= datax[x][0]) {
                        run = runtab[(datay[y] - data)/nobj];
                        attained[run]--;
                        if (!attained[run])
                            nattained--;
                    }
#if DEBUG > 1
                    for (k = 0; k < nruns; k++)
                        fprintf (stderr, "%d ", attained[k]);
                    fprintf (stderr, "\n");
#endif
                    y++;
                } while (y < ntotal && datay[y][1] == datay[y - 1][1]);
            } while (nattained >= level && y < ntotal);

            assert (nattained < level);

            eaf_store_point (eaf[l], datax[x][0], datay[y-1][1],
                             save_attained);

        } while (x < ntotal - 1 && y < ntotal);
    }
    free(save_attained);
    free(attained);
    free(runtab);
    free(datay);
    free(datax);

    return eaf;
}
