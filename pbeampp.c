/**************************************************************************
PBEAMPP.C of ZIB optimizer MCF, SPEC version

This software was developed at ZIB Berlin. Maintenance and revisions
solely on responsibility of Andreas Loebel

Dr. Andreas Loebel
Ortlerweg 29b, 12207 Berlin

Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)
Scientific Computing - Optimization
Takustr. 7, 14195 Berlin-Dahlem

Copyright (c) 1998-2000 ZIB.
Copyright (c) 2000-2002 ZIB & Loebel.
Copyright (c) 2003-2005 Andreas Loebel.
**************************************************************************/
/*  LAST EDIT: Sun Nov 21 16:22:04 2004 by Andreas Loebel (boss.local.de)  */
/*  $Id: pbeampp.c,v 1.10 2005/02/17 19:42:32 bzfloebe Exp $  */

#define K 300
#define B 50

#include "pbeampp.h"

#ifdef _PROTO_
int bea_is_dual_infeasible(arc_t *arc, cost_t red_cost)
#else
int bea_is_dual_infeasible(arc, red_cost)
arc_t *arc;
cost_t red_cost;
#endif
{
    return ((red_cost < 0 && arc->ident == AT_LOWER) || (red_cost > 0 && arc->ident == AT_UPPER));
}

typedef struct basket
{
    arc_t *a;
    cost_t cost;
    cost_t abs_cost;
} BASKET;

static long basket_size;
static BASKET basket[B + K + 1];
static BASKET *perm[B + K + 1];

#ifdef _PROTO_
void sort_basket(long min, long max)
#else
void sort_basket(min, max) long min, max;
#endif
{
    long l, r;
    cost_t cut;
    BASKET *xchange;

    l = min;
    r = max;

    cut = perm[(long)((l + r) / 2)]->abs_cost;

    do
    {
        while (perm[l]->abs_cost > cut)
            l++;
        while (cut > perm[r]->abs_cost)
            r--;

        if (l < r)
        {
            xchange = perm[l];
            perm[l] = perm[r];
            perm[r] = xchange;
        }
        if (l <= r)
        {
            l++;
            r--;
        }

    } while (l <= r);

    if (min < r)
        sort_basket(min, r);
    if (l < max && l <= B)
        sort_basket(l, max);
}

static long nr_group;
static long group_pos;

static long initialize = 1;

#ifdef _PROTO_
arc_t *primal_bea_mpp(long m, arc_t *arcs, arc_t *stop_arcs,
    cost_t *red_cost_of_bea)
#else
arc_t *primal_bea_mpp(m, arcs, stop_arcs, red_cost_of_bea) long m;
arc_t *arcs;
arc_t *stop_arcs;
cost_t *red_cost_of_bea;
#endif
{
    long old_group_pos;

    if (initialize)
    {
        for (long i = 1; i < K + B + 1; i++)
            perm[i] = &(basket[i]);
        nr_group = ((m - 1) / K) + 1;
        group_pos = 0;
        basket_size = 0;
        initialize = 0;
    }
    else
    {

        /*****************************************************************/
        /******************* BEGINNING FIRST FOR LOOP ********************/
        /*****************************************************************/

        int next_array[B + K + 1];
        int next_increased_array[B + K + 1];
        int next_increase_count[5];
        arc_t *arc_array[B + K + 1];
        cost_t red_cost_array[B + K + 1];
        int min = B < basket_size ? B : basket_size;
        int chunk_size = (min - 2 + 1) / 4;

#pragma omp parallel for
        for (long j = 0; j < 4; j++)
        {
            long chunk_start = 2 + j * chunk_size;
            long chunk_end = j == 3 ? min : 2 + j * chunk_size + chunk_size - 1;
            next_increase_count[j + 1] = 0;
            for (long i = chunk_start; i <= chunk_end; i++)
            {
                arc_t *arc = perm[i]->a;
                cost_t red_cost = arc->cost - arc->tail->potential + arc->head->potential;
                int ident = arc->ident;

                arc_array[i] = arc;
                red_cost_array[i] = red_cost;

                int predicate = (red_cost < 0 && ident == AT_LOWER) || (red_cost > 0 && ident == AT_UPPER);
                next_increased_array[i] = predicate ? 1 : 0;
                next_increase_count[j + 1] += predicate;
            }
        }

        next_increase_count[0] = 0;
        for (long j = 1; j <= 4; j++)
        {
            next_increase_count[j] += next_increase_count[j - 1];
        }

#pragma omp parallel for
        for (long j = 0; j < 4; j++)
        {
            long next = 0;
            long global_next;
            long chunk_start = 2 + j * chunk_size;
            long chunk_end = j == 3 ? min : 2 + j * chunk_size + chunk_size - 1;
            for (long i = chunk_start; i <= chunk_end; i++)
            {

                if (next_increased_array[i])
                {
                    next++;
                    global_next = next_increase_count[j] + next;
                    BASKET *current_prem = perm[global_next];
                    cost_t red_cost = red_cost_array[i];
                    current_prem->a = arc_array[i];
                    current_prem->cost = red_cost;
                    current_prem->abs_cost = ABS(red_cost);
                }
            }
        }

        basket_size = next_increase_count[4];
    }

    /*****************************************************************/
    /********************* END FIRST FOR LOOP ************************/
    /*****************************************************************/

    old_group_pos = group_pos;

    /*****************************************************************/
    /********************** BEGINNING GOTO LOOP **********************/
    /*****************************************************************/

    // int new_group_pos = group_pos;
    // int new_group_pos_set = 0;

    // int chunk_size = nr_group / 4;
    // int basket_size_increase_count[5];

    // int basket_size_increased_array_size = (((stop_arcs - arcs) / nr_group) + 1) * nr_group;

    // // int bla = stop_arcs - arcs;
    // // printf("bla: %d\n", basket_size_increased_array_size / 8);

    // int *basket_size_increased_array = malloc(basket_size_increased_array_size * sizeof(int));
    // for (long j = 0; j < 4; j++)
    // {
    //   // printf("HELLO!!!!!\n");
    //   long chunk_start = j * chunk_size;
    //   long chunk_end = j == 3 ? nr_group : j * chunk_size + chunk_size;
    //   long real_chunk_size = chunk_end - chunk_start + 1;
    //   basket_size_increase_count[j + 1] = 0;
    //   for (long group_pos_index = 0; group_pos_index < real_chunk_size; group_pos_index++)
    //   {
    //     long current_group_pos = (group_pos + (j * chunk_size) + group_pos_index) % nr_group;
    //     arc_t *arc = arcs + current_group_pos;

    //     int index = 0;
    //     for (; arc < stop_arcs; arc += nr_group)
    //     {
    //       long basket_size_increased_index = nr_group * current_group_pos + index;
    //       // if (basket_size_increased_index < B)
    //       // {
    //       // printf("HELsdfsdfsdfLO2!!!!!\n");
    //       if (arc->ident > BASIC)
    //       {
    //         cost_t red_cost = arc->cost - arc->tail->potential + arc->head->potential;
    //         int predicate = bea_is_dual_infeasible(arc, red_cost);

    //         basket_size_increased_array[basket_size_increased_index] = predicate;
    //         basket_size_increase_count[j + 1] += predicate;
    //       }
    //       // }
    //       // else if (new_group_pos_set == 0)
    //       // {
    //       //   new_group_pos_set = 1;
    //       //   new_group_pos = current_group_pos;
    //       // }

    //       index++;
    //     }
    //   }
    // }

    // basket_size_increase_count[0] = 0;
    // for (long j = 1; j <= 4; j++)
    // {
    //   basket_size_increase_count[j] += basket_size_increase_count[j - 1];
    // }

    // for (long j = 0; j < 4; j++)
    // {
    //   long global_basket_size = basket_size_increase_count[j];
    //   long chunk_start = j * chunk_size;
    //   long chunk_end = j == 3 ? nr_group : j * chunk_size + chunk_size;
    //   basket_size_increase_count[j + 1] = 0;

    //   for (long group_pos_index = 0; group_pos_index < chunk_size; group_pos_index++)
    //   {
    //     long current_group_pos = (group_pos + (j * chunk_size) + group_pos_index) % nr_group;
    //     arc_t *arc = arcs + current_group_pos;
    //     int index = 0;
    //     long local_basket_size = 0;
    //     for (; arc < stop_arcs; arc += nr_group)
    //     {
    //       long basket_size_increased_index = nr_group * current_group_pos + index;
    //       if (global_basket_size < B)
    //       {
    //         if (arc->ident > BASIC)
    //         {

    //           if (basket_size_increased_array[basket_size_increased_index])
    //           {
    //             local_basket_size++;
    //             global_basket_size = basket_size_increase_count[j] + local_basket_size;
    //             BASKET *current_prem = perm[global_basket_size];
    //             cost_t red_cost = arc->cost - arc->tail->potential + arc->head->potential;
    //             current_prem->a = arc;
    //             current_prem->cost = red_cost;
    //             current_prem->abs_cost = ABS(red_cost);
    //           }
    //         }
    //       }
    //       else if (new_group_pos_set == 0)
    //       {
    //         new_group_pos_set = 1;
    //         new_group_pos = current_group_pos;
    //       }

    //       index++;
    //     }
    //   }
    // }

    // group_pos = new_group_pos;
    // basket_size = basket_size_increase_count[4];

    /*****************************************************************/
    /************************* END GOTO LOOP *************************/
    /*****************************************************************/

    arc_t *arc;
    cost_t red_cost;

NEXT:
    /* price next group */
    arc = arcs + group_pos;
    for (; arc < stop_arcs; arc += nr_group)
    {
        if (arc->ident > BASIC)
        {
            /* red_cost = bea_compute_red_cost( arc ); */
            red_cost = arc->cost - arc->tail->potential + arc->head->potential;
            if (bea_is_dual_infeasible(arc, red_cost))
            {
                basket_size++;
                perm[basket_size]->a = arc;
                perm[basket_size]->cost = red_cost;
                perm[basket_size]->abs_cost = ABS(red_cost);
            }
        }
    }

    if (++group_pos == nr_group)
        group_pos = 0;

    if (basket_size < B && group_pos != old_group_pos)
        goto NEXT;

    if (basket_size == 0)
    {
        initialize = 1;
        *red_cost_of_bea = 0;
        return NULL;
    }

    sort_basket(1, basket_size);

    *red_cost_of_bea = perm[1]->cost;
    return (perm[1]->a);
}
