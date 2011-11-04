/*Copyright ESIEE (2009) 

Author :
Camille Couprie (c.couprie@esiee.fr)

Contributors : 
Hugues Talbot (h.talbot@esiee.fr)
Leo Grady (leo.grady@siemens.com)
Laurent Najman (l.najman@esiee.fr)

This software contains some image processing algorithms whose purpose is to be
used primarily for research.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

/* templated version
 *
 * Hugues Talbot	 2 Dec 2010
 */

#ifndef CCSORT_HXX
#define CCSORT_HXX

#include <assert.h>
#include <stdlib.h>
#include <stdint.h>
//#include <mtrand64.h> // 64-bit random number generator
#include <stdint.h>
//#include <powerwatsegm.h>
#include<stdint.h>
#include <stdio.h>

template <class wtype> wtype * linsortimageup(wtype *F,uint32_t  N, int maxi);
template <class wtype> wtype * linsortimagedown(wtype *F, uint32_t N, int maxi);
template <class wtype> void TriRapideStochastique_dec (wtype * A, uint32_t *I, long p, long r);
template <class wtype> void TriRapideStochastique_inc (wtype * A, uint32_t *I, long p, long r);

/* =============================================================== */
template <class wtype> long Partitionner_dec(wtype *A, uint32_t * I, long p, long r)
/* =============================================================== */
/*
  partitionne les elements de A entre l'indice p (compris) et l'indice r (compris)
  en deux groupes : ceux <= A[p] et les autres.
*/
{
    wtype t;
    int t1;
    wtype x = A[p];
    long i = p - 1;
    long j = r + 1;
    while (1)
    {
        do j--; while (A[j] < x);
        do i++; while (A[i] > x);
        if (i < j) 
        { 
            t = A[i];
            A[i] = A[j];
            A[j] = t; 
            t1 = I[i];
            I[i] = I[j];
            I[j] = t1; 
        }
        else return j;
    } /* while (1) */   
} /* Partitionner() */

/* =============================================================== */
template <class wtype> long PartitionStochastique_dec (wtype *A, uint32_t * I, long p, long r)
/* =============================================================== */
/*
  partitionne les elements de A entre l'indice p (compris) et l'indice r (compris)
  en deux groupes : ceux <= A[q] et les autres, avec q tire au hasard dans [p,r].
*/
{
  wtype t;
  int t1;
  long q;
  q = p + (rand() % (r - p + 1));
  // q = p + (genrand64_int64() % (r - p + 1)); /* rand must be 64-bit safe, should be OK now */
  //assert((q >= p) && (q <= r)) ; /* you'd be surprised */
  t = A[p];         /* echange A[p] et A[q] */
  A[p] = A[q]; 
  A[q] = t;
  
  t1 = I[p];         /* echange I[p] et I[q] */
  I[p] = I[q]; 
  I[q] = t1;

  return Partitionner_dec(A, I, p, r);
} /* PartitionStochastique() */


/* =============================================================== */
template <class wtype> void TriRapideStochastique_dec (wtype * A, uint32_t *I, long p, long r)
/* =============================================================== */
/* 
  trie les valeurs du tableau A de l'indice p (compris) a l'indice r (compris) 
  par ordre decroissant 
*/
{
  long q; 
  if (p < r)
  {
    q = PartitionStochastique_dec(A, I, p, r);
    TriRapideStochastique_dec (A, I, p, q) ;
    TriRapideStochastique_dec (A, I, q+1, r) ;
  }
} /* TriRapideStochastique() */


//=== increasing




/* =============================================================== */
template <class wtype> long Partitionner_inc(wtype *A, uint32_t * I, long p, long r)
/* =============================================================== */
/*
  partitionne les elements de A entre l'indice p (compris) et l'indice r (compris)
  en deux groupes : ceux <= A[p] et les autres.
*/
{
  wtype t;
  int t1;
  wtype x = A[p];
  long i = p - 1;
  long j = r + 1;
  while (1)
  {
    do j--; while (A[j] > x);
    do i++; while (A[i] < x);
    if (i < j) 
      { 
	t = A[i];
	A[i] = A[j];
	A[j] = t; 
	t1 = I[i];
	I[i] = I[j];
	I[j] = t1; 
      }
    else return j;
  } /* while (1) */   
} 

/* =============================================================== */
template <class wtype> long PartitionStochastique_inc (wtype *A, uint32_t * I, long p, long r)
/* =============================================================== */
/*
  partitionne les elements de A entre l'indice p (compris) et l'indice r (compris)
  en deux groupes : ceux <= A[q] et les autres, avec q tire au hasard dans [p,r].
*/
{
  wtype t;
  int t1;
  long q;
 q = p + (rand() % (r - p + 1));
 // q = p + (genrand64_int64() % (r - p + 1)); /* random number generator must be 64-bit safe */
 //assert((q >= p) && (q <= r)) ; /* you'd be surprised */
  t = A[p];         /* echange A[p] et A[q] */
  A[p] = A[q]; 
  A[q] = t;
  
  t1 = I[p];         /* echange I[p] et I[q] */
  I[p] = I[q]; 
  I[q] = t1;

  return Partitionner_inc(A, I, p, r);
} 

/* =============================================================== */
template <class wtype> void TriRapideStochastique_inc (wtype * A, uint32_t * I, long p, long r)
/* =============================================================== */
/* 
  trie les valeurs du tableau A de l'indice p (compris) a l'indice r (compris) 
  par ordre croissant 
*/
{
  long q; 
  if (p < r)
  {
    q = PartitionStochastique_inc(A, I, p, r);
    TriRapideStochastique_inc (A, I, p, q) ;
    TriRapideStochastique_inc (A, I, q+1, r) ;
  }
} /* TriRapideStochastique() */


/* ==================================== */
template <class wtype> wtype * linsortimageup(wtype *F, uint32_t N, int maxi)
/* ==================================== */
/*
  Tri par denombrement - cf. Cormen & al., "Introduction a l'algorithmique"
  version pour une image sur 32 bits
  F: l'image
  N: le nombre de pixels
  retourne: un tableau d'indices de pixels (taille N) 
            dans l'ordre croissant des valeurs
*/
#undef F_NAME
#define F_NAME "linsortimageup"
{
    int i, j, k, H[maxi+1];
    uint32_t *T = (uint32_t *)calloc(1,N * sizeof(uint32_t));
    if (T == NULL)
    {   fprintf(stderr, "%s() : calloc failed for T\n", F_NAME);
        return NULL;
    }
    for (i = 0; i < maxi+1; i++)
        H[i] = 0;   // initialise l'histogramme
    for (i = 0; i < (int)N; i++)
        H[F[i]] += 1; // calcule l'histogramme
    
    j = H[0];
    H[0] = 0;
    // calcule l'histogramme cumule
    for (i = 1; i < maxi+1; i++) {
        k = H[i];
        H[i] = j;
        j += k;
    }
    for (i = 0; i < (int)N; i++)               // tri lineaire
    {
        k = F[i];
	j = H[k];
	T[j] = i;
	H[k] += 1;
    }
    return T;
} /* linsortimageup() */

/* ==================================== */
template <class wtype> wtype * linsortimagedown(wtype *F, uint32_t N, int maxi)
/* ==================================== */
/*
  Tri par denombrement - cf. Cormen & al., "Introduction a l'algorithmique"
  version pour une image sur 8 bits
  F: l'image
  N: le nombre de pixels
  retourne: un tableau d'indices de pixels (taille N) 
            dans l'ordre decroissant des valeurs
*/
#undef F_NAME
#define F_NAME "linsortimagedown"
{
    int i, j, k, H[maxi+1];
    uint32_t *T = (uint32_t *)calloc(1,N* sizeof(uint32_t));
    if (T == NULL)
    {   fprintf(stderr, "%s() : calloc failed for T\n", F_NAME);
        return NULL;
    }
    for (i = 0; i < maxi+1; i++)
        H[i] = 0;   // initialise l'histogramme
  
    for (i = 0; i < N; i++) 
        H[F[i]] += 1; // calcule l'histogramme
  
    j = H[maxi]; 
    H[maxi] = 0;
    // calcule l'histogramme cumule
    for (i = maxi-1; i >= 0; i--) 
    { 
        k = H[i];
        H[i] = j; 
        j += k; 
    }
    for (i = 0; i < (int)N; i++)               // tri lineaire
    {
        k = F[i]; 
        j = H[k];
        T[j] = i;
        H[k] += 1;
    }
    return T;
} /* linsortimagedown() */



#endif // CCSORT_HXX
