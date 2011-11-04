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

/*! \file MSF_RW.c 
\brief Implements the power watershed algorithm with q=2 
\ingroup morpho

\author Camille Couprie
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
//#include <mccodimage.h>
//#include <mcimage.h>
//#include <mcutil.h>
//#include <cccodimage.h>
#include <union_find.h>
#include <MSF_RW.h>
#include <mclifo.h>
#include <sys/time.h>
#include <time.h>


/*========================================================================================================*/
void memory_allocation_PW(Lifo ** LIFO,       /* stack for plateau labeling */
			  Lifo ** LCP,        /* list of the edges belonging to a plateau */        
			  bool ** indic_E,    /* indicator for edges */
			  bool ** indic_P,    /* indicator for edges */
			  uint32_t ** indic_VP,    /* indicator for nodes */
			  uint32_t ** Rnk,         /* array needed for union-find efficiency */
			  uint32_t ** Fth,         /* array for storing roots of merged nodes trees */
			  uint32_t ** local_seeds, /* array for storing the index of seeded nodes of a plateau */
			  uint32_t ** LCVP ,       /* list of vertices of a plateau */
			  uint32_t ** NEs,         /* array of sorted edges according their weight */
			  int N,              /* number of vertices */  
			  int M)              /* number of edges*/
/*==========================================================================================================*/
{
#undef F_NAME
#define F_NAME "memory_allocation_PW"
 
  *LIFO = CreeLifoVide(M);
  if (*LIFO == NULL) { fprintf(stderr, "%s : CreeLifoVide failed\n", F_NAME); exit(0); }

  *LCP = CreeLifoVide(M);
  if (*LCP == NULL) { fprintf(stderr, "%s : CreeLifoVide failed\n", F_NAME); exit(0); }
 
  *indic_E = (bool*)calloc(M ,sizeof(bool)); 
  if (*indic_E == NULL) { fprintf(stderr, "%s: calloc indic_E failed\n", F_NAME); exit(0);}
  *indic_P = (bool*)calloc(M ,sizeof(bool)); 
  if (*indic_P == NULL) { fprintf(stderr, "%s: calloc indic_P failed\n", F_NAME); exit(0);}
  *indic_VP = (uint32_t*)calloc(N ,sizeof(uint32_t));
  if (*indic_VP == NULL) { fprintf(stderr, "%s: calloc indic_VP failed\n", F_NAME); exit(0);}
  *Rnk = (uint32_t*)calloc(N, sizeof(uint32_t));
  if (*Rnk == NULL) { fprintf(stderr, "%s: malloc Rnk failed\n", F_NAME); exit(0);}
  *Fth = (uint32_t*)malloc(N*sizeof(uint32_t));
  if (*Fth == NULL) { fprintf(stderr, "%s: malloc Fth failed\n", F_NAME); exit(0);}
  
  *local_seeds = (uint32_t*)malloc(N*sizeof(uint32_t));
  if (local_seeds == NULL) { fprintf(stderr, "%s: malloc local_seeds failed\n", F_NAME); exit(0);}
  
  *LCVP = (uint32_t*)malloc(N*sizeof(uint32_t)); // vertices of a plateau. 
  if (*LCVP == NULL) { fprintf(stderr, "%s: malloc LCVP failed\n", F_NAME); exit(0);}
 
  
  *NEs = ( uint32_t*)malloc(M*sizeof( uint32_t)); 
  if (*NEs == NULL) { fprintf(stderr, "%s: malloc NEs failed\n", F_NAME); exit(0);}
}



/*===================================================================================*/
void merge_node (int e1,            /* index of node 1 */
		 int e2,            /* index of node 2 */
		 uint32_t *Rnk,         /* array needed for union-find efficiency */
		 uint32_t *Fth,          /* array for storing roots of merged nodes trees */
		 DBL_TYPE ** proba, /* array for storing the result x */
		 int nb_labels)     /* nb of labels */
/*===================================================================================*/
/* update the result, Rnk and Fth arrays when 2 nodes are merged */
{
  int k,re1, re2;
  re1 = element_find(e1, Fth );
  re2 = element_find(e2, Fth );
 
  if ((re1 != re2) && (!(proba[0][re1]>=0 && proba[0][re2]>=0))) 
    {
      element_link(re1,re2, Rnk, Fth);
      //  printf("link %d %d \n", re1,re2); 
      if (proba[0][re2]>=0 && proba[0][re1]<0) 
	for(k=0;k<nb_labels;k++)
	  proba[k][re1]= proba[k][re2];
      else if (proba[0][re1]>=0 && proba[0][re2]<0)
	for(k=0;k<nb_labels;k++)
	  proba[k][re2]= proba[k][re1];
    }
}


