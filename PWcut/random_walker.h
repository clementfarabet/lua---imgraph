/*Copyright ESIEE (2009) 

Author :
 
Camille Couprie (c.couprie@esiee.fr)

Contributors : 
Leo Grady (leo.grady@siemens.com)
Hugues Talbot (h.talbot@esiee.fr)
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

//#include <powerwatsegm.h>
#ifdef __cplusplus
extern "C" {
#endif

#include <cs.h>
#include <stdint.h>
  typedef float DBL_TYPE;

/*******************************************************
Function RandomWalker computes the solution to the Dirichlet problem (RW potential function) 
on a general graph represented by an edge list, given boundary conditions (seeds, etc.)
*******************************************************/
extern bool RandomWalker(uint32_t** index_edges,           /* list of edges */
			 int M,                       /* number of edges */
			 uint32_t * index,                 /* list of vertices */
			 uint32_t * indic_vertex,          /* boolean array of vertices*/
			 int N,                       /* number of vertices */  
			 uint32_t* index_seeds,            /* list of nodes that are seeded*/
			 DBL_TYPE ** boundary_values,  /* associated values for seeds (labels)*/
			 int numb_boundary,           /* number of seeded nodes */
			 int nb_labels,               /* nb of different labels*/
			 DBL_TYPE ** proba) ;         /* solution to the Dirichlet problem*/

	       
  #ifdef __cplusplus
  }
#endif
