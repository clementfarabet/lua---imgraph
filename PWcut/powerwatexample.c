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
/*! \file powerwatexample.c
\defgroup graph Graph Example

\brief Example showing how to call the power watershed solver on a graph

<B>Usage:</B> ./ powerwatexample.exe 

\ingroup graph

\author Camille Couprie
*/


#include <graph_utils.h>
#include <sys/types.h>
#include <MSF_RW.h>


/* =============================================================== */
int main(int argc, char **argv) 
/* =============================================================== */
/*Example showing how to call the power watershed solver on the 
following weighted graph

               (Seeded) node 2, value=0  
		       /       \
		     1/         \6
		     /    2.5    \
		   node0 ----- node1
		     \            /
		     5\          /6
		       \        /
	        (Seeded) node 3, value = 1
*/
{
  struct graph <double> *G;
  int j;
  G = (struct graph <double>*)malloc(sizeof(struct graph<double>));
  G->weight_type=1; // 1:double 0:uint32_t

  Allocate_Graph(G, 
		 4,/*Number of nodes */ 
		 2, /*Max Number of seeded nodes */ 
		 5); /*Number of edges */ 
		
  G->P = 1; /*we are solving one problem (multiple problems solved for multilabel segmentation)*/ 
 
  /*Fill the seed values */
  G->S=2;/*there are 2 seeded nodes*/
  G->SeededNodes[0] = 2; // The node 2 is the first seeded node
  G->Labels[0]=0; // its value is 0
  G->SeededNodes[1] = 3; // The node 3 is the second seeded node
  G->Labels[1]=1; // its value is 1

  /*Add the weighted edges*/
  AddEdge<double>(G, 2/*node 2*/, 0/*node 0*/, 1/*weight value*/,0/*edge index*/);
  AddEdge<double>(G, 2, 1, 6, 1);
  AddEdge<double>(G, 0, 1, 2.5, 2);
  AddEdge<double>(G, 0, 3, 5, 3);
  AddEdge<double>(G, 1, 3, 6, 4);

  G->max_weight= 6; /*maximum weight value*/

  /*Solving problem with Power Watersheds*/
     
  PowerWatershed_q2<double>(G);
  
  // Writing results 
  printf("SOLUTION \n");
  for (j = 0; j < 4; j++)
    printf("%f \n", G->Solution[0][j]);
  
 
  Free_Graph(G, G->N);


return 0;
} 




