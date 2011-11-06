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

#ifndef GA_GEODESIC_H
#define GA_GEODESIC_H

#include <graph_utils.h>
#include <union_find.h>
#include <MSF_RW.h>
#include <ccsort.h>
#include <time.h>


/*===========================================================*/
template<class wtype> void element_link_geod_dilate( int n,
						     int p,
						     uint32_t *Fth,
						     wtype *G,
						     wtype *O, 
						     wtype MAX)
/*============================================================*/  
{ 
  int r = element_find(n, Fth);
  
  if (r != p)
    {
      if((G[r] == G[p])||(G[p]>=O[r]))
	{
	  Fth[r] = p;
	  //O[p] = max(O[r],O[p]);
	  if(O[r]>O[p])O[p] = O[r]; 
	}
      else O[p] = MAX;
    } 
}

/* ==================================================================================================== */
template<class wtype> void geodilate_union_find( wtype *H,  /* h : image weights */
						 wtype *O,  /* f : image seeds O : 
						 result of the reconstruction by dilation of h under f */
						 struct graph  <wtype>*G,
						 wtype max_weight) 
						 
						 
/* ===================================================================================================== */
/* reconstruction by dilation of h under f.  Union-find method described by Luc Vicent. */
{
  int k, p, n, e_jp;
  bool * Mrk = (bool*)calloc(G->M,sizeof(bool));
  uint32_t * Fth = (uint32_t*)malloc(G->M*sizeof(uint32_t));
  uint32_t * Es  = (uint32_t*)malloc(G->M*sizeof(uint32_t));
  wtype *W = (wtype*)malloc(G->M*sizeof(wtype));
  TypListechainee * lis; 

  for(k=0;k<G->M;k++) 
    {    
      Fth[k]=k;
      Es[k]=k;
      W[k]=G->Weights[k];
    }
int t1, t2,t3;
  t1=clock();   
  
  // Sort values of F 
  // if (G->weight_type==0)
  //  Es = linsortimageup(H, G->M, max_weight+1);
  //else  
    TriRapideStochastique_inc(H, Es, 0, G->M-1);
    for(k=0;k<G->M;k++) 
      G->Weights[k]=W[k];
    //int i;for(i=0;i<G->M;i++) fprintf(stderr,"%d  %d  %f\n" , i, Es[i], H[i] );
  int r;
  t2=clock();
  // fprintf(stderr, "sort : %.6lf seconds elapsed\n", ((double)t2-t1)/CLOCKS_PER_SEC);
  /* first pass */
  for(k=G->M-1;k>=0;k--)
    {
      p = Es[k];
     
      // parcourt les 4 premiers voisins de p 	
      e_jp = G->Edges[0][p];
      // printf("edge %d, node %d \n neigh. edges", p, e_jp);
       
	  for (lis=G->Neighbors[e_jp]; lis != NULL; lis = lis->suite) 
	    {
	    n= lis->elt;
	    //  printf(" %d ", n);
	      if(Mrk[n]==true) 
		{
		  /* element_link_geod_dilate(n,p, Fth,H, O, max_weight); */
		  r = element_find(n, Fth);
		  
		      if (r != p)
			{
			  if((abs(H[r] - H[p])<0.00001)||(H[p]>=O[r]))
			    {
			      Fth[r] = p;
			      if(O[r]>O[p])O[p] = O[r]; 
			    }
			  else O[p] = max_weight;
			} 
		}	
	      Mrk[p]=true;
	    }
	  // printf("\n");
	
	

      // parcourt les 4 derniers voisins de p
      e_jp = G->Edges[1][p];
     
	  for (lis=G->Neighbors[e_jp]; lis != NULL; lis = lis->suite) 
	    {
	      n= lis->elt;
	      if(Mrk[n]==true) 
		{
		  /* element_link_geod_dilate(n,p, Fth,H, O, max_weight); */
		  //r = element_find_iteratif(n, Fth );
		   r = element_find(n, Fth);
		  
		  if (r != p)
		    {
		      if((H[r] == H[p])||(H[p]>=O[r]))
			{
			  Fth[r] = p;
			  if(O[r]>O[p])O[p] = O[r]; 
			}
		      else O[p] = max_weight;
		    } 
		}	
	    }
	
    }
  free(Mrk);    
  t3=clock();
  //  fprintf(stderr, "Premiere passe reconstruction : %.6lf seconds elapsed\n", ((double)t3-t2)/CLOCKS_PER_SEC);
  /* second pass */
  
  for(k=0;k<G->M;k++)
    {
      p = Es[k];
      if (Fth[p]==(uint32_t)p) // p is root
	{
	  if (O[p]==max_weight) O[p]=H[p];
	}
      else O[p]= O[Fth[p]];
    }
  t1=clock();
  //  fprintf(stderr, "Deuxieme passe reconstruction : %.6lf seconds elapsed\n", ((double)t1-t3)/CLOCKS_PER_SEC);
  free(Es);
  free(Fth);

 free(W);
    
} 

#endif
