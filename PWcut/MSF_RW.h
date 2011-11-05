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

#ifndef MSF_RW_H
#define MSF_RW_H 


#include <math.h>
#include <union_find.h>
#include <ccsort.h>
#include <random_walker.h>
#include <time.h>
#include<graph_utils.h>
#include <memory.h>
#include <mclistechainee.h>
#include <mclifo.h>
#include <gageodesic.h>
//epsilon plateau
#define eps 0.000000001


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
			  int M)   ;           /* number of edges*/

extern void merge_node (int e1,            /* index of node 1 */
			int e2,            /* index of node 2 */
			uint32_t * Rnk,         /* array needed for union-find efficiency */
			uint32_t *Fth,          /* array for storing roots of merged nodes trees */
			DBL_TYPE ** proba, /* array for storing the result x */
			int nb_labels) ;    /* nb of labels */


template<class wtype>bool equal (wtype a, wtype b);

/*==================================================================================================================*/
template<class wtype>void PowerWatershed_q2(struct graph <wtype> *G) /* informations on the graph : nb of vertices, nodes, etc... */
/*==================================================================================================================*/
/* IN/OUT : potential/proba map x minimizing Epq, q=2, p-> infty*/
/*returns the result x (named here  G->Solution) of the energy minimization : 
min_x lim_p_inf sum_e_of_E w_{ij}^p |x_i-x_j|^2 
subject to data constraints */
{  
#undef F_NAME
#define F_NAME "PowerWatershed_q2"

  int i, j, k, x, y, e1, e2, re1,re2, p, xr, e_jx;
  int N;
  double val;
  int nb_vertices, e_max, Ne_max, nb_edges, Nnb_edges;
  int cpt_aretes, Ncpt_aretes;
  bool success, different_seeds;
  wtype wmax;
  Lifo *LIFO, *LCP;
  bool * indic_E, *indic_P; 
  uint32_t *indic_VP, *LCVP ;
  uint32_t *local_seeds;
  uint32_t *NEs, *Es, *Rnk, *Fth;
  uint32_t ** edgesLCP;
  DBL_TYPE ** local_labels;
  wtype * sorted_weights;
  TypListechainee * lis;
  wtype max_weight = G->max_weight;

 /* Prepare markers for reconstruction */
   for (j=0;j<G->S;j++)
     {
       i = G->SeededNodes[j];
       for (lis=G->Neighbors[i]; lis != NULL; lis = lis->suite) 
	 {
	   k=lis->elt;
	   G->RecWeights[k]= max_weight;//G->Weights[k];
	 }
     }
   //Print_Graph_Content(G);
   geodilate_union_find<wtype>( G->Weights, G->RecWeights, G, max_weight);

   
   // for(i=0;i<G->M;i++) fprintf(stderr,"%d  %f  %f\n" , i, G->Weights[i], G->RecWeights[i] );

  /* Initialisation of the solution with seeds values*/
  Allocate_Graph_Solution(G);
  for (i=0;i<G->S;i++)
    if (G->Labels[i]==1)
      G->Solution[0][G->SeededNodes[i]] = 1;
    else G->Solution[0][G->SeededNodes[i]] = 0;
  
 
 
  N = G->N; 

  memory_allocation_PW( &LIFO, &LCP, &indic_E,  &indic_P,  &indic_VP, &Rnk,  &Fth, &local_seeds, &LCVP, &NEs, G->N, G->M);
 
  edgesLCP = (uint32_t**)malloc(2*sizeof(uint32_t*));
  if (edgesLCP == NULL) { fprintf(stderr, "%s: malloc edgesLCP failed\n", F_NAME); exit(0);}
  for(k=0;k<2;k++) 
    {
      edgesLCP[k] = (uint32_t*)malloc(G->M*sizeof(uint32_t));
      if ( edgesLCP[k]== NULL) { fprintf(stderr, "%s: malloc edgesLCP failed\n", F_NAME); exit(0);}
    }

  for(k=0;k<G->N;k++) {Fth[k]=k;}

  local_labels = (DBL_TYPE**)malloc((G->P)*sizeof(DBL_TYPE*));
  if (local_labels == NULL) { fprintf(stderr, "%s: malloc local_labels failed\n", F_NAME); exit(0);}
  for (i=0; i<G->P; i++)
    {
      local_labels[i]= (DBL_TYPE *)malloc(G->N *sizeof(DBL_TYPE));
      if ( local_labels[i]== NULL)	
	{ fprintf(stderr, "%s: malloc local_labels failed\n", F_NAME); exit(0);}
    }
 
  // fprintf(stderr,"maxi = %d \n",max_weight);
  // if (G->weight_type==0)
  // Es = linsortimagedown<wtype>(G->RecWeights, G->M, max_weight+1);
  //else  
  Es = (uint32_t*)malloc(G->M*sizeof(uint32_t));
  wtype * TW = (wtype*)malloc(G->M*sizeof(wtype));
  for(k=0;k<G->M;k++) 
    {
      Es[k]=k;
      TW[k]=G->RecWeights[k];
    }
  TriRapideStochastique_dec<wtype>(G->RecWeights, Es, 0, G->M-1);
 for(k=0;k<G->M;k++) 
   G->RecWeights[k]=TW[k];
   
 free(TW);
  

 
  cpt_aretes = 0;
  Ncpt_aretes = 0;
 
 /* beginning of main loop */   
  while (cpt_aretes < G->M)
    {
      do 
	{
	  e_max=Es[cpt_aretes];
	  cpt_aretes=cpt_aretes+1;
	  if(cpt_aretes==G->M) break;
	}while(indic_E[e_max]==true);
     
      if(cpt_aretes==G->M) break;
      //fprintf(stderr, "apres\n");
      //1. Computing the edges of the plateau LCP linked to the edge e_max
      LifoPush(LIFO, e_max);
      indic_P[e_max]=true;
      indic_E[e_max]=true;
      LifoPush(LCP, e_max);
      nb_vertices=0;
      nb_edges = 0;
      wmax = G->RecWeights[e_max];
      //  printf("emax=%d, wmax=%f\n", e_max, wmax);
      //t1=clock();
      // 2. putting the edges and vertices of the plateau into arrays 
      while (! LifoVide(LIFO))
	{ 
	  x = LifoPop(LIFO);
	  e1 = G->Edges[0][x]; e2 = G->Edges[1][x];
	  re1 = element_find(e1, Fth );
	  re2 = element_find(e2, Fth );
	  if ( G->Solution[0][re1]<0 ||  G->Solution[0][re2]<0) 
	    {
	      if (indic_VP[e1]==0) 
		{
		  LCVP[nb_vertices]=e1;
		  nb_vertices++;
		  indic_VP[e1]=1;
		}
	      if (indic_VP[e2]==0) 
		{
		  LCVP[nb_vertices]=e2;
		  nb_vertices++;
		  indic_VP[e2]=1;
		}
	      edgesLCP[0][ nb_edges] = e1;
	      edgesLCP[1][ nb_edges] = e2;
	      NEs[nb_edges]=x;
	      nb_edges ++;
	    }
	  
	  for (j = 0; j <= 1; j += 1) 
	    {
	      e_jx = G->Edges[j][x];
	      if(e_jx< N)//si le sommet n'est pas marqué ATTENTION A VERIFIER 
		{
		  /* for (k = 0; k < G->Neighbors[e_jx]->Sp; k += 1) 
		    {
		      y = G->Neighbors[e_jx]->Pts[k]; 
		  */
		  for (lis=G->Neighbors[e_jx]; lis != NULL; lis = lis->suite)
		    { 
		    y=lis->elt;
		      if ((indic_P[y]==false) && (G->RecWeights[y]== wmax))
			{
			  indic_P[y]=true;
			  LifoPush(LIFO, y);
			  LifoPush(LCP, y);
			  indic_E[y]= true;
			} 
		    }
		}
	    }
	}
      //t2=clock();  printf("Computation time : %.6lf seconds elapsed\n", ((double)t2-t1)/CLOCKS_PER_SEC);
      for (j=0;j<nb_vertices;j++)
	indic_VP[LCVP[j]]=0;
      for (j=0;j<LCP->Sp;j++) 
	indic_P[LCP->Pts[j]]=false;

      // 3. If e_max belongs to a plateau
      if (nb_edges > 0)
	{
	  // 4. Evaluate if there are differents seeds on the plateau
	  p=0;  
	  different_seeds = false;
	  
	  for (i=0;i<G->P;i++)
	    { 
	      val = -0.5;
	      for (j=0;j<nb_vertices;j++)
		{
		  x = LCVP[j];
		  xr = element_find(x, Fth);
		  if(fabs( G->Solution[i][xr]-val)>epsilon &&  G->Solution[i][xr]>=0 ) 
		    {
		      p++; val =  G->Solution[i][xr]; 
		    }
		}
	      if (p>=2) 
		{
		  different_seeds = true;
		  break;
		}
	      else p=0;
	    }
	  
	  if (different_seeds == true)
	    {
	      
	      // 5. Sort the edges of the plateau according to their normal weight
	      sorted_weights = (wtype *)malloc(nb_edges*sizeof(wtype));
	      if (sorted_weights == NULL) { fprintf(stderr, "%s: malloc sorted_weights failed\n", F_NAME); exit(0);}

	       for(k=0;k<nb_edges;k++) 
		 {
		 sorted_weights[k]=G->Weights[NEs[k]]; 
		 //printf("sorted weights = %f \n", sorted_weights[k]) ; 
		 }
	       // NEs = linsortimagedown(sorted_weights, nb_edges, max_weight+1);
	       
	       TriRapideStochastique_dec(sorted_weights,NEs, 0, nb_edges-1); 
	       free(sorted_weights);

	      // Merge nodes for edges of real max weight
		nb_vertices=0;
		Nnb_edges = 0;
		for(Ncpt_aretes = 0; Ncpt_aretes< nb_edges; Ncpt_aretes++)
		  {
		    Ne_max=NEs[Ncpt_aretes];
		    e1 = G->Edges[0][ Ne_max];
		    e2 = G->Edges[1][ Ne_max];
		    if (G->Weights[Ne_max] != wmax)
		      merge_node (e1, e2,  Rnk, Fth,  G->Solution, G->P);
		    else 
		      {
			re1 = element_find(e1, Fth );
			re2 = element_find(e2, Fth );
			if ((re1 !=re2)&&(( G->Solution[0][re1]<0 || G->Solution[0][re2]<0)))
			  {
			    if (indic_VP[re1]==0) 
			      {
				LCVP[nb_vertices]=re1;
				nb_vertices++;
				indic_VP[re1]=1;
			      }
			    if (indic_VP[re2]==0) 
			      {
				LCVP[nb_vertices]=re2;
				nb_vertices++;
				indic_VP[re2]=1;
			      }
			    edgesLCP[0][ Nnb_edges] = re1;
			    edgesLCP[1][ Nnb_edges] = re2;
			    Nnb_edges ++;
			  }
		      }
		  }
	
	      for (i=0;i<G->P;i++)
		{ 
		  k=0;
		  for (j=0;j<nb_vertices;j++)
		    {
		      xr = LCVP[j];
		      if ( G->Solution[i][xr]>=0)
			{
			  local_labels[i][k] =  G->Solution[i][xr];
			  // printf("%f \n",local_labels[i][k] );
			  local_seeds[k] = xr;
			  k++;
			}
		    }
		}
		    
	      // 6. Execute Random Walker on plateaus
	      
	      if(nb_vertices<SIZE_MAX_PLATEAU)
		{
		  // printf("Plateau of (%d vertices,%d edges %d seeds) RW \n", nb_vertices, Nnb_edges, k);
		  // int t1=clock();
		success = RandomWalker(edgesLCP, Nnb_edges, LCVP, indic_VP, nb_vertices, local_seeds, local_labels, k, G->P+1,  G->Solution);
	        //int t2=clock();
	       	//printf("Random walker : %.6lf seconds elapsed\n", ((double)t2-t1)/CLOCKS_PER_SEC);
	      
		}
	      
	      //printf("Plateau of (%d vertices,%d edges) RW \n", nb_vertices, Nnb_edges);
	      if ((nb_vertices>=SIZE_MAX_PLATEAU)||(success==false))
		{
		  printf("Plateau of a big size (%d vertices,%d edges) the RW is not performed on it\n", nb_vertices, Nnb_edges);
		  for (j=0;j<Nnb_edges;j++)
		    {
		      e1 = edgesLCP[0][j];
		      e2 = edgesLCP[1][j];
		      merge_node (e1, e2,  Rnk, Fth,  G->Solution, G->P);
		    }
		}
	      for (j=0;j<nb_vertices;j++)
		indic_VP[LCVP[j]]=0;
	    }
	  else // if different seeds = false 
	    // 7. Merge nodes for edges of max weight
	    {
	     
	      for (j=0;j<nb_edges;j++)
		{
		  e1 = edgesLCP[0][j];
		  e2 = edgesLCP[1][j];
		  merge_node (e1, e2,  Rnk, Fth,  G->Solution, G->P);
		}
	      
	    }
	}
      LifoFlush(LCP);
    } // end main loop
  //fprintf(stderr, "juste avant\n");

  //building the final proba G->Solution map (find the root vertex of each tree)
  for(i=0; i<G->N; i++) 
    {
      j=i;
      xr = i;
      while(Fth[i] != (uint32_t)i)
	{ 
	  i = xr;
	  xr = Fth[i];
	}
      for(k=0; k< G->P;k++)  G->Solution[k][j] = G->Solution[k][i];
      i=j;
    }

  // free memory 
  LifoTermine(LCP);
  LifoTermine(LIFO);

  for (i=0;i<2;i++) 
    free(edgesLCP[i]); free(edgesLCP);
  
  free(Rnk);
  free(local_seeds);
  for (i=0; i<G->P; i++)  free(local_labels[i]);
  free(local_labels);
  
  free(LCVP);
  free(Es);
  free(NEs);
  free(indic_E);
  free(indic_VP);
  free(indic_P);
  free(Fth);
 
 
}

template <class wtype>
bool equal (wtype a, wtype b)
{
  if (a<b) return b-a<eps;
  else return a-b<eps;
}







#endif
