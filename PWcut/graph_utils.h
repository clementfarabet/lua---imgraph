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

#ifndef GRAPH_UTILS_H
#define GRAPH_UTILS_H

#define epsilon 0.00001
#define MAX_DEGREE 3000
#define SIZE_MAX_PLATEAU 1800000
typedef float DBL_TYPE;

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <mclistechainee.h>

/*!
  \file graph_utils.h
  \brief Generic graph structures

 Solution = argmin_x { sum_{e_{ij}=1...M} {w_{ij}^p|x_i-x_j|^q + sum_{v_i=1...N} {w_{ij}^p|x_i-y_i|^q} 
y_i 
 <B>
 Description: A graph structure is composed of 
 </B>
*/

template<class wtype> struct graph {
  int       N;           /*! nb_vertices */
  int       M;           /*! nb_edges */
  int       S;           /*! nb of seeded nodes */
  int       P;           /*! Number of different problems to solve on the graph, 
			    i.e. the number of labels in a multilabel segmentation */
  
  wtype max_weight;         /* maximum weight value */
  int weight_type;          /* 0: integer weights >0 double*/
  uint32_t** Edges;         /*! array of node indexes composing edges */
  TypListechainee **     Neighbors;     /*! array of lists of edge neighbors of each node */
 
  uint32_t*  SeededNodes;   /*! list of seeded nodes*/
  uint8_t*   Labels;        /*! list of labels corresponding to seeded nodes */
  wtype * Weights;          /*! edges weights */
  wtype * RecWeights;    /*! integer reconstructed weights relatively to the seeds */
 
  DBL_TYPE** Solution;      /*! proba value assigned to each node */
};	



/* =============================================================== */
template<class wtype>void Allocate_Minimum_Graph( struct graph <wtype> *G, /*Graph structure*/
			     int nb_nodes,  /*Number of nodes*/ 
			     int nb_max_seeded_nodes, /*Max Number of seeded nodes*/
			     int nb_edges) /*Number of edges*/ 
/* =============================================================== */
{
  G->N = nb_nodes;
  G->M = nb_edges;
  G->S = nb_max_seeded_nodes;
  G->SeededNodes = (uint32_t*)malloc(nb_max_seeded_nodes *sizeof(uint32_t));
  G->Labels = (uint8_t*)malloc(nb_max_seeded_nodes *sizeof(uint8_t));
}


/* =============================================================== */
template<class wtype> void Free_Minimum_Graph( struct graph  <wtype> *G, /*Graph structure*/
			 int nb_max_seeded_nodes) /*Max Number of seeded nodes*/
/* =============================================================== */
{
  free(G->SeededNodes);
  free(G->Labels);
}


/* =============================================================== */
template<class wtype> void Allocate_Graph( struct graph  <wtype> *G, /*Graph structure*/
		     int nb_nodes,  /*Number of nodes*/ 
		     int nb_max_seeded_nodes, /*Max Number of seeded nodes*/
		     int nb_edges) /*Number of edges*/ 
		    
/* =============================================================== */
{
  int i;
  Allocate_Minimum_Graph(G, nb_nodes,  nb_max_seeded_nodes, nb_edges);
  G->Neighbors =(TypListechainee **)malloc(G->N*sizeof(TypListechainee *));
  G->Edges =  (uint32_t**)malloc(2*sizeof(uint32_t*));
  for(i=0;i<2;i++) G->Edges[i] = (uint32_t*)malloc(G->M*sizeof(uint32_t));
  G->Weights = (wtype*)malloc(sizeof(wtype)*G->M);
 
  /*Initialise useful values for the power watershed algorithm*/
  G->RecWeights = (wtype*)calloc(sizeof(wtype),G->M);

  /*initializing the incidence array G->Neighbors*/
  for (i=0;i<G->N;i++)
    G->Neighbors[i] = ListechaineeVide();
}


/* =============================================================== */
template<class wtype> void Free_Graph( struct graph  <wtype> *G, /*Graph structure*/
				       int nb_max_seeded_nodes) /*Max Number of seeded nodes*/
		
/* =============================================================== */
{
  int i;
  Free_Minimum_Graph(G, nb_max_seeded_nodes);

  for(i=0;i<G->N;i++) 
  DetruitListechainee(G->Neighbors[i]); 
  free(G->Neighbors);

  free(G->RecWeights);

  for (i=0;i<2;i++) 
    free(G->Edges[i]);
  free(G->Edges);
  free(G->Weights);

 for (i=0; i<G->P; i++)  
    free(G->Solution[i]);
  free(G->Solution);
}


/* =============================================================== */
template<class wtype> void Free_Partial_Graph( struct graph  <wtype> *G, /*Graph structure*/
			 int nb_max_seeded_nodes,
			 int doubleweights) /*0 : weights are integers */
                                    /*1 : weights are DBL_TYPE (may be float) */
                                    /*2 : weights are double */
/* =============================================================== */
{
  int i;
 
  G->SeededNodes = (uint32_t*) realloc(G->SeededNodes ,sizeof(uint32_t)*G->S);
  G->Labels = (uint8_t*) realloc(G->Labels ,sizeof(uint8_t)*G->S);


   for(i=G->N+1; i<nb_max_seeded_nodes;i++) 
    G->Neighbors[i] = (TypListechainee *)realloc( G->Neighbors[i] ,0);
 
    G->Neighbors =(TypListechainee **)realloc(G->Neighbors, G->N*sizeof(TypListechainee *));
  
  for (i=0;i<2;i++) 
    G->Edges[i] = (uint32_t*)realloc(G->Edges[i] ,sizeof(uint32_t)* G->M);
  
 
  if (doubleweights==2) 
    G->dWeights = (double*) realloc(G->Weights ,sizeof(double)*G->M);
  else if (doubleweights==1) 
    G->DWeights = (DBL_TYPE*) realloc(G->Weights ,sizeof(DBL_TYPE)*G->M);
  else
   G->Weights = (uint32_t*) realloc(G->Weights ,sizeof(uint32_t)*G->M);
}





/* =============================================================== */
template<class wtype> void Allocate_Graph_Solution( struct graph  <wtype> *G) /*Graph structure*/
			    
/* =============================================================== */
{
  int i, j;
  G->Solution = (DBL_TYPE **)malloc((G->P) *sizeof(DBL_TYPE*));
  for (i=0; i<G->P; i++) 
    {
      G->Solution[i]= (DBL_TYPE *)malloc(G->N *sizeof(DBL_TYPE));
      for (j=0; j<G->N; j++) G->Solution[i][j]=-1;
    }
}


/* =============================================================== */
template<class wtype> void Print_Graph_Content( struct graph  <wtype> *G) /*Graph structure*/
			  
/* =============================================================== */
{
  int i;
  fprintf(stderr,"************************************\n");
  fprintf(stderr,"Nb of nodes %d \n", G->N);
  fprintf(stderr,"Nb of edges %d \n", G->M);
  fprintf(stderr,"Nb of seeded nodes %d \n", G->S);

  fprintf(stderr,"List of edges : \n");
  for(i=0;i<G->M;i++)
    fprintf(stderr, "edge %d : (node %d, node %d) \n", i, G->Edges[0][i], G->Edges[1][i] );
  fprintf(stderr,"************************************\n");

    fprintf(stderr,"List of edge neighbors for each node : \n");
    TypListechainee * lis;
  for(i=0;i<G->N;i++)
    {
      fprintf(stderr, "node %d has for neighbors edges ", i);
      for (lis=G->Neighbors[i]; lis != NULL; lis = lis->suite) 
	fprintf(stderr, " %d, ", lis->elt);
      fprintf(stderr, "\n ");
    }
  fprintf(stderr,"************************************\n");

  fprintf(stderr,"Edge weights : \n");
  if (G->weight_type>0) 
    for(i=0;i<G->M;i++)
      fprintf(stderr, "weight of float edge %d = %f  \n", i, G->Weights[i]);
  else
    for(i=0;i<G->M;i++)
      fprintf(stderr, "weight of integer edge %d = %d  \n",i,  G->Weights[i]);
  fprintf(stderr,"************************************\n");

 fprintf(stderr,"List of seeded nodes : \n");
  for(i=0;i<G->S;i++)
    fprintf(stderr, "label ( node %d) = %d \n", G->SeededNodes[i], G->Labels[i]);
  fprintf(stderr,"************************************\n");

  
}



/* ======================================================================================================== */
template<class wtype> void AddEdge( struct graph  <wtype> * G, int node1, int node2, wtype weight, int index_edge)
/* ======================================================================================================== */
{
  G->Edges[0][index_edge]= node1;
  G->Edges[1][index_edge]= node2;
  G->Neighbors[node1] = Cons(index_edge, G->Neighbors[node1]);
  G->Neighbors[node2] = Cons(index_edge, G->Neighbors[node2]);
  G->Weights[index_edge] = weight;
}






#endif
