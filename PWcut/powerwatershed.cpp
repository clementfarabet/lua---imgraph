/*
  Power watershed algorithm for Maximum Spanning Forest (MSF) computation
  implemented to compute an MSF cut in a tree (hierarchy)
  author: Camille Couprie
  31 oct. 2011
*/

#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
/*#include <stdlib.h>

#include <mccodimage.h>
#include <mcimage.h>
#include <mclifo.h>
#include <mcindic.h>
#include <mcutil.h>
#include <jcgraphes.h>
#include <jccomptree.h>*/


#include <MTree_utils.h>
#include <graph_utils.h>
#include <sys/types.h>
#include <MSF_RW.h>

#ifndef _MERGETREESTRUCT_
#define _MERGETREESTRUCT_
typedef struct {
  mtree *tree;
  RAG *rag;
  struct xvimage *labels;
  int32_t *altitudes;
  float *weights;
  int cs;
  int rs;
} MergeTree;
#endif


#ifndef _LISTSTRUCT_
#define _LISTSTRUCT_
typedef struct list
{
  int index;
  struct list *next;
} list ;
#endif


/*================================================*/
void Insert(list **sl, int index)
/*================================================*/
{
  list *tmp = NULL;
  list *csl = *sl;
  list *elem = (list*) malloc(sizeof(list));
  if(!elem) exit(EXIT_FAILURE);
  elem->index = index;
  while(csl)
    {
      tmp = csl;
      csl = csl->next;
    }
  elem->next = csl;
  if(tmp) tmp->next = elem;
  else *sl = elem;
}

//#define false 0
//#define true 1


/*=====================================================================================*/
list * Powerwatershed(MergeTree * MT)
/*=====================================================================================*/
/*Segment a tree into two components.
  Returns a list of nodes correspunding to the Max Spanning Forest cut,
  computed using Power watershed's algorithm */
{
  int y, i, j, tmp;
  int nb_markers; int nb_leafs;
  int N, M; 
  float val=1; //weight parameter for leafs.

 // -------- Gathering usefull input graph (MT) informations -----------

  mtree * T= MT->tree;
  float * W = MT->weights;
  JCctree *CT = T->CT;
  int root_node = CT->root;
  JCsoncell *s;

  //nb nodes
   M = CT->nbnodes;
  
   // nb_edges
   nb_leafs = 0;
   for (i = 0; i < M; i++)
     if (CT->tabnodes[i].nbsons == 0)
       nb_leafs++;
   
   nb_markers = nb_leafs+1;
   N=M+nb_markers;
   M=N-1;
   printf("Nb nodes:%d Nb edges: %d Nb leafs :%d \n", N, M, nb_leafs);


  struct graph <double> *G;
  G = (struct graph <double>*)malloc(sizeof(struct graph<double>));
  G->weight_type=1; // 1:double 0:uint32_t

  Allocate_Graph(G, 
		 N,/*Number of nodes */ 
		 nb_markers, /*Max Number of seeded nodes */ 
		 M); /*Number of edges */ 
		
  G->P = 1; /*we are solving one problem (multiple problems solved for multilabel segmentation)*/ 
 
  /*Fill the seeds values */
  G->S=nb_markers;/*there are nb_markers seeded nodes*/
  G->SeededNodes[0]= M;
  G->Labels[0]=0;
  j=1;
  for (i = 0; i < CT->nbnodes; i++)
    if (CT->tabnodes[i].nbsons == 0)
      {
	G->SeededNodes[j]= i+CT->nbnodes;
	G->Labels[j]=1;
	j++;
      }
  // weights
   
   float * weights = (float *)malloc(N*sizeof(float));
   for(i=0;i<CT->nbnodes;i++)
     weights[i]=W[i];
   
   for(i=0;i<nb_markers;i++)
     weights[CT->nbnodes+i]=val;
     


  /*Add the weighted edges*/
  for  (i=0;i<CT->nbnodes;i++)
    {
      tmp = CT->tabnodes[i].nbsons;
      if (tmp==0) //leaf
	{
	  y= G->SeededNodes[i+1]; // edge index
	  // fprintf(stderr,"edge (%d %d) %d \n", i,y,  (int)(weights[y]*CTE_WEIGHTS ));
	  AddEdge<double>(G, i, G->SeededNodes[i], weights[y],y/*edge index*/);
	 
	} 
      else
	{
	  for ( s = CT->tabnodes[i].sonlist;s!=NULL;s = s->next)  
	    {
	      y=s->son;
	      //  fprintf(stderr,"edge (%d %d) %d \n", i,y,  (int)(weights[y]*CTE_WEIGHTS ));
	      AddEdge<double>(G, i, y, weights[y],y/*edge index*/);
	    }
	}
    }
  AddEdge<double>(G, root_node, M, weights[y],y/*edge index*/);

  G->max_weight= val; /*maximum weight value*/

  /*Solving problem with Power Watersheds*/
     
  PowerWatershed_q2<double>(G);
  
  // Writing results 
  printf("SOLUTION \n");
  for (j = 0; j < G->N; j++)
    printf("%f \n", G->Solution[0][j]);

 // ------------------ Process the tree to find the cut ----------------------
  list * cut = NULL;
 for (i = 0; i < CT->nbnodes; i++)
    {
      // nodes having a different value than their father are in the cut
      if ((CT->tabnodes[i].father != -1) && (G->Solution[0][CT->tabnodes[i].father] != G->Solution[0][i]))
        Insert(&cut, i);
      // leafs having the same label as the root are in the cut
      if ((CT->tabnodes[i].nbsons == 0) && (G->Solution[0][i]==0))
        Insert(&cut, i);
    }

 
 if (cut == NULL)  Insert(&cut, root_node);
 //PrintList(cut);

 
  Free_Graph(G, G->N);

  return cut;



}

