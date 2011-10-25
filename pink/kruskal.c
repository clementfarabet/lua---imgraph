/*
  Kruskal algorithm for maximum spanning forest computation
  author: Camille Couprie
  21 oct. 2011
*/

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <sys/types.h>
#include <mccodimage.h>
#include <mcimage.h>
#include <mclifo.h>
#include <mcindic.h>
#include <mcutil.h>
#include <jcgraphes.h>
#include <jccomptree.h>
#include <mtrand64.h>
#include "kruskal.h"

#define false 0
#define true 1

void Insert(list **sl, int index)
{
  list *tmp = NULL;
  list *csl = *sl;
  list *elem = malloc(sizeof(list));
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

/*=====================================================================================*/
list * MSF_Kruskal(MergeTree * MT)  // max_weight    /* maximum weight value */ )
/*=====================================================================================*/
/*Segment a tree into two components.
  Returns a list of nodes correspunding to the Max Spanning Forest cut,
  computed using Kruskal's algorithm */
{
  int i, j, k, x, y, e1, e2;
  int nb_markers; int nb_leafs;
  long N, M;
  float val=1; //weight parameter for leafs.
  mtree * T= MT->tree;
  float * W = MT->weights;
  // mergeTreePrint(T);
  JCctree *CT = T->CT;
  int root_node = CT->root;

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
  //printf("Nb nodes:%d Nb edges: %d Nb leafs :%d \n", N, M, nb_leafs);

  // indexes of edges : son's nodes indexes
  //Memory allocation of temporary arrays for Krukal's algorithm
  Lifo * LIFO;
  LIFO = CreeLifoVide(M);
  if (LIFO == NULL) { fprintf(stderr, "kruskal : CreeLifoVide failed\n"); exit(0); }

  int * Mrk = (int*)calloc(N ,sizeof(int));
  if (Mrk == NULL) { fprintf(stderr, "kruskal : malloc failed\n"); exit(0); }

  uint32_t * SeededNodes = (uint32_t*)malloc(nb_markers*sizeof(uint32_t));
  if (SeededNodes == NULL) { fprintf(stderr, "kruskal : malloc failed\n"); exit(0); }

  // markers
  SeededNodes[0]= M;
  j=1;
  for (i = 0; i < CT->nbnodes; i++)
    if (CT->tabnodes[i].nbsons == 0)
      {
        SeededNodes[j]= i+CT->nbnodes;
        Mrk[SeededNodes[j]] = 1;
	j++;
      }
  Mrk[M] = 1;
 
  uint32_t * Rnk = (uint32_t*)calloc(N, sizeof(uint32_t));
  if (Rnk == NULL) { fprintf(stderr, "kruskal : malloc failed\n"); exit(0); }
  uint32_t * Fth = (uint32_t*)malloc(N*sizeof(uint32_t));
  if (Fth == NULL) { fprintf(stderr, "kruskal : malloc failed\n"); exit(0); }
  for(k=0;k<N;k++) { Fth[k]=k; }

  // Es : E sorted by decreasing weights
  uint32_t * Es = (uint32_t*)malloc(M*sizeof(uint32_t));
  if (Es == NULL) { fprintf(stderr, "kruskal : malloc failed\n"); exit(0); }
  for(k=0;k<M;k++) Es[k]=k;

  float * sorted_weights = (float *)malloc(M*sizeof(float));
  for(k=0;k<M;k++)
    sorted_weights[k]=W[k];

  for(k=0;k<nb_leafs;k++)
    sorted_weights[CT->nbnodes+k]=val;
    
  TriRapideStochastique_dec(sorted_weights,Es, 0, M-1);
  free(sorted_weights);

  long nb_arete = 0;
  int e_max, root;
  long cpt_aretes = 0;

  // beginning of main loop
  while (nb_arete < N-nb_markers)
    {
      e_max=Es[cpt_aretes];
      cpt_aretes=cpt_aretes+1;
      e1= e_max;  // e1 = Edges[0][e_max];
      if (e_max<CT->nbnodes) e2= CT->tabnodes[e_max].father;
      else e2= e_max-CT->nbnodes;
      if (e2==-1)e2=M;  //e2 = Edges[1][e_max];
      //printf("(%d %d)\n", e1,e2);
      x = element_find(e1, Fth );
      y = element_find(e2, Fth );
      if ((x != y) && (!(Mrk[x]>=1 && Mrk[y]>=1)))
        {
          root = element_link( x,y, Rnk, Fth);
	  //printf("link\n");
          nb_arete=nb_arete+1;
          if ( Mrk[x]>=1) Mrk[root]= Mrk[x];
          else if ( Mrk[y]>=1) Mrk[root]= Mrk[y];
        }
    }

  //building the labeling for each individual markers in map
  // (find the root vertex of each tree)
  int * Map2 = (int *)malloc(N*sizeof(int));
  int * Map = (int *)malloc(N*sizeof(int));
   for (i=0; i<N; i++)
    Map2[i] = element_find(i, Fth);
    
  // Compute the binary labeling in Map
  for (i = 1; i < nb_markers; i++)
    Map[SeededNodes[i]] = 1;
  Map[M]=0;

  for (i=0;i<N;i++) Mrk[i] = false;
  for (i=0;i<nb_markers; i++)
    {
      LifoPush(LIFO, SeededNodes[i]);
      while (!LifoVide(LIFO))
        {
          x = LifoPop(LIFO);
          Mrk[x]=true;
          j= nb_neighbors(x, CT, nb_leafs);
	  for (k=0;k<j;k++)
            {
	    
              y = neighbor(x, k, CT, nb_leafs, SeededNodes);
	   
              if (y==-1)y=M;
	      if (Map2[y]==Map2[SeededNodes[i]]  && Mrk[y]==false)
                {
                  LifoPush(LIFO, y);
                  if (i==0) Map[y]= 0;
                  else  Map[y]= 1;
                  Mrk[y]=true;
                }
            }
	}
      LifoFlush(LIFO);
    }
  for (i = 1; i < nb_markers; i++)
    Map[SeededNodes[i]] = 1;
  Map[M]=0;

  /*  for (i=0; i<N; i++) {
    printf("Map[%d]=%d \n",i,Map[i]);
    }*/

  // Process the tree to find the cut
  list * cut = NULL;
  for (i = 0; i < CT->nbnodes; i++)
    {
      // nodes having a different value than their father are in the cut
      if ((CT->tabnodes[i].father != -1) && (Map[CT->tabnodes[i].father] != Map[i]))
        Insert(&cut, i);
      // leafs having the same label as the root are in the cut
      if ((CT->tabnodes[i].nbsons == 0) && (Map[i]==0))
        Insert(&cut, i);
    }

  //PrintList(cut);

  LifoTermine(LIFO);
  free(Mrk);
  free(SeededNodes);
  free(Rnk);
  free(Fth);
  free(Es);
  free(Map);
  free(Map2);
  return cut;
}
/*================================================*/
void PrintList(list *sl)
/*================================================*/
{
  fprintf(stderr, "Nodes of the cut:\n");
  while(sl)
    {
      printf("%d\n",sl->index);
      sl = sl->next;
    }
}

/*================================================*/
int nb_neighbors(int x, JCctree *CT, int nb_leafs)
/*================================================*/
{
  int tmp;
  if (x<CT->nbnodes)
    {
      tmp = CT->tabnodes[x].nbsons;
      if (tmp==0) tmp=1;
      return tmp+1;
    }
  else if (x<CT->nbnodes+nb_leafs)
    return 1; //CT->tabnodes[CT->root].nbsons;
  else return 1;
}


/*================================================*/
int neighbor(int x, int k, JCctree *CT, int nb_leafs, int * SeededNodes)
/*================================================*/
{

  JCsoncell *s;int i, tmp;
  if (x<CT->nbnodes)
    {
      tmp = CT->tabnodes[x].nbsons;
      if (tmp==0)
        { if (k==0) return CT->tabnodes[x].father;
          return SeededNodes[x+1];
        }
      else if (k<=tmp)
        {
          if (k==tmp) return CT->tabnodes[x].father;
          s = CT->tabnodes[x].sonlist;
          for (i=0;i<k;i++) s = s->next; // INEFFICACE A REFAIRE
	  //fprintf(stderr," ici ");
          return s->son;
        }
    }
  else if (x<CT->nbnodes+nb_leafs)
    return x-CT->nbnodes;
  else
    {
      return CT->root;
    }
}

/*================================================*/
int element_link( int x,int y, uint32_t *Rnk, uint32_t *Fth)
/*================================================*/
{
  if( Rnk[x] > Rnk[y])
    {
      int t;
      t=x;
      x=y;
      y=t;
    }
  if( Rnk[x] == Rnk[y])
    {
      Rnk[y]=Rnk[y]+1;
    }
  Fth[x] = y;
  return y;
}

/*===============================*/
int element_find(int x, uint32_t *Fth )
/*===============================*/
{
  if (Fth[x] != x)
    Fth[x] = element_find(Fth[x], Fth);
  return Fth[x];
}

/* =============================================================== */
long Partitionner_dec(float *A, uint32_t * I, long p, long r)
/* =============================================================== */
/*
  partitionne les elements de A entre l'indice p (compris) et l'indice r (compris)
  en deux groupes : ceux <= A[p] et les autres.
*/
{
  float t;
  int t1;
  float x = A[p];
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
long PartitionStochastique_dec (float *A, uint32_t * I, long p, long r)
/* =============================================================== */
/*
  partitionne les elements de A entre l'indice p (compris) et l'indice r (compris)
  en deux groupes : ceux <= A[q] et les autres, avec q tire au hasard dans [p,r].
*/
{
  float t;
  int t1;
  long q;

  q = p + (genrand64_int64() % (r - p + 1)); /* rand must be 64-bit safe, should be OK now */
  t = A[p];         /* echange A[p] et A[q] */
  A[p] = A[q];
  A[q] = t;

  t1 = I[p];         /* echange I[p] et I[q] */
  I[p] = I[q];
  I[q] = t1;

  return Partitionner_dec(A, I, p, r);
} /* PartitionStochastique() */

/* =============================================================== */
void TriRapideStochastique_dec (float * A, uint32_t *I, long p, long r)
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
