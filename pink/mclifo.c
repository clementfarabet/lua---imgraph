/*
Copyright ESIEE (2009) 

m.couprie@esiee.fr

This software is an image processing library whose purpose is to be
used primarily for research and teaching.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software. You can  use, 
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
/* $Id: mclifo.c,v 1.6 2006/02/28 07:49:16 michel Exp $ */
/* 
   Librairie mclifo :

   fonctions pour la gestion d'une liste lifo

   Michel Couprie 1996
*/

/* #define TESTLifo */
#include <stdio.h>
#include <stdint.h>
#include <sys/types.h>
#include <stdlib.h>
#include <mclifo.h>

/* ==================================== */
Lifo * CreeLifoVide(
  int32_t taillemax)
/* ==================================== */
{
  Lifo * L = (Lifo *)calloc(1,sizeof(Lifo) + sizeof(int32_t) * (taillemax-1));
  if (L == NULL)
  {   
    fprintf(stderr, "CreeLifoVide() : malloc failed : %d bytes\n", 
            sizeof(Lifo) + sizeof(int32_t) * (taillemax-1));
    return NULL;
  }
  L->Max = taillemax;
  L->Sp = 0;
  return L;
}

/* ==================================== */
void LifoFlush(
  Lifo * L)
/* ==================================== */
{
  L->Sp = 0;
}

/* ==================================== */
int32_t LifoVide(
  Lifo * L)
/* ==================================== */
{
  return (L->Sp == 0);
}

/* ==================================== */
int32_t LifoPop(
  Lifo * L)
/* ==================================== */
{
  if (L->Sp == 0)
  {
    fprintf(stderr, "erreur Lifo vide\n");
    exit(1);
  }
  L->Sp -= 1;
  return L->Pts[L->Sp];
}

/* ==================================== */
int32_t LifoHead(
  Lifo * L)
/* ==================================== */
{
  if (L->Sp == 0)
  {
    fprintf(stderr, "erreur Lifo vide\n");
    exit(1);
  }
  return L->Pts[L->Sp-1];
}
  
/* ==================================== */
void LifoPush(Lifo * L, int32_t V)
/* ==================================== */
{
  if (L->Sp > L->Max - 1)
  {
    fprintf(stderr, "erreur Lifo pleine\n");
    exit(1);
  }
  L->Pts[L->Sp] = V;
  L->Sp += 1;
}

/* ==================================== */
void LifoPrint(Lifo * L)
/* ==================================== */
{
  int32_t i;
  if (LifoVide(L)) {printf("[]"); return;}
  printf("[ ");
  for (i = 0; i < L->Sp; i++)
    printf("%d ", L->Pts[i]);
  printf("]");
}

/* ==================================== */
void LifoPrintLine(Lifo * L)
/* ==================================== */
{
  int32_t i;
  if (LifoVide(L)) {printf("[]\n"); return;}
/*
  printf("Max = %d ; Sp = %d \n", L->Max, L->Sp);
*/
  printf("[ ");
  for (i = 0; i < L->Sp; i++)
    printf("%d ", L->Pts[i]);
  printf("]\n");
}

/* ==================================== */
void LifoTermine(
  Lifo * L)
/* ==================================== */
{
  free(L);
}

#ifdef TESTLifo
void main()
{
  Lifo * L = CreeLifoVide(3);
  LifoPrint(L);
  if (LifoVide(L)) printf("LifoVide OUI\n");
  LifoPush(L,1);
  LifoPrint(L);
  if (!LifoVide(L)) printf("LifoVide NON\n");
  LifoPush(L,2);
  LifoPrint(L);
  LifoPush(L,3);
  LifoPrint(L);
  printf("LifoPop %d attendu 3\n", LifoPop(L));
  LifoPrint(L);
  LifoPush(L,4);
  LifoPrint(L);
  printf("LifoPop %d attendu 4\n", LifoPop(L));
  LifoPrint(L);
  printf("LifoPop %d attendu 2\n", LifoPop(L));
  LifoPrint(L);
  printf("LifoPop %d attendu 1\n", LifoPop(L));
  LifoPrint(L);
  if (LifoVide(L)) printf("LifoVide OUI\n");
  printf("maintenant sortie attendue sur lifo pleine :\n");
  LifoPush(L,3);
  LifoPush(L,3);
  LifoPush(L,3);
  LifoPush(L,3);  
}
#endif

