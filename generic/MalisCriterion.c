#include <iostream>
using namespace std;

#ifndef TH_GENERIC_FILE
#define TH_GENERIC_FILE "generic/MalisCriterion.c"
#else

/**********************************************************************
 * Compute the MALIS loss function and its derivative wrt an 
 * affinity graph (a maximum spanning tree is computed)
 *
 * Original Author: Srini Turaga (sturaga@mit.edu)
 * All rights reserved
 *
 * Clement Farabet: reformatted Srini's code to fit Torch's
 * criterion's system, readapted the code to 2D graphs, and got rid
 * of a couple of hard-coded elements
 * 
 * Srini Turaga (Aug 10, 2012): fixed some bugs, changed the elementary
 * loss from the sq-sq to the hinge loss, added ability to compute only
 * positive examples, only negative examples, or both. Re-re-adapted to
 * bring it back to 3d graphs (sorry Clement) and got rid of hard-coded
 * elements like the threshold and margin (yay!)
 **********************************************************************/


class nn_(AffinityGraphCompare){
private:
  const real * mEdgeWeightArray;
public:
  nn_(AffinityGraphCompare)(const real * EdgeWeightArray){
    mEdgeWeightArray = EdgeWeightArray;
  }
  bool operator() (const long& ind1, const long& ind2) const {
    return (mEdgeWeightArray[ind1] > mEdgeWeightArray[ind2]);
  }
};

static int nn_(MalisCriterion_forward)(lua_State *L)
{
  // 4d connectivity graph [#edges * depth * height * width]
  THTensor *conn = (THTensor *)luaT_checkudata(L, 2, torch_Tensor);  
  conn = THTensor_(newContiguous)(conn);
  long conn_ndims = THTensor_(nDimension)(conn);
  long conn_nelts = THTensor_(nElement)(conn);
  long conn_nedges = conn->size[0];
  real *conn_data = THTensor_(data)(conn);

  // target segmentation [depth * height * width]
  THTensor *target = (THTensor *)luaT_checkudata(L, 3, torch_Tensor);
  target = THTensor_(newContiguous)(target);
  long target_ndims = THTensor_(nDimension)(target);
  long target_nelts = THTensor_(nElement)(target);
  real *target_data = THTensor_(data)(target);

  // graph neighborhood descriptor [#edges * 3]
  THTensor *nhood = (THTensor *)luaT_getfieldcheckudata(L, 1, "nhood", torch_Tensor);
  long nhood_ndims = THTensor_(nDimension)(nhood);
  long nhood_nelts = THTensor_(nElement)(nhood);
  long nhood_nedges = nhood->size[0];
  long nhood_geom = nhood->size[1];
  real *nhood_data = THTensor_(data)(nhood);

  // hinge loss threshold and margin
  const real threshold = luaT_getfieldchecknumber(L, 1, "threshold");
  const real margin = luaT_getfieldchecknumber(L, 1, "margin");

  // is this a positive example (true) or negative (false)
  int pos = luaT_getfieldcheckboolean(L, 1, "posexample");
  int neg = luaT_getfieldcheckboolean(L, 1, "negexample");

  // output: gradient wrt connectivity graph
  THTensor *dloss = (THTensor *)luaT_getfieldcheckudata(L, 1, "gradInput", torch_Tensor);
  THTensor_(resizeAs)(dloss, conn);
  real *dloss_data = THTensor_(data)(dloss);
  THTensor_(zero)(dloss);

  // checks
  if (nhood_ndims != 2)
    THError("nhood should be 2d");
  if (nhood_geom != (conn_ndims-1))
    THError("nhood should be #edges * 3");

  // nb of vertices
  long nVert = target_nelts;

  // convert n-d offset vectors (nhood) into linear array offset scalars
  vector<long> nHood(nhood_nedges);
  for (long i=0; i<nhood_nedges; ++i) {
    nHood[i] = 0;
    for (long j=0; j<nhood_geom; ++j) {
      nHood[i] += (long)nhood_data[i*nhood_geom + j] * conn->stride[j+1];
    }
// cout << "nHood[" << i << "]=" << nHood[i] << endl;
  }

  // disjoint sets and sparse overlap vectors
  long nLabeledVert=0;
  long nPairPos=0;
  vector<map<long,long> > overlap(nVert);
  vector<long> rank(nVert);
  vector<long> parent(nVert);
  map<long,long> segSizes;
  boost::disjoint_sets<long*, long*> dsets(&rank[0],&parent[0]);
  for (long i=0; i<nVert; ++i) {
    dsets.make_set(i);
    if (target_data[i] != 0) {
      overlap[i].insert(pair<long,long>(target_data[i],1));
      ++nLabeledVert;
      ++segSizes[target_data[i]];
      nPairPos += (segSizes[target_data[i]]-1);
    }
  }
  long nPairTot = (nLabeledVert*(nLabeledVert-1))/2;
  long nPairNeg = nPairTot - nPairPos;
  long nPairNorm = 0;
  if (pos) nPairNorm += nPairPos;
  if (neg) nPairNorm += nPairNeg;
  if (nPairNorm == 0) {
    if (!pos && !neg) {
      THError("found 0 pairs, aborting... atleast one of 'posexample' or 'negexample' should be true!");
    } else {
      cout << "MalisCriterion.c: Zero pixel pair error!" << endl;
      for (map<long,long>::iterator it = segSizes.begin(); it != segSizes.end(); ++it) {
        cout << it->first << "," << it->second << endl;
      }
      THError("found 0 pairs, aborting... (this is a bug, please fix me) \n\
        pos = %d, neg = %d \n\
        nPairPos = %d, nPairNeg = %d \n\
        nPairTot = %d, nVert = %d \n",
        pos, neg,
        nPairPos, nPairNeg,
        nPairTot, nVert);
    }
  }

  // Sort all the edges in increasing order of weight
  vector<long> pqueue( conn_nelts );
  long j = 0, idx[3], nextidx;
  bool oob;
  for (long e = 0, i = 0; e < conn->size[0]; ++e)
    for (idx[0] = 0; idx[0] < conn->size[1]; ++idx[0])
      for (idx[1] = 0; idx[1] < conn->size[2]; ++idx[1])
        for (idx[2] = 0; idx[2] < conn->size[3]; ++idx[2], ++i) {
          // check bounds to make sure the other vertex isn't outside the image bounds
          oob = false;
          for (int coord=0; coord < 3; ++coord) {
            nextidx = idx[coord]+nhood_data[e*nhood_geom + coord];
            oob |= (nextidx < 0) || (nextidx >= target->size[coord]);
          }
          if ( !oob )
            pqueue[ j++ ] = i;
        }
// cout << "pqueue size: " << pqueue.size() << endl;
  pqueue.resize(j);
// cout << "pqueue size: " << pqueue.size() << endl;
  sort(pqueue.begin(), pqueue.end(), nn_(AffinityGraphCompare)(conn_data));

  // start MST (Kruskal)
  long minEdge;
  long e, v1, v2;
  long set1, set2, tmp;
  long nPair = 0;
  long nVert1, nVert2;
  long seg1, seg2;
  real loss=0, l=0;
  long nPairFalseNeg = 0, nPairFalsePos = 0, nPairTrueNeg = 0, nPairTruePos = 0;
  map<long,long>::iterator it1, it2;

  for ( long i = 0; i < pqueue.size(); ++i ) {
    minEdge = pqueue[i];
    e = minEdge/nVert; v1 = minEdge%nVert; v2 = v1+nHood[e];

    set1 = dsets.find_set(v1);
    set2 = dsets.find_set(v2);

    if (set1!=set2){
      dsets.link(set1, set2);

      // compute the dloss for this MST edge
      for (it1 = overlap[set1].begin();
           it1 != overlap[set1].end(); ++it1) {
        for (it2 = overlap[set2].begin();
             it2 != overlap[set2].end(); ++it2) {

          seg1 = it1->first;
          seg2 = it2->first;
          nVert1 = it1->second;
          nVert2 = it2->second;
          nPair = nVert1*nVert2;


          // +ve example pairs
          if (pos && (seg1 == seg2)) {
            if (conn_data[minEdge] <= threshold) { // an error
                nPairFalseNeg += nPair;
            } else {
              nPairTruePos += nPair;
            }
            // hinge loss is used here
            l = max(0.0,margin-(conn_data[minEdge]-threshold));
            loss += nPair * l;
            dloss_data[minEdge] -= nPair * (l > 0);
          }
          // -ve example pairs
          if (neg && (seg1 != seg2)) {
            if (conn_data[minEdge] > threshold) { // an error
                nPairFalsePos += nPair;
            } else {
              nPairTrueNeg += nPair;
            }
            // hinge loss is used here
            l = max(0.0,margin+(conn_data[minEdge]-threshold));
            loss += nPair * l;
            dloss_data[minEdge] += nPair * (l > 0);
          }
        }
      }
      dloss_data[minEdge] /= nPairNorm;

      // move the pixel bags of the non-representative to the representative
      if (dsets.find_set(set1) == set2) // make set1 the rep to keep and set2 the rep to empty
        swap(set1,set2);

      for (it2 = overlap[set2].begin();
          it2 != overlap[set2].end(); ++it2) {
        it1 = overlap[set1].find(it2->first);
        if (it1 == overlap[set1].end()) {
          overlap[set1].insert(pair<long,long>(it2->first,it2->second));
        } else {
          it1->second += it2->second;
        }
      }
      overlap[set2].clear();
    } // end link

  } // end while

  // cleanup
  THTensor_(free)(conn);
  THTensor_(free)(target);

  // return loss + class error + rand
  loss /= nPairNorm;
  real classerr = (real) (nPairFalsePos+nPairFalseNeg) / (real)nPairNorm;
  real rand = 1-classerr;
  lua_pushnumber(L, loss);
  lua_pushnumber(L, classerr);
  lua_pushnumber(L, rand);
  return 3;
}

static int nn_(MalisCriterion_backward)(lua_State *L)
{
  return 0;
}

static const struct luaL_Reg nn_(MalisCriterion__) [] = {
  {"MalisCriterion_forward", nn_(MalisCriterion_forward)},
  {"MalisCriterion_backward", nn_(MalisCriterion_backward)},
  {NULL, NULL}
};

static void nn_(MalisCriterion_init)(lua_State *L)
{
  luaT_pushmetatable(L, torch_Tensor);
  luaT_registeratname(L, nn_(MalisCriterion__), "nn");
  lua_pop(L,1);
}

#endif
