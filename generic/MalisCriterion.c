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
 **********************************************************************/

using namespace std;

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
  // 3d connectivity graph [#edges * height * width]
  THTensor *conn = (THTensor *)luaT_checkudata(L, 2, torch_Tensor);  
  conn = THTensor_(newContiguous)(conn);
  long conn_ndims = THTensor_(nDimension)(conn);
  long conn_nelts = THTensor_(nElement)(conn);
  long conn_nedges = conn->size[0];
  long conn_height = conn->size[1];
  long conn_width = conn->size[2];
  real *conn_data = THTensor_(data)(conn);

  // target segmentation [height * width]
  THTensor *target = (THTensor *)luaT_checkudata(L, 3, torch_Tensor);  
  target = THTensor_(newContiguous)(target);
  long target_ndims = THTensor_(nDimension)(target);
  long target_nelts = THTensor_(nElement)(target);
  long target_height = target->size[0];
  long target_width = target->size[1];
  real *target_data = THTensor_(data)(target);

  // graph neighborhood descriptor [#edges * 2]
  THTensor *nhood = (THTensor *)luaT_getfieldcheckudata(L, 1, "nhood", torch_Tensor);
  long nhood_ndims = THTensor_(nDimension)(nhood);
  long nhood_nelts = THTensor_(nElement)(nhood);
  long nhood_nedges = nhood->size[0];
  long nhood_geom = nhood->size[1];
  real *nhood_data = THTensor_(data)(nhood);

  // sq-sq loss margin
  const real margin = luaT_getfieldchecknumber(L, 1, "margin");

  // is this a positive example (true) or negative (false)
  int pos = luaT_getfieldcheckboolean(L, 1, "posexample");

  // output: gradient wrt connectivity graph
  THTensor *dloss = (THTensor *)luaT_getfieldcheckudata(L, 1, "gradInput", torch_Tensor);
  THTensor_(resizeAs)(dloss, conn);
  real *dloss_data = THTensor_(data)(dloss);
  THTensor_(zero)(dloss);

  // checks
  if (nhood_ndims != 2)
    THError("nhood should be 2d");
  if (nhood_geom != (conn_ndims-1))
    THError("nhood should be #edges * 2");

  // nb of vertices
  long nVert = conn_nelts / conn_nedges;

  // convert n-d offset vectors (nhood) into linear array offset scalars
  vector<long> nHood(nhood_nedges);
  for (long i=0; i<nhood_nedges; ++i) {
    nHood[i] = 0;
    for (long j=0; j<nhood_geom; ++j) {
      nHood[i] += (long)nhood_data[i*nhood_geom + j] * conn->stride[j+1];
    }
  }

  // disjoint sets and sparse overlap vectors
  vector<map<long,long> > overlap(nVert);
  vector<long> rank(nVert);
  vector<long> parent(nVert);
  map<long,long> segSizes;
  long nLabeledVert=0;
  long nPairPos=0;
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
  long nPairNorm;
  if (pos) {
    if (nPairPos == 0) {
      pos = 0;
      nPairNorm = nPairNeg;
    } else {
      nPairNorm = nPairPos;
    }
  } else {
    if (nPairNeg == 0) {
      pos = 1;
      nPairNorm = nPairPos;
    } else {
      nPairNorm = nPairNeg;
    }
  }
  if (nPairNorm == 0) {
    THError("found 0 pairs, aborting... (this is a bug, please fix me)");
  }

  // Sort all the edges in increasing order of weight
  vector<long> pqueue( conn->size[0] * (conn->size[1]-1) * (conn->size[2]-1) );
  long j = 0;
  for (long d = 0, i = 0; d < conn->size[0]; ++d)
    for (long y = 0; y < conn->size[1]; ++y)
      for (long x = 0; x < conn->size[2]; ++x, ++i)
        if (x < (conn->size[2]-1) && y < (conn->size[1]-1))
          pqueue[ j++ ] = i;
  sort(pqueue.begin(), pqueue.end(), nn_(AffinityGraphCompare)(conn_data));

  // start MST (Kruskal)
  long minEdge;
  long e, v1, v2;
  long set1, set2, tmp;
  long nPair = 0;
  real loss=0, dl=0;
  long nPairIncorrect = 0;
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

          nPair = it1->second * it2->second;

          if (pos && (it1->first == it2->first)) {
            // +ve example pairs
            // Sq-Sq loss is used here
            dl = max(0.0,0.5+margin-conn_data[minEdge]);
            loss += 0.5*dl*dl*nPair;
            dloss_data[minEdge] += dl*nPair;
            if (conn_data[minEdge] <= 0.5) { // an error
              nPairIncorrect += nPair;
            }

          } else if ((!pos) && (it1->first != it2->first)) {
            // -ve example pairs
            // Sq-Sq loss is used here
            dl = max(0.0,conn_data[minEdge]-0.5+margin);
            loss += 0.5*dl*dl*nPair;
            dloss_data[minEdge] += dl*nPair;
            if (conn_data[minEdge] > 0.5) { // an error
              nPairIncorrect += nPair;
            }
          }
        }
      }
      dloss_data[minEdge] /= nPairNorm;

      // move the pixel bags of the non-representative to the representative
      if (dsets.find_set(set1) == set2) // make set1 the rep to keep and set2 the rep to empty
        swap(set1,set2);

      it2 = overlap[set2].begin();
      while (it2 != overlap[set2].end()) {
        it1 = overlap[set1].find(it2->first);
        if (it1 == overlap[set1].end()) {
          overlap[set1].insert(pair<long,long>(it2->first,it2->second));
        } else {
          it1->second += it2->second;
        }
        overlap[set2].erase(it2++);
      }
    } // end link

  } // end while

  // cleanup
  THTensor_(free)(conn);
  THTensor_(free)(target);

  // return loss + class error + rand index
  loss /= nPairNorm;
  real classerr = (real)nPairIncorrect / (real)nPairNorm;
  real randIndex = 1.0 - ((real)nPairIncorrect / (real)nPairNorm);
  lua_pushnumber(L, loss);
  lua_pushnumber(L, classerr);
  lua_pushnumber(L, randIndex);
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
