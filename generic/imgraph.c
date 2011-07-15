#ifndef TH_GENERIC_FILE
#define TH_GENERIC_FILE "generic/imgraph.c"
#else

#ifdef square
#undef square
#endif
#define square(x) ((x)*(x))

#ifdef epsilon
#undef epsilon
#endif epsilon
#define epsilon 1e-8

static inline real imgraph_(ndiff)(THTensor *img,
                                   int x1, int y1, int x2, int y2, char dt) {
  int nfeats = img->size[0];
  real dist  = 0;
  real dot   = 0;
  real normx = 0;
  real normy = 0;
  real res = 0;
  int i;
  for (i=0; i<nfeats; i++) {
    if (dt == 'e') {
      dist  += square(THTensor_(get3d)(img, i, y1, x1)-THTensor_(get3d)(img, i, y2, x2));
    } else if (dt == 'a') {
      dot   += THTensor_(get3d)(img, i, y1, x1) * THTensor_(get3d)(img, i, y2, x2);
      normx += square(THTensor_(get3d)(img, i, y1, x1));
      normy += square(THTensor_(get3d)(img, i, y2, x2));
    }
  }
  if (dt == 'e') res = sqrt(dist);
  else if (dt == 'a') res = acos(dot/(sqrt(normx)*sqrt(normy) + epsilon));
  return res;
}

static int imgraph_(graph)(lua_State *L) {
  // get args
  THTensor *dst = luaT_checkudata(L, 1, torch_(Tensor_id));
  THTensor *src = luaT_checkudata(L, 2, torch_(Tensor_id));
  int connex = lua_tonumber(L, 3);
  const char *dist = lua_tostring(L, 4);
  char dt = dist[0];

  // compute all edge weights
  if (connex == 4) {

    // get input dims
    THTensor *src3d = NULL;
    long channels, height, width;
    if (src->nDimension == 3) {
      channels = src->size[0];
      height = src->size[1];
      width = src->size[2];
      src3d = src;
    } else if (src->nDimension == 2) {
      channels = 1;
      height = src->size[0];
      width = src->size[1];
      src3d = THTensor_(newWithTensor)(src);
      THTensor_(resize3d)(src3d, channels, height, width);
    }

    // resize output, and fill it with -1 (which means non-valid edge)
    THTensor_(resize3d)(dst, 2, height, width);
    THTensor_(fill)(dst, -1);

    // build graph with 4-connexity
    long num = 0;
    long x,y;
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x < width-1) {
          THTensor_(set3d)(dst, 0, y, x, imgraph_(ndiff)(src3d, x, y, x+1, y, dt));
          num++;
        }
        if (y < height-1) {
          THTensor_(set3d)(dst, 1, y, x, imgraph_(ndiff)(src3d, x, y, x, y+1, dt));
          num++;
        }
      }
    }

    // cleanup
    if (src->nDimension == 2) {
      THTensor_(free)(src3d);
    }

  } else if (connex == 8) {

    // get input dims
    THTensor *src3d = NULL;
    long channels, height, width;
    if (src->nDimension == 3) {
      channels = src->size[0];
      height = src->size[1];
      width = src->size[2];
      src3d = src;
    } else if (src->nDimension == 2) {
      channels = 1;
      height = src->size[0];
      width = src->size[1];
      src3d = THTensor_(newWithTensor)(src);
      THTensor_(resize3d)(src3d, channels, height, width);
    }

    // resize output, and fill it with -1 (which means non-valid edge)
    THTensor_(resize3d)(dst, 4, height, width);
    THTensor_(fill)(dst, -1);

    // build graph with 8-connexity
    long num = 0;
    long x,y;
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x < width-1) {
          THTensor_(set3d)(dst, 0, y, x, imgraph_(ndiff)(src3d, x, y, x+1, y, dt));
          num++;
        }
        if (y < height-1) {
          THTensor_(set3d)(dst, 1, y, x, imgraph_(ndiff)(src3d, x, y, x, y+1, dt));
          num++;
        }
        if ((x < width-1) && (y < height-1)) {
          THTensor_(set3d)(dst, 2, y, x, imgraph_(ndiff)(src3d, x, y, x+1, y+1, dt));
          num++;
        }
        if ((x < width-1) && (y > 0)) {
          THTensor_(set3d)(dst, 3, y, x, imgraph_(ndiff)(src3d, x, y, x+1, y-1, dt));
          num++;
        }
      }
    }

    // cleanup
    if (src->nDimension == 2) {
      THTensor_(free)(src3d);
    }

  }

  return 0;
}

static int imgraph_(connectedcomponents)(lua_State *L) {
  // get args
  THTensor *dst = luaT_checkudata(L, 1, torch_(Tensor_id));
  THTensor *src = luaT_checkudata(L, 2, torch_(Tensor_id));
  real threshold = lua_tonumber(L, 3);
  int color = lua_toboolean(L, 4);

  // dims
  long edges = src->size[0];
  long height = src->size[1];
  long width = src->size[2];

  // make a disjoint-set forest
  Set *set = set_new(width*height);

  // process in one pass
  int x,y;
  for (y = 0; y < height-1; y++) {
    for (x = 0; x < width-1; x++) {
      int a = set_find(set, y*width + x);
      int b = set_find(set, y*width + x+1);
      int c = set_find(set, (y+1)*width + x);
      if ((a != b) && (THTensor_(get3d)(src, 0, y, x) < threshold)) set_join(set, a, b);
      if ((a != c) && (THTensor_(get3d)(src, 1, y, x) < threshold)) set_join(set, a, c);
    }
  }

  // generate output
  if (color) {
    THTensor *colormap = THTensor_(newWithSize2d)(width*height, 3);
    THTensor_(fill)(colormap, -1);
    THTensor_(resize3d)(dst, 3, height, width);
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        int comp = set_find(set, y * width + x);
        real check = THTensor_(get2d)(colormap, comp, 0);
        if (check == -1) {
          THTensor_(set2d)(colormap, comp, 0, random());
          THTensor_(set2d)(colormap, comp, 1, random());
          THTensor_(set2d)(colormap, comp, 2, random());
        }
        real r = THTensor_(get2d)(colormap, comp, 0);
        real g = THTensor_(get2d)(colormap, comp, 1);
        real b = THTensor_(get2d)(colormap, comp, 2);
        THTensor_(set3d)(dst, 0, y, x, r);
        THTensor_(set3d)(dst, 1, y, x, g);
        THTensor_(set3d)(dst, 2, y, x, b);
      }
    }
  } else {
    THTensor_(resize2d)(dst, height, width);
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        int comp = set_find(set, y * width + x);
        THTensor_(set2d)(dst, y, x, comp);
      }
    }
  }

  // push number of components
  lua_pushnumber(L, set->nelts);

  // cleanup
  set_free(set);

  // return
  return 1;
}

#ifndef _EDGE_STRUCT_
#define _EDGE_STRUCT_
typedef struct {
  float w;
  int a, b;
} Edge;

void sort_edges(Edge *data, int N)
{
  int i, j;
  real v;
  Edge t;

  if(N<=1) return;

  // Partition elements
  v = data[0].w;
  i = 0;
  j = N;
  for(;;)
    {
      while(data[++i].w < v && i < N) { }
      while(data[--j].w > v) { }
      if(i >= j) break;
      t = data[i]; data[i] = data[j]; data[j] = t;
    }
  t = data[i-1]; data[i-1] = data[0]; data[0] = t;
  sort_edges(data, i-1);
  sort_edges(data+i, N-i);
}
#endif

static int imgraph_(segmentmst)(lua_State *L) {
  // get args
  THTensor *dst = luaT_checkudata(L, 1, torch_(Tensor_id));
  THTensor *src = luaT_checkudata(L, 2, torch_(Tensor_id));
  real thres = lua_tonumber(L, 3);
  int minsize = lua_tonumber(L, 4);
  int color = lua_toboolean(L, 5);

  // dims
  long nmaps = src->size[0];
  long height = src->size[1];
  long width = src->size[2];

  // create edge list from graph (src)
  Edge *edges = NULL; int nedges = 0;
  edges = calloc(width*height*nmaps, sizeof(Edge));
  int x,y;
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x++) {
      if (x < width-1) {
        edges[nedges].a = y*width+x;
        edges[nedges].b = y*width+(x+1);
        edges[nedges].w = THTensor_(get3d)(src, 0, y, x);
        nedges++;
      }
      if (y < height-1) {
        edges[nedges].a = y*width+x;
        edges[nedges].b = (y+1)*width+x;
        edges[nedges].w = THTensor_(get3d)(src, 1, y, x);
        nedges++;
      }
      if (nmaps >= 4) {
        if ((x < width-1) && (y < height-1)) {
          edges[nedges].a = y * width + x;
          edges[nedges].b = (y+1) * width + (x+1);
          edges[nedges].w = THTensor_(get3d)(src, 2, y, x);
          nedges++;
        }
        if ((x < width-1) && (y > 0)) {
          edges[nedges].a = y * width + x;
          edges[nedges].b = (y-1) * width + (x+1);
          edges[nedges].w = THTensor_(get3d)(src, 3, y, x);
          nedges++;
        }
      }
    }
  }

  // sort edges by weight
  sort_edges(edges, nedges);

  // make a disjoint-set forest
  Set *set = set_new(width*height);

  // init thresholds
  real *threshold = calloc(width*height, sizeof(real));
  int i;
  for (i = 0; i < width*height; i++) threshold[i] = thres;

  // for each edge, in non-decreasing weight order,
  // decide to merge or not, depending on current threshold
  for (i = 0; i < nedges; i++) {
    // components conected by this edge
    int a = set_find(set, edges[i].a);
    int b = set_find(set, edges[i].b);
    if (a != b) {
      if ((edges[i].w <= threshold[a]) && (edges[i].w <= threshold[b])) {
        set_join(set, a, b);
        a = set_find(set, a);
        threshold[a] = edges[i].w + thres/set->elts[a].surface;
      }
    }
  }

  // post process small components
  for (i = 0; i < nedges; i++) {
    int a = set_find(set, edges[i].a);
    int b = set_find(set, edges[i].b);
    if ((a != b) && ((set->elts[a].surface < minsize) || (set->elts[b].surface < minsize)))
      set_join(set, a, b);
  }

  // generate output
  if (color) {
    THTensor *colormap = THTensor_(newWithSize2d)(width*height, 3);
    THTensor_(fill)(colormap, -1);
    THTensor_(resize3d)(dst, 3, height, width);
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        int comp = set_find(set, y * width + x);
        real check = THTensor_(get2d)(colormap, comp, 0);
        if (check == -1) {
          THTensor_(set2d)(colormap, comp, 0, random());
          THTensor_(set2d)(colormap, comp, 1, random());
          THTensor_(set2d)(colormap, comp, 2, random());
        }
        real r = THTensor_(get2d)(colormap, comp, 0);
        real g = THTensor_(get2d)(colormap, comp, 1);
        real b = THTensor_(get2d)(colormap, comp, 2);
        THTensor_(set3d)(dst, 0, y, x, r);
        THTensor_(set3d)(dst, 1, y, x, g);
        THTensor_(set3d)(dst, 2, y, x, b);
      }
    }
  } else {
    THTensor_(resize2d)(dst, height, width);
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        int comp = set_find(set, y * width + x);
        THTensor_(set2d)(dst, y, x, comp);
      }
    }
  }

  // push number of components
  lua_pushnumber(L, set->nelts);

  // cleanup
  set_free(set);
  free(edges);
  free(threshold);

  // return
  return 1;
}

real imgraph_(max)(real *a, int n) {
  int i;
  real max = -1;
  for (i = 0; i < n; i++)
    if (a[i] > max) max = a[i];
  return max;
}

static int imgraph_(gradient)(lua_State *L) {
  // get args
  THTensor *output = luaT_checkudata(L, 1, torch_(Tensor_id));
  THTensor *input = luaT_checkudata(L, 2, torch_(Tensor_id));

  // dims
  long nmaps = input->size[0];
  long height = input->size[1];
  long width = input->size[2];

  // resize output
  THTensor_(resize2d)(output, height, width);

  // compute max gradient
  int x,y;
  if (nmaps == 2) { // 4-connex
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        real edges[] = {-1,-1,-1,-1};
        if (x > 0) edges[0] = THTensor_(get3d)(input, 0, y, x-1);
        if (x < (width-1)) edges[1] = THTensor_(get3d)(input, 0, y, x);
        if (y > 0) edges[2] = THTensor_(get3d)(input, 1, y-1, x);
        if (y < (height-1)) edges[3] = THTensor_(get3d)(input, 1, y, x);
        real max = imgraph_(max)(edges, 4);
        THTensor_(set2d)(output, y, x, max);
      }
    }
  } else if (nmaps == 4) { // 8-connex
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        real edges[] = {-1,-1,-1,-1,-1,-1,-1,-1};
        if (x > 0) edges[0] = THTensor_(get3d)(input, 0, y, x-1);
        if (x < (width-1)) edges[1] = THTensor_(get3d)(input, 0, y, x);
        if (y > 0) edges[2] = THTensor_(get3d)(input, 1, y-1, x);
        if (y < (height-1)) edges[3] = THTensor_(get3d)(input, 1, y, x);

        if ((x > 0) && (y > 0)) edges[4] = THTensor_(get3d)(input, 0, y-1, x-1);
        if ((x < (width-1)) && (y > 0)) edges[5] = THTensor_(get3d)(input, 0, y-1, x);
        if ((y < (height-1)) && (x > 0)) edges[6] = THTensor_(get3d)(input, 1, y, x-1);
        if ((x < (width-1)) && (y < (height-1))) edges[7] = THTensor_(get3d)(input, 1, y, x);

        real max = imgraph_(max)(edges, 8);
        THTensor_(set2d)(output, y, x, max);
      }
    }
  }

  // return number of components
  lua_pushnumber(L, 0);
  return 1;
}

static int imgraph_(watershed)(lua_State *L) {
  // get args
  THTensor *output = luaT_checkudata(L, 1, torch_(Tensor_id));
  THTensor *input = luaT_checkudata(L, 2, torch_(Tensor_id));
  int color = lua_toboolean(L, 3);

  // dims
  long nedges = input->size[0];
  long height = input->size[1];
  long width = input->size[2];

  // return number of components
  lua_pushnumber(L, 0);
  return 1;
}

int imgraph_(colorize)(lua_State *L) {
  // get args
  THTensor *output = luaT_checkudata(L, 1, torch_(Tensor_id));
  THTensor *input = luaT_checkudata(L, 2, torch_(Tensor_id));

  // dims
  long height = input->size[0];
  long width = input->size[1];

  // generate output
  THTensor *colormap = THTensor_(newWithSize2d)(width*height, 3);
  THTensor_(fill)(colormap, -1);
  THTensor_(resize3d)(output, 3, height, width);
  int x,y;
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x++) {
      int id = THTensor_(get2d)(input, y, x);
      real check = THTensor_(get2d)(colormap, id, 0);
      if (check == -1) {
        THTensor_(set2d)(colormap, id, 0, random());
        THTensor_(set2d)(colormap, id, 1, random());
        THTensor_(set2d)(colormap, id, 2, random());
      }
      real r = THTensor_(get2d)(colormap, id, 0);
      real g = THTensor_(get2d)(colormap, id, 1);
      real b = THTensor_(get2d)(colormap, id, 2);
      THTensor_(set3d)(output, 0, y, x, r);
      THTensor_(set3d)(output, 1, y, x, g);
      THTensor_(set3d)(output, 2, y, x, b);
    }
  }

  // dont return anything
  return  0;
}

int imgraph_(histpooling)(lua_State *L) {
  // get args
  THTensor *vectors = luaT_checkudata(L, 1, torch_(Tensor_id));
  THTensor *segm = luaT_checkudata(L, 2, torch_(Tensor_id));
  int computeLists = lua_toboolean(L, 3);
  int histmax = lua_toboolean(L, 4);
  real minConfidence = lua_tonumber(L, 5);

  // check dims
  if ((vectors->nDimension != 3) || (segm->nDimension != 2))
    THError("<imgraph.histpooling> vectors must be KxHxW and segm HxW");

  // get dims
  int depth = vectors->size[0];
  int height = vectors->size[1];
  int width = vectors->size[2];

  // (0) create all necessary tables
  // final cluster list
  lua_newtable(L);  // a = {}
  int table_clean = lua_gettop(L);

  // temporary geometry list
  lua_newtable(L);  // g = {}
  int table_geometry = lua_gettop(L);

  // temporary histogram list
  lua_newtable(L);  // c = {}
  int table_hists = lua_gettop(L);

  // temporary sizes list
  lua_newtable(L);  // s = {}
  int table_sizes = lua_gettop(L);

  // optional confidence map
  THTensor *confidence = THTensor_(newWithSize2d)(width, height);
  THTensor *helper = THTensor_(newWithSize1d)(depth);

  // (1) loop over segm, and accumulate histograms of vectors pixels
  int x,y,k;
  int size;
  THTensor *histo = NULL;
  THTensor *select1 = THTensor_(new)();
  THTensor *select2 = THTensor_(new)();
  for (y=0; y<height; y++) {
    for (x=0; x<width; x++) {
      // compute hash codes for vectors and segm
      int segm_id = THTensor_(get2d)(segm,y,x);
      // is this hash already registered ?
      lua_rawgeti(L,table_hists,segm_id);   // c[segm_id]
      if (lua_isnil(L,-1)) {    // c[segm_id] == nil ?
        lua_pop(L,1);
        // then create a vector to accumulate an histogram of classes,
        // for this cluster
        histo = THTensor_(newWithSize1d)(depth);
        THTensor_(zero)(histo);
        luaT_pushudata(L, histo, torch_(Tensor_id));
        lua_rawseti(L,table_hists,segm_id); // c[segm_id] = histo
        lua_pushnumber(L, 1);
        lua_rawseti(L,table_sizes,segm_id); // s[segm_id] = 1 (initial size = 1)
        size = 1;
      } else {
        // retrieve histo
        histo = luaT_toudata(L, -1, torch_(Tensor_id));
        lua_pop(L,1);
        // retrieve size
        lua_rawgeti(L,table_sizes,segm_id);   // s[segm_id]
        size = lua_tonumber(L,-1);
        lua_pop(L,1);
      }

      // slice the class vector
      THTensor_(select)(select1, vectors, 2, x);
      THTensor_(select)(select2, select1, 1, y);

      // measure confidence
      real local_conf = 1;
      if (minConfidence > 0 ) {
        THTensor_(copy)(helper, select2);
        real max = -1000;
        real idx = 0;
        for (k=0; k<depth; k++) {
          real val = THTensor_(get1d)(helper, k);
          if (val > max) {
            max = val; idx = k;
          }
        }
        THTensor_(set1d)(helper, idx, -1000);
        real max2 = -1000;
        for (k=0; k<depth; k++) {
          real val = THTensor_(get1d)(helper, k);
          if (val > max2) {
            max2 = val;
          }
        }
        local_conf = max-max2;
        if (local_conf < 0) THError("assert error : max < 2nd max");

        // store confidence
        THTensor_(set2d)(confidence, y, x, local_conf);
      }

      // accumulate current vector into histo
      if (local_conf >= minConfidence) {
        THTensor_(cadd)(histo, 1, select2);
        lua_pushnumber(L, ++size);
        lua_rawseti(L,table_sizes,segm_id); // s[segm_id] = ++size
      }
    }
  }

  // (2) then merge vectors into segm, based on the histogram's winners
  THTensor_(zero)(vectors);
  for (y=0; y<height; y++) {
    for (x=0; x<width; x++) {
      // compute hash codes for vectors and segm
      int segm_id = THTensor_(get2d)(segm,y,x);
      // get max
      int argmax = 0, max = -1;
      // retrieve histogram
      lua_rawgeti(L,table_hists,segm_id);   // c[segm_id]  (= histo)
      histo = luaT_toudata(L, -1, torch_(Tensor_id));
      lua_pop(L,1);
      // get geometry entry
      lua_pushinteger(L,segm_id);
      lua_rawget(L,table_geometry);
      if (lua_isnil(L,-1)) {    // g[segm_id] == nil ?
        lua_pop(L,1);
        int i;
        // retrieve size to normalize histogram
        lua_rawgeti(L,table_sizes,segm_id);   // s[segm_id]  (= size)
        size = lua_tonumber(L, -1);
        lua_pop(L,1);
        for (i=0; i<depth; i++) {
          THTensor_(set1d)(histo, i, THTensor_(get1d)(histo, i) / size);
        }
        // compute max
        for (i=0; i<depth; i++) {
          if (max <= THTensor_(get1d)(histo,i)) {
            argmax = i; max = THTensor_(get1d)(histo,i);
          }
        }
        // then create a table to store geometry of component:
        // x,y,size,class,hash
        lua_pushinteger(L,segm_id);
        lua_newtable(L);
        int entry = lua_gettop(L);
        lua_pushnumber(L, x+1);
        lua_rawseti(L,entry,1); // entry[1] = x
        lua_pushnumber(L, y+1);
        lua_rawseti(L,entry,2); // entry[2] = y
        lua_pushnumber(L, 1);
        lua_rawseti(L,entry,3); // entry[3] = size (=1)
        lua_pushnumber(L, argmax+1);
        lua_rawseti(L,entry,4); // entry[4] = class (=argmax+1)
        lua_pushnumber(L, segm_id);
        lua_rawseti(L,entry,5); // entry[5] = hash
        // store entry
        lua_rawset(L,table_geometry); // g[segm_id] = entry
      } else {
        // retrieve entry
        int entry = lua_gettop(L);
        lua_rawgeti(L, entry, 1);
        long cx = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_pushnumber(L, cx+x+1);
        lua_rawseti(L, entry, 1); // entry[1] = cx + x + 1
        lua_rawgeti(L, entry, 2);
        long cy = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_pushnumber(L, cy+y+1);
        lua_rawseti(L, entry, 2); // entry[2] = cy + y + 1
        lua_rawgeti(L, entry, 3);
        long size = lua_tonumber(L, -1) + 1; lua_pop(L, 1);
        lua_pushnumber(L, size);
        lua_rawseti(L, entry, 3); // entry[3] = size + 1
        lua_rawgeti(L, entry, 4);
        argmax = lua_tonumber(L, -1) - 1; lua_pop(L, 1);
        // and clear entry
        lua_pop(L,1);
      }
      // set argmax (winning class) to 1
      if (histmax) {
        THTensor_(set3d)(vectors, argmax, y, x, 1);
      } else {
        int i;
        for (i=0; i<depth; i++) {
          THTensor_(set3d)(vectors, i, y, x, THTensor_(get1d)(histo,i));
        }
      }
    }
  }

  // (3) traverse geometry table to produce final component list
  lua_pushnil(L);
  int cur = 1;
  while (lua_next(L, table_geometry) != 0) {
    // uses 'key' (at index -2) and 'value' (at index -1)

    // normalize cx and cy, by component's size
    int entry = lua_gettop(L);
    lua_rawgeti(L, entry, 3);
    long size = lua_tonumber(L, -1) + 1; lua_pop(L, 1);
    lua_rawgeti(L, entry, 1);
    long cx = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_pushnumber(L, cx/size);
    lua_rawseti(L, entry, 1); // entry[1] = cx/size
    lua_rawgeti(L, entry, 2);
    long cy = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_pushnumber(L, cy/size);
    lua_rawseti(L, entry, 2); // entry[2] = cy/size

    // store entry table into clean table
    lua_rawseti(L, table_clean, cur++);
  }

  // pop/remove histograms + sizes
  lua_pop(L, 2);

  // cleanup
  THTensor_(free)(select1);
  THTensor_(free)(select2);
  THTensor_(free)(helper);

  // return two tables: indexed and hashed, plus the confidence map
  luaT_pushudata(L, confidence, torch_(Tensor_id));
  return 3;
}

static const struct luaL_Reg imgraph_(methods__) [] = {
  {"graph", imgraph_(graph)},
  {"gradient", imgraph_(gradient)},
  {"connectedcomponents", imgraph_(connectedcomponents)},
  {"segmentmst", imgraph_(segmentmst)},
  {"watershed", imgraph_(watershed)},
  {"histpooling", imgraph_(histpooling)},
  {"colorize", imgraph_(colorize)},
  {NULL, NULL}
};

static void imgraph_(Init)(lua_State *L)
{
  luaT_pushmetaclass(L, torch_(Tensor_id));
  luaT_registeratname(L, imgraph_(methods__), "imgraph");
  lua_pop(L,1);
}

#endif
