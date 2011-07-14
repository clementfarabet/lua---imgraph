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

static int imgraph_(tensor2graph)(lua_State *L) {
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

int imgraph_(spatialhistpooling)(lua_State *L) {
  // get args
  THTensor *vectors = luaT_checkudata(L, 1, torch_(Tensor_id));
  THTensor *segm = luaT_checkudata(L, 2, torch_(Tensor_id));
  int minConfidence = 0;
  if (lua_isnumber(L,3)) minConfidence = lua_tonumber(L, 3);

  // check dims
  if ((vectors->nDimension != 3) || (segm->nDimension != 2))
    THError("<imgraph.spatialhist> vectors must be KxHxW and segm HxW");

  // get dims
  int nbClasses = vectors->size[0];
  int height = vectors->size[1];
  int width = vectors->size[2];

  // final cluster list
  lua_newtable(L);  // f = {}
  int table_clean = lua_gettop(L);

  // temporary geometry list
  lua_newtable(L);  // g = {}
  int table_geometry = lua_gettop(L);

  // temporary histogram list
  lua_newtable(L);  // c = {}
  int table_hists = lua_gettop(L);

  // optional confidence map
  THTensor *confidence = THTensor_(newWithSize2d)(width, height);
  THTensor *helper = THTensor_(newWithSize1d)(nbClasses);

  // loop over segm, and accumulate histograms of vectors pixels
  int x,y,k;
  THTensor *histo = NULL;
  THTensor *select1 = THTensor_(new)();
  THTensor *select2 = THTensor_(new)();
  for (y=0; y<height; y++) {
    for (x=0; x<width; x++) {
      // compute hash codes for vectors and segm
      int segm_id = THTensor_(get2d)(segm,y,x);
      // is this hash already registered ?
      lua_pushinteger(L,segm_id);
      lua_rawget(L,table_hists);   // c[segm_id]
      if (lua_isnil(L,-1)) {    // c[segm_id] == nil ?
        lua_pop(L,1);
        // then create a vector to accumulate an histogram of classes,
        // for this cluster
        histo = THTensor_(newWithSize1d)(nbClasses);
        THTensor_(zero)(histo);
        lua_pushinteger(L,segm_id);
        luaT_pushudata(L, histo, torch_(Tensor_id));
        lua_rawset(L,table_hists); // c[segm_id] = histo
      } else {
        // retrieve histo
        histo = luaT_toudata(L, -1, torch_(Tensor_id));
        lua_pop(L,1);
      }

      // slice the class vector
      THTensor_(select)(select1, vectors, 2, x);
      THTensor_(select)(select2, select1, 1, y);

      // measure confidence
      THTensor_(copy)(helper, select2);
      real max = -1000;
      real idx = 0;
      for (k=0; k<nbClasses; k++) {
        real val = THTensor_(get1d)(helper, k);
        if (val > max) {
          max = val; idx = k;
        }
      }
      THTensor_(set1d)(helper, idx, -1000);
      real max2 = -1000;
      for (k=0; k<nbClasses; k++) {
        real val = THTensor_(get1d)(helper, k);
        if (val > max2) {
          max2 = val;
        }
      }
      real local_conf = max-max2;
      if (local_conf < 0) THError("assert error : max < 2nd max");

      // accumulate current vector into histo
      if (local_conf >= minConfidence)
        THTensor_(cadd)(histo, 1, select2);

      // store confidence
      THTensor_(set2d)(confidence, y, x, local_conf);
    }
  }

  // then merge vectors into segm, based on the histogram's winners
  THTensor_(zero)(vectors);
  for (y=0; y<height; y++) {
    for (x=0; x<width; x++) {
      // compute hash codes for vectors and segm
      int segm_id = THTensor_(get2d)(segm,y,x);
      // get max
      int argmax = 0, max = -1;
      // get geometry entry
      lua_pushinteger(L,segm_id);
      lua_rawget(L,table_geometry);
      if (lua_isnil(L,-1)) {    // g[segm_id] == nil ?
        lua_pop(L,1);
        // retrieve histogram
        lua_pushinteger(L,segm_id);
        lua_rawget(L,table_hists);   // c[segm_id]  (= histo)
        histo = luaT_toudata(L, -1, torch_(Tensor_id));
        lua_pop(L,1);
        // compute max
        int i;
        for (i=0; i<nbClasses; i++) {
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
      THTensor_(set3d)(vectors, argmax, y, x, 1);
    }
  }

  // traverse geometry table to produce final component list
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

  // pop/remove histograms
  lua_pop(L, 1);

  // cleanup
  THTensor_(free)(select1);
  THTensor_(free)(select2);
  THTensor_(free)(helper);

  // return two tables: indexed and hashed, plus the confidence map
  luaT_pushudata(L, confidence, torch_(Tensor_id));
  return 3;
}

static const struct luaL_Reg imgraph_(methods__) [] = {
  {"tensor2graph", imgraph_(tensor2graph)},
  {"connectedcomponents", imgraph_(connectedcomponents)},
  {"watershed", imgraph_(watershed)},
  {"spatialhistpooling", imgraph_(spatialhistpooling)},
  {NULL, NULL}
};

static void imgraph_(Init)(lua_State *L)
{
  luaT_pushmetaclass(L, torch_(Tensor_id));
  luaT_registeratname(L, imgraph_(methods__), "imgraph");
  lua_pop(L,1);
}

#endif
