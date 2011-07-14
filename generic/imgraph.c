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

static int imgraph_(graph2tensor)(lua_State *L) {
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

static const struct luaL_Reg imgraph_(methods__) [] = {
  {"tensor2graph", imgraph_(tensor2graph)},
  {"graph2tensor", imgraph_(graph2tensor)},
  {"watershed", imgraph_(watershed)},
  {NULL, NULL}
};

static void imgraph_(Init)(lua_State *L)
{
  luaT_pushmetaclass(L, torch_(Tensor_id));
  luaT_registeratname(L, imgraph_(methods__), "imgraph");
  lua_pop(L,1);
}

#endif
