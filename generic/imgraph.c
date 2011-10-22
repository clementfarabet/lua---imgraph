#ifndef TH_GENERIC_FILE
#define TH_GENERIC_FILE "generic/imgraph.c"
#else

#ifdef square
#undef square
#endif
#define square(x) ((x)*(x))

#ifdef rand0to1
#undef rand0to1
#endif
#define rand0to1() ((float)rand()/(float)RAND_MAX)

#ifdef epsilon
#undef epsilon
#endif epsilon
#define epsilon 1e-8

static inline real imgraph_(ndiff)(real *img,
                                   int nfeats, int height, int width,
                                   int x1, int y1, int x2, int y2, char dt) {
  real dist  = 0;
  real dot   = 0;
  real normx = 0;
  real normy = 0;
  real res = 0;
  int i;
  for (i=0; i<nfeats; i++) {
    if (dt == 'e') {
      dist  += square( img[(i*height+y1)*width+x1] - img[(i*height+y2)*width+x2] );
    } else if (dt == 'a') {
      dot   += img[(i*height+y1)*width+x1] * img[(i*height+y2)*width+x2];
      normx += square(img[(i*height+y1)*width+x1]);
      normy += square(img[(i*height+y2)*width+x2]);
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

  // make sure input is contiguous
  src = THTensor_(newContiguous)(src);

  // compute all edge weights
  if (connex == 4) {

    // get input dims
    long channels, height, width;
    if (src->nDimension == 3) {
      channels = src->size[0];
      height = src->size[1];
      width = src->size[2];
    } else if (src->nDimension == 2) {
      channels = 1;
      height = src->size[0];
      width = src->size[1];
    }

    // resize output, and fill it with -1 (which means non-valid edge)
    THTensor_(resize3d)(dst, 2, height, width);
    THTensor_(fill)(dst, 0);

    // get raw pointers
    real *src_data = THTensor_(data)(src);
    real *dst_data = THTensor_(data)(dst);

    // build graph with 4-connexity
    long num = 0;
    long x,y;
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x < width-1) {
          dst_data[(0*height+y)*width+x] = imgraph_(ndiff)(src_data, channels, height, width,
                                                           x, y, x+1, y, dt);
          num++;
        }
        if (y < height-1) {
          dst_data[(1*height+y)*width+x] = imgraph_(ndiff)(src_data, channels, height, width,
                                                           x, y, x, y+1, dt);
          num++;
        }
      }
    }

  } else if (connex == 8) {

    // get input dims
    long channels, height, width;
    if (src->nDimension == 3) {
      channels = src->size[0];
      height = src->size[1];
      width = src->size[2];
    } else if (src->nDimension == 2) {
      channels = 1;
      height = src->size[0];
      width = src->size[1];
    }

    // resize output, and fill it with -1 (which means non-valid edge)
    THTensor_(resize3d)(dst, 4, height, width);
    THTensor_(fill)(dst, 0);

    // get raw pointers
    real *src_data = THTensor_(data)(src);
    real *dst_data = THTensor_(data)(dst);

    // build graph with 8-connexity
    long num = 0;
    long x,y;
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        if (x < width-1) {
          dst_data[(0*height+y)*width+x] = imgraph_(ndiff)(src_data, channels, height, width,
                                                           x, y, x+1, y, dt);
          num++;
        }
        if (y < height-1) {
          dst_data[(1*height+y)*width+x] = imgraph_(ndiff)(src_data, channels, height, width,
                                                           x, y, x, y+1, dt);
          num++;
        }
        if ((x < width-1) && (y < height-1)) {
          dst_data[(2*height+y)*width+x] = imgraph_(ndiff)(src_data, channels, height, width,
                                                           x, y, x+1, y+1, dt);
          num++;
        }
        if ((x < width-1) && (y > 0)) {
          dst_data[(3*height+y)*width+x] = imgraph_(ndiff)(src_data, channels, height, width,
                                                           x, y, x+1, y-1, dt);
          num++;
        }
      }
    }

  }

  // cleanup
  THTensor_(free)(src);

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
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x++) {
      int a = set_find(set, y*width + x);
      // 4-connex:
      if (x < width-1) {
        int b = set_find(set, y*width + x+1);
        if ((a != b) && (THTensor_(get3d)(src, 0, y, x) < threshold)) set_join(set, a, b);
      }
      if (y < height-1) {
        int c = set_find(set, (y+1)*width + x);
        if ((a != c) && (THTensor_(get3d)(src, 1, y, x) < threshold)) set_join(set, a, c);
      }
      // 8-connex:
      if (edges >= 4) {
        if ((x < width-1) && (y < height-1)) {
          int d = set_find(set, (y+1)*width + x+1);
          if ((a != d) && (THTensor_(get3d)(src, 2, y, x) < threshold)) set_join(set, a, d);
        }
        if ((y > 0) && (x < width-1)) {
          int e = set_find(set, (y-1)*width + x+1);
          if ((a != e) && (THTensor_(get3d)(src, 3, y, x) < threshold)) set_join(set, a, e);
        }
      }
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
          THTensor_(set2d)(colormap, comp, 0, rand0to1());
          THTensor_(set2d)(colormap, comp, 1, rand0to1());
          THTensor_(set2d)(colormap, comp, 2, rand0to1());
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

  // make sure input is contiguous
  src = THTensor_(newContiguous)(src);
  real *src_data = THTensor_(data)(src);

  // create edge list from graph (src)
  Edge *edges = NULL; int nedges = 0;
  edges = calloc(width*height*nmaps, sizeof(Edge));
  int x,y;
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x++) {
      if (x < width-1) {
        edges[nedges].a = y*width+x;
        edges[nedges].b = y*width+(x+1);
        edges[nedges].w = src_data[(0*height+y)*width+x];
        nedges++;
      }
      if (y < height-1) {
        edges[nedges].a = y*width+x;
        edges[nedges].b = (y+1)*width+x;
        edges[nedges].w = src_data[(1*height+y)*width+x];
        nedges++;
      }
      if (nmaps >= 4) {
        if ((x < width-1) && (y < height-1)) {
          edges[nedges].a = y * width + x;
          edges[nedges].b = (y+1) * width + (x+1);
          edges[nedges].w = src_data[(2*height+y)*width+x];
          nedges++;
        }
        if ((x < width-1) && (y > 0)) {
          edges[nedges].a = y * width + x;
          edges[nedges].b = (y-1) * width + (x+1);
          edges[nedges].w = src_data[(3*height+y)*width+x];
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
          THTensor_(set2d)(colormap, comp, 0, rand0to1());
          THTensor_(set2d)(colormap, comp, 1, rand0to1());
          THTensor_(set2d)(colormap, comp, 2, rand0to1());
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
    real *dst_data = THTensor_(data)(dst);
    for (y = 0; y < height; y++) {
      for (x = 0; x < width; x++) {
        dst_data[y*width+x] = set_find(set, y * width + x);
      }
    }
  }

  // push number of components
  lua_pushnumber(L, set->nelts);

  // cleanup
  set_free(set);
  free(edges);
  free(threshold);
  THTensor_(free)(src);

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

// conversion functions
static struct xvimage * imgraph_(tensor2xv)(THTensor *src, struct xvimage *dst) {
  // create output
  if (dst == NULL)
    dst = allocimage((char *)NULL, src->size[1], src->size[0], 1, VFF_TYP_1_BYTE);

  // get pointer to dst
  uint8_t *dst_data = UCHARDATA(dst);

  // get pointer to src
  THTensor *src_c = THTensor_(newContiguous)(src);
  real *src_data = THTensor_(data)(src_c);

  // copy data
  long i;
  for (i=0; i<(src->size[0]*src->size[1]); i++) {
    real val = ((*src_data++) * 255);
    val = (val < 0) ? 0 : val;
    val = (val > 255) ? 255 : val;
    *dst_data++ = (uint8_t)val;
  }

  // cleanup
  THTensor_(free)(src_c);

  // return copy
  return dst;
}

static THTensor * imgraph_(xv2tensor)(struct xvimage *src, THTensor *dst) {
  // create/resize output
  if (dst == NULL) dst = THTensor_(new)();
  THTensor_(resize2d)(dst, src->col_size, src->row_size);
  real *dst_data = THTensor_(data)(dst);

  // get pointer
  uint8_t *src_data = UCHARDATA(src);

  // copy data
  long i;
  for (i=0; i<(src->col_size*src->row_size); i++) {
    dst_data[i] = ((real)src_data[i]) / 255.0;
  }

  // return copy
  return dst;
}

static struct xvimage * imgraph_(tensor2xvg)(THTensor *src, struct xvimage *dst) {
  // create output
  if (dst == NULL)
    dst = allocGAimage((char *)NULL, src->size[1], src->size[0]/2, 1, VFF_TYP_GABYTE);

  // get pointer to dst
  uint8_t *dst_data = UCHARDATA(dst);

  // get pointer to src
  THTensor *src_c = THTensor_(newContiguous)(src);
  real *src_data = THTensor_(data)(src_c);

  // copy data
  long i;
  for (i=0; i<(src->size[0]*src->size[1]); i++) {
    real val = ((*src_data++) * 255);
    val = (val < 0) ? 0 : val;
    val = (val > 255) ? 255 : val;
    *dst_data++ = (uint8_t)val;
  }

  // cleanup
  THTensor_(free)(src_c);

  // return copy
  return dst;
}

static THTensor * imgraph_(xvg2tensor)(struct xvimage *src, THTensor *dst) {
  // create/resize output
  if (dst == NULL) dst = THTensor_(new)();
  THTensor_(resize2d)(dst, src->col_size*2, src->row_size);
  real *dst_data = THTensor_(data)(dst);

  // get pointer
  uint8_t *src_data = UCHARDATA(src);

  // copy data
  long i;
  for (i=0; i<(dst->size[0]*dst->size[1]); i++) {
    dst_data[i] = ((real)src_data[i]) / 255.0;
  }

  // return copy
  return dst;
}

#ifndef _XV_INVERSE_
#define _XV_INVERSE_
static void inverse(struct xvimage * image) {
  int32_t i, N = rowsize(image) * colsize(image) * depth(image);
  uint8_t *pt;
  for (pt = UCHARDATA(image), i = 0; i < N; i++, pt++)
    *pt = NDG_MAX - *pt;
}
#endif

static int imgraph_(watershed)(lua_State *L) {
  // get args
  THTensor *watershed = luaT_checkudata(L, 1, torch_(Tensor_id));
  THTensor *input = luaT_checkudata(L, 2, torch_(Tensor_id));
  int minHeight = (lua_tonumber(L, 3)) * 255;
  int connex = lua_tonumber(L, 4);

  // dims
  long height = input->size[0];
  long width = input->size[1];

  // make input contiguous
  THTensor *inputc = THTensor_(newContiguous)(input);

  // (1) filter out noisy minimas
  struct xvimage *filtered_xv = imgraph_(tensor2xv)(inputc, NULL);
  inverse(filtered_xv);
  lheightmaxima(filtered_xv, connex, minHeight);
  inverse(filtered_xv);

  // (2) compute minimas of the gradient map
  int nblabels;
  struct xvimage *labels_xv = allocimage(NULL, rowsize(filtered_xv), colsize(filtered_xv), 1, VFF_TYP_4_BYTE);
  llabelextrema(filtered_xv, connex, LABMIN, labels_xv, &nblabels);

  uint8_t *filtered_data = UCHARDATA(filtered_xv);
  int32_t *labels_data = SLONGDATA(labels_xv);
  int i;
  for (i = 0; i < (rowsize(filtered_xv)*colsize(filtered_xv)); i++)
    if (labels_data[i]) filtered_data[i] = NDG_MAX; 
    else filtered_data[i] = NDG_MIN;

  imgraph_(xv2tensor)(filtered_xv, input);

  // (3) compute watershed of input, using minimas
  struct xvimage *input_xv = imgraph_(tensor2xv)(inputc, NULL);
  lwshedtopobin(input_xv, filtered_xv, connex);
  imgraph_(xv2tensor)(input_xv, watershed);

  // cleanup
  freeimage(input_xv);
  THTensor_(free)(inputc);

  // return number of components
  lua_pushnumber(L, nblabels);
  return 1;
}

static int imgraph_(mergetree)(lua_State *L) {
  // get args
  THTensor *graph = luaT_checkudata(L, 1, torch_(Tensor_id));
  graph = THTensor_(newContiguous)(graph);

  // convert
  struct xvimage *graph_xv = imgraph_(tensor2xvg)(graph, NULL);

  // dims
  int32_t rs = rowsize(graph_xv);
  int32_t cs = colsize(graph_xv);

  // compute labels
  struct xvimage *labels = allocimage(NULL,rs,cs,1,VFF_TYP_4_BYTE);
  int32_t *labels_data = SLONGDATA(labels);
  int i; 
  for (i=0; i<rs*cs; i++) labels_data[i] = i;
  //flowMapping(graph_xv, labels_data);

  // construct adjacency graph
  RAG *rag = construitRAG(graph_xv, labels, NULL);

  // compute merge tree
  mtree *mt;
  mergeTree(rag, &mt);

  // compute altitudes
  int32_t *altitudes = (int32_t *)malloc(sizeof(int32_t) * 2 * rag->g->nsom);
  for (i = 0; i < mt->CT->nbnodes; i++) {
    altitudes[i] = mt->CT->tabnodes[i].data;
  }

  // cleanup
  freeimage(graph_xv);
  THTensor_(free)(graph);

  // return tree
  MergeTree *t = lua_pushMergeTree(L);
  t->tree = mt;
  t->labels = labels;
  t->rag = rag;
  t->cs = cs;
  t->rs = rs;
  t->altitudes = altitudes;
  return 1;
}

static int imgraph_(filtertree)(lua_State *L) {
  // get args
  MergeTree *t = lua_toMergeTree(L, 1);
  int mode = 0;
  if (lua_isnumber(L, 2)) mode = lua_tonumber(L, 2);

  // compute merge attributes
  int32_t *attribute;
  switch(mode){
  case 0:
    attribute = surfaceMergeTree(t->tree->CT,t->rag);
    break;
  case 1:
    attribute = dynaMergeTree(t->tree->CT,t->rag); 
    break;
  case 2:
    attribute = volumeMergeTree(t->tree->CT,t->rag);
    break;
  case 3:
    attribute = omegaMergeTree(t->tree->CT,t->rag);
    break;
  }

  // compute spanning tree with attributes
  int32_t *mst = (int32_t *)malloc(sizeof(int32_t) * t->rag->g->nsom - 1);
  int32_t *value = (int32_t *)malloc(sizeof(int32_t) * t->rag->g->nsom - 1);
  int32_t *staltitude = (int32_t *)malloc(sizeof(int32_t) * 2 * t->rag->g->nsom);
  JCctree *st;
  mstCompute(t->tree, mst, value, attribute);
  jcSaliencyTree_b(&st, mst, value, t->rag, staltitude);

  // store new comp tree and attributes
  componentTreeFree(t->tree->CT);
  t->tree->CT = st;
  if (t->altitudes) free(t->altitudes);
  t->altitudes = staltitude;

  // cleanup
  free(mst);
  free(value);
  free(attribute);

  // done
  return 1;
}

static int imgraph_(tree2graph)(lua_State *L) {
  // get args
  MergeTree *t = lua_toMergeTree(L, 1);
  THTensor *graph = luaT_checkudata(L, 2, torch_(Tensor_id));

  // alloc graph
  struct xvimage *graph_xv = allocGAimage((char *)NULL, t->rs, t->cs, 1, VFF_TYP_GABYTE);

  // generate graph
  computeSaliencyMap(t->tree->CT, graph_xv, SLONGDATA(t->labels), t->altitudes);

  // export
  imgraph_(xvg2tensor)(graph_xv, graph);

  // done
  return 1;
}

static int imgraph_(graph2map)(lua_State *L) {
  // get args
  THTensor *rendered = luaT_checkudata(L, 1, torch_(Tensor_id));
  THTensor *graph = luaT_checkudata(L, 2, torch_(Tensor_id));

  // render graph using khalimsky representation
  struct xvimage *graph_xv = imgraph_(tensor2xvg)(graph, NULL);
  struct xvimage *rendered_xv = allocimage(NULL, rowsize(graph_xv)*2 , colsize(graph_xv)*2, 1, VFF_TYP_1_BYTE);
  lga2khalimsky(graph_xv, rendered_xv, 0);
  imgraph_(xv2tensor)(rendered_xv, rendered);

  // cleanup
  freeimage(graph_xv);

  // return
  return 0;
}

int imgraph_(colorize)(lua_State *L) {
  // get args
  THTensor *output = luaT_checkudata(L, 1, torch_(Tensor_id));
  THTensor *input = luaT_checkudata(L, 2, torch_(Tensor_id));
  THTensor *colormap = luaT_checkudata(L, 3, torch_(Tensor_id));

  // dims
  long height = input->size[0];
  long width = input->size[1];

  // generate color map if not given
  if (THTensor_(nElement)(colormap) == 0) {
    THTensor_(resize2d)(colormap, width*height, 3);
    THTensor_(fill)(colormap, -1);
  }

  // colormap channels
  int channels = colormap->size[1];

  // generate output
  THTensor_(resize3d)(output, channels, height, width);
  int x,y,k;
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x++) {
      int id = THTensor_(get2d)(input, y, x);
      real check = THTensor_(get2d)(colormap, id, 0);
      if (check == -1) {
        for (k = 0; k < channels; k++) {
          THTensor_(set2d)(colormap, id, k, rand0to1());
        }
      }
      for (k = 0; k < channels; k++) {
        real color = THTensor_(get2d)(colormap, id, k);
        THTensor_(set3d)(output, k, y, x, color);
      }
    }
  }

  // return nothing
  return 0;
}

#ifndef __setneighbor__
#define __setneighbor__
static inline void setneighbor(lua_State *L, long matrix, long id, long idn) {
  // retrieve or create table at index 'id'
  lua_rawgeti(L, matrix, id);
  if (lua_isnil(L, -1)) {
    lua_pop(L, 1);
    lua_createtable(L, 32, 32); // pre-alloc for 32 neighbors
  }
  // append idn
  lua_pushboolean(L, 1);
  lua_rawseti(L, -2, idn);
  // write table back
  lua_rawseti(L, matrix, id);
}
#endif

int imgraph_(adjacency)(lua_State *L) {
  // get args
  THTensor *input = THTensor_(newContiguous)(luaT_checkudata(L, 1, torch_(Tensor_id)));
  long matrix = 2;

  // dims
  long height = input->size[0];
  long width = input->size[1];

  // raw pointers
  real *input_data = THTensor_(data)(input);

  // generate output
  int x,y;
  for (y = 0; y < height; y++) {
    for (x = 0; x < width; x++) {
      long id = input_data[width*y+x];
      if (x < (width-1)) {
        long id_east = input_data[width*y+x+1];
        if (id != id_east) {
          setneighbor(L, matrix, id, id_east);
          setneighbor(L, matrix, id_east, id);
        }
      }
      if (y < (height-1)) {
        long id_south = input_data[width*(y+1)+x];
        if (id != id_south) {
          setneighbor(L, matrix, id, id_south);
          setneighbor(L, matrix, id_south, id);
        }
      }
    }
  }

  // cleanup
  THTensor_(free)(input);

  // return matrix
  return 1;
}

int imgraph_(histpooling)(lua_State *L) {
  // get args
  THTensor *vectors = luaT_checkudata(L, 1, torch_(Tensor_id));
  THTensor *segm = luaT_checkudata(L, 2, torch_(Tensor_id));
  int computeLists = lua_toboolean(L, 3);
  int histmax = lua_toboolean(L, 4);
  real minConfidence = lua_tonumber(L, 5);

  // check dims
  if ((vectors->nDimension != 3) || (segm->nDimension != 2) 
      || (segm->size[0] != vectors->size[1]) || (segm->size[1] != vectors->size[2]))
    THError("<imgraph.histpooling> vectors must be KxHxW and segm HxW");

  // get dims
  int depth = vectors->size[0];
  int height = vectors->size[1];
  int width = vectors->size[2];

  // get raw pointers to tensors
  segm = THTensor_(newContiguous)(segm);
  real *segm_data = THTensor_(data)(segm);

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
  THTensor *confidence = THTensor_(newWithSize2d)(height, width);
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
      int segm_id = segm_data[y*width+x];
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
        real max = -THInf;
        real idx = 0;
        for (k=0; k<depth; k++) {
          real val = THTensor_(get1d)(helper, k);
          if (val > max) {
            max = val; idx = k;
          }
        }
        THTensor_(set1d)(helper, idx, -THInf);
        real max2 = -THInf;
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
      int segm_id = segm_data[y*width+x];
      // get max
      int argmax = 0;
      real max = -THInf;
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
        THTensor *entry = THTensor_(newWithSize1d)(10);
        real *data = THTensor_(data)(entry);
        data[0] = x+1;       // x
        data[1] = y+1;       // y
        data[2] = 1;         // size
        data[3] = argmax+1;  // class
        data[4] = segm_id;   // hash

        // store entry
        lua_pushinteger(L,segm_id);
        luaT_pushudata(L, entry, torch_(Tensor_id));
        lua_rawset(L,table_geometry); // g[segm_id] = entry
      } else {
        // retrieve entry
        THTensor *entry = luaT_toudata(L, -1, torch_(Tensor_id));
        real *data = THTensor_(data)(entry);
        data[0] += x+1;       // cx + x + 1
        data[1] += y+1;       // cy + y + 1
        data[2] += 1;         // csize + size + 1
        argmax = data[3]-1;   // retrieve argmax

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

    // normalize cx and cy, by component's size */
    THTensor *entry = luaT_toudata(L, -1, torch_(Tensor_id));
    real *data = THTensor_(data)(entry);
    long size = data[2];
    data[0] /= size;  // cx/size
    data[1] /= size;  // cy/size

    // store entry table into clean table
    lua_rawseti(L, table_clean, cur++);
  }

  // pop/remove histograms + sizes
  lua_pop(L, 2);

  // cleanup
  THTensor_(free)(select1);
  THTensor_(free)(select2);
  THTensor_(free)(helper);
  THTensor_(free)(segm);

  // return two tables: indexed and hashed, plus the confidence map
  luaT_pushudata(L, confidence, torch_(Tensor_id));
  return 3;
}

int imgraph_(segm2components)(lua_State *L) {
  // get args
  THTensor *segm = luaT_checkudata(L, 1, torch_(Tensor_id));
  real *segm_data = THTensor_(data)(segm);

  // check dims
  if ((segm->nDimension != 2))
    THError("<imgraph.segm2components> segm must be HxW");

  // get dims
  int height = segm->size[0];
  int width = segm->size[1];

  // (0) create two tables, one indexable, the other one hash-indexable
  lua_newtable(L);
  int table_hash = lua_gettop(L);

  // (1) extra components' info
  long x,y;
  for (y=0; y<height; y++) {
    for (x=0; x<width; x++) {
      // get component ID
      int segm_id = segm_data[width*y+x];

      // get geometry entry
      lua_pushinteger(L,segm_id);
      lua_rawget(L,table_hash);
      if (lua_isnil(L,-1)) {
        // g[segm_id] = nil
        lua_pop(L,1);

        // then create a table to store geometry of component:
        // x,y,size,class,hash
        THTensor *entry = THTensor_(newWithSize1d)(13);
        real *data = THTensor_(data)(entry);
        data[0] = x+1;       // x
        data[1] = y+1;       // y
        data[2] = 1;         // size
        data[3] = 0;         // compat with 'histpooling' method
        data[4] = segm_id;   // hash
        data[5] = x+1;       // left_x
        data[6] = x+1;       // right_x
        data[7] = y+1;       // top_y
        data[8] = y+1;       // bottom_y

        // store entry
        lua_pushinteger(L,segm_id);
        luaT_pushudata(L, entry, torch_(Tensor_id));
        lua_rawset(L,table_hash); // g[segm_id] = entry

      } else {
        // retrieve entry
        THTensor *entry = luaT_toudata(L, -1, torch_(Tensor_id)); 
        lua_pop(L,1);

        // update content
        real *data = THTensor_(data)(entry);
        data[0] += x+1;       // x += x + 1
        data[1] += y+1;       // y += y + 1
        data[2] += 1;         // size += 1
        data[5] = (x+1)<data[5] ? x+1 : data[5];   // left_x
        data[6] = (x+1)>data[6] ? x+1 : data[6];   // right_x
        data[7] = (y+1)<data[7] ? y+1 : data[7];   // top_y
        data[8] = (y+1)>data[8] ? y+1 : data[8];   // bottom_y
      }
    }
  }

  // (2) traverse geometry table to produce final component list
  lua_pushnil(L);
  int cur = 1;
  while (lua_next(L, table_hash) != 0) {
    // retrieve entry
    THTensor *entry = luaT_toudata(L, -1, torch_(Tensor_id)); lua_pop(L,1);
    real *data = THTensor_(data)(entry);

    // normalize cx and cy, by component's size
    long size = data[2];
    data[0] /= size;  // cx/size
    data[1] /= size;  // cy/size

    // extra info
    data[9] = data[6] - data[5] + 1;     // box width
    data[10] = data[8] - data[7] + 1;    // box height
    data[11] = (data[6] + data[5]) / 2;  // box center x
    data[12] = (data[8] + data[7]) / 1;  // box center y
  }

  // return component table
  return 1;
}

static const struct luaL_Reg imgraph_(methods__) [] = {
  {"graph", imgraph_(graph)},
  {"gradient", imgraph_(gradient)},
  {"connectedcomponents", imgraph_(connectedcomponents)},
  {"segmentmst", imgraph_(segmentmst)},
  {"watershed", imgraph_(watershed)},
  {"histpooling", imgraph_(histpooling)},
  {"colorize", imgraph_(colorize)},
  {"graph2map", imgraph_(graph2map)},
  {"mergetree", imgraph_(mergetree)},
  {"filtertree", imgraph_(filtertree)},
  {"tree2graph", imgraph_(tree2graph)},
  {"adjacency", imgraph_(adjacency)},
  {"segm2components", imgraph_(segm2components)},
  {NULL, NULL}
};

static void imgraph_(Init)(lua_State *L)
{
  luaT_pushmetaclass(L, torch_(Tensor_id));
  luaT_registeratname(L, imgraph_(methods__), "imgraph");
  lua_pop(L,1);

  MergeTree_register(L);
}

#endif
