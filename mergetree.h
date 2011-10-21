#ifndef _MERGETREE_
#define _MERGETREE_

// create a proper Lua class to represent a merge tree
#define MT "imgraph.MergeTree"

typedef struct {
  mtree *tree;
  RAG *rag;
  struct xvimage *labels;
  int32_t *altitudes;
  int cs;
  int rs;
} MergeTree;

static MergeTree *lua_toMergeTree (lua_State *L, int index)
{
  MergeTree *mt = (MergeTree *)lua_touserdata(L, index);
  if (mt == NULL) luaL_typerror(L, index, MT);
  return mt;
}

static MergeTree *lua_checkMergeTree (lua_State *L, int index)
{
  MergeTree *mt;
  luaL_checktype(L, index, LUA_TUSERDATA);
  mt = (MergeTree *)luaL_checkudata(L, index, MT);
  if (mt == NULL) luaL_typerror(L, index, MT);
  return mt;
}

static MergeTree *lua_pushMergeTree (lua_State *L)
{
  MergeTree *mt = (MergeTree *)lua_newuserdata(L, sizeof(MergeTree));
  mt->tree = NULL;
  mt->labels = NULL;
  mt->rag = NULL;
  mt->altitudes = NULL;
  mt->cs = 0;
  mt->rs = 0;
  luaL_getmetatable(L, MT);
  lua_setmetatable(L, -2);
  return mt;
}

static int MergeTree_gc (lua_State *L)
{
  MergeTree *t = lua_toMergeTree(L, 1);
  if (t->tree) mergeTreeFree(t->tree);
  if (t->labels) freeimage(t->labels);
  if (t->rag) termineRAG(t->rag);
  if (t->altitudes) free(t->altitudes);
  return 0;
}

static int MergeTree_tostring (lua_State *L)
{
  MergeTree *t = lua_toMergeTree(L, 1);
  char *cstr = malloc(10*1024);
  char *str = cstr;
  str += sprintf(str, "<%s>\n", MT);

  int32_t i;
  JCsoncell *s;
  JCctree *CT = t->tree->CT;
  str += sprintf(str, " + root: %d ;  nbnodes: %d ; nbsoncells: %d", CT->root, CT->nbnodes, CT->nbsoncells);
  for (i = 0; i < CT->nbnodes; i++) 
  {
    str += sprintf(str, "\n");
    str += sprintf(str, " - node: %d ; level %d ; nbsons: %d ; father: %d ; ", 
            i, CT->tabnodes[i].data, CT->tabnodes[i].nbsons, CT->tabnodes[i].father);
    if (CT->tabnodes[i].nbsons > 0)
    {
      str += sprintf(str, "sons: ");
      for (s = CT->tabnodes[i].sonlist; s != NULL; s = s->next)
        str += sprintf(str, "%d  ", s->son);
    }
  }

  lua_pushfstring(L, "%s", cstr);
  free(cstr);
  return 1;
}

static const luaL_reg MergeTree_meta[] = {
  {"__gc",       MergeTree_gc},
  {"__tostring", MergeTree_tostring},
  {0, 0}
};

static int MergeTree_register (lua_State *L)
{
  luaL_newmetatable(L, MT);          /* create metatable for Graph,
                                        and add it to the Lua registry */
  luaL_openlib(L, 0, MergeTree_meta, 0);    /* fill metatable */
  lua_pushliteral(L, "__metatable");
  lua_pushvalue(L, -3);               /* dup methods table*/
  lua_rawset(L, -3);                  /* hide metatable:
                                         metatable.__metatable = methods */
  lua_pop(L, 1);                      /* drop metatable */
}

#endif
