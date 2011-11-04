#ifndef TH_GENERIC_FILE
#define TH_GENERIC_FILE "generic/MalisCriterion.c"
#else

using namespace std;

static int nn_(MalisCriterion_forward)(lua_State *L)
{
  THTensor *input = (THTensor *)luaT_checkudata(L, 2, torch_(Tensor_id));  

  return 0;
}

static int nn_(MalisCriterion_backward)(lua_State *L)
{
  THTensor *input = (THTensor *)luaT_checkudata(L, 2, torch_(Tensor_id));
  THTensor *gradInput = (THTensor *)luaT_checkudata(L, 3, torch_(Tensor_id));

  return 0;
}

static const struct luaL_Reg nn_(MalisCriterion__) [] = {
  {"MalisCriterion_forward", nn_(MalisCriterion_forward)},
  {"MalisCriterion_backward", nn_(MalisCriterion_backward)},
  {NULL, NULL}
};

static void nn_(MalisCriterion_init)(lua_State *L)
{
  luaT_pushmetaclass(L, torch_(Tensor_id));
  luaT_registeratname(L, nn_(MalisCriterion__), "nn");
  lua_pop(L,1);
}

#endif
