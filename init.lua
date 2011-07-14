----------------------------------------------------------------------
--
-- Copyright (c) 2011 Clement Farabet
-- 
-- Permission is hereby granted, free of charge, to any person obtaining
-- a copy of this software and associated documentation files (the
-- "Software"), to deal in the Software without restriction, including
-- without limitation the rights to use, copy, modify, merge, publish,
-- distribute, sublicense, and/or sell copies of the Software, and to
-- permit persons to whom the Software is furnished to do so, subject to
-- the following conditions:
-- 
-- The above copyright notice and this permission notice shall be
-- included in all copies or substantial portions of the Software.
-- 
-- THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
-- EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
-- MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
-- NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
-- LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
-- OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
-- WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
-- 
----------------------------------------------------------------------
-- description:
--     imgraph - a graph package for images.
--
-- history: 
--     July 13, 2011, 6:16PM - first draft - Clement Farabet
----------------------------------------------------------------------

require 'torch'
require 'xlua'
require 'image'

-- create global nnx table:
imgraph = {}

-- c lib:
require 'libimgraph'

----------------------------------------------------------------------
-- convert image <-> graph
--
function imgraph.tensor2graph(...)
   -- get args
   local args = {...}
   local dest, img, connex, distance
   local arg2 = torch.typename(args[2])
   if arg2 and arg2:find('Tensor') then
      dest = args[1]
      img = args[2]
      connex = args[3]
      distance = args[4]
   else
      dest = torch.Tensor()
      img = args[1]
      connex = args[2]
      distance = args[3]
   end

   -- defaults
   connex = connex or 4
   distance = ((not distance) and 'e') or ((distance == 'euclid') and 'e') or ((distance == 'angle') and 'a')

   -- usage
   if not img or (connex ~= 4 and connex ~= 8) or (distance ~= 'e' and distance ~= 'a') then
      print(xlua.usage('imgraph.tensor2graph',
                       'computes an edge-weighted graph of an image', nil,
                       {type='torch.Tensor', help='input tensor (for now KxHxW or HxW)', req=true},
                       {type='number', help='connexity (edges per vertex): 4 | 8', default=4},
                       {type='string', help='distance metric: euclid | angle', req='euclid'},
                       "",
                       {type='torch.Tensor', help='destination: existing graph', req=true},
                       {type='torch.Tensor', help='input tensor (for now KxHxW or HxW)', req=true},
                       {type='number', help='connexity (edges per vertex): 4 | 8', default=4},
                       {type='string', help='distance metric: euclid | angle', req='euclid'}))
      xlua.error('incorrect arguments', 'imgraph.tensor2graph')
   end

   -- compute graph
   img.imgraph.tensor2graph(dest, img, connex, distance)

   -- return result
   return dest
end

function imgraph.graph2tensor(...)
   --get args
   local args = {...}
   local dest, graph, threshold, colorize, watershed
   if arg2 and arg2:find('Tensor') then
      dest = args[1]
      graph = args[2]
      threshold = args[3]
      colorize = args[4]
   else
      dest = torch.Tensor()
      graph = args[1]
      threshold = args[2]
      colorize = args[3]
   end

   -- defaults
   if threshold == 'watershed' then watershed = true end
   threshold = threshold or 0.5

   -- usage
   if not graph then
      print(xlua.usage('imgraph.graph2tensor',
                       'generates an image from an edge-weighted graph (connected components)', nil,
                       {type='torch.Tensor', help='input graph', req=true},
                       {type='number | string', help='threshold for connected components, or "watershed"', default=0.5},
                       {type='boolean', help='replace components id by random colors', default=false},
                       "",
                       {type='torch.Tensor', help='destination tensor', req=true},
                       {type='torch.Tensor', help='input graph', req=true},
                       {type='number | string', help='threshold for connected components, or "watershed"', default=0.5},
                       {type='boolean', help='replace components id by random colors', default=false}))
      xlua.error('incorrect arguments', 'imgraph.graph2tensor')
   end

   -- compute image
   local nelts
   if watershed then
      nelts = graph.imgraph.watershed(dest, graph, colorize)
   else
      nelts = graph.imgraph.connectedcomponents(dest, graph, threshold, colorize)
   end

   -- return image
   return dest, nelts
end

----------------------------------------------------------------------
-- a simple test me function
--
function imgraph.testme()
   -- load an image, compute its graph, and convert back to an image
   local lena = image.lena()
   lena = image.convolve(lena, image.gaussian(3))
   local graph = imgraph.tensor2graph(lena)
   local imaged = imgraph.graph2tensor(graph,0.1,true)
   image.display{image=graph, legend='left: horizontal edges, right: vertical edges'}
   image.display{image=imaged, legend='thresholded graph'}
end
