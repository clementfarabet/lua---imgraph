----------------------------------------------------------------------
--
-- Copyright (c) 2011 Clement Farabet
--               2006 Pedro Felzenszwalb
-- 
-- This program is free software; you can redistribute it and/or modify
-- it under the terms of the GNU General Public License as published by
-- the Free Software Foundation; either version 2 of the License, or
-- (at your option) any later version.
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
-- 
-- You should have received a copy of the GNU General Public License
-- along with this program; if not, write to the Free Software
-- Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
-- 
----------------------------------------------------------------------
-- description:
--     imgraph - a graph package for images:
--               this package contains several routines to build graphs
--               on images, and then compute their connected components,
--               watershed, min-spanning trees, and so on.
--               
--               The min-spanning three based segmentation (segmentmst)
--               originates from Felzenszwalb's 2004 paper:
--               "Efficient Graph-Based Image Segmentation".
-- 
-- history: 
--     July 14, 2011, 12:49PM - MST + HistPooling    - Clement Farabet
--     July 13, 2011, 6:16PM  - first draft          - Clement Farabet
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
function imgraph.graph(...)
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
      print(xlua.usage('imgraph.graph',
                       'computes an edge-weighted graph of an image', nil,
                       {type='torch.Tensor', help='input tensor (for now KxHxW or HxW)', req=true},
                       {type='number', help='connexity (edges per vertex): 4 | 8', default=4},
                       {type='string', help='distance metric: euclid | angle', req='euclid'},
                       "",
                       {type='torch.Tensor', help='destination: existing graph', req=true},
                       {type='torch.Tensor', help='input tensor (for now KxHxW or HxW)', req=true},
                       {type='number', help='connexity (edges per vertex): 4 | 8', default=4},
                       {type='string', help='distance metric: euclid | angle', req='euclid'}))
      xlua.error('incorrect arguments', 'imgraph.graph')
   end

   -- compute graph
   img.imgraph.tensor2graph(dest, img, connex, distance)

   -- return result
   return dest
end

function imgraph.connectcomponents(...)
   --get args
   local args = {...}
   local dest, graph, threshold, colorize
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
   threshold = threshold or 0.5

   -- usage
   if not graph then
      print(xlua.usage('imgraph.connectcomponents',
                       'simple connected components on an edge-weighted graph', nil,
                       {type='torch.Tensor', help='input graph', req=true},
                       {type='number', help='threshold for connecting components', default=0.5},
                       {type='boolean', help='replace components id by random colors', default=false},
                       "",
                       {type='torch.Tensor', help='destination tensor', req=true},
                       {type='torch.Tensor', help='input graph', req=true},
                       {type='number', help='threshold for connecting components', default=0.5},
                       {type='boolean', help='replace components id by random colors', default=false}))
      xlua.error('incorrect arguments', 'imgraph.connectcomponents')
   end

   -- compute image
   local nelts = graph.imgraph.connectedcomponents(dest, graph, threshold, colorize)

   -- return image
   return dest, nelts
end

function imgraph.watershed(...)
   --get args
   local args = {...}
   local dest, graph, colorize
   if arg2 and arg2:find('Tensor') then
      dest = args[1]
      graph = args[2]
      colorize = args[3]
   else
      dest = torch.Tensor()
      graph = args[1]
      colorize = args[2]
   end

   -- usage
   if not graph then
      print(xlua.usage('imgraph.watershed',
                       'copmutes the watershed of an edge-weighted graph', nil,
                       {type='torch.Tensor', help='input graph', req=true},
                       {type='boolean', help='replace components id by random colors', default=false},
                       "",
                       {type='torch.Tensor', help='destination tensor', req=true},
                       {type='torch.Tensor', help='input graph', req=true},
                       {type='boolean', help='replace components id by random colors', default=false}))
      xlua.error('incorrect arguments', 'imgraph.watershed')
   end

   -- compute image
   local nelts = graph.imgraph.watershed(dest, graph, colorize)

   -- return image
   return dest, nelts
end

function imgraph.segmentmst(...)
   --get args
   local args = {...}
   local dest, graph, thres, minsize, colorize
   if arg2 and arg2:find('Tensor') then
      dest = args[1]
      graph = args[2]
      thres = args[3]
      minsize = args[4]
      colorize = args[5]
   else
      dest = torch.Tensor()
      graph = args[1]
      thres = args[2]
      minsize = args[3]
      colorize = args[4]
   end

   -- usage
   if not graph then
      print(xlua.usage('imgraph.segmentmst',
                       'segments an edge-weighted graph, using a surface adaptive criterion\n'
                          .. 'on the min-spanning tree of the graph (see Felzenszwalb et al. 2004)',
                       nil,
                       {type='torch.Tensor', help='input graph', req=true},
                       {type='boolean', help='replace components id by random colors', default=false},
                       {type='number', help='base threshold for merging', default=3},
                       {type='number', help='min size: merge components of smaller size', default=20},
                       "",
                       {type='torch.Tensor', help='destination tensor', req=true},
                       {type='torch.Tensor', help='input graph', req=true},
                       {type='number', help='base threshold for merging', default=3},
                       {type='number', help='min size: merge components of smaller size', default=20},
                       {type='boolean', help='replace components id by random colors', default=false}))
      xlua.error('incorrect arguments', 'imgraph.segmentmst')
   end

   -- compute image
   local nelts = graph.imgraph.segmentmst(dest, graph, thres, minsize, colorize)

   -- return image
   return dest, nelts
end

function imgraph.histpooling(...)
   --get args
   local args = {...}
   local srcdest, segmentation, histmax, minconfidence
   srcdest = args[1]
   segmentation = args[2]
   histmax = args[3]
   minconfidence = args[4]

   -- defaults
   histmax = histmax or false
   minconfidence = minconfidence or 0

   -- usage
   if not srcdest or not segmentation
      or not torch.typename(srcdest):find('Tensor')
      or not torch.typename(segmentation):find('Tensor')
   then
      print(xlua.usage('imgraph.histpooling',
                       'pools the features (or pixels) of an image into a segmentation map,\n'
                          .. 'using histogram accumulation. this is useful for colorazing a\n'
                          .. 'segmentation with the original pixel colors, or for cleaning up\n'
                          .. 'a dense prediction map.\n'
                          .. 'the pooling is done in place (the input is replaced), and two\n'
                          .. 'extra lists are returned, that contain the geometry of all components\n'
                          .. 'found in the process.',
                       nil,
                       {type='torch.Tensor', help='input image/map/matrix to pool (must be KxHxW)', req=true},
                       {type='torch.Tensor', help='segmentation to guide pooling (must be HxW)', req=true},
                       {type='number', help='hist max: replace histograms by their max bin', default=false},
                       {type='number', help='min confidence: vectors with a low confidence are not accumulated', default=0}))
      xlua.error('incorrect arguments', 'imgraph.histpooling')
   end

   -- compute image
   local iresults, hresults = srcdest.imgraph.histpooling(srcdest, segmentation, histmax, minconfidence)

   -- return image
   return srcdest, iresults, hresults
end

----------------------------------------------------------------------
-- a simple test me function
--
function imgraph.testme()
   -- (1) load an image & compute its graph
   local lena = image.lena()
   local lenag = image.convolve(lena, image.gaussian(3), 'full')
   local graph = imgraph.graph(lenag)

   -- (2) compute its connected components, and mst segmentation
   local cc = imgraph.connectcomponents(graph, 0.1, true)
   local mstsegm = imgraph.segmentmst(graph, 3, 20, true)

   -- (3) do a histogram pooling of the original image:
   local mst = imgraph.segmentmst(graph, 3, 20)
   local pool = imgraph.histpooling(lena, mst)

   -- (4) display results
   image.display{image=image.lena(), legend='input image'}
   image.display{image=graph, legend='left: horizontal edges, right: vertical edges'}
   image.display{image=cc, legend='thresholded graph'}
   image.display{image=mstsegm, legend='segmented graph, using min-spanning tree'}
   image.display{image=pool, legend='original imaged hist-pooled by segmentation'}
end
