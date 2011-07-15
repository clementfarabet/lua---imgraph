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
--     July 14, 2011, 5:52PM  - colorizing function  - Clement Farabet
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
-- computes a graph from an 2D or 3D image
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
                       'compute an edge-weighted graph on an image', nil,
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
   img.imgraph.graph(dest, img, connex, distance)

   -- return result
   return dest
end

----------------------------------------------------------------------
-- compute the connected components of a graph
--
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
                       'compute the connected components of an edge-weighted graph', nil,
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

----------------------------------------------------------------------
-- compute the gradient amplitude of a graph
--
function imgraph.gradient(...)
   -- get args
   local args = {...}
   local gradient, graph
   local arg2 = torch.typename(args[2])
   if arg2 and arg2:find('Tensor') then
      gradient = args[1]
      graph = args[2]
   else
      gradient = torch.Tensor()
      graph = args[1]
   end

   -- usage
   if not graph then
      print(xlua.usage('imgraph.gradient',
                       'compute an approximated gradient map of a graph\n'
                          .. 'the map has the size of the original image on which the graph\n'
                          .. 'was constructed',
                       nil,
                       {type='torch.Tensor', help='input graph', req=true},
                       "",
                       {type='torch.Tensor', help='destination tensor', req=true},
                       {type='torch.Tensor', help='input graph', req=true}))
      xlua.error('incorrect arguments', 'imgraph.gradient')
   end

   -- compute graph
   graph.imgraph.gradient(gradient, graph)

   -- return result
   return gradient
end

----------------------------------------------------------------------
-- compute the watershed of a graph
--
function imgraph.watershed(...)
   --get args
   local args = {...}
   local dest, gradient, minHeight, connex, colorize
   if arg2 and arg2:find('Tensor') then
      dest = args[1]
      gradient = args[2]
      minHeight = args[3]
      connex = args[4]
      colorize = args[5]
   else
      dest = torch.Tensor()
      gradient = args[1]
      minHeight = args[2]
      connex = args[3]
      colorize = args[4]
   end

   -- defaults
   minHeight = minHeight or 0.05
   connex = connex or 4

   -- usage
   if not gradient or (gradient:nDimension() ~= 2) then
      print(xlua.usage('imgraph.watershed',
                       'compute the watershed of a gradient map (or arbitrary grayscale image)', nil,
                       {type='torch.Tensor', help='input gradient map (HxW tensor)', req=true},
                       {type='number', help='filter minimas by imposing a minimum height', default=0.05},
                       {type='number', help='connexity: 4 | 8', default=4},
                       {type='boolean', help='replace components id by random colors', default=false},
                       "",
                       {type='torch.Tensor', help='destination tensor', req=true},
                       {type='torch.Tensor', help='input gradient map (HxW tensor)', req=true},
                       {type='number', help='filter minimas by imposing a minimum height', default=0.05},
                       {type='number', help='connexity: 4 | 8', default=4},
                       {type='boolean', help='replace components id by random colors', default=false}))
      xlua.error('incorrect arguments', 'imgraph.watershed')
   end

   -- compute image
   local nelts = gradient.imgraph.watershed(dest, gradient, minHeight, connex, colorize)

   -- return image
   return dest, nelts
end

----------------------------------------------------------------------
-- segment a graph, by computing its min-spanning tree and
-- merging vertices based on a dynamic threshold
--
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

   -- defaults
   thres = thres or 3
   minsize = minsize or 20
   colorize = colorize or false

   -- usage
   if not graph then
      print(xlua.usage('imgraph.segmentmst',
                       'segment an edge-weighted graph, using a surface adaptive criterion\n'
                          .. 'on the min-spanning tree of the graph (see Felzenszwalb et al. 2004)',
                       nil,
                       {type='torch.Tensor', help='input graph', req=true},
                       {type='number', help='base threshold for merging', default=3},
                       {type='number', help='min size: merge components of smaller size', default=20},
                       {type='boolean', help='replace components id by random colors', default=false},
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

----------------------------------------------------------------------
-- pool the features (or pixels) of an image into a segmentation map
--
function imgraph.histpooling(...)
   --get args
   local args = {...}
   local srcdest, segmentation, lists, histmax, minconfidence
   srcdest = args[1]
   segmentation = args[2]
   histmax = args[3]
   minconfidence = args[4]

   -- defaults
   histmax = histmax or false
   minconfidence = minconfidence or 0

   -- usage
   if not srcdest or not segmentation or not torch.typename(srcdest):find('Tensor') or not torch.typename(segmentation):find('Tensor') then
      print(xlua.usage('imgraph.histpooling',
                       'pool the features (or pixels) of an image into a segmentation map,\n'
                          .. 'using histogram accumulation. this is useful for colorazing a\n'
                          .. 'segmentation with the original pixel colors, or for cleaning up\n'
                          .. 'a dense prediction map.\n\n'
                          .. 'the pooling is done in place (the input is replaced)\n\n'
                          .. 'two extra lists of components are optionally generated:\n'
                          .. 'the first list is an array of these components, \n'
                          .. 'while the second is a hash table; each entry has this format:\n'
                          .. 'entry = {centroid_x, centroid_y, surface, hist_max, id}',
                       nil,
                       {type='torch.Tensor', help='input image/map/matrix to pool (must be KxHxW)', req=true},
                       {type='torch.Tensor', help='segmentation to guide pooling (must be HxW)', req=true},
                       {type='number', help='hist max: replace histograms by their max bin', default=false},
                       {type='number', help='min confidence: vectors with a low confidence are not accumulated', default=0}))
      xlua.error('incorrect arguments', 'imgraph.histpooling')
   end

   -- compute image
   local iresults, hresults = srcdest.imgraph.histpooling(srcdest, segmentation, 
                                                          true, histmax, minconfidence)

   -- return image
   return srcdest, iresults, hresults
end

----------------------------------------------------------------------
-- colorize a segmentation map
--
function imgraph.colorize(...)
   -- get args
   local args = {...}
   local segmentation, colorized
   if select('#',...) == 2 then
      colorized = args[1]
      segmentation = args[2]
   else
      colorized = torch.Tensor()
      segmentation = args[1]
   end

   -- usage
   if not segmentation then
      print(xlua.usage('imgraph.colorize',
                       'colorize a segmentation map',
                       'graph = imgraph.graph(image.lena())\n'
                          .. 'segm = imgraph.segmentmst(graph)\n'
                          .. 'colored = imgraph.colorize(segm)',
                       {type='torch.Tensor', help='input segmentation map (must be HxW)', req=true}))
      xlua.error('incorrect arguments', 'imgraph.colorize')
   end

   -- colorize !
   segmentation.imgraph.colorize(colorized, segmentation)

   -- return colorized segmentation
   return colorized
end

----------------------------------------------------------------------
-- a simple test me function
--
imgraph._example = [[
      -- (1) load an image & compute its graph
      local lena = image.lena()
      local lenag = image.convolve(lena, image.gaussian(3), 'full')
      local graph = imgraph.graph(lenag)

      -- (2) compute its connected components, and mst segmentation
      local cc = imgraph.connectcomponents(graph, 0.1, true)
      local mstsegm = imgraph.segmentmst(graph, 3, 20)
      local mstsegmcolor = imgraph.colorize(mstsegm)

      -- (3) do a histogram pooling of the original image:
      local pool = imgraph.histpooling(lena, mstsegm)

      -- (4) display results
      image.display{image=image.lena(), legend='input image'}
      image.display{image=graph, legend='left: horizontal edges, right: vertical edges'}
      image.display{image=cc, legend='thresholded graph'}
      image.display{image=mstsegmcolor, legend='segmented graph, using min-spanning tree'}
      image.display{image=pool, legend='original imaged hist-pooled by segmentation'}
]]
function imgraph.testme()
   local example = loadstring(imgraph._example)
   print 'imgraph sample code {\n'
   print (imgraph._example)
   print '}'
   example()
end
