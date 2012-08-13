require 'nn'
local MalisCriterion, parent = torch.class('nn.MalisCriterion', 'nn.Criterion')

function MalisCriterion:__init(tree, metric, connex, threshold, margin, alternateposneg)
   parent.__init(self)
   connex = connex or 4
   -- 2d connectivities
   if connex == 4 then
      self.twoD = true
      self.nhood = torch.Tensor(2, 3)
      self.nhood[1][1] = 0  -- edge1: z offset
      self.nhood[1][2] = 0  -- edge1: y offset
      self.nhood[1][3] = 1  -- edge1: x offset
      self.nhood[2][1] = 0  -- edge1: z offset
      self.nhood[2][2] = 1  -- edge2: y offset
      self.nhood[2][3] = 0  -- edge2: x offset
   elseif connex == 8 then
      self.twoD = true
      self.nhood = torch.Tensor(4, 3)
      self.nhood[1][1] = 0  -- edge1: z offset
      self.nhood[1][2] = 0  -- edge1: y offset
      self.nhood[1][3] = 1  -- edge1: x offset
      self.nhood[2][1] = 0  -- edge1: z offset
      self.nhood[2][2] = 1  -- edge2: y offset
      self.nhood[2][3] = 0  -- edge2: x offset
      self.nhood[3][1] = 0  -- edge1: z offset
      self.nhood[3][2] = 1  -- edge3: y offset
      self.nhood[3][3] = 1  -- edge3: x offset
      self.nhood[4][1] = 0  -- edge1: z offset
      self.nhood[4][2] = 1  -- edge4: y offset
      self.nhood[4][3] = -1 -- edge4: x offset
   -- 3d connectivities
   elseif connex == 6 then
      self.twoD = false
      self.nhood = torch.Tensor(3, 3)
      self.nhood[1][1] = 0  -- edge1: z offset
      self.nhood[1][2] = 0  -- edge1: y offset
      self.nhood[1][3] = -1 -- edge1: x offset
      self.nhood[2][1] = 0  -- edge1: z offset
      self.nhood[2][2] = -1 -- edge2: y offset
      self.nhood[2][3] = 0  -- edge2: x offset
      self.nhood[3][1] = -1 -- edge1: z offset
      self.nhood[3][2] = 0  -- edge2: y offset
      self.nhood[3][3] = 0  -- edge2: x offset
   else
      error('connexity must be one of 4 | 8 for 2d images and 6 for 3d images')
   end
   self.threshold = threshold or 0.5
   self.margin = margin or 0.3
   self.posexample = true
   self.negexample = true
   self.alternateposneg = alternateposneg or false
   if self.alternateposneg then
      self.negexample = false -- negative example first
   end
   self.metric = metric or 'loss' -- loss | classified | rand
   self.tree = tree or 'max'
   if self.tree == 'min' then
      self.threshold = -self.threshold
   end
   self.inputtemp = torch.Tensor()
   self.targettemp = torch.Tensor()
end

function MalisCriterion:forward(input, target, posexample, negexample)
   self._i = input
   self._t = target

   -- pos? neg?
   if self.alternateposneg then
      self.posexample = not self.posexample
      self.negexample = not self.negexample
   end
   if (posexample ~= nil) then
      self.posexample = posexample
   end
   if (negexample ~= nil) then
      self.negexample = negexample
   end

   -- ensure contiguity by copying (contiguous() seems not to work properly...? or my scoping is bad...?)
   self.inputtemp:resizeAs(input):copy(input)
   self.targettemp:resizeAs(target):copy(target)

   -- shapes and signs
   if self.tree == 'min' then
      self.inputtemp:mul(-1)
   end
   if self.twoD then -- convert to 3d by inserting a dimension
      self.inputtemp:resize(input:size(1),1,input:size(2),input:size(3))
      self.targettemp:resize(1,target:size(1),target:size(2))
   end

   -- run the MST and compute the loss and its derivative
   local e = {}
   e.loss,e.classified,e.rand = input.nn.MalisCriterion_forward(self, self.inputtemp, self.targettemp)

   -- shapes and signs
   if self.tree == 'min' then
      self.gradInput:mul(-1)
   end
   if self.twoD then
      self.gradInput:resize(input:size())
   end

   -- return the goodness score
   self.output = e[self.metric]
   return self.output
end

function MalisCriterion:backward(input, target)
   if self._i ~= input or self._t ~= target then
      error('forward must be called once before backward()')
   end
   return self.gradInput
end
