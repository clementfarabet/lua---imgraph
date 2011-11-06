require 'nn'
local MalisCriterion, parent = torch.class('nn.MalisCriterion', 'nn.Criterion')

function MalisCriterion:__init(tree, connex, margin)
   parent.__init(self)
   connex = connex or 4
   if connex == 4 then
      self.nhood = torch.Tensor(2, 2)
      self.nhood[1][1] = 0  -- edge1: y offset
      self.nhood[1][2] = 1  -- edge1: x offset
      self.nhood[2][1] = 1  -- edge2: y offset
      self.nhood[2][2] = 0  -- edge2: x offset
   elseif connex == 8 then
      self.nhood = torch.Tensor(4, 2)
      self.nhood[1][1] = 0  -- edge1: y offset
      self.nhood[1][2] = 1  -- edge1: x offset
      self.nhood[2][1] = 1  -- edge2: y offset
      self.nhood[2][2] = 0  -- edge2: x offset
      self.nhood[3][1] = 1  -- edge3: y offset
      self.nhood[3][2] = 1  -- edge3: x offset
      self.nhood[4][1] = 1  -- edge4: y offset
      self.nhood[4][2] = -1 -- edge4: x offset
   else
      error('connexity must be one of: 4 | 8')
   end
   self.margin = margin or 0.3
   self.posexample = false
   self.tree = tree or 'max'
end

function MalisCriterion:forward(input, target)
   local inputt = input
   if self.tree == 'min' then
      inputt = -input+1
   end
   local loss,classerr,rand = input.nn.MalisCriterion_forward(self, inputt, target)
   self.posexample = not self.posexample
   self.output = loss
   self._i = input
   self._t = target
   return loss,classerr,rand
end

function MalisCriterion:backward(input, target)
   if self._i ~= input or self._t ~= target then
      input.nn.MalisCriterion_forward(self, input, target)
      self._i = input
      self._t = target
   end
   return self.gradInput
end
