require 'nn'
local MalisCriterion, parent = torch.class('nn.MalisCriterion', 'nn.Criterion')

function MalisCriterion:__init()
   parent.__init(self)
end

function MalisCriterion:forward(input, target)
   self.output = input.nn.MalisCriterion_forward(input, target)
   return self.output
end

function MalisCriterion:backward(input, target)
   self.gradInput = input.nn.MalisCriterion_backward(input, target)
   return self.gradInput
end
