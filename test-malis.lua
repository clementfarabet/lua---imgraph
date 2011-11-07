
require 'imgraph'

m = nn.MalisCriterion('min','rand')

i = torch.Tensor(8,8):fill(1)
t = torch.Tensor(8,8):zero()

i:narrow(1,5,4):fill(1.8)
i[5][4] = 1.6

print(i)

ig = imgraph.graph(i)

t:narrow(1,1,4):fill(1)
t:narrow(1,5,4):fill(2)

print(ig)
print(t)

loss = m:forward(ig,t)
grad = m:backward(ig,t):clone()

m:forward(ig,t)
grad = grad + m:backward(ig,t):clone()

print('rand index = ' .. loss)
print(grad)
