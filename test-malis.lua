
require 'imgraph'

m = nn.MalisCriterion('min')

i = torch.Tensor(8,8):fill(1)
t = torch.Tensor(8,8):zero()

i:narrow(1,5,4):fill(1.8)
i[5][4] = 1.5

print(i)

ig = imgraph.graph(i)

t:narrow(1,1,4):fill(1)
t:narrow(1,5,4):fill(2)

print(ig)
print(t)

loss,classerr,rand = m:forward(ig,t)
grad = m:backward(ig,t)

print('loss = ' .. loss)
print('class error = ' .. classerr)
print('rand = ' .. rand)
print(grad)
