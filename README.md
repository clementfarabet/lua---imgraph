# imgraph: a package to create/manipulate graphs on images

This package provides standard functions to
create and manipulate edge-weighted graphs 
of images: create a graph, segment it, 
compute its watershed, or its connected
components...

## Install dependencies 

1/ third-party libraries:

On Linux (Ubuntu > 9.04):

``` sh
$ apt-get install gcc g++ git libreadline5-dev cmake wget libqt4-core libqt4-gui libqt4-dev
```

On Mac OS (Leopard, or more), using [Homebrew](http://mxcl.github.com/homebrew/):

``` sh
$ brew install git readline cmake wget qt
```

2/ Lua 5.1 + Luarocks + xLua:

``` sh
$ git clone https://github.com/clementfarabet/lua4torch
$ cd lua4torch
$ make install PREFIX=/usr/local
```

3/ imgraph:

``` sh
$ luarocks install imgraph
```

## Use the library

First run xlua, and load imgraph:

``` sh
$ xlua
``` 

``` lua
> require 'imgraph'
```

Once loaded, tab-completion will help you navigate through the
library:

``` lua
> imgraph. + TAB
imgraph.colorize(           imgraph.connectcomponents(  
imgraph.graph(              imgraph.histpooling(        
imgraph.segmentmst(         imgraph.testme(             
imgraph.watershed(          imgraph.gradient(
```

To get quickly started, run the testme() function:

``` lua
> imgraph.testme()
```

which computes a few things on the famous image of Lena:

![results](http://data.neuflow.org/share/imgraph-testme.png)
