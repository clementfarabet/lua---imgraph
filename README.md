# imgraph: a package to create/manipulate graphs on images

This package provides standard functions to
create and manipulate edge-weighted graphs 
of images: create a graph, segment it, 
compute its watershed, or its connected
components...

## Install dependencies 

On Linux (Ubuntu > 9.04):

``` sh
$ apt-get install gcc g++ git libreadline5-dev cmake wget libqt4-core libqt4-gui libqt4-dev
```

On Mac OS (Leopard, or more), using [Homebrew](http://mxcl.github.com/homebrew/):

``` sh
$ brew install git readline cmake wget qt
```

## Install Lua, Luarocks and Torch7

1/ Lua 5.1 + Luarocks + xLua 

``` sh
$ git clone https://github.com/clementfarabet/lua4torch
$ cd lua4torch
$ make install PREFIX=/usr/local
```

2/ Torch7 (a numeric package for Lua)

``` sh
$ luarocks install torch
```

## Install imgraph

``` sh
$ luarocks install imgraph
```
