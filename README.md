Barnes-Hut Algorithm
==========

Before go to detail of the algorithm, it is necessary to understand the structure of a quad-tree. In fact, it is similar to a binary tree, except each node has 4 childrens and can be demonstrated by a 2D plane
![screenshot](https://github.com/ntta/barnes-hut-algorithm/blob/master/quadtreeEx.gif)
The algorithm to compute the force on a particle was first proposed in 1986 by Barnes and Hut:
+ Step 1: Build the quad-tree by create the root node. Then, for each particle, insert it in the quad-tree
+ 
