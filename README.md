## Quadtree

Before go to detail of the algorithm, it is necessary to understand the structure of a quad-tree. In fact, it is similar to a binary tree, except each node has 4 childrens and can be demonstrated by a 2D plane (as shown below).

![screenshot](/quadtreeEx.gif)

## Barnes-Hut Algorithm
The algorithm to compute the force on a particle was first proposed in 1986 by Barnes and Hut:
+ Step 1: Build the quad-tree by create the root node. Then, for each particle, insert it in the quad-tree.
+ Step 2: For each sub-square in the quad-tree, compute the total mass it contains and the center of mass.
+ Step 3: Traverse the quad-tree to compute the force on each particle by using approximation
