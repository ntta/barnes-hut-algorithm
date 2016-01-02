from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random
import time

class Node:
    children = None
    mass = None
    center_of_mass = None
    bbox = None

theta = 0.5
G = 6.67408e-11 #m^3 kg^-1 s^-2
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')

def quad_insert(root, x, y, m):
    if root.mass is None:   #when the root is empty, add the first particle
        root.mass = m
        root.center_of_mass = [x,y]
        return
    elif root.children is None:
        root.children = [None,None,None,None]
        old_quadrant = quadrant_of_particle(root.bbox, root.center_of_mass[0], root.center_of_mass[1])
        if root.children[old_quadrant] is None:
            root.children[old_quadrant] = Node()
            root.children[old_quadrant].bbox = quadrants(root.bbox,old_quadrant)
        quad_insert(root.children[old_quadrant], root.center_of_mass[0], root.center_of_mass[1], root.mass)
        new_quadrant = quadrant_of_particle(root.bbox, x, y)
        if root.children[new_quadrant] is None:
            root.children[new_quadrant] = Node()
            root.children[new_quadrant].bbox = quadrants(root.bbox,new_quadrant)
        quad_insert(root.children[new_quadrant], x, y, m)
        root.center_of_mass[0] = (root.center_of_mass[0]*root.mass + x*m) / (root.mass + m)
        root.center_of_mass[1] = (root.center_of_mass[1]*root.mass + y*m) / (root.mass + m)
        root.mass = root.mass + m
    else:
        new_quadrant = quadrant_of_particle(root.bbox, x, y)
        if root.children[new_quadrant] is None:
            root.children[new_quadrant] = Node()
            root.children[new_quadrant].bbox = quadrants(root.bbox, new_quadrant)
        quad_insert(root.children[new_quadrant], x, y, m)
        root.center_of_mass[0] = (root.center_of_mass[0]*root.mass + x*m) / (root.mass + m)
        root.center_of_mass[1] = (root.center_of_mass[1]*root.mass + y*m) / (root.mass + m)
        root.mass = root.mass + m

def display(root):
    if root.mass is None:
        return
    if root.children is not None:
        x = (root.bbox[0] + root.bbox[1]) / 2
        y = (root.bbox[2] + root.bbox[3]) / 2
        width = x-root.bbox[0]
        plt_node(root.bbox[0], root.bbox[2], width)
        plt_node(root.bbox[0], y, width)
        plt_node(x, root.bbox[2], width)
        plt_node(x, y, width)
        for i in xrange(4):
            if root.children[i] is not None:
                display(root.children[i])
    else:
        plt.scatter(root.center_of_mass[0], root.center_of_mass[1])

def compute_force(root,x,y,m):
    f = 0
    r = distance(x, y, root.center_of_mass[0], root.center_of_mass[1])
    d = root.bbox[1]-root.bbox[0]
    if root.mass is None or r == 0:
        return 0
    elif d/r <= theta or root.children is None:
        f = f + force(x, y, m, root.center_of_mass[0], root.center_of_mass[1], root.mass)
    else:
        for i in xrange(4):
            if root.children[i] is not None:
                f = f + compute_force(root.children[i], x, y, m)
    return f

################################################# SUPPORTING FUNCTION ##############################################################

def force(x, y, m, xcm, ycm, mcm):
    return (G*m*mcm)/(distance(x,y,xcm,ycm)**2)

def distance(x1, y1, x2, y2):
    return math.sqrt((x2-x1)**2+(y2-y1)**2)

def plt_node(x, y, width):
    ax.add_patch(patches.Rectangle((x, y), width, width, fill = False))

def find_root_bbox(array):
    """ Create a suitable square boundary box for the input particles
    """
    if len(array) == 0 or len(array) == 1:
        return None
    bmin = np.min(array, axis = 0)
    bmax = np.max(array, axis = 0)
    xmin, xmax, ymin, ymax = bmin[0], bmax[0], bmin[1], bmax[1]
    if xmax - xmin == ymax - ymin:
        return xmin, xmax, ymin, ymax
    elif xmax - xmin > ymax - ymin:
        return xmin, xmax, ymin, ymax+(xmax-xmin-ymax+ymin)
    else:
        return xmin, xmax+(ymax-ymin-xmax+xmin), ymin, ymax

def quadrant_of_particle(bbox, x, y):
    """Return position of quadrant of the particle (x,y)
    """
    if y >= (bbox[3] + bbox[2])/2:
        if x <= (bbox[1] + bbox[0])/2:
            return 0
        else:
            return 1
    else:
        if x >= (bbox[1] + bbox[0])/2:
            return 2
        else:
            return 3

def quadrants(bbox,quadrant):
    """Return the coordinate of the quadrant
    """
    x = (bbox[0] + bbox[1]) / 2
    y = (bbox[2] + bbox[3]) / 2
    #Quadrant 0: (xmin, x, y, ymax)
    if quadrant == 0:
        return bbox[0], x, y, bbox[3]
    #Quadrant 1: (x, xmax, y, ymax)
    elif quadrant == 1:
        return x, bbox[1], y, bbox[3]
    #Quadrant 2: (x, xmax, ymin, y)
    elif quadrant == 2:
        return x, bbox[1], bbox[2], y
    #Quadrant 3: (xmin, x, ymin, y)
    elif quadrant == 3:
        return bbox[0], x, bbox[2], y

def data_from_file(filename, array):
    with open(filename) as f:
        for line in f:
            if line[0] == '#':
                continue
            else:
                x,y,m = line.split(',')
                array.append([float(x),float(y),float(m)])

def data_from_random(quantity, array):
    for i in xrange(quantity):
        x = random.randrange(0, 10000)
        y = random.randrange(0, 10000)
        m = random.randrange(100, 10000)
        array.append([x, y, m])

if __name__ == '__main__':
    filename = ('testcase.txt')
    particles = []
    #data_from_random(10000, particles)
    data_from_file(filename, particles)
    root = Node()
    root.center_of_mass = []
    root.bbox = find_root_bbox(particles)
    for i in xrange(len(particles)):
        quad_insert(root, particles[i][0], particles[i][1], particles[i][2])
    print 'Boundary box: ',root.bbox
    print 'Total mass: ',root.mass
    print 'Coordinate of center of mass: ',root.center_of_mass
    plt.scatter(root.center_of_mass[0], root.center_of_mass[1], c='r', marker='x', s=50)
    display(root)
    print 'Theta: ', theta
    print 'Force interact on particle [0,0,5.972e20]:', compute_force(root, 0, 0, 5.972e20)
    plt.show()
