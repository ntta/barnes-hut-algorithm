from __future__ import division
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class Node:
    children = None
    mass = None
    center_of_mass = None
    bbox = None

G = 6.67408e-11 #m^3 kg^-1 s^-2
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')

def quad_insert(root, x, y, m):
    new_particle = Node()
    new_particle.mass = m
    new_particle.center_of_mass = [x,y]
    if root.mass is None:   #when the root is empty, add the first particle
        root.mass = m
        root.center_of_mass = [x,y]
    elif root.children is None:
        root.children = [None,None,None,None]
        #add the old particle in a children
        old_quadrant = quadrant_of_particle(root.bbox, root.center_of_mass[0], root.center_of_mass[1])
        old_particle = Node()
        old_particle.mass = root.mass
        old_particle.center_of_mass = [root.center_of_mass[0], root.center_of_mass[1]]
        old_particle.bbox = quadrants(root.bbox,old_quadrant)
        root.children[old_quadrant] = old_particle
        #Put the new particle in the children which is the same level with the old one
        new_quadrant = quadrant_of_particle(root.bbox, x, y)
        if root.children[new_quadrant] is None:
            new_particle.bbox = quadrants(root.bbox,new_quadrant)
            root.children[new_quadrant] = new_particle
        else:
            quad_insert(root.children[new_quadrant], x, y, m)
        root.mass = root.mass + new_particle.mass
        root.center_of_mass[0] = (old_particle.center_of_mass[0]*old_particle.mass + new_particle.center_of_mass[0]*new_particle.mass) / root.mass
        root.center_of_mass[1] = (old_particle.center_of_mass[1]*old_particle.mass + new_particle.center_of_mass[1]*new_particle.mass) / root.mass
    else:
        new_quadrant = quadrant_of_particle(root.bbox, x, y)
        if root.children[new_quadrant] is None:
            new_particle.bbox = quadrants(root.bbox, new_quadrant)
            root.children[new_quadrant] = new_particle
        else:
            quad_insert(root.children[new_quadrant], x, y, m)
        root.center_of_mass[0] = (root.center_of_mass[0] * root.mass + new_particle.center_of_mass[0] * new_particle.mass) / (root.mass + new_particle.mass)
        root.center_of_mass[1] = (root.center_of_mass[1] * root.mass + new_particle.center_of_mass[1] * new_particle.mass) / (root.mass + new_particle.mass)
        root.mass = root.mass + new_particle.mass

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

################################################# SUPPORTING FUNCTION ##############################################################

def force(x, y, m, xcm, ycm, mcm):
    r = math.sqrt((xcm-x)**2+(ycm-y)**2)
    xf = G*m*mcm*(xcm-x)/r**3
    yf = G*m*mcm*(ycm-y)/r**3
    return [xf, yf]

def plt_node(x, y, width):
    ax.add_patch(patches.Rectangle((x, y), width, width, fill = False))

def find_root_bbox(array):
    """ Create a suitable square boundary box for the input particles
    """
    xmin = array[0][0]
    xmax = array[0][0]
    ymin = array[0][1]
    ymax = array[0][1]
    if len(array) == 0 or len(array) == 1:
        return None
    for i in xrange(len(array)):
        if array[i][0] < xmin:
            xmin = array[i][0]
        if array[i][0] > xmax:
            xmax = array[i][0]
        if array[i][1] < ymin:
            ymin = array[i][1]
        if array[i][1] > ymax:
            ymax = array[i][1]
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

if __name__ == '__main__':
    filename = ('testcase.txt')
    particles = []
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
    plt.show()
