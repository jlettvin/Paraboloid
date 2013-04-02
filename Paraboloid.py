#!/usr/bin/env python
"""Paraboloid.py"""

###############################################################################
__date__       = "20130131"
__author__     = "jlettvin"
__maintainer__ = "jlettvin"
__email__      = "jlettvin@gmail.com"
__copyright__  = "Copyright(c) 2013 Jonathan D. Lettvin, All Rights Reserved"
__license__    = "Trade Secret"
__status__     = "Production"
__version__    = "0.0.1"

###############################################################################
from scipy      import array, sqrt, arange, rot90, fabs, where, around
from scipy      import set_printoptions
from itertools  import product
from optparse   import OptionParser
from subprocess import call

import pylab as p
#import mpl_toolkits.mplot3d.axes3d as p3

set_printoptions(precision=2, linewidth=500)

class Paraboloid(object):
    def __init__(self, **kw):
        self.kw = kw
        self.E  = kw.get('E', 12.0) # Edge halflength of the containing cube.
        self.dE = kw.get('dE', 0.2) # Edge increment length.
        self.dR = kw.get('dR', 0.1) # Radial increment for coincidence.
        #print self.E, self.dE, self.dR
        rng     = list(arange(-self.E,self.E+self.dE, self.dE)) # Edge values.
        # Force rebalance
        rng    -= (rng[-1] + rng[0])/2.0
        self.X  = array([rng,] * len(rng))      # x values plane segment
        self.Y  = rot90(self.X)                 # y values plane segment
        self.r  = sqrt(self.X**2 + self.Y**2)   # r values plane segment
        #print self.kw
        #print self.X
        #print self.Y
        #print self.r
        # persistent information about nodes, rings of nodes, and edges between
        self.ring   = {}            # radius: ring dictionary (duplicates)
        self.unique = {}            # radius: ring dictionary (unique)
        self.edges  = []            # pairs of vertices
        self.point  = {}            # unique points in paraboloid
        self.keyed  = {}            # keyed access to rings of unique points
        last        = None
        # Build a ring stack comprising a discrete paraboloid of revolution.
        z           = -self.dE
        zdiv        = 30
        self.zmax   = 0
        while True:
            z          += self.dE
            self.zmax   = z / zdiv
            r = around(sqrt(z / self.E), decimals=1)
            if r > self.E: break
            Q = fabs(self.r-r)
            this = [(x,y,z)
                    for x,y,c in zip(
                        around(self.X[Q <= self.dR], decimals=2).flatten(),
                        around(self.Y[Q <= self.dR], decimals=2).flatten(),
                        around(     Q[Q <= self.dR], decimals=2).flatten())
                    if c]
            if last == None or this != last:
                if this != last:
                    # connect vertices from last to this
                    pass
                last = this

                self.ring[r] = this
                self.unique[r] = this
                self.keyed[r] = []

                for i,j,k in this:
                    xykey = '%+3.2e %+3.2e' % (i,j)
                    # Keep a dictionary of unique points.
                    self.point[xykey] = (i,j)
                    n = sqrt(i*i+j*j)
                    if not self.point.has_key(xykey):
                        self.keyed[r] += [
                                {'key':xykey, 'xy':(i,j), 'tip':(i/n,j/n) },]

                if kw.get('verbose', False):
                    print r, array(self.ring[r])
                    print where(Q <= self.dR, '*', ' ')
            else:
                # Don't eliminate duplicates.
                self.ring[r] = last

        filename = 'graph/E%3.1fe%3.1fr%3.2f' % (self.E,self.dE, self.dR)
        #terminal = "png nocrop enhanced font verdana 12 size 640,480"
        terminal = "png nocrop size 640,480"

        print filename, 'to', filename+'.png'
        with open(filename+'.obj', 'w') as objfile:
            with open(filename+'.dat', 'w') as datfile:
                print>>datfile, "# %s (this file)" % (filename)
                print>>datfile, "# %f %f %f" % (self.E, self.dE, self.dR)
                print>>datfile, "# %d points" % (len(self.point))
                print>>datfile, "# %d z values" % (len(self.ring.keys()))
                print>>datfile, "# %d (x,y) sets" % (len(self.unique.keys()))
                print>>datfile, "set terminal %s" % (terminal)
                print>>datfile, "set output '%s.png'" % (filename)
                print>>datfile, "set zrange[0:%d]" % (self.zmax+1)
                print>>datfile, "set isosamples 100"
                print>>datfile, "splot '-' with vectors title 'Paraboloid'"
                print>>datfile
                for r, layer in self.unique.iteritems():
                    for i,j,kl in layer:
                        n = sqrt(i*i+j*j)
                        n += float(n==0.0)
                        k = float(kl)/zdiv
                        print>>datfile, '%+3.2e '*6 % (i, j, k, i/n, j/n, 0)
                        print>>objfile, ('v'+' %+3.2e'*3) % (i, j, k)
        call(['gnuplot', filename+'.dat'])

    def save(self, filename):
        pass

if __name__ == "__main__":

    parser = OptionParser()
    parser.add_option( "--E", default=11.70, type=float, help="max edge length")
    parser.add_option("--dE", default= 1.00, type=float, help="edge increment")
    parser.add_option("--dR", default= 0.05, type=float, help="radial increment")

    parser.add_option("-p", "--plain",
            action="store_true", default=False,
            help="make gnuplot file plain, not fancy")

    parser.add_option( "-v", "--verbose",
            action="store_true", default=False,
            help="announce actions and sizes")
    (opts, args) = parser.parse_args()
    kw = vars(opts)

    P = Paraboloid(**kw)
    #print len(P.point), 'points'
    #print len(P.ring.keys()), 'z values'
    #print len(P.unique.keys()), '(x,y) sets'
    #print P.ring