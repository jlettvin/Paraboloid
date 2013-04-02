#!/usr/bin/env python

###############################################################################
"""
Paraboloid.py

A utility for constructing truncated paraboloid shells in discrete meshes.
Copyright(c) 2013 Jonathan D. Lettvin, All Rights Reserved"

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__module__     = "Paraboloid.py"
__author__     = "Jonathan D. Lettvin"
__copyright__  = """\
Copyright(C) 2011 Jonathan D. Lettvin, All Rights Reserved"""
__credits__    = [ "Jonathan D. Lettvin" ]
__license__    = "GPLv3"
__version__    = "0.0.3"
__maintainer__ = "Jonathan D. Lettvin"
__email__      = "jlettvin@gmail.com"
__contact__    = "jlettvin@gmail.com"
__status__     = "Demonstration"
__date__       = "20111110"

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
