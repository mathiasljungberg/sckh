from __future__ import print_function
import numpy as np
from ase import Atoms
from mc import Mc, Mc_observer
from h2o_pjt2_calculator import H2o_pjt2_calculator, get_internal_coords_h2o, get_cart_coords_h2o
from h2o_pjt2_calculator import get_G_matrix

# some utility classes

class Cov_xy:
    def __init__(self):
        self.x_sum =0.0
        self.y_sum =0.0
        self.xy_sum =0.0
        self.n = 0.0
    
    def update(self,x,y):
        self.x_sum += x
        self.y_sum += y
        self.xy_sum += x*y
        self.n += 1.0

    def comp_av_cov(self):
        x_av = self.x_sum / self.n
        y_av = self.y_sum / self.n
        cov = (self.xy_sum - self.x_sum * self.y_sum / self.n ) / self.n
        return x_av, y_av, cov


class Cov:
    def __init__(self, dim=1):
        self.x_sum = np.zeros(dim)
        self.y_sum = np.zeros(dim)
        self.xy_sum =np.zeros((dim,dim))
        self.n = 0.0
    
    def update(self,x,y):
        self.x_sum += x
        self.y_sum += y
        self.xy_sum += np.outer(x,y)
        self.n += 1.0

    def comp_av_cov(self):
        x_av = self.x_sum / self.n
        y_av = self.y_sum / self.n
        cov = (self.xy_sum - np.outer(self.x_sum, self.y_sum) / self.n ) / self.n
        return x_av, y_av, cov

    

    
class Histogram:
    
    def __init__(self, nbins=10, limits=[-1.0,1.0]):
        self.nbins = nbins
        

        self.limits = limits
        self.dx = (limits[1]-limits[0]) / (nbins-1.0)
        self.bin_edges = np.zeros(nbins)
        self.hist = np.zeros(nbins)

        # points defined at left of interval
        for i in range(nbins):
            self.bin_edges[i] = self.limits[0] + i * self.dx

            
    def add(self, x, y):
        rx = int((x -self.limits[0]) / self.dx) 

        # force into a bin
        if rx < 0:
            rx =0
        if rx > self.nbins-1:
            rx = self.nbins-1

        self.hist[rx] += y

    def write(self, outfile):
        f = open(outfile, 'w')
        
        for i in range(self.nbins):
            f.write('{} {} \n'.format(self.bin_edges[i], self.hist[i]))
            
        f.close()


    
class Mc_coord_observer(Mc_observer):

    def __init__(self, dyn, atoms):

        Mc_observer.__init__(self, dyn, atoms)
        
        # covariances 
        self.cov = []
        for i in range(3):
            tmp = []
            for j in range(3):
                tmp.append(Cov_xy())
            self.cov.append(tmp)
            
        self.cov_cart = Cov(9)
        self.cov_one_d = Cov_xy()
            
        # histograms 
        self.histograms =[]
        self.histograms.append(Histogram(nbins=100, limits=[0.6, 1.4]))
        self.histograms.append(Histogram(nbins=100, limits=[0.6, 1.4]))
        self.histograms.append(Histogram(nbins=100, limits=[0.0, 6.3]))

        # G matrix
        self.G= np.zeros((3,3))
        
    def __call__(self):

        Mc_observer.__call__(self)
        
        # average internal coordinates
        pos=self.atoms.get_positions() 
        q1, q2, theta = get_internal_coords_h2o(pos[0], pos[1], pos[2])
        #q =[q1, q2, theta *180.0 /np.pi]
        q =[q1, q2, theta]

        for i in range(3):
            for j in range(3):
                self.cov[i][j].update(q[i],q[j])
                
        for i in range(3):        
            self.histograms[i].add(q[i],1.0)

        # mass weighting cartesian coordiates
        #print(atoms.get_masses())
        #print(np.sqrt(np.array(atoms.get_masses())))
        #print(pos)
        masses = self.atoms.get_masses()

        pos_sqrt_mass = []
        for i in range(len(pos)):
            pos_sqrt_mass.append(np.array(pos[i]) * np.sqrt(masses[i]))

        pos_sqrt_mass = np.array(pos_sqrt_mass)
        #print(pos_sqrt_mass)
        pos_sqrt_mass = pos_sqrt_mass.flatten()
        #print(pos_sqrt_mass)
        
        self.cov_one_d.update(pos[1][0],pos[1][0])
        
        self.cov_cart.update(pos_sqrt_mass, pos_sqrt_mass)

        #  G matrix
        self.G += get_G_matrix(q1, q2, theta, masses)

        
