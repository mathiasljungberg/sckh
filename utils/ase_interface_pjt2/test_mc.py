from __future__ import print_function
import numpy as np
from ase import Atoms
from h2o_pjt2_calculator import H2o_pjt2_calculator, get_internal_coords_h2o, get_cart_coords_h2o
from ase.optimize import BFGS
from ase.optimize.bfgslinesearch import BFGSLineSearch
from ase.vibrations import Vibrations
from mc import Mc, Mc_observer

r_OH1 = 0.95792059
r_OH2 = r_OH1
theta = 104.4996469200000 * np.pi / 180.0
#r_OH1 = 0.99792059
#r_OH2 = 1.90792059
#theta = 114.4996469200000 * np.pi / 180.0

r_o, r_h1, r_h2 = get_cart_coords_h2o(r_OH1, r_OH2, theta, [0.0,0.0, 0.0])

atoms = Atoms('OHH', positions=[r_o, r_h1, r_h2])

calc  = H2o_pjt2_calculator()
atoms.set_calculator(calc)

E = atoms.get_potential_energy()

print(E)

# Monte Carlo
monte_carlo = Mc(atoms, max_step=0.15, T=298)


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


class Mc_coord_observer(Mc_observer):

    def __init__(self, dyn, atoms):

        Mc_observer.__init__(self, dyn, atoms)
        
        #self.q1_sum = 0.0
        #self.q2_sum = 0.0
        #self.theta_sum = 0.0

        # test to do things with covariances instead
        self.cov = []
        for i in range(3):
            tmp = []
            for j in range(3):
                tmp.append(Cov_xy())
            self.cov.append(tmp)

        #print(self.cov)
            
    def __call__(self):

        Mc_observer.__call__(self)
        
        # average internal coordinates
        pos=atoms.get_positions() 
        q1, q2, theta = get_internal_coords_h2o(pos[0], pos[1], pos[2])
        q =[q1, q2, theta]
        
        #self.q1_sum += q1
        #self.q2_sum += q2
        #self.theta_sum += theta * 180.0/ np.pi

        for i in range(3):
            for j in range(3):
                self.cov[i][j].update(q[i],q[j])
                
        
        
mc_obs1 = Mc_observer(monte_carlo, atoms)
mc_obs2 = Mc_coord_observer(monte_carlo, atoms)

monte_carlo.attach(mc_obs1, interval=1)
monte_carlo.attach(mc_obs2, interval=100)

monte_carlo.run(5000)

print('mc_obs1:')
print('Number of moves {}, number of accepted moves {}, acceptance ratio {}'.format(mc_obs1.n_steps,
                                                                                    mc_obs1.n_acc,
                                                                                    mc_obs1.n_acc/ mc_obs1.n_steps))    
print()
print('mc_obs2:')
print('Number of moves {}, number of accepted moves {}, acceptance ratio {}'.format(mc_obs2.n_steps,
                                                                                    mc_obs2.n_acc,
                                                                                    mc_obs2.n_acc/ mc_obs2.n_steps))    
#print('<q1> = {} Ang, <q2> = {} Ang, <theta> ={} Ang '.format(mc_obs2.q1_sum / mc_obs2.n_steps,
#                                                              mc_obs2.q2_sum / mc_obs2.n_steps,
#                                                              mc_obs2.theta_sum / mc_obs2.n_steps))

q1_avs =[]
q2_avs =[]
q_covs =[]
for i in range(3):
    #for j in range(3):
    q1_av, q2_av, q1_q2_cov =mc_obs2.cov[i][i].comp_av_cov()
    q1_avs.append(q1_av)
    q2_avs.append(q2_av)
    q_covs.append(q1_q2_cov)
    

print('<q1> = {} Ang, cov(q1,q1) = {} Ang**2, sigma(q1,q1) = {} Ang '.format(q1_avs[0], q_covs[0], np.sqrt(q_covs[0]) ))
print('<q2> = {} Ang, cov(q2,q2) = {} Ang**2, sigma(q2,q2) = {} Ang '.format(q1_avs[1], q_covs[1], np.sqrt(q_covs[1])))
print('<theta> = {} deg, cov(theta,theta) = {} deg**2, sigma(theta,theta) = {} deg '.format(q1_avs[2]* 180.0 / np.pi,
                                                                                            q_covs[2]* (180.0 / np.pi)**2,
                                                                                            np.sqrt(q_covs[2]) * 180.0 / np.pi ))

