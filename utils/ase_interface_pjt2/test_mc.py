from __future__ import print_function
import numpy as np
from ase import Atoms
from h2o_pjt2_calculator import H2o_pjt2_calculator
from h2o_pjt2_calculator import get_internal_coords_h2o
from h2o_pjt2_calculator import  get_cart_coords_h2o
from h2o_pjt2_calculator import  get_dq_dr
from h2o_pjt2_calculator import  get_G_matrix
from ase.optimize import BFGS
from ase.optimize.bfgslinesearch import BFGSLineSearch
from ase.vibrations import Vibrations
from mc import Mc, Mc_observer, Fix_centre_of_mass
from ase.constraints import FixAtoms, FixInternals
from mc_utils import Cov_xy, Cov, Histogram, Mc_coord_observer
import scipy 

#r_OH1 = 0.95792059
#r_OH2 = r_OH1
#theta = 104.4996469200000 * np.pi / 180.0
r_OH1 = 0.99792059
r_OH2 = 1.90792059
theta = 114.4996469200000 * np.pi / 180.0

r_o, r_h1, r_h2 = get_cart_coords_h2o(r_OH1, r_OH2, theta, [0.0,0.0, 0.0])


#r_h1 =[r_OH1, 0.0, 0.0]
#r_h2 =[r_OH2 * np.cos(theta), r_OH2 * np.sin(theta), 0.0]
#r_o= [0.0, 0.0, 0.0]

atoms = Atoms('OHH', positions=[r_o, r_h1, r_h2])

# set new positions from MC simulation 1000K
atoms.set_positions([[-17.3806585,   -0.06235601, -29.08992769],
                      [-18.35480593,  -0.35371981, -29.22545189],
                      [-17.24796855,   0.76312459, -29.51363892]])


#    [[-26.02658496, -12.33583177, -21.62541   ],
#                      [-26.30546605, -13.12107032, -21.26245388],
#                      [-25.17565037, -12.09760938, -21.28311626]])
                    


#    [[-17.58224226, -10.32889305,   1.8931163 ],
#                      [-18.05133387, -10.84621147,   2.52740663],
#                      [-18.41062381,  -9.74530678,   1.64014181]])
                    

#    [[ -4.7605108,  -14.09472455,  -0.88578233],
#                      [ -4.24582249, -13.38882193,  -0.48102673],
#                      [ -5.61440173, -14.02165564,  -0.49473533]])

## test G matrix
#dq_dr = get_dq_dr(r_OH1, r_OH2, theta, [0.0,0.0,0.0], h=0.00001)
#
#masses= atoms.get_masses()
#print('masses',masses)
#G = get_G_matrix(r_OH1, r_OH2, theta, masses, h=0.001)
#
#print(G)
#
#eigval, eigvec = np.linalg.eig(G)
#print(eigval)
#print(eigvec)
#asd



#bond1=[0.95792059, [0,1]]
#angle_indices1=[1,0,2]
#angle1=[atoms.get_angle(angle_indices1), angle_indices1]
#constr = FixInternals(bonds=[bond1], angles=[angle1])
#constr = FixAtoms(indices=[0,2])
#print('constr',constr)
#print(dir(FixInternals))
#asd

constr = Fix_centre_of_mass()
atoms.set_constraint(constr)

calc  = H2o_pjt2_calculator()
atoms.set_calculator(calc)

E = atoms.get_potential_energy()

print(E)

# Monte Carlo
T=1000
monte_carlo = Mc(atoms, max_step=0.025 * np.sqrt(T/10), T=T) #, one_d_coord=[1,0])
        
mc_obs1 = Mc_observer(monte_carlo, atoms)
mc_obs2 = Mc_coord_observer(monte_carlo, atoms)

monte_carlo.attach(mc_obs1, interval=1)
monte_carlo.attach(mc_obs2, interval=10)

monte_carlo.run(10000)

print(atoms.get_positions())

print('mc_obs1:')
print('Number of moves {}, number of accepted moves {}, acceptance ratio {}'.format(mc_obs1.n_steps,
                                                                                    mc_obs1.n_acc,
                                                                                    mc_obs1.n_acc/ mc_obs1.n_steps))    
print()

# write histograms
for i in range(3):
    mc_obs2.histograms[i].write('hist_q_' + str(i) + '.txt')

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

#q1_avs =[]
#q2_avs =[]
cov_qq =np.zeros((3,3))
for i in range(3):
    for j in range(3):
        q1_av, q2_av, cov_qq[i,j] =mc_obs2.cov[i][j].comp_av_cov()


    

print('<q1> = {} Ang, cov(q1,q1) = {} Ang**2, sigma(q1,q1) = {} Ang '.format(q1_avs[0], q_covs[0], np.sqrt(q_covs[0]) ))
print('<q2> = {} Ang, cov(q2,q2) = {} Ang**2, sigma(q2,q2) = {} Ang '.format(q1_avs[1], q_covs[1], np.sqrt(q_covs[1])))
print('<theta> = {} deg, cov(theta,theta) = {} deg**2, sigma(theta,theta) = {} deg '.format(q1_avs[2] * 180.0 / np.pi,
                                                                                            q_covs[2]* (180.0 / np.pi)**2,
                                                                                            np.sqrt(q_covs[2])* 180.0 / np.pi ))
from ase.units import Bohr, Ha, kB
amu = 1822.88839

J_to_cm = 5.03411759319722e22
Ha_to_J= 4.35974418e-18
Ha_to_cm = Ha_to_J * J_to_cm
eV_to_cm = Ha_to_cm / Ha 

# detta stmmer ju ej. maste ta mer hansyn...
print('internal 1d')
print(np.sqrt(T*kB / Ha / (q_covs[0] * amu / Bohr**2 )) * Ha_to_cm)
#print(np.sqrt(T*kB / Ha / (q_covs[1] * amu / Bohr**2 )) * Ha_to_cm)

q1_cart, q2_cart, cov_cart = mc_obs2.cov_cart.comp_av_cov()
q1_one_d, q2_one_d, cov_one_d = mc_obs2.cov_one_d.comp_av_cov()
#cov_cart = 0.5*(cov_cart+ np.transpose(cov_cart))

print('cartesian 1d')
print(np.sqrt( T * kB / Ha / (cov_one_d * amu / Bohr**2)) * Ha_to_cm)

print('cartesian 1d from mass scaled covariance')
eigval = cov_cart[3,3]
print(np.sqrt( T * kB / Ha / (eigval * amu / Bohr**2)) * Ha_to_cm)


eigval = np.linalg.eigvalsh(cov_cart)
print('eigval', eigval)
print(np.sqrt( T * kB / Ha / (eigval[3:10] * amu / Bohr**2)) * Ha_to_cm)

# solve problem in internal coordinates
G = mc_obs2.G / mc_obs2.n_steps
cov_qq= 0.5*(cov_qq + np.transpose(cov_qq))
G= 0.5*(G + np.transpose(G))

eigval, eigvec = scipy.linalg.eig(G,cov_qq)

print('eigval, internal', eigval)
print('eigvec', eigvec)
print(np.sqrt( T * kB / Ha *eigval / (amu / Bohr**2)) * Ha_to_cm)


# symetrize the stretch coordinates
tmp = 0.5 * (G[0,0] + G[1,1])
G[0,0] = tmp
G[1,1] = tmp

tmp = 0.5 * (cov_qq[0,0] + cov_qq[1,1])
cov_qq[0,0] = tmp
cov_qq[1,1] = tmp

eigval, eigvec = scipy.linalg.eig(G,cov_qq)

print('Symmetrized eigval, internal', eigval)
print('Symmetrized eigvec', eigvec)
print(np.sqrt( T * kB / Ha *eigval / (amu / Bohr**2)) * Ha_to_cm)

#print(q1_cart)
#print(q2_cart)
#print(cov_cart)
