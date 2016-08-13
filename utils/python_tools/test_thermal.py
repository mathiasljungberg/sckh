from __future__ import print_function
import eigenfunctions
import eigenfunctions_fortran
import numpy as np
import ase.units as units 
import sample_funcs

amu = 1822.888486192
T_au = 3.15777504e5

n=0
#omega=44.1e-3 / units.Ha # 443.1 meV = 3574.0 cm^{-1}
omega=443.1e-3 / units.Ha # 443.1 meV = 3574.0 cm^{-1}
hartree_to_cm = eigenfunctions_fortran.parameters.cm * eigenfunctions_fortran.parameters.hartree # this is correct 
omega_cm= omega * hartree_to_cm
print('omega', omega)
print('omega_cm', omega_cm)
mu=1.0

# ground state gaussian exponent
gauss_exp = omega*mu*amu / 2.0
gauss2_exp = 2.0 * gauss_exp
fwhm = 2.0 * np.sqrt(np.log(2)/gauss2_exp)
print('ground state squared gaussian exponent', gauss2_exp, 'fwhm', fwhm)
gauss_exp_fft = 1.0/(4.0*gauss_exp)
gauss2_exp_fft = 2.0 * gauss_exp_fft
fwhm_fft = 2.0 * np.sqrt(np.log(2)/gauss2_exp_fft)
print('ground state squared FFT gaussian exponent', gauss2_exp_fft, 'fwhm', fwhm_fft)

n_x=1000
x=np.linspace(-10.0,10.0, n_x)
dx=x[1]-x[0]
l_x = x[n_x-1]-x[0]
dx_Ang=dx * units.Bohr 

n_f = 20
fun = np.zeros((n_x, n_f))
fun_fft = np.zeros((n_x, n_f), dtype=np.complex)

# create eigenfunctions
for i in range(n_x):
    for j in range(n_f):
        fun[i,j] = eigenfunctions.harm_eigenfun(j, omega, mu, x[i] , units='Bohr' )


# get the momenta by fft
for j in range(n_f):
    fun_fft[:,j] = np.fft.fft(fun[:,j])

fun_fft *= dx / (np.sqrt(2.0 *np.pi)) 
# frequencies in the fft
freq_fft = np.fft.fftfreq(n_x) * 2.0 *np.pi / dx
print('freq_fft', freq_fft)

# reorder things (scipy already contains this reorder function that I struggled with a lot)
freq_fft= np.fft.fftshift(freq_fft)
fun_fft= np.fft.fftshift(fun_fft,axes=0)
print('freq_fft after shift', freq_fft)

# check frequency spacing and integral
dx_fft = 2*np.pi / l_x
print('dx_fft old', dx_fft)
dx_fft = freq_fft[1]-freq_fft[0]
print('dx_fft new', dx_fft)

print(np.sum(abs(fun_fft[:,0]**2))*dx_fft)
print(np.sum(abs(fun_fft[:,2]**2))*dx_fft)
print(np.sum(abs(fun_fft[:,3]**2))*dx_fft)


# write eigenfunciton in position and momentum to files
for j in range(n_f):
    f = open('fft_eig_' + str(j) + '.txt', 'w')
    for i in range(n_x):
        f.write(str(freq_fft[i])+ " ")
        f.write(str(fun_fft[i,j].real) + " ")
        f.write(str(fun_fft[i,j].imag) + " ")
        f.write(str(abs(fun_fft[i,j])) + " ")
        f.write(str(abs(fun_fft[i,j])**2) + " ")
        f.write("\n")
    f.close()

    f = open('harm_eig_' + str(j) + '.txt', 'w')
    for i in range(n_x):
        f.write(str(x[i])+ " ")
        f.write(str(fun[i,j]) + " ")
        f.write(str(abs(fun[i,j])) + " ")
        f.write(str(abs(fun[i,j])**2) + " ")
        f.write("\n")
    f.close()

    
print('Energy, omega * (0.5 +n) [au]', omega * (0.5 +n) )
print('Energy, omega * (0.5 +n) [cm-1]', omega * (0.5 +n) * hartree_to_cm)

# test sample_even
x_sample = sample_funcs.sample_even(x, fun[:,0]**2, 10)
print('x_sample', x_sample)

mom_sample = sample_funcs.sample_even(freq_fft, abs(fun_fft[:,0])**2, 10)
print('mom_sample', mom_sample)

# test sample_random
#x_sample = sample_funcs.sample_random(x, fun[:,0]**2, np.amin(x), np.amax(x), 10000)
x_sample = sample_funcs.sample_random(x, fun[:,0]**2, -1.0, 1.0, 100000)
print('x_sample, random', x_sample)
hist_x, hist_edges_x = np.histogram(x_sample, bins='auto')
edges_mean_x = 0.5 *(hist_edges_x[0:len(hist_edges_x)-1] + hist_edges_x[1:len(hist_edges_x)])   

f = open('hist_x.txt', 'w')
for i in range(len(hist_x)):
    f.write(str(edges_mean_x[i])+ " ")
    f.write(str(hist_x[i]))
    f.write("\n")
f.close()


# for several temperatures, compute thermal averages of delta(x-x_0) and delta(p-p_0)
for T_K in [10,100,300,1000,5000]:
    # compute harmonic density matrix.
    #T_K= 1000 
    T= T_K / T_au
    beta=1.0/T
    
    nn = np.array(range(n_f))
    rho = (1.0-np.exp(-beta*omega))*np.exp(-beta*omega*nn)
    print(rho[:])
    print('np.sum(rho)', np.sum(rho))
    
    # compute thermal distribution
    dist = np.zeros(n_x)
    #for i in range(n_x):
    #    dist[i] = eigenfunctions.thermal_pos_dist(fun, rho, i)
    dist = np.dot(np.abs(fun[:,:]) ** 2, rho[:])
    
    f = open('thermal_dist_' + str(T_K) + '.txt', 'w')
    for i in range(len(x)):
        f.write(str(x[i])+ " ")
        f.write(str(dist[i]))
        f.write("\n")
    f.close()

    # thermal distribution of momenta
    dist_mom = np.zeros(n_x)
    #for i in range(n_x):
    #    dist[i] = eigenfunctions.thermal_pos_dist(fun, rho, i)
    dist_mom = np.dot(np.abs(fun_fft[:,:]) ** 2, rho[:])
    
    f = open('thermal_dist_mom_' + str(T_K) + '.txt', 'w')
    for i in range(len(x)):
        f.write(str(freq_fft[i])+ " ")
        f.write(str(dist_mom[i]))
        f.write("\n")
    f.close()
    

    
    
