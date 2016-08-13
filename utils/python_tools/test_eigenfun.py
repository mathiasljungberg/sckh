from __future__ import print_function
import eigenfunctions
import eigenfunctions_fortran
import numpy as np
import ase.units as units 

amu = 1822.888486192

print('Ha', units.Ha)
print('Bohr', units.Bohr)
print('eV', units.eV)
n=0
#omega=443.1e-3 / units.Ha # 443.1 meV = 3574.0 cm^{-1}
omega=4.1e-3 / units.Ha # 443.1 meV = 3574.0 cm^{-1}
hartree_to_cm = eigenfunctions_fortran.parameters.cm * eigenfunctions_fortran.parameters.hartree # this is correct 
omega_cm= omega * hartree_to_cm
print('omega', omega)
print('omega_cm', omega_cm)
mu=1.0
D=0.3
a=np.sqrt(mu * amu * omega ** 2 / (2.0 * D) ) # this gives the harmonic frequency omega to the morse potential

print('a', a)

n_x=10000
x=np.linspace(-2.0,2.0, n_x)
dx=x[1]-x[0]
dx_Ang=dx * units.Bohr 
fun = np.zeros(n_x)
fun_f = np.zeros(n_x)
fun_morse = np.zeros(n_x)
fun_morse_f = np.zeros(n_x)

for i in range(len(x)):
    fun[i] = eigenfunctions.harm_eigenfun(n, omega, mu, x[i] , units='Bohr' )

    # fortran routine works with omega [cm^{-1}], x [Ang], mu [amu]
    fun_f[i] = eigenfunctions_fortran.harm_eigenfun(n, omega_cm , mu, x[i] * units.Bohr )

    fun_morse[i] = eigenfunctions.morse_eigenfun(n, D, a, mu, x[i], units='Bohr')

    # fortran routine works with a [Ang^{-1}], x [Ang], D Hartree, mu [amu]
    fun_morse_f[i] = eigenfunctions_fortran.morse_eigenfun(n, D, a / units.Bohr, mu, x[i] * units.Bohr)

# test
d2fun_dx2 = np.gradient(np.gradient(fun, dx),dx)
omega_eval = np.sum(fun**2 * 0.5 * (mu*amu) * omega **2 * x**2 -0.5*fun*d2fun_dx2 / (mu*amu))*dx

print('Energy, omega * (0.5 +n) [au]', omega * (0.5 +n) )
print('Energy, omega * (0.5 +n) [cm-1]', omega * (0.5 +n) * hartree_to_cm)

print('Energy, evaluated [au]', omega_eval)
print('Energy, evaluated [cm-1]', omega_eval * hartree_to_cm)

E_morse = omega*(0.5 + n) - omega**2/(4*D)*(0.5 +n)**2
E_morse_next = omega*(0.5 + n+1) - omega**2/(4*D)*(0.5 +n+1)**2

print('Morse,  energy, [au]', E_morse)
print('Morse,  energy, [cm-1]', E_morse * hartree_to_cm)

d2fun_dx2 = np.gradient(np.gradient(fun_morse, dx),dx)
omega_eval = np.sum(fun_morse**2 * D * (1.0 -np.exp(-a*x))**2 -0.5 * fun_morse*d2fun_dx2 / (mu*amu))*dx

print('Morse,  energy, evaluated [au]', omega_eval)
print('Morse,  energy, evaluated [cm-1]', omega_eval* hartree_to_cm)

print('Morse, diff energy, [au]', E_morse_next-E_morse)
print('Morse, diff energy, [cm-1]', (E_morse_next -E_morse) * hartree_to_cm)

# compute harmonic density matrix. Temperature in what units?
T_au = 3.15777504e5 
print('300 K in au',300 / T_au)

nn =np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20])
for i in range(10):
    T_K=(1 +100 * i) 
    T= T_K / T_au
    beta=1.0/T
    rho = (1.0-np.exp(-beta*omega))*np.exp(-beta*omega*nn)
    print(T_K, " ", rho)
    print('np.sum(rho)', np.sum(rho))
#print('Energy, evaluated [cm-1]', omega_eval * hartree_to_cm)



f = open('eigenfun.txt', 'w')
for i in range(len(x)):
    f.write(str(x[i])+ " ")
    f.write(str(fun[i]))
    f.write("\n")
f.close()

f = open('eigenfun_f.txt', 'w')
for i in range(len(x)):
    f.write(str(x[i])+ " ")
    f.write(str(fun_f[i]))
    #f.write( " ")
    #f.write(str(fun_f[i]))
    f.write("\n")
f.close()

f = open('eigenfun_morse.txt', 'w')
for i in range(len(x)):
    f.write(str(x[i])+ " ")
    f.write(str(fun_morse[i]))
    f.write("\n")
f.close()

f = open('eigenfun_morse_f.txt', 'w')
for i in range(len(x)):
    f.write(str(x[i])+ " ")
    f.write(str(fun_morse_f[i]))
    f.write("\n")
f.close()
    
# fun normalized by integration in atomic units
#print(np.sum(fun**2)*dx)
print(np.sum(fun**2)*dx)
# fun_f normalized by integration in Angstroms
print(np.sum(fun_f**2)*dx_Ang)

print(np.sum(fun_morse**2)*dx)
print(np.sum(fun_morse_f**2)*dx_Ang)

