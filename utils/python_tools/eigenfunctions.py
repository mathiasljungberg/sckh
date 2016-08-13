import numpy as np
import scipy.special as spec
import scipy.misc as misc 
import ase.units 

#
# everything in atomic units
#

# potential V(x) = 1/2 * mu * omega^2 * x^2
# if units==Bohr, then normalization of integral in bohr
# if units==Ang, then normalaization of integral in Ang
def harm_eigenfun(n, omega, mu, x, units='Bohr'):
    if units == 'Ang':
        x = x / ase.units.Bohr
        norm = 1.0 / np.sqrt(ase.units.Bohr)
    elif units == 'Bohr':
        norm = 1.0
    else:
        raise RuntimeError("units must be either 'Bohr' or 'Ang'")
    
    amu = 1822.888486192
    s = mu * amu *omega
    sx = np.sqrt(s) * x
    nn = 1.0 / np.sqrt(2.0 ** n * misc.factorial(n))
    fun = norm * nn * (s/np.pi)**(1.0/4.0) * np.exp(-(s/2.0)*x ** 2.0) * spec.eval_hermite(n,sx)
    return fun

# potential V(x) = D(1-e^{-ax})^2  
def morse_eigenfun(n, D, a, mu, x, units='Bohr'):
    if units == 'Ang':
        x = x #/ ase.units.Bohr
        a = a #/ ase.units.Bohr
        norm = 1.0 / np.sqrt(ase.units.Bohr)
    elif units == 'Bohr':
        norm = 1.0
    else:
        raise RuntimeError("units must be either 'Bohr' or 'Ang'")
    
    amu = 1822.888486192
    k =2.0 * np.sqrt(2.0 * mu * amu * D) / a
    z =k * np.exp(-a*x)
    b=k -2.0*n -1.0
    assert b>0, 'b = {} is less than 0, there is bound eigenfunction for n={} in this potential'.format(b, n)
    nn = np.sqrt(a * b * misc.factorial(n) / spec.gamma(k-n)) 
    fun = norm * nn * np.exp(-z/2.0) * (z**(b/2.0)) * spec.eval_genlaguerre(n, b, z)
    return fun

