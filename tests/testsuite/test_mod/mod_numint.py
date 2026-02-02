from __future__ import print_function

#
#
#
def numint_trapz(*args) :
  import sys
 
  try : 
    import scipy.integrate
    if len(args)==1 :
      return scipy.integrate.trapezoid(args[0])
    elif  len(args)==2 :
      return scipy.integrate.trapezoid(args[0], args[1])
    else :
      print('numint_trapz: first 2 arguments will be used with scipy.integrate.trapezoid!')
      return scipy.integrate.trapezoid(args[0], args[1])
  
  except  ImportError : 
    if len(args)==1 :
      return numint_own_trapz(args[0])
    elif  len(args)==2 :
      return numint_own_trapz(args[0], args[1])
    else :
      print('numint_trapz: first 2 arguments will be used with numint_own_trapz!')
      return numint_own_trapz(args[0], args[1])

#
#
#
def numint_own_trapz(*args) :
  la = len(args)
  if la<1 : 
    print('numint_own_trapz: len(args)<1', la)
    sys.exit(1)
  elif la==2 : 
    f = args[0]
    x = args[1]
  elif la==1 : 
    f = args[0]
    x = range(len(f))
  else:
    print('numint_own_trapz: len(args)>2 ?', la)
    sys.exit(1)
    
  lenf = len(f)
  if not lenf==len(x) : 
    print('numint_own_trapz: len(f), len(x) ', len(f), len(x))
    print('numint_trapez: len(f)/=len(x)! ')
    sys.exit(1)

  lmin = 0
  lmax = lenf-1
  integral = -x[lmin]*(f[lmin]+f[lmin+1])
  for i in range(lmin+1,lmax) :
    integral = integral + (f[i-1]-f[i+1])*x[i]
  integral = integral + x[lmax]*(f[lmax-1]+f[lmax])
  integral = integral/2
  return integral
