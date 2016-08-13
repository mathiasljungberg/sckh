from __future__ import print_function
import numpy as np
import ase.units as units 
import scipy.interpolate as interpolate

# rewritten from the fortran routines in SCKH_PES
def sample_even(x, fun, n_sample):

    # find integral of fun
    dx=x[1]-x[0]
    n_x=len(x)
    x_sample = np.zeros(n_sample)
    integral = np.cumsum(fun) * dx

    # divide the interval [0, integral[n_x-1]] into n_sample equal parts
    # then the sample points are placed at the center of each interval
    i_sample = np.array(range(n_sample)) 
    points_sample = ( (i_sample+0.5) * integral[n_x-1]) / n_sample
    print('points_sample', points_sample )

    # 
    for i in range(n_sample):
        for j in range(n_x):
            if(integral[j] < points_sample[i]):
                x_sample[i] = x[j]
            else:
                break

    return x_sample
            
def sample_random(x, fun, x_low, x_high, n_sample):
    
    n_x = len(x)
    x_sample = np.zeros(n_sample)
    f_max = np.amax(fun)
    
    # spline function
    spl_coeff = interpolate.splrep(x, fun, s=0)

 
    
    for i in range(n_sample):
        condition = False 
        counter = 0
        while not condition:
            x_sample_i = x_low + np.random.random() * (x_high-x_low)
            y_sample_i = np.random.random() * f_max
            f = interpolate.splev(x_sample_i, spl_coeff, der=0)
            #x_sample_i =  sample_min + np.random.random() * (sample_max-sample_min)
            #y_sample_j =  np.random.random() * f_max
            #f = interpolate.splev(x_sample_i, spl_coeff, der=0)
            counter += 1
            
            if counter > 10000:
                raise RuntimeError('sample could not be done in 10000 tries') 
            
            if (f> y_sample_i) : 
                condition = True
                x_sample[i] = x_sample_i
                
    return x_sample

    
#
#    
#subroutine sample_random(x, funct, x_sampl)
#use m_precision, only: wp
#use m_splines, only: spline, splint
#!use m_constants, only: const
#
#! passed variables
#real(kind=wp), dimension(:), intent(in)::x, funct
#real(kind=wp), dimension(:), intent(out)::x_sampl
#! local variables
#integer,parameter:: npoints_new=10000
#integer:: npoints, npoints_sampl
#real(kind=wp), dimension(:),allocatable::yder
#real(kind=wp):: funct_x_tmp, x_tmp, y_tmp, tol, x_min, x_max,y_min, y_max, x_l, y_l,ran1,ran2
#integer:: i,j
#
#npoints = size(x)
#npoints_sampl = size(x_sampl)
#
#allocate(yder(npoints))
#
#do i=1,npoints
#if(funct(i) .lt. 0) then
#write(6,*) "function has a negative value in sample_random!", i, funct(i)
#stop
#end if
#end do
#
#
#tol = maxval(funct) * 5.0d-3
#
#!find the interval of sampling
#y_min = 0
#y_max = maxval(funct)
#y_l = y_max - y_min
#
#do i=1,npoints
#if(funct(i) .le. tol) then
#x_min = x(i)
#else
#exit
#end if
#end do
#
#do i=npoints, 1, -1
#if(funct(i) .le. tol) then
#x_max = x(i)
#else
#exit
#end if
#end do
#
#x_l = x_max - x_min
#
#write(6,*) "x_min, x_max, x_l", x_min, x_max, x_l
#
#! sample!
#
#j=1
#call random_seed
#call spline(x,funct,npoints,1.0d30,1.0d30,yder)
#
#do while(j .le. npoints_sampl )
#call random_number(ran1)
#call random_number(ran2)
#
#x_tmp = x_l * ran1 + x_min
#y_tmp = y_l * ran2 + y_min
#
#!call spline_one(x, funct, npoints, x_tmp, funct_x_tmp)
#call splint(x,funct,yder,npoints,x_tmp, funct_x_tmp)
#
#if(y_tmp .le. funct_x_tmp) then
#x_sampl(j) = x_tmp
#j=j+1
#end if
#end do
#
#deallocate(yder)
#
#end subroutine sample_random
#
