import inspect

curdir = ''
executable_fullpath = ''

#
#
#
def lineno():
    """Returns the current line number in our program."""
    return inspect.currentframe().f_back.f_lineno

#
#
#
def report_passed(msg, fname='unknown', line=-1) :
  from my_bash import bcolors as bc
  if(fname=='unknown' and line==-1) :
    print bc.OKGREEN, msg, bc.ENDC
  else :
    print bc.OKGREEN, msg, ' in ', fname,'(',line,')', bc.ENDC

#
#
#
def report_failed(msg, fname='unknown', line=-1) :
  from my_bash import bcolors as bc
  import sys
  print bc.FAIL, msg, ' in ', fname,'(',line,')', bc.ENDC
  print >>sys.stderr, bc.FAIL, msg, ' in ', fname,'(',line,')', bc.ENDC

#
#
#
def run_tests_py(prefix, bindir, rdirs, iv) : 
  import os,sys,time
  from my_bash import bcolors as bc, my_cd
  import my_bash

  global curdir
  global executable_fullpath
  curdir = ''
  executable_fullpath = ''
  n = 0
  for subdir in rdirs :
    if(len(subdir.strip())==0) : continue
    curdir=prefix+'/'+ subdir
    os.chdir(curdir)
    if(iv>0) : print 
    if(iv>0) : print '  cd '+curdir, time.asctime()
    t1 = time.time()
    f = open(curdir+'/tests.py', 'r')
    exec(f)
    f.close()
    t2 = time.time()
    print ' Finished '+curdir+'.', 'ET: ', bc.WARNING, t2-t1, bc.ENDC, time.asctime()
  return n
# END of run_tests

#
#
#
def my_diff(fname, iv=1, error_measure='REL', threshold=1e-4) :
  import sys
  from mod_numint import numint_trapz
  import numpy as np
  from my_bash import bcolors as bc
  import math
  
  try :
    a = np.genfromtxt(fname, usecols=(0,1) )
  except :
    report_failed('! genfromtxt() '+fname+ ' failed.', __file__, lineno())
    return 1
  try :
    b = np.genfromtxt(fname+'-ref', usecols=(0,1))
  except :
    report_failed('! genfromtxt() '+fname+ '-ref failed.', __file__, lineno())
    return 1

  try :
    err_x = numint_trapz(abs(a[:,0]-b[:,0]))
  except :
    report_failed('shape error?', __file__, lineno())
    err_x = 999.0

  x_ok = True
  if err_x>threshold or math.isnan(err_x):
    x_ok = False
    report_failed('! '+fname+ ' err_x>'+str(threshold)+ ' ' +str(err_x), __file__, lineno())

  try :
    area_1 = numint_trapz(a[:,1], a[:,0])
  except :
    report_failed('shape error?', __file__, lineno())
    area_1 = 0.0

  try :
    area_d = numint_trapz(abs(a[:,1]-b[:,1]), a[:,0])
  except :
    report_failed('shape error?', __file__, lineno())
    area_d = 999.0

  if(error_measure=='ABS') :
    err_y = abs(area_d)
  elif (error_measure=='REL') :
    if (area_1==0):
      err_y = abs(area_d)
    else:
      err_y = abs(area_d/area_1)
  else :
    report_failed("error_measure must be 'ABS' or 'REL' ")
    return 1
  
  if err_y>threshold or not x_ok or math.isnan(err_y):
    report_failed('! '+curdir+'/'+fname+' '+executable_fullpath+ 
      ' err_y>'+str(threshold)+': '+str(err_y), __file__, lineno())
    return 1
  else :
    report_passed('Ok '+fname + ' '+str(err_y)+' vs '+str(threshold))
    return 0


#
#
#
def my_diff_sum_ref(fname, refname, iv=0, maxdim=None, threshold=1e-5, threshold_max=1e-4) :
  import sys
  from mod_numint import numint_trapz
  import numpy as np
  from my_bash import bcolors as bc
  try :
    a = np.genfromtxt(fname)
  except IOError :
    report_failed('! '+fname+ ' does not exists.', __file__, lineno())
    return 1
  try :
    b = np.genfromtxt(refname)
  except IOError :
    report_failed('! '+refname+' does not exists.', __file__, lineno())
    return 1

  if maxdim is not None:
    err_sq_sum = np.sum((a[:,:maxdim]-b[:,:maxdim]) ** 2)
    err_max_sq = np.max((a[:,:maxdim]-b[:,:maxdim]) ** 2)
  else:
    err_sq_sum = np.sum((a-b) ** 2)
    err_max_sq = np.max((a-b) ** 2)

  err_sq_sum_ok = True
  if err_sq_sum > threshold : #1e-5 :
    err_sq_sum_ok = False
    report_failed('! '+fname+ ' err_sq_sum>'+str(threshold) + ' ' +str(err_sq_sum), __file__, lineno())
  
  if err_max_sq> threshold_max or not err_sq_sum_ok: #1e-4 
    report_failed('! '+fname+ ' err_max_sq>' + str(threshold_max) + ' ' +str(err_max_sq), __file__, lineno())
    return 1
  else : 
    report_passed('Ok '+fname + ' ' + str(err_sq_sum) + ' ' +str(err_max_sq), __file__, lineno())
    return 0

#
#
#
def my_get_value_file(fname, string, pos_of_value):
  import sys
  import numpy as np

  f = open(fname, 'r')
  
  for line in f:
    if line.find(string):
      print 'line', line
      line2 = line.split()
      value = line2[pos_of_value]
      print 'value', value
      break
    
  f.close()

  return  float(value)
  


def my_search_value_file(fname, string, pos_of_value, tol=1.0e-6):  
  import sys
  import numpy as np
  from my_bash import bcolors as bc

  value = my_get_value_file(fname, string, pos_of_value)
  value_ref = my_get_value_file(fname + '-ref', string, pos_of_value)

  error =  np.abs(value-value_ref)
  if(error > tol ):
    report_failed('! '+fname+ ' error =' + str(error) +' > '+str(tol), __file__, lineno())
  else:
    report_passed('Ok '+fname+ ' error =' + str(error) +' < '+str(tol), __file__, lineno())
    return 0


#
#
#
def run_test_no_comp(cmd, iv) :
  import os, my_bash
  bname = os.path.basename(cmd)
  n = 0
  #for fname in fnames_out : 
  #my_bash.my_exec('rm -f '+fname, 'rm.out', 'rm.err', iv-1)
  my_bash.my_exec(cmd, bname+'.out', bname+'.err', iv, False)
  #for fname in fnames_out : n = n + comp_fn(fname, iv)  #my_diff_arithmetic_ref(fname, iv)
  #return n    

#
#
#
def run_tests(prefix, bindir, rdirs, ls_progs_fname, iv) : 
  import os,sys,my_bash
  from my_bash import my_cd,my_exec
  from my_bash import bcolors as bc
  n = 0
  for subdir in rdirs :
    curdir=prefix+'/'+ subdir
    my_cd(curdir, iv)
    progs_map = my_bash.get_prog2results_map(ls_progs_fname, iv-1)
    if iv>1 : print '  progs_map: ', progs_map
    rules = list()
    for entry in progs_map :
      if len(entry)>1 :
        rules.append(entry)
      else :
        print bc.WARNING, 'run_tests: len(entry)<2', entry, bc.ENDC
    if iv>0 : print '  rules:',rules

    for rule in rules :
      cmd = bindir+'/'+rule[0]
      n = n + run_test(cmd, rule, iv)
  # END for subdir
  
  return n
# END of run_tests

#
#
#
def run_test(cmd, fnames_out, iv, comp_fn=my_diff) :
  import os, my_bash
  global executable_fullpath

  bname = os.path.basename(cmd)
  n = 0
  for fname in fnames_out : my_bash.my_exec('rm -f '+fname, 'rm.out', 'rm.err', iv-1)
  executable_fullpath = cmd
  err_code = my_bash.my_exec(cmd, bname+'.out', bname+'.err', iv, False)
  if(err_code==0) : 
    for fname in fnames_out : n = n + comp_fn(fname, iv)  #my_diff_arithmetic_ref(fname, iv)
  else :
    n = n + 1
  return n    
