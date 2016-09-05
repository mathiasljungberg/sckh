class bcolors:
    import os
    T = os.getenv('TERM')
    if ( T=='cygwin' or T=='mingw' ) :
        HEADER = '\033[01;35m'
        OKBLUE = '\033[01;34m'
        OKGREEN = '\033[01;32m'
        WARNING = '\033[01;33m'
        FAIL = '\033[01;31m'
        ENDC = '\033[0m'
    else :
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'

    def disable(self):
        self.HEADER = ''
        self.OKBLUE = ''
        self.OKGREEN = ''
        self.WARNING = ''
        self.FAIL = ''
        self.ENDC = ''

#
#
#
def my_exec(cmd, fname_out, fname_err, iv=0, stop_if_CalledProcessError=True) :
  import os,sys
  from subprocess import check_call,CalledProcessError
  cmd_out_err = cmd + ' 1>'+fname_out+ ' 2>'+fname_err
  if iv>0 : print('   '+cmd);
  err_code = 0
  try:
    check_call(cmd_out_err, shell=True);

  except CalledProcessError: 
    if stop_if_CalledProcessError :
      print 'cmd: '+cmd
      print 'cwd: '+os.getcwd() 
      print 'fname_out: '+fname_out
      print 'fname_err: '+fname_err
      print bcolors.FAIL, 'my_exec: CalledProcessError occured, look into logs, please. EXIT', bcolors.ENDC
      err_code = 1
      sys.exit(1)
    else :
      print bcolors.FAIL, 'my_exec: CalledProcessError occured, look into logs, please. CONTINUE', bcolors.ENDC
      err_code = 1
      return err_code
  
  return err_code

#
#
#
def my_cd(dir, iv) :
  import os
  if iv>0 : print('  cd '+dir);
  os.chdir(dir)

#
#
#
def my_diff_ref(fname, iv) : 
  import os,filecmp
  if not filecmp.cmp(fname, fname+'-ref') :
    if iv>0 : print bcolors.WARNING, '! differs '+os.getcwd()+'/'+fname, bcolors.ENDC
    return 1
  else :
    if iv>0 : print bcolors.OKGREEN, 'Ok '+os.getcwd()+'/'+fname, bcolors.ENDC
    return 0
#
#
#
def report_nerr(n, iv) :
  if n>0 : 
    if iv>0 : print bcolors.FAIL, '! Some errors occured...', bcolors.ENDC
    if iv>0 : print bcolors.FAIL, '! Number of errors ', n, bcolors.ENDC
  else :
    if iv>0 : print bcolors.OKBLUE, 'Everything seems to work fine !', bcolors.ENDC

#
#
#
def get_list(fname) :
  ls = []
  f = open(fname, 'r')
  for line in f: ls.append( line.strip() )
  f.close()
  return ls

#
#
#
def get_list_from_file(fname, iv) :
  import sys
  try :
    f = open(fname, 'r')
    list = [];
    for line in f: list.append( line.strip() )
    if iv>0 : 
      print '  '+fname+' defines: ', list
    f.close()
    return list
    
  except IOError :
    print bcolors.FAIL, '  '+fname+' I/O error with this file.', bcolors.ENDC
    print '  Please, create this file first.'
    if fname=='ls_rundirs' : print 'ls|head -n 2 >'+fname; sys.exit(1)
    if fname=='ls_archs' : print 'echo arch.inc.gfortran.standard >'+fname; sys.exit(1)
    sys.exit(1)

#
#
#
def get_list_of_progs(rdirs, fname_ls_progs, iv) : 
  set_of_progs = set()
  for rdir in rdirs :
    if iv>1 : print 'get_list_of_progs: ', rdir+'/'+fname_ls_progs
    map1 = get_prog2results_map(rdir+'/'+fname_ls_progs, iv-1)
    exe_list = list()
    for entry in map1 :
      exe_list.append(entry[0]);
    set_of_progs |= set(exe_list)
  list_of_progs = list(set_of_progs)
  return list_of_progs

#
#
#
def get_prog2results_map(fname_ls_prog2results, iv) : 
  prog2results = get_list_from_file(fname_ls_prog2results, iv)
  list_prog2results = list()
  for string in prog2results:
    p_tokens = string.split()
    list_prog2results.append(p_tokens) 
  if iv>0 : print 'list_prog2results from file '+fname_ls_prog2results
  if iv>0 : print list_prog2results
  return(list_prog2results)

