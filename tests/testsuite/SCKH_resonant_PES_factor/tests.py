from __future__ import print_function
import os, sys
import test_mod.my_bash as my_bash
import test_mod.mod_tests as mod_tests

def run_test():

    # get env variable
    sckh_path = os.environ.get('SCKH_PATH')

    # prepare rundir and go to it
    old_dir = os.getcwd()
    os.system('mkdir -p rundir')
    os.system('rm rundir/*')
    os.system('cp input/* rundir')
    os.chdir('rundir')

    # run calcultion
    err_code = my_bash.my_exec(sckh_path +'/sckh_main', 'test.out', 'test.err')

    # check spectrum
    err_code += mod_tests.my_diff_sum_ref('spectrum_resonant_sigma_532.04.dat',
                                          '../ref/spectrum_resonant_sigma_532.04.dat')
        
    # clean up
    os.chdir(old_dir)
    
    # return 0 if tests succeded
    return err_code
                                                     
