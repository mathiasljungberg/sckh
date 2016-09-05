from __future__ import print_function
import sys, os

# check if we have the SCKH_PATH environment variable set
sckh_path = os.environ.get('SCKH_PATH')
if sckh_path is None:
    raise RuntimeError('Please set' +
                       ' the SCKH_PATH' +
                       ' environment variable')

dirs = ['KH_simple', 'KH_resonant_simple']

old_dir = os.getcwd()

for d in dirs:

    print('Entering directory '+ d)
    os.chdir(d)
    sys.path.insert(0,'.')

    import tests 
    
    err_code = tests.run_test()
    
    if err_code > 0:
        print('test ' +d + ' failed!')
        
    if 'tests' in sys.modules:
        del sys.modules['tests']
        sys.path.pop(0)

    os.chdir(old_dir)
    print('Exiting')

