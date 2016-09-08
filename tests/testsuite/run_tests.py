from __future__ import print_function
import sys, os

# check if we have the SCKH_PATH environment variable set
sckh_path = os.environ.get('SCKH_PATH')
if sckh_path is None:
    raise RuntimeError('Please set' +
                       ' the SCKH_PATH' +
                       ' environment variable')


if (len(sys.argv) > 1):
    dirs=sys.argv[1:]
else:
    #dirs = ['XAS', 'SCKH_PES', 'SCKH','vib_finite_diff', 'KH', 'KH_resonant', 'KH_resonant_el']
    dirs= ['SCKH_PES', 'SCKH', 'KH', 'KH_resonant', 'KH_resonant_el']

print(dirs)
    
old_dir = os.getcwd()
err_tot=0

for d in dirs:

    print('Entering directory '+ d)
    os.chdir(d)
    sys.path.insert(0,'.')

    import tests 
    
    err_code = tests.run_test()

    err_tot += err_code 
    
    if err_code > 0:
        print('test ' +d + ' failed!')
        
    if 'tests' in sys.modules:
        del sys.modules['tests']
        sys.path.pop(0)

    os.chdir(old_dir)
    print('Exiting')

if err_tot ==0:
    print('All tests passed')
else:
    print('{} test failed'.format(err_tot))
