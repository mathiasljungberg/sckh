from __future__ import print_function
import sys, os
import importlib.util
import subprocess

# check if we have the SCKH_PATH environment variable set
sckh_path = os.environ.get('SCKH_PATH')
if sckh_path is None:
    raise RuntimeError('Please set' +
                       ' the SCKH_PATH' +
                       ' environment variable')


if (len(sys.argv) > 1):
    dirs=sys.argv[1:]
else:
    dirs= ['KH',
           'KH_resonant',
           'KH_resonant_el',
           'SCKH',
           'SCKH_PES',
           'SCKH_resonant_PES',
           'XAS',
           'SCXAS_PES',
           'vib_finite_diff',
           ]

print(dirs)
    
old_dir = os.getcwd()
err_tot=0

for d in dirs:

    print('Entering directory '+ d)
    test_dir = os.path.join(old_dir, d)
    test_file = os.path.join(test_dir, 'tests.py')

    # Load the module from the explicit file path
    spec = importlib.util.spec_from_file_location("tests", test_file)
    tests = importlib.util.module_from_spec(spec)

    # Temporarily change to test directory in case tests.py uses relative paths
    os.chdir(test_dir)
    spec.loader.exec_module(tests)

    err_code = tests.run_test()

    err_tot += err_code

    if err_code > 0:
        print('test ' +d + ' failed!')

    os.chdir(old_dir)
    print('Exiting')

if err_tot ==0:
    print('All tests passed')
else:
    print('{} test failed'.format(err_tot))
