"""
Compare the results for different files of std2 and report errors if the info they contain differs
"""

import argparse
import math

def get_tda_results(f):
    """
    Extract the list of excitations with corresponding fL
    """
    results = {
        'excitation energy': [],
        'fL': [],
    }
    
    lines = f.readlines()
    datxy_found = -1
    for i, line in enumerate(lines):
        if 'DATXY' in line:
            datxy_found = i
            break
    
    if datxy_found < 0:
        raise Exception('DATXY not found')
    
    for line in lines[datxy_found + 1:]:
        chunks = line.split()
        
        results['excitation energy'].append(float(chunks[1]))
        results['fL'].append(float(chunks[2]))
    
    return results

def get_beta_HRS_results(f):
    """
    Extract the list of beta HRS
    """
    
    results = {
        'beta HRS': [],
    }
    
    lines = f.readlines()
    
    for line in lines:
        chunks = line.split()
        results['beta HRS'].append(float(chunks[1]))
    
    return results

def get_2PA_results(f):
    """
    Extract the list of 2PA quantities
    """
    
    results = {
        'excitation energy': [],
        'Delta_2PA_//': [],
        'Delta_2PA__|_': [],
        'Delta_2PA_circ': [],
        'rho': [],
    }
    lines = f.readlines()
    
    for line in lines:
        chunks = line.split()
        results['excitation energy'].append(float(chunks[1]))
        results['Delta_2PA_//'].append(float(chunks[2]))
        results['Delta_2PA__|_'].append(float(chunks[3]))
        results['Delta_2PA_circ'].append(float(chunks[4]))
        results['rho'].append(float(chunks[5]))
    
    return results
    

def check_equal(expected, actual, label, delta=1e-3):
    """
    Check if two arrays contain the same data, with a difference < `delta`
    """
    
    if len(expected) != len(actual):
        print('! ERROR: the size of `expected` and `actual` differs')
        return abs(len(expected) - len(actual))
    
    n_errors = 0
    for left, right in zip(expected[label], actual[label]):
        if math.fabs(left - right) > delta:
            n_errors += 1
    
    if n_errors > 0:
        print('! ERROR: differences are larger than {} for {}'.format(delta, label))
        print(' expected ::', expected[label])
        print(' actual   ::', actual[label])
    
    return n_errors
            
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--type', choices=['tda', 'beta_HRS', '2PA'])
parser.add_argument('expected', type=argparse.FileType('r'))
parser.add_argument('actual', type=argparse.FileType('r'))

args = parser.parse_args()

# count errors
n_errors = 0

if args.type == 'tda':
    expected_data = get_tda_results(args.expected)
    actual_data = get_tda_results(args.actual)

    n_errors += check_equal(expected_data, actual_data, 'excitation energy')
    n_errors += check_equal(expected_data, actual_data, 'fL')

if args.type == 'beta_HRS':
    expected_data = get_beta_HRS_results(args.expected)
    actual_data = get_beta_HRS_results(args.actual)

    n_errors += check_equal(expected_data, actual_data, 'beta HRS')

if args.type == '2PA':
    expected_data = get_2PA_results(args.expected)
    actual_data = get_2PA_results(args.actual)

    n_errors += check_equal(expected_data, actual_data, 'excitation energy')
    n_errors += check_equal(expected_data, actual_data, 'Delta_2PA_//', delta=1e-2)
    n_errors += check_equal(expected_data, actual_data, 'Delta_2PA__|_', delta=1e-2)
    n_errors += check_equal(expected_data, actual_data, 'Delta_2PA_circ', delta=1e-2)
    n_errors += check_equal(expected_data, actual_data, 'rho', delta=1e-2)

# exit
if n_errors == 0:
    print('OK')

exit(n_errors)
