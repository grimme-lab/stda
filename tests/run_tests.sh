#!/bin/bash

# check parameters
if [ "$#" -ne 1 ]; then
    >&2 echo "Illegal number of parameters"
    >&2 echo "usage: run_test.sh STD2_EXE"
    exit 1
fi

# EXE
STD2_EXE=${1}

# directory for the results
EXPECTED_RESULTS_DIR="./expected_results"

# total number of errors in the whole suite
TOT_NERRORS=0

# -- test suite for std2 (excitations energies)
declare -A TESTS_STD2_EXCI

TESTS_STD2_EXCI[EXCI_TEST1]='-f water_sto3g.molden -sty 3 -ax 1.0 -e 20'
TESTS_STD2_EXCI[EXCI_TEST2]='-f water_sto3g.molden -sty 3 -ax 0.5 -e 30'
TESTS_STD2_EXCI[EXCI_TEST3]='-f water_sto3g.molden -sty 3 -ax 1.0 -e 30 -be 1.0 -al 3.0'

# run test suite
for TEST in "${!TESTS_STD2_EXCI[@]}"; do
  # fetch parameters and run test
  TEST_PARAMS=${TESTS_STD2_EXCI[${TEST}]}
  echo "** Running test ${TEST} using '${TEST_PARAMS}'"
  ${STD2_EXE} ${TEST_PARAMS} > /dev/null
  
  # compare with expected results
  python3 compare.py -t tda ${EXPECTED_RESULTS_DIR}/${TEST}.dat tda.dat
  
  # count errors
  TEST_NERR=$?
  TOT_NERRORS=$(( ${TOT_NERRORS} + ${TEST_NERR} ))
done

# -- test suite for std2 (beta)
declare -A TESTS_STD2_BETA

TESTS_STD2_BETA[BETA_TEST1]='-f water_sto3g.molden -sty 3 -ax 1.0 -e 20 -resp 2'
TESTS_STD2_BETA[BETA_TEST2]='-f water_sto3g.molden -sty 3 -ax 0.5 -e 30 -resp 2'
TESTS_STD2_BETA[BETA_TEST3]='-f water_sto3g.molden -sty 3 -ax 1.0 -e 30 -be 1.0 -al 3.0 -resp 2'

# run test suite
for TEST in "${!TESTS_STD2_BETA[@]}"; do
  # fetch parameters and run test
  TEST_PARAMS=${TESTS_STD2_BETA[${TEST}]}
  echo "** Running test ${TEST} using '${TEST_PARAMS}'"
  ${STD2_EXE} ${TEST_PARAMS} > /dev/null
  
  # compare with expected results
  python3 compare.py -t beta_HRS ${EXPECTED_RESULTS_DIR}/${TEST}.txt beta_HRS
  
  # count errors
  TEST_NERR=$?
  TOT_NERRORS=$(( ${TOT_NERRORS} + ${TEST_NERR} ))
done

# -- test suite for std2 (beta)
declare -A TESTS_STD2_2PA

TESTS_STD2_2PA[2PA_TEST1]='-f water_sto3g.molden -sty 3 -ax 1.0 -e 20 -2PA 3'
TESTS_STD2_2PA[2PA_TEST2]='-f water_sto3g.molden -sty 3 -ax 0.5 -e 30 -2PA 3'
TESTS_STD2_2PA[2PA_TEST3]='-f water_sto3g.molden -sty 3 -ax 1.0 -e 30 -be 1.0 -al 3.0 -2PA 3'

# run test suite
for TEST in "${!TESTS_STD2_2PA[@]}"; do
  # fetch parameters and run test
  TEST_PARAMS=${TESTS_STD2_2PA[${TEST}]}
  echo "** Running test ${TEST} using '${TEST_PARAMS}'"
  ${STD2_EXE} ${TEST_PARAMS} > /dev/null
  
  # compare with expected results
  python3 compare.py -t 2PA  ${EXPECTED_RESULTS_DIR}/${TEST}.txt 2PA-abs
  
  # count errors
  TEST_NERR=$?
  TOT_NERRORS=$(( ${TOT_NERRORS} + ${TEST_NERR} ))
done

# exit
exit ${TOT_NERRORS}
