#!/usr/local/python2.6.2/bin/python
# The above line uses a specific python installed on cluster
# this version has many of the packages people will need
import sys
import os, glob, adni_xml2
import multiprocessing as mp
from datetime import datetime
from optparse import OptionParser
from subprocess import call

'''
This script will run a matlab function/script on # processes in parallel
When this python script is specified in the SGE batch submission bash script,
min(NHOSTS, (# of SGE Jobs)) x (PROCESSES_PER_MACHINE) matlab scripts will be run simultaneously
'''

MATLAB_TO_RUN = 'parallel_test'
INPUT_VAR_NAME = 'subj' # used if running a matlab script with input, ignored otherwise
PROCESSES_PER_MACHINE = 2

TASK_ID = os.environ.get('SGE_TASK_ID')
print task_id

def splitInput(input_vals):
    '''
    Split the input values by task id
    Task IDs start indexing at 1
    Input bandwidth for the task is specified by PROCESSES_PER_MACHINE
    '''
    start_index = PROCESSES_PER_MACHINE * (int(TASK_ID) - 1)
    return input_vals[start_index:start_index+PROCESSES_PER_MACHINE]

def runMatlabScript(in_val):
    script_str = "startup;"
    script_str += "%s='%s';" % (INPUT_VAR_NAME, in_val)
    script_str += "run('%s');exit" % (MATLAB_TO_RUN,)
    cmds = ['matlab','-nosplash','-nojvm','-nodisplay','-r',script_str]
    call(cmds)

def runMatlabFunction(in_val):
    script_str = "startup;"
    script_str += "%s('%s');exit" % (MATLAB_TO_RUN, in_val)
    cmds = ['matlab','-nosplash','-nojvm','-nodisplay','-r',script_str]
    call(cmds)

def startProcessing(input_vals, matlab_script=False):
    '''
    Each each input value, start a subprocess, and pass the input value to the matlab script
    '''
    procs = []
    for val in input_val:
        if matlab_function is True:
            new_proc = mp.Process(target=runMatlabFunction, args=(val,))
        else:
            new_proc = mp.Process(target=runMatlabScript, args=(val,))
        new_procs.start()
        procs.append(new_proc)
    for p in procs:
        p.join()

if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-f", "--file", dest="input_file", default=None, help="File containing matlab script input values, one per line")
    parser.add_option("-i", "--input", dest="input_csv", default=None, help="Command line input in csv format")
    parser.add_option("-s", "--script", dest="script", default=False, action="store_true", help="Treat matlab task as script instead of function")
    (options, args) = parser.parse_args()
    in_file = options.input_file
    in_csv = options.input_csv
    is_script = options.script

    if in_file is not None:
        inputs = open(in_file).readlines()
    elif in_csv is not None:
        inputs = [_.strip() for _ in in_csv.split(',') if _.strip() != '']
    else:
        inputs = []

    task_inputs = splitInput(inputs)
    startProcessing(task_inputs, matlab_script=is_script)

