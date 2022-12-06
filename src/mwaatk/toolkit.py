from mwaatk.absorbancecalculator import AbsorbanceCalculator
from mwaatk.alphaminimizer import AlphaMinimizer
from mwaatk.compoundpreprocessor import CompoundPreprocessor
from mwaatk.datasetabsorbancecalculator import DatasetAbsorbanceCalculator
from mwaatk.dataset import Dataset
from mwaatk.filtersample import Filter
from mwaatk.preprocessor import Preprocessor
from mwaatk.angstrom import AngstromCalculator
import json
import sys
from pprint import pprint

def prepro(options):
    """
    Wrapper function for the integrals calculations from BLAnCA raw data.
    """
    with open(options, 'r') as f:
        opt = json.load(f)

    cp = CompoundPreprocessor(opt)
    cp.confirm_start()
    input()
    cp.raw_file_open()
    cp.set_gains_dark()
    cp.calculate_avg()
    cp.normalize_gains()
    cp.calculate_alpha()
    cp.integrate()
    bfile, wfile = cp.write()
    return bfile, wfile

def absorb(options, bfile=None, wfile=None):
    """
    Wrapper function for the ABS minimization routine.
    
    Parameters:
        options: string, 
            a path to the configuration file
        wfile: string (optional),
            a path to the input JSON file containing the white integrals
        bfile: string (optional),
            a path to the input JSON file containing the black integrals
    """
    with open(options, 'r') as f:
        opt = json.load(f)

    infile_b = opt['b_integrals'] if bfile is None else bfile
    infile_w = opt['w_integrals'] if wfile is None else wfile
    
    ds = Dataset(infile_b, infile_w, opt)
    ds.create_set()
    s = ds.get_full_set()
    dac = DatasetAbsorbanceCalculator(s, opt)
    dac.do_minimization()
    dac.calc_abs()
    abs_file = dac.write()
    return abs_file

def angstr(options, fname=None):
    """
    Wrapper function for the AAE routine.

    Parameters:
        options: string,
            a path to the configuration file
        fname: string (optional),
            a path to the input JSON file containing absorbance data
    """
    with open(options, 'r') as f:
        opt = json.load(f)
   
    infile = opt['abs_file'] if fname is None else fname
    
    ac = AngstromCalculator(infile, opt)
    ac.fit_angstrom()
    aae_file = ac.writeout()
    return aae_file


def run(options):
    """
    Wraps the entire functionality of the toolkit so far in a 
    sequential manner
    """
    integrals = prepro(options)
    absorbance = absorb(options, bfile=integrals[0], wfile=integrals[1])
    angstr(options, fname=absorbance)
 

def cli_script():
    """Command Line script to use as entry point"""
    try:
        if '-p' in sys.argv:
            prepro(sys.argv[-1])
        if '-a' in sys.argv:
            absorb(sys.argv[-1])
        if '-f' in sys.argv:
            angstr(sys.argv[-1])
        elif '-p' not in sys.argv and '-a' not in sys.argv and '-f' not in sys.argv:
            run(sys.argv[1])
    except IndexError:
        print("ERROR: unrecognized options and/or missing configuration file.\nUsage:\n\
    mwaatk [-p] [-a] [-f] <path_to_config_file>")
