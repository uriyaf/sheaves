# command line main file
# type -h for details.
# A generic command line looks like:
# main.py <sheaf name> -i <no. of steps> -s <random seed (optional)> -o <output file>
# e.g.
# main.py LSVq1_sheafdim17_1_3_275_s1 -i 3 -s 1 -o my_output.txt

import argparse
import sheaf_examples
import pickle
from modification_process import ModificationProcess


parser = argparse.ArgumentParser(description='Run an experiment using one of the sheaves from sheaf_examples.py.')
parser.add_argument('sheaf', type=str, 
                    help='name of the sheaf (e.g. LSVq1_sheafdim17_1_3_4643, LSVq1_sheafdim33_1_3_9011, LSVq1_sheafdim49_1_3_13379, etc.)')
parser.add_argument('-r', '--randseed', type=int, default=None, \
                    help='seed for the random function (default: None)')
parser.add_argument('-i', '--iters', type=int, default=3, \
                    help='number of iterations (steps) performed by the experiment (default: 3)')
parser.add_argument('-o', '--output', type=str, default=None, \
                    help='name of output file (default name is ./outputs/<sheaf>_<mode><dimE0>i<iters>r<randseed>.pkl)')
parser.add_argument('-l', '--log', type=str, default=None, \
                    help='name of log file (default name is ./outputs/<sheaf>_<mode><dimE0>i<iters>r<randseed>.log)')
parser.add_argument('-cy', '--cocycles', dest='mode', action='store_const', const='cy', default='cb', \
                    help='specify in order to use cocycles rather than coboundaries')
parser.add_argument('-n', dest='dimE0', type=int, default=1, \
                    help='initial dimension of E0 (defaul: one)')
parser.add_argument('--nolog', action='store_const', const=True, default=False,
                    help='if specified, no log file will be written.')
parser.add_argument('--nooutput', action='store_const', const=True, default=False,
                    help='if specified, no output file will be written.')
parser.add_argument('--loadsheaf', action='store_const', const=True, default=False,
                    help='if specified, SHEAF is interpreted as file and the sheaf is being loaded from that file using \
the pickle module. Load only files which you trust as pickle may run code on your processor.')

def no_log(*args, **kwargs):
    pass

def print_and_save(filename):
    """
    Returns a function that recieves arguments (same format as print), prints it and saves it into the given filename.
    """
    f = open(filename, mode='w')
    f.close() # deletes the file
    def res(*args, **kwargs):
        print(*args, **kwargs)
        f = open(filename, mode='a')
        print(*args, **kwargs, file=f)
        f.close()
    return res

def no_save(obj, *args, **kwargs):
    pass

def pickle_save(filename):
    """
    Returns a function that recieves arguments (same format as in pickle.dump) and saves it into the given filename.
    """
    def res(obj, *args, **kwargs):
        f = open(filename, mode='wb')
        pickle.dump(obj, f, *args, **kwargs)
        f.close()
    return res

def experiment(sheaf_name, mode, dimE0, steps, rand_seed, log, save, loadsheaf):
    if loadsheaf:
        log("Loading sheaf from "+sheaf_name+"... ", end="", flush=True)
        f = open(sheaf_name, mode='rb')
        sheaf = pickle.load(f)
        f.close()
        log("done.", flush=True)
    else:
        log("Constructing sheaf "+sheaf_name+"... ", end="")
        sheaf = getattr(sheaf_examples, sheaf_name)()
        log("done.")
    log("Sheaf dimension =", sheaf.face_dim(sheaf.base_complex().rand_face(0)))
    log("Initializing modification process... ")
    mp = ModificationProcess(sheaf, rand_seed)
    save(mp)
    if mode == 'cb':
        mp.init_from_1coboundary(dimE0)
    if mode == 'cy':
        mp.init_from_1cocycle_basis(dimE0)
    log("done.")
    log("dim(E0) =", len(mp.get_E_basis()))
    for i in range(steps):
        log("--Stage", i, "begins.")
        mp.step(log=log)
        log("dim(E'%d) =" % i, len(mp.get_E_prime_basis(i)), flush=True)
        save(mp)
    return mp

args = parser.parse_args()
## for non-command line main use, e.g.,
##args = parser.parse_args(['./sheaves/LSVq1_sheafdim65_1_3_17747.pkl', '-i', '3', '-r', '12', \
##                          '--loadsheaf', \
##                          '-o', './outputs/LSVq1_sheafdim65_1_3_17747_i3r12.pkl', \
##                          '-l', './outputs/LSVq1_sheafdim65_1_3_17747_i3r12.log'])

sheaf_name = args.sheaf
mode = args.mode
dimE0 = args.dimE0
steps = args.iters
rseed = args.randseed
loadsheaf = args.loadsheaf

all_is_fine = True

if args.nooutput:
    save = no_save
elif args.output is None:
    if loadsheaf:
        print("When loading a sheaf from a file, the output file must be specified (use the -o or --nooutput options).")
        all_is_fine = False
    else:
        save = pickle_save('./outputs/' + sheaf_name + "_" + mode + str(dimE0) + "i" + str(steps) + "r" + str(rseed) + ".pkl") 
else:
    save = pickle_save(args.output)

if args.nolog:
    log = no_log
if args.log is None:
    if loadsheaf:
        print("When loading a sheaf from a file, the log file must be specified (use the -l or --nolog options).")
        all_is_fine = False
    else:
        log = print_and_save('./outputs/' + sheaf_name + "_" + mode + str(dimE0) + "i" + str(steps) + "r" + str(rseed) + ".log")
else:
    log = print_and_save(args.log)

if all_is_fine:
    log('Parameters:', args, flush=True)
    log('')
    mp = experiment(sheaf_name=sheaf_name, \
                    mode=mode, \
                    dimE0=dimE0, \
                    steps=steps, \
                    rand_seed=rseed, \
                    log=log, \
                    save=save, \
                    loadsheaf=args.loadsheaf)
