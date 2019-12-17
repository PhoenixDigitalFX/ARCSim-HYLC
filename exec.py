#python exec.py -o resumeoffline OUTDIR FRAMENUM

import os
import sys
import subprocess
import argparse


def build(sourcedir, builddir, debug=False):
    cfg = 'Debug' if debug else 'Release'
    cmake_args = ['-DCMAKE_BUILD_TYPE=' + cfg]
    build_args = ['--config', cfg]
    build_args += ['--', '-j']

    os.makedirs(builddir, exist_ok=True)

    subprocess.check_call(['cmake', sourcedir] +
                          cmake_args, cwd=builddir)
    subprocess.check_call(['cmake', '--build', '.'] +
                          build_args, cwd=builddir)  # basically make



# get command line arguments
ap = argparse.ArgumentParser()
ap.add_argument("-b", "--build", default="1",
                help="...")
ap.add_argument("-r", "--run", default="1",
                help="...")
ap.add_argument("-d", "--debug", action='store_true',
                help="...")
# ap.add_argument("-c", "--conf",
#                 help="...")
# ap.add_argument("-C", "--conf",
#                 help="...")
ap.add_argument("-o", "--op", default="simulate",
                help="...")
# TODO add a  no compilation flag
args, unknownargs = ap.parse_known_args()
args = vars(args)
args['build'] = args['build'] != "0"
args['run'] = args['run'] != "0"

# BUILD
sourcedir = os.getcwd()  # where the main CMakeLists.txt is
builddir = os.path.join(
    os.getcwd(), "build-Debug" if args['debug'] else "build-Release")


if args['build']:
    build(sourcedir=sourcedir, builddir=builddir, debug=args['debug'])

if not args['run']:
    exit()
# RUN
workdir = os.getcwd()
executable = os.path.join(builddir, "bin", "arcsim_0.2.1")
# from workdir call
# ./build-Release/bin/arcsim_0.2.1 $op $configfile ${3}

# operation: simulate, simulateoffline, replay etc.
op = args['op']

simargs = unknownargs # remaining args
try:
    if op in ["replay","resume", "resumeoffline"]:
        print("Executing:", [executable, op] + simargs)
        subprocess.check_call([executable, op] + simargs, cwd=workdir)
    else:
        print("Executing:", [executable, op] + simargs)
        subprocess.check_call([executable, op] + simargs, cwd=workdir)
except KeyboardInterrupt:
    print("PY: Aborting execution")
