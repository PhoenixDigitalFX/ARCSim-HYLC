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
ap.add_argument("-d", "--debug", action='store_true',
                help="...")
ap.add_argument("-c", "--conf", default="0",
                help="...")
ap.add_argument("-o", "--op", default="simulate",
                help="...")
# TODO add a  no compilation flag
args, unknownargs = ap.parse_known_args()
args = vars(args)

# BUILD
sourcedir = os.getcwd()  # where the main CMakeLists.txt is
builddir = os.path.join(
    os.getcwd(), "build-Debug" if args['debug'] else "build-Release")

build(sourcedir=sourcedir, builddir=builddir, debug=args['debug'])

# RUN
workdir = os.getcwd()
executable = os.path.join(builddir, "bin", "arcsim_0.2.1")
# from workdir call
# ./build-Release/bin/arcsim_0.2.1 $op $configfile ${3}

# operation: simulate, simulateoffline etc.
op = args['op']
# configfile
if args['conf'].endswith(".json"):
    conf = args['conf']
else:
    conf = [
        "conf/hylc_sphere_noremesh.json",
        "conf/hylc_sphere.json",
        "conf/hylc_stretchx.json",
        "conf/hylc_stretchy.json",
        "conf/hylc_shearx.json"
    ][int(args['conf'])]

simargs = unknownargs # remaining args
print("Executing:", executable)
try:
    subprocess.check_call([executable, op, conf] + simargs, cwd=workdir)
except KeyboardInterrupt:
    print("PY: Aborting execution")
