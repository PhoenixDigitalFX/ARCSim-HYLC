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
# ap = argparse.ArgumentParser()
# ap.add_argument("-d", "--debug", action='store_true',
#                 help="...")
# ap.add_argument("-c", "--conf", default="0",
#                 help="...")
# ap.add_argument("-o", "--op", default="simulate",
#                 help="...")
# args, unknownargs = ap.parse_known_args()
# args = vars(args)

# BUILD
sourcedir = os.path.join(os.getcwd(),"..") 
builddir = os.path.join(
    os.getcwd(), "..", "build-Release")

build(sourcedir=sourcedir, builddir=builddir, debug=False)

# RUN
workdir = os.path.join(os.getcwd(),"..") 
executable = os.path.join(builddir, "bin", "arcsim_0.2.1")
# from workdir call
# ./build-Release/bin/arcsim_0.2.1 $op $configfile ${3}

op = "simulateoffline"

conffolder = os.path.join(os.getcwd(),"conf") 

try:
    for conf in os.listdir("conf"):
        if not conf.endswith(".json"):
            continue

        simname = os.path.splitext(os.path.basename(conf))[0]
        confpath = os.path.join(os.getcwd(),"conf",conf) 
        outputfolder = os.path.join(os.getcwd(),simname) 

        subprocess.check_call([executable, op, confpath, outputfolder], cwd=workdir)
except KeyboardInterrupt:
    print("PY: Aborting execution")
