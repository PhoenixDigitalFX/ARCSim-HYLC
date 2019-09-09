import os
import sys
import subprocess
import argparse



outputfolder = "sims"

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
args, unknownargs = ap.parse_known_args()
args = vars(args)
args['build'] = args['build'] != "0"
args['run'] = args['run'] != "0"

# BUILD
sourcedir = os.path.join(os.getcwd(),"..") 
builddir = os.path.join(
    os.getcwd(), "..", "build-Release")



if args['build']:
    build(sourcedir=sourcedir, builddir=builddir, debug=False)

if not args['run']:
    exit()

# RUN
workdir = os.path.join(os.getcwd(),"..") 
executable = os.path.join(builddir, "bin", "arcsim_0.2.1")
# from workdir call
# ./build-Release/bin/arcsim_0.2.1 $op $configfile ${3}

op = "simulateoffline"

conffolder = os.path.join(os.getcwd(),"conf") 
os.makedirs(os.path.join(os.getcwd(),outputfolder), exist_ok=True)

try:
    for conf in sorted(os.listdir("conf")):
        if not conf.endswith(".json"):
            continue
        if "tmp" in conf:
            continue
        # if not ("stock" in conf or "satin_stretch" in conf):
        #     continue
        # if not "basket_stretch" in conf:
        #     continue
        simname = os.path.splitext(os.path.basename(conf))[0]
        confpath = os.path.join(os.getcwd(),"conf",conf) 
        outputdir = os.path.join(os.getcwd(),outputfolder,simname) 
        if os.path.isdir(outputdir):
            continue

        subprocess.check_call([executable, op, confpath, outputdir], cwd=workdir)
except KeyboardInterrupt:
    print("PY: Aborting execution")


# python3 run.py && gnome-session-quit --power-off --force --no-prompt
