import os, sys, subprocess, argparse, platform
assert(platform.python_version().startswith("3")) # cluster safety
from multiprocessing import Pool
# import tempfile, re
import getpass, smtplib
from email.mime.text import MIMEText

outputfolder = "sims"

def try_get_password():
    # timeoutable getpass hack ...
    import signal, time
    signal.signal(signal.SIGALRM, lambda signum, frame: 1/0)
    signal.alarm(60)
    try:
        pw = getpass.getpass()
        signal.alarm(0)
    except ZeroDivisionError:
        print("Timed out getting PW.")
        pw = ""
        signal.alarm(0)
    except:
        print("Aborted getting PW.")
        pw = ""
        signal.alarm(0)
    return pw

def try_send_mail(header, content, pw):
    usermail = 'gsperl@ist.ac.at'
    user = 'gsperl'
    msg = MIMEText(content)
    msg['Subject'] = header
    msg['From'] = usermail
    msg['To'] = usermail

    try:
        server = smtplib.SMTP('owa.ist.ac.at', 587)
        server.ehlo()
        server.starttls()
        server.login(user, pw)
        server.send_message(msg)
        print("Sent email.")
        server.close()
    except:
        print("Email error.")

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
ap.add_argument("-p", "--processes", default="1",
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

pw = try_get_password()[::-1]

# RUN
workdir = os.path.join(os.getcwd(),"..") 
executable = os.path.join(builddir, "bin", "arcsim_0.2.1")
# from workdir call
# ./build-Release/bin/arcsim_0.2.1 $op $configfile ${3}

op = "simulateoffline"

confs = sorted(os.listdir("conf"))
confs = ["shirt_%s.json" % s for s in [
    "rib",
    "basket",
    "satin",
    "honey",
    "stock",
]] # NOTE own order
conffolder = os.path.join(os.getcwd(),"conf") 
os.makedirs(os.path.join(os.getcwd(),outputfolder), exist_ok=True)

def tasks():
    for conf in confs:
        if not conf.endswith(".json"):
            continue
        simname = os.path.splitext(os.path.basename(conf))[0]
        confpath = os.path.join(os.getcwd(),"conf",conf) 
        outputdir = os.path.join(os.getcwd(),outputfolder,simname) 
        if os.path.isdir(outputdir):
            continue

        yield [executable, op, confpath, outputdir]
def execute_task(task):
    subprocess.check_call(task, cwd=workdir)

try:
    n_processes = int(args['processes'])
    if n_processes > 1:
        pool = Pool(n_processes)
        iterator = tasks()
        pool_it = pool.imap_unordered(execute_task, iterator, chunksize=1)
        # wait and iterate pool results, keyboard abortable
        try:
            n_finished = 0
            for _ in pool_it:
                n_finished += 1
                print("Finished simulation #%d." % n_finished)
        except KeyboardInterrupt:
            print("Caught KeyboardInterrupt, terminating workers")
            pool.terminate() # terminating 
        else:
            pool.close() # normal termination
    else:
        for task in tasks():
            execute_task(task)
    
except KeyboardInterrupt:
    print("PY: Aborting execution")
    pw = "" # avoid sending mail

if not pw == "" and not pw is None:
    try_send_mail("HYLC Finished", "The simulations have finished.", pw[::-1])