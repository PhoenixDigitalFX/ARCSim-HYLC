import os
import sys
import subprocess
import argparse

# workdir = os.path.join(os.getcwd(),"..") 

folder1 = "/media/gsperl/Elements/HYLC/final_sim_copies/vids"
# folder1 = "/media/gsperl/Elements/HYLC/2019-06-24-hwrib/vids"
folder2 = "/home/gsperl/Projects/arcsim-hylc/run_tests/vids"
output = "vids_merged"
os.makedirs(output, exist_ok=True)

files1 = [os.path.join(folder1, f) for f in sorted(os.listdir(folder1))]
files2 = [os.path.join(folder2, f) for f in sorted(os.listdir(folder2))]


def namingfun(f):
    if "drapex" in f.lower():
        sim = "drapeX"
    elif "drapey" in f.lower():
        sim = "drapeY"
    elif "stretchx" in f.lower():
        sim = "stretchX"
    elif "stretchy" in f.lower():
        sim = "stretchY"
    else:
        sim = "other"

    if "stock" in f.lower():
        mat = "stock"
    elif "rib" in f.lower():
        mat = "rib"
    elif "basket" in f.lower():
        mat = "basket"
    else:
        mat = "other"

    if "cam_b" in f.lower():
        cam = "1"
    else:
        cam = "0"
    return "%s_%s_%s" % (mat, sim, cam)

try:
    for i,(f1,f2) in enumerate(zip(files1,files2)):

        name = namingfun(f1)
        # if not "rib" in name:
        #     continue

        print(name)
        outfile = os.path.join(output,"%s.mp4" % name)

        if os.path.isfile(outfile):
            continue

        # f2 = f1 # DEBUG
        command = " ".join([
            "ffmpeg",
            "-i %s -i %s" % (f1, f2),
            "-filter_complex",
            "'[0:v][1:v]hstack=inputs=2[vid]'",
            "-map [vid]",
            "-crf 20",
            outfile
        ])


        # # command = "ffmpeg -r 25 -start_number 0 -i " + inputfmt + " -vcodec libx264 -pix_fmt yuv420p" + output
        # command = "ffmpeg -r 25 -i " + inputfmt + " " + output # YLC
        # command = "ffmpeg -r 30 -i " + inputfmt + " " + output # HYLC
        # #-y

        subprocess.check_call(command, shell=True)#, cwd=workdir)
except KeyboardInterrupt:
    print("PY: Aborting execution")
