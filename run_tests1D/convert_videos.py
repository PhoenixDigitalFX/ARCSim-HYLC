import os
import sys
import subprocess
import argparse

# workdir = os.path.join(os.getcwd(),"..") 

os.makedirs("vids", exist_ok=True)
try:
    for anim in sorted(os.listdir("anims")):
        print("Anim", anim)
        inputfmt = "anims/%s/" % anim + "%4d.jpg"
        output = "vids/%s.mp4" % anim
        #os.makedirs("vids")
        # if not "basket_stretch" in anim:
        #     continue

        if os.path.isfile(output):
            continue
        # command = "ffmpeg -r 25 -start_number 0 -i " + inputfmt + " -vcodec libx264 -pix_fmt yuv420p" + output
        # command = "ffmpeg -r 25 -i " + inputfmt + " " + output # YLC
        command = "ffmpeg -r 30 -i " + inputfmt + " " + output # HYLC
        subprocess.check_call(command, shell=True)#, cwd=workdir)
except KeyboardInterrupt:
    print("PY: Aborting execution")
