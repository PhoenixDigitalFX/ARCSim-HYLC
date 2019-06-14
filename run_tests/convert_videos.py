import os
import sys
import subprocess
import argparse

# workdir = os.path.join(os.getcwd(),"..") 

try:
    for anim in os.listdir("anims"):
        print("Anim", anim)
        inputfmt = "anims/%s/" % anim + "%4d.jpg"
        output = "vids/%s.mp4" % anim

        # command = "ffmpeg -r 25 -start_number 0 -i " + inputfmt + " -vcodec libx264 -pix_fmt yuv420p" + output
        command = "ffmpeg -r 25 -i " + inputfmt + " " + output
        # command = "ffmpeg -r 50 -i " + inputfmt + " " + output
        subprocess.check_call(command, shell=True)#, cwd=workdir)
except KeyboardInterrupt:
    print("PY: Aborting execution")
