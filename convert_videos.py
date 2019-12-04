import os, sys, subprocess, argparse

ap = argparse.ArgumentParser()
ap.add_argument('folder', help='input folder with anims/ subfolder')

ap.add_argument("-i", "--input", default="anims",
                help="input subfolder")
ap.add_argument("-o", "--output", default="vids",
                help="output subfolder")
ap.add_argument("-r", "--fps", default=30, type=int,
                help="fps") # NOTE: 30 for HYLC and 25 for YLC
ap.add_argument("-f", "--filter", default=[], nargs="*",
                help="filter to include")
ap.add_argument("-F", "--filterex", default=[], nargs="*",
                help="filter to exclude")
args, unknownargs = ap.parse_known_args()
args = vars(args)

infolder = os.path.join(args['folder'],args['input'])
outfolder = os.path.join(args['folder'],args['output'])
skip_existing = True

os.makedirs(outfolder, exist_ok=True)
try:
    for anim in sorted(os.listdir(infolder)):
        print("Anim", anim)
        inputfmt = os.path.join(os.path.join(infolder,'%s' % anim), "%4d.jpg")
        output = os.path.join(outfolder,'%s.mp4' % anim)

        # exclusive filter
        if any([fil in anim for fil in args['filterex']]):
            continue
        # inclusive filter
        if not all([fil in anim for fil in args['filter']]):
            continue

        if skip_existing and os.path.isfile(output):
            continue

        # command = "ffmpeg -r 25 -start_number 0 -i " + inputfmt + " -vcodec libx264 -pix_fmt yuv420p" + output
        # command = "ffmpeg -r 25 -i " + inputfmt + " " + output # YLC
        # command = "ffmpeg -r 30 -i " + inputfmt + " " + output # HYLC
        command = "ffmpeg -i %s -vcodec libx264 -pix_fmt yuv420p -vf colormatrix=bt601:bt709 -r %d %s" % (inputfmt, args['fps'], output)
        subprocess.check_call(command, shell=True)#, cwd=workdir)
except KeyboardInterrupt:
    print("PY: Aborting execution")

# command = 'ffmpeg -i %s -filter_complex "[0]split=2[bg][fg];[bg]drawbox=c=white@1:replace=1:t=fill[bg]; [bg][fg]overlay=format=auto" -vcodec libx264 -pix_fmt yuv420p  -r %d %s' % (inputfmt, fps, output)
