import tempfile
import re

replacedict = {"name": "basket1"}


# open file conf
with open("conf/hylc_dress.json","r") as basefile:
    filestr = basefile.read()

# find and replace all expressions ${*}
exprs = set(re.findall("\$\{(.*?)\}",filestr))
print([e for e in exprs])
for e in exprs:
    if e in replacedict:
        filestr = filestr.replace("${%s}" % e,replacedict[e])

# create a new replaced file
with tempfile.NamedTemporaryFile(mode="w+") as tmp:
    tmp.write(filestr)
    tmp.flush()

    # use it
    with open(tmp.name,"r") as asdf:
        print(asdf.read())
