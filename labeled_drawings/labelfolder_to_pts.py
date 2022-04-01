import json
import argparse
import os
import shutil

parser = argparse.ArgumentParser(description="Transform .json labels from folder to our .pts format")
parser.add_argument("--input", "-i", default='mergedataset/')
args = parser.parse_args()

inpdir = args.input
assert(os.path.isdir(inpdir))
if not inpdir.endswith("/"):
    inpdir += "/"

outdir = inpdir + "labels/"
if not (os.path.isdir(outdir)):
    os.mkdir(outdir)

files_in_folder = os.listdir(inpdir)
jsonfiles_in_folder = [f for f in files_in_folder if f.endswith('.json')]
pngfiles_in_folder = [f for f in files_in_folder if f.endswith('.png')]

for jsonfile in jsonfiles_in_folder:
    outfile = outdir + jsonfile.replace(".json", ".png.pts")

    with open(inpdir + jsonfile) as f:
        load_dict = json.load(f)

    points_extracted = list()

    for shape in load_dict['shapes']:
        points_extracted.append(shape['points'][0])

    with open(outfile, "w") as f:
        f.write(str(len(points_extracted)) + "\n")
        for a in points_extracted:
            f.write("{:2.1f} {:2.1f}\n".format(a[1], a[0]))
    print(f"done with {jsonfile}")

    pngfile = jsonfile.replace('.json', '.png')
    if pngfile in pngfiles_in_folder:
        shutil.copy2(inpdir+pngfile, outdir)
    else:
        print(f" - No png was found for {jsonfile}, copy it manually!")