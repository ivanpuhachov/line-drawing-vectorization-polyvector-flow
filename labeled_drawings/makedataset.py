import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from PIL import Image
import json
import shutil
import pandas as pd


def png_image_to_npy(path_to_image, normalize=False):
    np_frame = 1 - np.array(Image.open(path_to_image).convert('L')).astype(np.float32) / 256.0
    if normalize:
        np_frame = np_frame / max(np.max(np_frame), 0.2)
    # np_frame = resize(np_frame, output_shape=(np_frame.shape[0] // 2, np_frame.shape[1] // 2))
    return np_frame


def json_to_npy_sample(path_to_json):
    sample = dict()
    load_dict = json.load(open(path_to_json))
    # print(load_dict.keys())
    endpoints = list()
    intersectionpoints = list()
    sharppoints = list()
    for shape in load_dict['shapes']:
        # print(shape['label'])
        # print(shape.keys())
        coords = shape['points'][0]
        if shape['label'] == '1':
            endpoints.append(coords)
        if shape['label'] == '3':
            intersectionpoints.append(coords)
        if shape['label'] == '5':
            sharppoints.append(coords)
    sample['endpoints'] = np.array(endpoints)
    sample['intersectionpoints'] = np.array(intersectionpoints)
    sample['sharppoints'] = np.array(sharppoints)
    # print(sample)
    return sample


def test():
    img = png_image_to_npy("part2/w0001.png")
    sample = json_to_npy_sample("part2json/w0001.json")
    plt.figure(figsize=(6,6))
    plt.imshow(img, cmap='gray_r')
    for key in sample.keys():
        if key.endswith('points'):
            plt.scatter(sample[key][:,0], sample[key][:,1], label=key, marker='x')
    plt.axis("off")
    plt.show()
    plt.close()


def viewsample(n):
    files = os.listdir("mergedataset/")
    files = [f for f in files if f.endswith(".png")]
    selectedfile = files[n]
    selectedfile_json = selectedfile.replace("png", "json")
    img = png_image_to_npy(f"mergedataset/{selectedfile}")
    sample = json_to_npy_sample(f"mergedataset/{selectedfile_json}")
    plt.figure(figsize=(6,6))
    plt.imshow(img, cmap='gray_r')
    for key in sample.keys():
        if key.endswith('points'):
            if len(sample[key])==0:
                continue
            plt.scatter(sample[key][:,0], sample[key][:,1], label=key, marker='x')
    plt.axis("off")
    plt.title(selectedfile)
    plt.show()
    plt.close()


def copy_labeled_files(olddir, newdir):
    oldfiles = os.listdir(olddir)
    for pngfile in [f for f in oldfiles if f.endswith(".png")]:
        jsonfile = pngfile.replace(".png", ".json")
        if jsonfile in oldfiles:
            shutil.copy(olddir + pngfile, newdir)
            shutil.copy(olddir + jsonfile, newdir)


def load_to_dataset_item(pngfilepath, outdir):
    jsonfilepath = pngfilepath.replace(".png", ".json")
    pngfilename = pngfilepath[pngfilepath.rfind("/")+1:]
    npzfilename = pngfilename.replace(".png", ".npz")
    npzfilepath = outdir + npzfilename
    assert(os.path.exists(jsonfilepath))
    img = png_image_to_npy(pngfilepath)
    sample = json_to_npy_sample(jsonfilepath)
    sample['image'] = img
    with open(npzfilepath, "wb") as fl:
        np.savez_compressed(fl, **sample)


def view_dataset_item(itempath):
    with open(itempath, 'rb') as f:
        loaded = np.load(f)
        sample = dict()
        for key in loaded.keys():
            if len(loaded[key] > 0):
                sample[key] = loaded[key].astype(np.float32)

        plt.imshow(sample['image'], cmap='gray_r')
        for key in sample.keys():
            if key.endswith('points'):
                if len(sample[key])==0:
                    continue
                plt.scatter(sample[key][:,0], sample[key][:,1], label=key, marker='x')
        plt.axis("off")
        plt.title(itempath)
        plt.show()
        plt.close()


def create_dataset_csv(datset_folder):
    coll = [{"filename": f, "category": "random"}
           for f in os.listdir(datset_folder)
           if f.endswith(".npz")]
    df = pd.DataFrame(coll)
    print(df.head())
    df.to_csv(datset_folder + "dataset.csv", index=False)


# copy_labeled_files("part1/", "mergedataset/")
# copy_labeled_files('part2/', 'mergedataset/')
# copy_labeled_files('part3/', 'mergedataset/')

# for i in range(120):
#     viewsample(i)

# load_to_dataset_item("mergedataset/0002.png", "dataset/")

# for f in os.listdir("mergedataset/"):
#     if not (f.endswith("png")):
#         continue
#     print(f)
#     load_to_dataset_item("mergedataset/"+f, "dataset/")

# view_dataset_item("dataset/0001.npz")

create_dataset_csv("dataset/")