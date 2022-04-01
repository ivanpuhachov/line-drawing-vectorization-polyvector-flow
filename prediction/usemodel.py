import torch
import numpy as np
from layers import BaseHeatmapModel, CenterNet, CenterNet_classification, PyramidalNet
from PIL import Image
from skimage.feature import peak_local_max
import argparse
import os
import matplotlib.pyplot as plt
from datetime import datetime
from skimage.transform import resize
from skimage import exposure
from sketchyimage import SketchyImage


def load_checkpoint(filepath):
    """
    Loads model from a checkpoint
    """
    checkpoint = torch.load(filepath)
    loadmodel = PyramidalNet()
    loadmodel.load_state_dict(checkpoint['model_state_dict'])
    return loadmodel


def process_image(mymodel: BaseHeatmapModel, sketchyimage):
    """
    Loads object of SketchyImage and updates its fields by predicting from model
    """
    mymodel.eval()
    maxtile = 720
    if max(sketchyimage.np_image.shape) > maxtile:
        print("Big image, using tiling!")
        hm = np.zeros_like(sketchyimage.np_image)
        classes_hms = [np.zeros_like(sketchyimage.np_image) for _ in range(3)]
        for i in range(int(np.ceil(sketchyimage.np_image.shape[0] / maxtile))):
            for j in range(int(np.ceil(sketchyimage.np_image.shape[1] / maxtile))):
                cropped_image = sketchyimage.np_image[maxtile*i:maxtile*(i+1), maxtile*j:maxtile*(j+1)]
                crop_hm, crop_cls = model_predict(mymodel=mymodel, img=cropped_image)
                hm[maxtile*i:maxtile*(i+1), maxtile*j:maxtile*(j+1)] = crop_hm
                for c in range(3):
                    classes_hms[c][maxtile*i:maxtile*(i+1), maxtile*j:maxtile*(j+1)] = crop_cls[c]
    else:
        hm, classes_hms = model_predict(mymodel=mymodel, img=sketchyimage.np_image)
    sketchyimage.set_pred_hm(pred_heatmap=hm, pred_class_heatmaps=classes_hms)
    return sketchyimage


def model_predict(mymodel: BaseHeatmapModel, img):
    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    mymodel.eval()
    a = torch.from_numpy(img).unsqueeze(0).unsqueeze(0).to(device)
    if mymodel.can_classify:
        with torch.no_grad():
            outsideheads, classes = mymodel.forward_class_heatmaps(a)
            return outsideheads[-1].detach().cpu().squeeze(0).squeeze(0).numpy(), \
                   [x.detach().cpu().squeeze(0).squeeze(0).numpy() for x in classes]
    with torch.no_grad():
        pred = mymodel.forward_point_heatmap(a)
    heatmap = pred.detach().cpu().squeeze(0).squeeze(0).numpy()
    return heatmap


def predict_image(image_path: str, model_centernet, outfile: str, threshold: float, do_vis: bool):
    assert (os.path.exists(image_path))

    model_centernet.eval()

    sketchyimage = SketchyImage.from_image(image_path, normalize=True)

    process_image(mymodel=model_centernet, sketchyimage=sketchyimage)

    if do_vis:
        # save predicted points
        if not os.path.isdir("predicted_points/"):
            os.mkdir("predicted_points/")
        filename = image_path[image_path.rfind("/") + 1:-4]

        plt.figure(figsize=(10, 10))
        sketchyimage.plot_result()
        plt.savefig("predicted_points/" + filename + "_pred.png", bbox_inches='tight', dpi=600)
        plt.show()
        plt.close()

    sketchyimage.save_pred_points_and_labels(outfile)



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run a model from checkpoint")
    parser.add_argument("--input", "-i", default='example.png')
    parser.add_argument("--model", "-m", default="best_model_checkpoint.pth")
    parser.add_argument("--output", "-o", default="example.pts")
    parser.add_argument("--thr", "-t", default=0.5, help="threshold")
    parser.add_argument("--vis", action="store_true", default=False)
    args = parser.parse_args()

    assert (os.path.exists(args.model))

    device = 'cuda' if torch.cuda.is_available() else 'cpu'
    model = load_checkpoint(args.model).to(device)

    predict_image(image_path=args.input, model_centernet=model, outfile=args.output,
                  threshold=args.thr, do_vis=args.vis)
    # run_folder(folder_path=args.inputdir, model_centernet=model, outfolder=args.outputdir, threshold=args.thr)
    # process_labeled_folder(model, path_to_folder='benchmark_selected/', output_dir='testdir/')
    # stats = compute_stats_on_labeled(model, '/home/ivan/projects/svg_detector/my_labeled/labels/')
    # stats = compute_stats_on_labeled(model, '/home/ivan/projects/sketch-dataset/mergedataset/labels/')
    # print(f"Mean F1: \t\t {np.mean(stats[:, 0])}")
