import numpy as np
from PIL import Image
from skimage.feature import peak_local_max
import matplotlib.pyplot as plt
from skimage.transform import resize
from skimage import exposure
import json


def png_to_numpy(image_path, normalize=False):
    np_frame = 1 - np.array(Image.open(image_path).convert('L')).astype(np.float32) / 256.0
    if normalize:
        np_frame = np_frame / max(np.max(np_frame), 0.2)
        p2, p98 = np.percentile(np_frame, (2, 98))
        np_frame = exposure.rescale_intensity(np_frame, in_range=(p2, p98))
    return np_frame


def pad_image(np_image):
    pad_ax0 = 8 - (np_image.shape[0] % 8)
    pad_ax1 = 8 - (np_image.shape[1] % 8)
    np_image = np.pad(np_image, pad_width=[(0, pad_ax0), (0, pad_ax1)], constant_values=0)
    return np_image


def get_coords_from_heatmap(heatmap, thr=0.5):
    coords = peak_local_max(heatmap, min_distance=1, threshold_abs=thr)
    # transform coordinates to be typical (x,y) coords
    coords = np.array([
        [point[1], point[0]] for point in coords
    ])
    if len(coords)==0:
        coords = np.array([[0,0],[1,1]])
    return coords


def findclosestdistance(point, otherpoints):
    distances = [np.sqrt(np.sum((point - p) ** 2)) for p in otherpoints]
    return np.min(distances)


class SketchyImage(object):
    def __init__(self, np_image,
                 gt_heatmap=None,
                 gt_points_coords=None,
                 gt_labels=None,
                 gt_class_heatmaps=None,
                 gt_class_points=None,
                 pred_heatmap=None,
                 pred_points=None,
                 pred_labels=None,
                 pred_class_heatmaps=None,
                 pred_class_points=None,
                 heatmap_threshold=0.5
                 ):
        self.np_image = np_image  # stores np array representing the image

        self.gt_heatmap = gt_heatmap  # stores ground truth heatmap of labeled points
        self.gt_points = gt_points_coords  # stores ground truth coordinates of labeled points
        self.gt_labels = gt_labels  # stores gt labels of points in 1-3-5 format, 1 - endpoint, 3 - junction, 5 - sharp
        self.gt_class_heatmaps = gt_class_heatmaps  # stores 3 gt heatmaps of point classification
        self.gt_class_points = gt_class_points  # stores 3  gt point coordinates extracted from classification heatmaps

        self.pred_heatmap = pred_heatmap  # stores predicted heatmap for point positions
        self.pred_points = pred_points  # stores point coordinates extracted from predicted heatmap
        self.pred_labels = pred_labels  # stores labels of predicted points in 1-3-5 format
        self.pred_class_heatmaps = pred_class_heatmaps  # stores 3 predicted heatmaps for point classification
        self.pred_class_points = pred_class_points  # stores 3 point coordinates extracted from classification heatmaps

        self.heatmap_threshold = heatmap_threshold  # threshold for heatmap-to-coordinates extraction

    @classmethod
    def from_image(cls, image_path, normalize=False):
        np_frame = png_to_numpy(image_path, normalize=normalize)
        return cls(np_image=pad_image(np_frame))

    @classmethod
    def from_labeled_sample_pts(cls,
                                path_to_image="davinci_inputs/da5_gt.png",
                                path_to_points="davinci_inputs/da5_gt.png.pts",
                                has_class_info=False):
        image = png_to_numpy(path_to_image, normalize=True)
        with open(path_to_points) as fpts:
            n = int(next(fpts))
            lines = [[float(x) for x in line.split()] for line in fpts]
            points = [[line[1], line[0]] for line in lines]
            if len(points) != n:
                print("Something went wrong: number of points != indicated; check your file")
            points = np.array(points)
            if has_class_info:
                point_types = [line[2] for line in lines]
                point_types = np.array(point_types)
        limit_size = 1024
        if np.max(image.shape) > limit_size:
            # print("File may be too big for CUDA, decrease size")
            maxdim = np.max(image.shape)
            scaling_coeff = limit_size / maxdim
            newdim0, newdim1 = int(scaling_coeff * image.shape[0]), int(scaling_coeff * image.shape[1])
            points *= scaling_coeff
            image = resize(image, (newdim0, newdim1), order=3)
        pad_ax0 = 8 - (image.shape[0] % 8)
        pad_ax1 = 8 - (image.shape[1] % 8)
        image = np.pad(image, pad_width=[(0, pad_ax0), (0, pad_ax1)], constant_values=0)
        if has_class_info:
            return cls(np_image=image, gt_points_coords=points, gt_labels=point_types)
        else:
            return cls(np_image=image, gt_points_coords=points)

    @classmethod
    def from_labeled_sample_json(cls, imagefile, jsonfile):
        image = pad_image(png_to_numpy(imagefile, normalize=True))
        points_extracted = list()
        point_classes = list()
        with open(jsonfile) as f:
            load_dict = json.load(f)

        for shape in load_dict['shapes']:
            points_extracted.append(shape['points'][0])
            point_classes.append(int(shape['label']))

        points_extracted = np.array(points_extracted)
        point_classes = np.array(point_classes)
        return cls(np_image=image, gt_points_coords=points_extracted, gt_labels=point_classes)

    @classmethod
    def from_training_heatmaps(cls, imagetensor, gt_heatmap_tensor, gt_class_heatmaps_tensor,
                               pred_heatmap_tensor, pred_class_heatmaps_tensor):
        image, gt_heatmap, gt_class_heatmaps, pred_heatmap, pred_class_heatmaps = [x.detach().cpu().numpy() for x in
                                                                                  [imagetensor, gt_heatmap_tensor, gt_class_heatmaps_tensor, pred_heatmap_tensor, pred_class_heatmaps_tensor]]
        return cls(np_image=image, gt_heatmap=gt_heatmap, gt_class_heatmaps=gt_class_heatmaps, pred_heatmap=pred_heatmap, pred_class_heatmaps=pred_class_heatmaps)

    def set_pred_hm(self, pred_heatmap, pred_class_heatmaps):
        self.pred_heatmap = pred_heatmap
        self.pred_class_heatmaps = pred_class_heatmaps
        self.compute_coords_from_pred_heatmaps()
        self.compute_point_labels_from_heatmap()

    def compute_coords_from_pred_heatmaps(self):
        if self.pred_heatmap is not None:
            self.pred_points = get_coords_from_heatmap(self.pred_heatmap, thr=self.heatmap_threshold)

        if self.pred_class_heatmaps is not None:
            self.pred_class_points = list()
            for hm in self.pred_class_heatmaps:
                self.pred_class_points.append(get_coords_from_heatmap(heatmap=hm, thr=self.heatmap_threshold))

    def compute_coords_from_gt_heatmaps(self):
        if self.gt_heatmap is not None:
            self.gt_points = get_coords_from_heatmap(self.gt_heatmap, thr=self.heatmap_threshold)

        if self.gt_class_heatmaps is not None:
            self.gt_class_points = list()
            for hm in self.gt_class_heatmaps:
                self.gt_class_points.append(get_coords_from_heatmap(heatmap=hm, thr=self.heatmap_threshold))

    def compute_gt_point_labels_from_heatmaps(self):
        self.gt_labels = list()
        for point in self.gt_points:
            point_type = np.argmax([self.gt_class_heatmaps[i][point[1], point[0]] for i in range(3)])
            if point_type == 0:
                # This is endpoint
                self.gt_labels.append(1)
            if point_type == 1:
                # This is sharp
                self.gt_labels.append(5)
            if point_type == 2:
                # This is intersection
                self.gt_labels.append(3)
        self.gt_labels = np.array(self.gt_labels)

    def compute_point_labels_from_heatmap(self):
        """
        extract labels associated with self.pred_points by taking argmax of 3 predicted classification heatmaps
        """
        assert (self.pred_points is not None) and (self.pred_class_heatmaps is not None)
        self.pred_labels = list()
        for point in self.pred_points:
            point_type = np.argmax([self.pred_class_heatmaps[i][point[1], point[0]] for i in range(3)])
            if point_type == 0:
                # This is endpoint
                self.pred_labels.append(1)
            if point_type == 1:
                # This is sharp
                self.pred_labels.append(5)
            if point_type == 2:
                # This is intersection
                self.pred_labels.append(3)
        self.pred_labels = np.array(self.pred_labels)

    def plot_gt(self, legend=False):
        plt.imshow(self.np_image, cmap='gray_r')
        plt.clim(0, 1)
        if self.gt_points is not None:
            if self.gt_labels is not None:
                endpoints = self.gt_points[self.gt_labels == 1]
                sharppoints = self.gt_points[self.gt_labels == 5]
                junctions = self.gt_points[self.gt_labels == 3]
                plt.scatter(endpoints[:, 0], endpoints[:, 1], marker='+', c='lime', label='gt end', alpha=0.7, linewidths=1)
                plt.scatter(sharppoints[:, 0], sharppoints[:, 1], marker='+', c='magenta', label='gt sharp', alpha=0.7, linewidths=1)
                plt.scatter(junctions[:, 0], junctions[:, 1], marker='+', c='cyan', label='gt inter', alpha=0.7, linewidths=1)
            else:
                plt.scatter(self.gt_points[:, 0], self.gt_points[:, 1], c='red', marker='+', alpha=0.7, linewidths=0.5,
                            label='gt')
        if legend:
            plt.legend()
        plt.axis('off')

    def plot_result(self):
        plt.imshow(self.np_image, cmap='gray_r')
        plt.clim(0, 1)
        if not (self.pred_heatmap is None):
            hm = self.pred_heatmap.copy()
            hm[hm < self.heatmap_threshold] = np.nan
            plt.imshow(hm, cmap='hot_r', alpha=0.5)
            plt.clim(0, 1)
        if self.gt_points is not None:
            plt.scatter(self.gt_points[:, 0], self.gt_points[:, 1], c='red', marker='+', alpha=0.7, linewidths=0.5,
                        label='gt')
        if self.pred_points is not None:
            endpoints = self.pred_points[self.pred_labels == 1]
            sharppoints = self.pred_points[self.pred_labels == 5]
            junctions = self.pred_points[self.pred_labels == 3]
            plt.scatter(endpoints[:, 0], endpoints[:, 1], marker='+', c='lime', label='end', linewidths=0.3)
            plt.scatter(sharppoints[:, 0], sharppoints[:, 1], marker='+', c='magenta', label='sharp', linewidths=0.3)
            plt.scatter(junctions[:, 0], junctions[:, 1], marker='+', c='cyan', label='inter', linewidths=0.3)
        plt.legend()
        plt.axis('off')

    def plot_class_result(self, legend=False, heatmap_only=False):
        if not heatmap_only:
            plt.imshow(self.np_image, cmap='gray_r', interpolation="nearest",)
            plt.clim(0, 1)
        if not (self.pred_class_heatmaps[0] is None):
            hm = self.pred_class_heatmaps[0].copy()
            hm[hm < self.heatmap_threshold] = np.nan
            plt.imshow(hm, cmap='hot_r', alpha=0.5, interpolation="nearest",)
            # plt.clim(0, 1)
            # plt.colorbar()
            # print(hm)
        else:
            print("is none")
        if not heatmap_only:
            if self.gt_points is not None:
                plt.scatter(self.gt_points[:, 0], self.gt_points[:, 1], c='red', marker='+', alpha=0.7, linewidths=0.5, label='gt')
            if self.pred_class_points is not None:
                if len(self.pred_class_points[0]) > 0:
                    plt.scatter(self.pred_class_points[0][:, 0], self.pred_class_points[0][:, 1], marker='+', c='lime', alpha=0.7, label='end', linewidths=1)
                if len(self.pred_class_points[1]) > 0:
                    plt.scatter(self.pred_class_points[1][:, 0], self.pred_class_points[1][:, 1], marker='+', c='magenta', alpha=0.7, label='sharp', linewidths=1)
                if len(self.pred_class_points[2]) > 0:
                    plt.scatter(self.pred_class_points[2][:, 0], self.pred_class_points[2][:, 1], marker='+', c='cyan', alpha=0.7, label='inter', linewidths=1)
            else:
                if self.pred_points is not None:
                    plt.scatter(self.pred_points[:, 0], self.pred_points[:, 1], marker='+', c='blue', alpha=0.7,
                                linewidths=0.5, label='preds')
        if legend:
            plt.legend()
        plt.axis('off')

    def save_pred_points(self, filepath):
        """
        saves point coords in .pts format
        """
        with open(filepath, "w") as f:
            f.write(f"{len(self.pred_points)} \n")
            for coo in self.pred_points:
                f.write(f"{coo[1]} {coo[0]}\n")

    def save_pred_points_and_labels(self, filepath):
        """
        saves point coords and class data in .pts format
        """
        assert len(self.pred_points) == len(self.pred_labels)
        with open(filepath, "w") as f:
            f.write(f"{len(self.pred_points)} \n")
            for i in range(len(self.pred_points)):
                coo = self.pred_points[i]
                f.write(f"{coo[1]} {coo[0]} {self.pred_labels[i]}\n")

    def save_gt_points_and_labels(self, filepath):
        assert len(self.gt_points) == len(self.gt_labels)
        with open(filepath, "w") as f:
            f.write(f"{len(self.gt_points)} \n")
            for i in range(len(self.gt_points)):
                coo = self.gt_points[i]
                f.write(f"{coo[1]} {coo[0]} {self.gt_labels[i]}\n")

    def compute_f1_score(self, ground_truth_points_coords, predicted_points_coords,
                         closestpoint_relative_threshold=0.01, return_all_stats=False):
        closestpoint_THRESHOLD = self.np_image.shape[0] * closestpoint_relative_threshold
        if len(predicted_points_coords) > 600:
            print("Too many points predicted, score is 0")
            if return_all_stats:
                return 0, 0, 0
            return 0
        if len(predicted_points_coords) == 0:
            print("No points predicted, score is 0")
            if return_all_stats:
                return 0, 0, 0
            return 0
        if len(ground_truth_points_coords) == 0:
            print("No GT points here, score is 0.5")
            if return_all_stats:
                return 0.5, 0.5, 0.5
            return 0.5
        true_point_distances = [findclosestdistance(point, predicted_points_coords) for point in ground_truth_points_coords]
        true_point_found = [x < closestpoint_THRESHOLD for x in true_point_distances]
        pred_point_distances = [findclosestdistance(point, ground_truth_points_coords) for point in predicted_points_coords]
        pred_point_good = [x < closestpoint_THRESHOLD for x in pred_point_distances]
        tp = np.sum(true_point_found)
        fn = len(ground_truth_points_coords) - tp
        fp = len(predicted_points_coords) - np.sum(pred_point_good)
        if tp == 0:
            if return_all_stats:
                return 0, 0, 0
            return 0
        precision = tp / (tp + fp)
        recall = tp / (tp + fn)
        f1 = 2 * precision * recall / (precision + recall)
        if return_all_stats:
            return f1, precision, recall
        return f1

    def compute_score(self, closestpoint_relative_threshold=0.01, return_all_stats=False):
        return self.compute_f1_score(ground_truth_points_coords=self.gt_points,
                                     predicted_points_coords=self.pred_points,
                                     closestpoint_relative_threshold=closestpoint_relative_threshold,
                                     return_all_stats=return_all_stats)

    def compute_classification_score(self, closestpoint_relative_threshold=0.01):
        """
        Computes 3 f1 scores for each point type
        """
        point_types = [1,3,5]
        f1scores = list()
        for mytype in point_types:
            gt_typepoints = self.gt_points[self.gt_labels==mytype]
            pred_typepoints = self.pred_points[self.pred_labels==mytype]
            f1 = self.compute_f1_score(ground_truth_points_coords=gt_typepoints,
                                       predicted_points_coords=pred_typepoints,
                                       closestpoint_relative_threshold=closestpoint_relative_threshold,
                                       return_all_stats=False)
            f1scores.append(f1)
        return f1scores
    
    @staticmethod
    def save_heatmap(heatmap, note="test"):
        np.savetxt(fname=f"logs/{note}.txt", X=heatmap)
        plt.imsave(fname=f"logs/{note}.png", arr=heatmap, cmap="Reds")


if __name__=="__main__":
    si = SketchyImage.from_labeled_sample_pts(path_to_image="my_labeled/labels_with_class/msfry.png",
                                              path_to_points="my_labeled/labels_with_class/msfry.png.pts",
                                              has_class_info=True)
    si.plot_gt()
    plt.show()
