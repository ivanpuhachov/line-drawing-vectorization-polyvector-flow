# Sketch keypoint handlabeled dataset
We are publishing public URLs to images we labeled manually. The labels are available using these links: 
- [part 1](https://drive.google.com/file/d/1KdfigmwfwmumP0hJaUMHMm1XpmbBbzms/view?usp=sharing)
- [part 2](https://drive.google.com/file/d/1ygrdA9JzYH-xiDu1HOPetNGhMXgvb5fG/view?usp=sharing)
- [part 3](https://drive.google.com/file/d/1sFQF33nJSwlCxrhsrbSK20oO-49OfrOw/view?usp=sharing)


Run 
```bash
bash download1.sh
```


## Convention
Class `1` - endpoints
Class `3` - junctions and intersections
Class `5` - sharp points


****

# Labelme - labeling tool

Link: https://github.com/wkentaro/labelme

Conda installation:
```bash
conda create --name=labelme python=3.6
conda activate labelme
pip install labelme
```

## Usage instructions
0. To run GUI: `labelme`
1. Open image
2. Create points (top bar `Edit` -> `Create Point`)
3. Click on point, put label according to convention
4. Save `.json` file here
5. Run `labelme_to_pts.py` to transform labelme json file to our format


### Change shortcuts

After the first usage, `~/.labelmerc` is created. Go there, line 91, set `create_point: Ctrl+[`