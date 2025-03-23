# Keypoint-Driven Line Drawing Vectorization via PolyVector Flow

This repository contains: 
 - C++ code with vectorization algorithm
 - Matlab code for debugging polyvector flow
 - Python code for keypoint inference

 * Dataset for training was released in a separate repository: [line-drawing-vectorization-polyvector-flow-dataset](https://github.com/ivanpuhachov/line-drawing-vectorization-polyvector-flow-dataset)

 * Paper webpage [link](https://puhachov.xyz/publications/keypoint-driven-polyvector-flow/)

Contact me if you need code for speed benchmarks. This release was tuned for better vectorization accuracy.

##### Citation
```
@article{Puhachov2021KeypointPolyvector,
    author = {Ivan Puhachov and William Neveu and Edward Chien and Mikhail Bessmeltsev},
    title = {Keypoint-Driven Line Drawing Vectorization via PolyVector Flow},
    journal = {ACM Transactions on Graphics (Proceedings of SIGGRAPH Asia)},
    volume = {40}, number = {6}, year = {2021}, month = dec,
    doi = {10.1145/3478513.3480529}
}
```

## Installation

### Cmake build
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

### Using Docker (alternative)
Requires [Docker](https://docs.docker.com/engine/install/) and [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html) (for GPU acceleration).

Retrieve a [Gurobi WLS license](https://www.gurobi.com/features/web-license-service/). This can be done by running:

```
grbgetkey
```

This will produce a `gurobi.lic` file. Place this into the repository directory.

To build the Docker image using docker compose:

```
docker compose build
```

Alternatively, build it using `docker build`:

```
docker build -t line-drawing-vectorization-polyvector-flow-gpu .
```

Then, run a Docker container based on the built image using docker compose:

```
docker compose run gpu inputs/example.png
```

Alternatively, run it using `docker run`:

```
docker run --rm -v ./gurobi.lic:/opt/gurobi/gurobi.lic:ro -v ./inputs:/workspace/line-drawing-vectorization-polyvector-flow/inputs --gpus all line-drawing-vectorization-polyvector-flow-gpu inputs/example.png
```

## Usage

Running `build/vectorize` will start vectorization on default input parameters: `inputs/example.png` and `inputs/example.pts`

It will output a debug window:

![](window.png)

To suppress debug window, provide output destination as third parameter:

> `build/vectorize image.png image.pts image.svg`

This bash script runs both keypoint prediciton and vectorization give a single image as input:

```
bash run_vectorization.sh inputs/example.png
```

***
## Keypoints prediction
Uses PyTorch (tested on 1.10) and CUDA. Can work with CPU, but expect longer inference time. 

To predict keypoints run:
```bash
 python prediction/usemodel.py --model prediction/best_model_checkpoint.pth --input image.png --output image.pts
```
***

### Illustrations Copyrights
* 'Dracolion', 'Mouse', 'Muten', 'Sheriff' are from [Noris et al. 2013]. Please ask the authors for the .pngs.

* 'Banana-Tree', 'Elephant', 'Hippo', 'Kitten', 'Puppy', 'Leaf', 'Dog06', 'Dog14': (c) Ivan Huska. https://www.easy-drawings-and-sketches.com/


***

## Build Instructions

### Gurobi
Tested with Gurobi 9.1.1, 9.0.3, 8.1.1. Install Gurobi and activate license, then update `GUROBI_HOME` in `CMakeLists.txt`:
```
set(GUROBI_HOME "/opt/gurobi911/linux64/")
```
Alternativelly, you can set your own path in `cmake/FindGUROBI.cmake`

### Qt
Instructions: [ubuntu](https://wiki.qt.io/Install_Qt_5_on_Ubuntu)

You may want to:
```bash
sudo apt-get install libz-dev
sudo apt-get install libbz2-dev
sudo apt install libatlas-base-dev
```

### Boost
https://www.boost.org/
```bash
sudo apt-get install libboost-all-dev
```

### OpenMP
```bash
sudo apt-get install libomp-dev
```

### Paal
http://paal.mimuw.edu.pl/

Then move `paal/include/paal/` folder here to `paal/`

### Gurobi
We use environment variable, don't forget to setup `GUROBI_HOME` and 

Linux: you may need to recompile (see https://support.gurobi.com/hc/en-us/articles/360039093112-How-do-I-resolve-undefined-reference-errors-while-linking-Gurobi-in-C-)

```bash
cd /opt/gurobi903/linux64/src/build
make
cp libgurobi_c++.a ../../lib/
```

### OpenCV

This one can be finicky. Here's what worked for me (William, Ubuntu 18.04).

Install the OpenCV library [this way](https://docs.opencv.org/master/d2/de6/tutorial_py_setup_in_ubuntu.html). Make sure to do ```sudo make install``` at the very end.

Make sure that you have a FindOpenCV.cmake file in polyvector_flow/cpp/cmake.

Make sure that the CMakelists.txt file in polyvector_flow/cpp specifies the path to the OpenCV static libraries (the .a ones). 
Example : ```find_package(OpenCV REQUIRED PATHS "/usr/lib/x86_64-linux-gnu")```

Building will take at least 30 minutes


### Eigen
`sudo apt install libeigen3-dev`

### Miscellaneous
* ["How do I resolve "undefined reference" errors while linking Gurobi in C++?"](https://support.gurobi.com/hc/en-us/articles/360039093112-How-do-I-resolve-undefined-reference-errors-while-linking-Gurobi-in-C-)
* `cmake ../` should finish with `-- Generating done`. It also outputs log info with paths, for example:
```
-- Qt found
-- Found OpenMP_C: -fopenmp (found version "4.5") 
-- Found OpenMP_CXX: -fopenmp (found version "4.5") 
-- Found OpenMP: TRUE (found version "4.5")  
-- Found OpenCV: /usr/local (found version "4.5.0") 
-- Found Boost: /usr/lib/x86_64-linux-gnu/cmake/Boost-1.71.0/BoostConfig.cmake (found suitable version "1.71.0", minimum required is "1.48")  
-- Boost_INCLUDE_DIRS = /usr/include
-- GUROBI_HOME = /opt/gurobi911/linux64/
-- Found GUROBI: /opt/gurobi911/linux64/lib/libgurobi_c++.a  
-- Gurobi found
-- Paal = /home/ivan/projects/polyvector_cpp/paal
-- Configuring done

```
