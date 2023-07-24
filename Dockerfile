FROM gurobi/optimizer:9.1.2

WORKDIR /opt/build

# Set opencv version and disable gui
ENV OPENCV_VERSION="4.7.0" \
    QT_QPA_PLATFORM=offscreen

# Install dependencies (boost, openmp, qt, opencv, eigen, numpy)
RUN apt-get -qq update \
    && apt-get -qq install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        wget \
        unzip \
        yasm \
        pkg-config \
        python3-dev \
        libswscale-dev \
        libtbb-dev \
        libjpeg-dev \
        libpng-dev \
        libtiff-dev \
        libopenjp2-7-dev \
        libavformat-dev \
        libpq-dev \
        libz-dev \
        libbz2-dev \
        libatlas-base-dev \
        qtbase5-dev \
        libboost-all-dev \
        libeigen3-dev \
        libomp-dev \
	time \
    && pip install numpy \
    && wget -q https://github.com/opencv/opencv/archive/refs/tags/${OPENCV_VERSION}.zip -O opencv.zip \
    && unzip -qq opencv.zip -d /opt \
    && rm -rf opencv.zip \
    && cmake \
        -D BUILD_TIFF=ON \
        -D BUILD_opencv_java=OFF \
        -D WITH_CUDA=ON \
        -D WITH_OPENGL=ON \
        -D WITH_OPENCL=ON \
        -D WITH_IPP=ON \
        -D WITH_TBB=ON \
        -D WITH_EIGEN=ON \
        -D WITH_V4L=ON \
        -D BUILD_TESTS=OFF \
        -D BUILD_PERF_TESTS=OFF \
        -D CMAKE_BUILD_TYPE=RELEASE \
        -D CMAKE_INSTALL_PREFIX=$(python3.11 -c "import sys; print(sys.prefix)") \
        -D PYTHON_EXECUTABLE=$(which python3.11) \
        -D PYTHON_INCLUDE_DIR=$(python3.11 -c "from distutils.sysconfig import get_python_inc; print(get_python_inc())") \
        -D PYTHON_PACKAGES_PATH=$(python3.11 -c "from distutils.sysconfig import get_python_lib; print(get_python_lib())") \
        /opt/opencv-${OPENCV_VERSION} \
    && make -j$(nproc) \
    && make install \
    && rm -rf /opt/build/* \
    && rm -rf /opt/opencv-${OPENCV_VERSION} \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get -qq autoremove \
    && apt-get -qq clean


COPY . /workspace/line-drawing-vectorization-polyvector-flow
WORKDIR /workspace/line-drawing-vectorization-polyvector-flow
# Install python dependencies
RUN pip install -U pip wheel \
    && pip install -r requirements.txt \
    && pip cache purge

# Install paal
RUN wget http://paal.mimuw.edu.pl/paal.zip && unzip paal.zip && mv home/paal/sources/paal/include/paal . && rm -rf home

# Fix for: https://support.gurobi.com/hc/en-us/articles/360039093112-How-do-I-resolve-undefined-reference-errors-while-linking-Gurobi-in-C-
RUN ln -s /opt/gurobi /opt/gurobi911 \
    && cd /opt/gurobi/linux64/src/build \
    && make \
    && cp /opt/gurobi/linux64/src/build/libgurobi_c++.a /opt/gurobi/linux64/lib/

# Build the vectorization algorithm
RUN mkdir build 
WORKDIR /workspace/line-drawing-vectorization-polyvector-flow/build
RUN cmake -DCMAKE_BUILD_TYPE=Release .. && make -j$(nproc)
WORKDIR /workspace/line-drawing-vectorization-polyvector-flow

STOPSIGNAL SIGINT
ENTRYPOINT ["/bin/bash", "run_vectorization.sh"]
