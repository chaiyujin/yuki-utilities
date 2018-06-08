# Yuki utilities
Following modules are included:
- math
    - types extened from basic types: int2, int3, float2, float3, double2, double3, rgba
    - useful types: range, shape
    - 3 dim tensor
- image
    - color map, for feature/data visualization
    - draw functions (by **OpenCV**)
- audio
    - simple wav file reader
    - features:
        - power spectrum, mel-filterbank, mfcc
        - lpc [reference](http://cs.haifa.ac.il/~nimrod/Compression/Speech/S4LinearPredictionCoding2009.pdf)
- ffmpeg wrapper (to-do)
    - media reader
    - media writer
    - media transcode
    - custom codec and extented format ".mkv" for **depth stream**.

- statistic 3d face model (to-do)

# Dependence
- Eigen, **REQUIRED**
- FFTW3, **REQUIRED** to calculate audio features
- nanogui, **OPTIONAL** for gui
