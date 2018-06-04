#pragma once

#include "common.h"
#include "unsupported/Eigen/FFT"

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(audio)

class Features
{
public:
    static Eigen::VectorXf mfcc(const Eigen::VectorXf &signal)
    {
        Eigen::FFT<float> fft;
        Eigen::VectorXcf freq;
        fft.fwd(freq, signal);
    }
};

NAMESPACE_END(audio)
NAMESPACE_END(yuki)