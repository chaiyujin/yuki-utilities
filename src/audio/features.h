#pragma once

#include "common.h"
#include "unsupported/Eigen/FFT"
#include <stdio.h>
#include <cfloat>

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(audio)

class Features
{
public:
    class Window
    {
    public:
        static AudioMask ones(int length);
    };

    static AudioSamples pre_emphasis(const AudioSamples &signal, double preemph);

    static std::vector<AudioFeature> power_spectrum(
        const AudioSamples &signal, int samplerate,
        double winlen=0.025, double winstep=0.01, int nfft=0,
        AudioMask(*func)(int length)=Window::ones);

    static std::vector<AudioFeature> mfcc(
        const AudioSamples &signal, int samplerate,
        double winlen=0.025, double winstep=0.01,
        int nfft=0, int numcep=13, int nfilt=26,
        double lowfreq=0, double highfreq=DBL_MAX,
        double preemph=0.97, int ceplifter=22,
        bool append_energy=true, AudioMask(*func)(int length)=Window::ones);

};

NAMESPACE_END(audio)
NAMESPACE_END(yuki)