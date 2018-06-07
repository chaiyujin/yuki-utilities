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

    static double hz2mel(double hz) { return 2595.0 * std::log10(1 + hz / 700.0); }
    static double mel2hz(double mel) { return 700.0 * (std::pow(10, mel / 2595.0) - 1); }
    // return (nfilt, nfft / 2 + 1),  which are nfilt triangle filters.
    static Eigen::MatrixXd mel_filters(int nfilt, int nfft, int samplerate, double lowfreq, double highfreq);

    static AudioSamples pre_emphasis(const AudioSamples &signal, double preemph);

    static AudioFeatureList power_spectrum(
        const AudioSamples &signal, int samplerate,
        double winlen=0.025, double winstep=0.01, int nfft=0,
        AudioMask(*func)(int length)=Window::ones);

    static AudioFeatureList filter_bank(
        const AudioSamples &signal, int samplerate,
        double winlen=0.025, double winstep=0.01, int nfft=0,
        int nfilt=40, double lowfreq=0, double highfreq=DBL_MAX,
        double preemph=0.63, AudioMask(*func)(int length)=Window::ones);

    static std::pair<AudioFeatureList, AudioFeatureList> filter_bank_energy(
        const AudioSamples &signal, int samplerate,
        double winlen=0.025, double winstep=0.01, int nfft=0,
        int nfilt=40, double lowfreq=0, double highfreq=DBL_MAX,
        double preemph=0.63, AudioMask(*func)(int length)=Window::ones);

    static AudioFeatureList mfcc(
        const AudioSamples &signal, int samplerate,
        double winlen=0.025, double winstep=0.01,
        int nfft=0, int numcep=13, int nfilt=40,
        double lowfreq=0, double highfreq=DBL_MAX,
        double preemph=0.63, int ceplifter=22,
        bool append_energy=true, AudioMask(*func)(int length)=Window::ones);

};

NAMESPACE_END(audio)
NAMESPACE_END(yuki)