#pragma once
#include "common.h"

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(audio)

class DSP
{
public:
    static double energy2dB(double x) { return std::log10(x) * 20; }
    static double hz2mel(double hz) { return 2595.0 * std::log10(1 + hz / 700.0); }
    static double mel2hz(double mel) { return 700.0 * (std::pow(10, mel / 2595.0) - 1); }
    static Eigen::MatrixXd dct(const Eigen::MatrixXd &in, bool normalization=true);
    /*
     * return (nfilt, nfft / 2 + 1),  which are nfilt triangle filters.
    */
    static Eigen::MatrixXd mel_filters(
        int nfilt, int nfft, int samplerate,
        double lowfreq, double highfreq);
    /*
     * pre emphasis audio signal, if preemph <= 0, original signal will be returned.
    */
    static AudioSamples pre_emphasis(
        const AudioSamples &signal,
        double preemph);
    /*
     * input: a vector of signal frames, each frame has length `L`
     * return (`L` + 1, frames)
     * auto correlation, using ifft(fft(signal) ** 2) to calculate
    */
    static Eigen::MatrixXd autocorrelation(
        const std::vector<AudioSamples> &signal_list,
        bool biased_acorre_estimator=false);

    static std::pair<Eigen::VectorXd, Eigen::VectorXcd> freqz(
        const Eigen::VectorXd &b,
        const Eigen::VectorXd &a,
        int wor_n   = 512,
        bool whole  = false);
    
    static std::pair<Eigen::VectorXd, Eigen::VectorXcd> freqz(
        double b,
        const Eigen::VectorXd &a,
        int wor_n   = 512,
        bool whole  = false);

    static Eigen::MatrixX2d WH_to_FreqdB(
        const std::pair<Eigen::VectorXd, Eigen::VectorXcd> &WH,
        double frequency);
};

NAMESPACE_END(audio)
NAMESPACE_END(yuki)