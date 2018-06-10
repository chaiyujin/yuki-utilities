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

    static double   default_winlen;
    static double   default_winstep;
    static int      default_nfft;
    static int      default_numcep_dct;
    static int      default_nfilt_mel;
    static int      default_ceplifter_mfcc;
    static double   default_preemph;
    static double   dafault_lowfreq;
    static double   dafault_highfreq;
    static AudioMask(*default_winfunc)(int);
    

    static double hz2mel(double hz) { return 2595.0 * std::log10(1 + hz / 700.0); }
    static double mel2hz(double mel) { return 700.0 * (std::pow(10, mel / 2595.0) - 1); }
    static Eigen::MatrixXd dct(const Eigen::MatrixXd &in, bool normalization=true);

    // return (nfilt, nfft / 2 + 1),  which are nfilt triangle filters.
    static Eigen::MatrixXd mel_filters(int nfilt, int nfft, int samplerate, double lowfreq, double highfreq);

    static AudioSamples pre_emphasis(const AudioSamples &signal, double preemph);

    static AudioFeatureList power_spectrum(
        const AudioSamples &signal, int samplerate,
        double  winlen      = default_winlen,
        double  winstep     = default_winstep,
        int     nfft        = default_nfft,
        AudioMask(*winfunc)(int length)=default_winfunc);

    static AudioFeatureList filter_bank(
        const AudioSamples &signal, int samplerate,
        double  winlen      = default_winlen,
        double  winstep     = default_winstep,
        int     nfft        = default_nfft,
        int     nfilt       = default_nfilt_mel,
        double  lowfreq     = dafault_lowfreq,
        double  highfreq    = dafault_highfreq,
        double  preemph     = default_preemph,
        AudioMask(*winfunc)(int length)=default_winfunc);

    static std::pair<AudioFeatureList, AudioFeatureList> filter_bank_energy(
        const AudioSamples &signal, int samplerate,
        double  winlen      = default_winlen,
        double  winstep     = default_winstep,
        int     nfft        = default_nfft,
        int     nfilt       = default_nfilt_mel,
        double  lowfreq     = dafault_lowfreq,
        double  highfreq    = dafault_highfreq,
        double  preemph     = default_preemph,
        AudioMask(*winfunc)(int length)=default_winfunc);

    static AudioFeatureList mfcc(
        const AudioSamples &signal, int samplerate,
        double  winlen      = default_winlen,
        double  winstep     = default_winstep,
        int     nfft        = default_nfft,
        int     numcep      = default_numcep_dct,
        int     nfilt       = default_nfilt_mel,
        double  lowfreq     = dafault_lowfreq,
        double  highfreq    = dafault_highfreq,
        double  preemph     = default_preemph,
        int     ceplifter   = default_ceplifter_mfcc,
        bool    append_energy=true,
        AudioMask(*winfunc)(int length)=default_winfunc);

    static AudioFeatureList autocorrelation(
        const std::vector<AudioSamples> &signal_list,
        bool biased_acorre_estimator=false);

    static AudioFeatureList lpc(
        const AudioSamples &signal, int samplerate, int order,
        double  winlen      = default_winlen,
        double  winstep     = default_winstep,
        double  preemph     = default_preemph,
        AudioMask(*winfunc)(int length)=default_winfunc);

    static std::vector<std::vector<double>> formants(
        const AudioFeatureList &lpc_A, int samplerate); 
};

NAMESPACE_END(audio)
NAMESPACE_END(yuki)