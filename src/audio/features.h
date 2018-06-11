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

    static int num_samples(
        int samplerate,
        int winnum,
        double  winlen      = default_winlen,
        double  winstep     = default_winstep)
    {
        int samples = winlen * samplerate;
        int step = winstep * samplerate;
        return (winnum - 1) * step + samples;
    }

    static int samples_aligned_to_feature(
        int samplerate, int feature_cols,
        double winstep     = default_winstep)
    { return feature_cols * samplerate * winstep; }

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
    
    static AudioFeatureList log_filter_bank(
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

    static AudioFeatureList lpc(
        const AudioSamples &signal, int samplerate, int order,
        double  winlen      = default_winlen,
        double  winstep     = default_winstep,
        double  preemph     = default_preemph,
        AudioMask(*winfunc)(int length)=default_winfunc);

    /*
     * formants estimation,
     * reference: https://www.mathworks.com/help/signal/ug/formant-estimation-with-lpc-coefficients.html
    */
    static std::vector<std::vector<double>> formants(
        const AudioFeatureList &lpc_A, int samplerate); 

};

NAMESPACE_END(audio)
NAMESPACE_END(yuki)