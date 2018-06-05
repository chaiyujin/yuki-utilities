#include "features.h"
#include "fftw3.h"

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(audio)

AudioMask Features::Window::ones(int length)
{
    return AudioMask(length, 1.0);
}

AudioSamples Features::pre_emphasis(const AudioSamples &signal, double preemph)
{
    AudioSamples ret(signal.size());
    ret[0] = signal[0];
    for (size_t i = 1; i < signal.size(); ++i)
    {
        ret[i] = signal[i] - preemph * signal[i - 1];
    }
    return ret;
}

std::vector<AudioFeature> Features::power_spectrum(
    const AudioSamples &signal, int samplerate,
    double winlen, double winstep, int nfft,
    AudioMask(*winfunc)(int length))
{
    std::vector<AudioFeature> ret;
    int N = 0;
    int samples = winlen * samplerate;
    int step = winstep * samplerate;
    if (nfft > samples)
        N = nfft;
    else
        N = samples;

    auto win_mask = winfunc(samples);
    double *fftw_in = fftw_alloc_real(N);
    fftw_complex *fftw_out = fftw_alloc_complex(N / 2 + 1);
    auto p = fftw_plan_dft_r2c_1d(N, fftw_in, fftw_out, FFTW_ESTIMATE);
    memset(fftw_in, 0, sizeof(double) * N);

    for (int s = 0; s + samples < signal.size(); s += step)
    {
        // fullfill data
        for (int i = 0; i < samples; ++i)
            fftw_in[i] = win_mask[i] * signal[s + i];
        fftw_execute(p);

        // append new feature
        ret.emplace_back(N / 2 + 1);
        auto &feature = ret.back();
        
        feature[0] = fftw_out[0][0] * fftw_out[0][0];
        for (int i = 1; i < (N + 1) / 2; ++i)
            feature[i] = fftw_out[i][0] * fftw_out[i][0] + fftw_out[i][1] * fftw_out[i][1];
        if (N % 2 == 0)
            feature[N / 2] = fftw_out[N / 2][0] * fftw_out[N / 2][0];
        // scale according to NFFT.
        for (size_t i = 0; i < feature.size(); ++i)
            feature[i] *= 1.0 / N;
    }

    fftw_destroy_plan(p);
    return ret;
}

std::vector<AudioFeature> Features::mfcc(
    const AudioSamples &signal, int samplerate,
    double winlen, double winstep,
    int nfft, int numcep, int nfilt,
    double lowfreq, double highfreq,
    double preemph, int ceplifter,
    bool append_energy, AudioMask(*winfunc)(int length))
{

    auto calc = [&](int n, double *in, fftw_complex *out, AudioFeature &feature) -> void
    {
        auto p = fftw_plan_dft_r2c_1d(n, in, out, FFTW_ESTIMATE);
        fftw_execute(p);
        feature.resize(n / 2 + 1);
        feature[0] = out[0][0] * out[0][0];
        for (int i = 1; i < (n + 1) / 2; ++i)
        {
            feature[i] = out[i][0] * out[i][0] + out[i][1] * out[i][1];
        }
        if (n % 2 == 0)
            feature[n / 2] = out[n / 2][0] * out[n / 2][0];

        // power_spectrum[0] = out[0]*out[0];  /* DC component */
        // for (k = 1; k < (N+1)/2; ++k)  /* (k < N/2 rounded up) */
        //     power_spectrum[k] = out[k]*out[k] + out[N-k]*out[N-k];
        // if (N % 2 == 0) /* N is even */
        //     power_spectrum[N/2] = out[N/2]*out[N/2];  /* Nyquist freq. */

        fftw_destroy_plan(p);
    };


    return std::vector<AudioFeature>();
}


NAMESPACE_END(audio)
NAMESPACE_END(yuki)