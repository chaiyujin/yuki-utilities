#include "features.h"
#include "fftw3.h"
#include <iostream>

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(audio)

AudioMask Features::Window::ones(int length)
{
    return AudioMask(length, 1.0);
}

AudioSamples Features::pre_emphasis(const AudioSamples &signal, double preemph)
{
    if (preemph > 0)
    {
        AudioSamples ret(signal.size());
        ret[0] = signal[0];
        for (size_t i = 1; i < signal.size(); ++i)
        {
            ret[i] = signal[i] - preemph * signal[i - 1];
        }
        return ret;
    }
    else
        return signal;
}

Eigen::MatrixXd Features::mel_filters(int nfilt, int nfft, int samplerate, double lowfreq, double highfreq)
{
    Eigen::MatrixXd filters(nfilt, nfft / 2 + 1);
    if (highfreq > samplerate / 2)
        highfreq = samplerate / 2.0;
    double low_mel = hz2mel(lowfreq);
    double high_mel = hz2mel(highfreq);
    double mel_step = (high_mel - low_mel) / (nfilt + 1);  // triangles, need one more point
    std::vector<int> bin(nfilt + 2);
    for (double mel = low_mel, i = 0; mel < high_mel; mel += mel_step, i++)
        bin[i] = (int)std::floor(mel2hz(mel) * (nfft + 1.0) / samplerate);
    bin.back() = (int)std::floor(highfreq * (nfft + 1.0) / samplerate);
    // clear zero
    filters.setZero();
    for (int i = 0; i < nfilt; ++i)
    {
        for (int j = bin[i]; j < bin[i + 1]; ++j)
            filters(i, j) = (double)(j - bin[i]) / (double)(bin[i + 1] - bin[i]);
        for (int j = bin[i + 1]; j < bin[i + 2]; ++j)
            filters(i, j) = (double)(j - bin[i + 2]) / (double)(bin[i + 1] - bin[i + 2]);
    }
    return filters;
}

AudioFeatureList Features::power_spectrum(
    const AudioSamples &signal, int samplerate,
    double winlen, double winstep, int nfft,
    AudioMask(*winfunc)(int length))
{
    int samples = winlen * samplerate;
    int step = winstep * samplerate;
    int N = (nfft > samples) ? nfft : samples;

    auto win_mask = winfunc(samples);
    double *fftw_in = fftw_alloc_real(N);
    fftw_complex *fftw_out = fftw_alloc_complex(N / 2 + 1);
    auto p = fftw_plan_dft_r2c_1d(N, fftw_in, fftw_out, FFTW_ESTIMATE);
    memset(fftw_in, 0, sizeof(double) * N);

    int frames = 0;
    if ((int)signal.size() >= samples)
        frames = ((int)signal.size() - samples) / step + 1;

    AudioFeatureList ret(N / 2 + 1, frames);

    for (int s = 0; s + samples < signal.size(); s += step)
    {
        // fullfill data
        for (int i = 0; i < samples; ++i)
            fftw_in[i] = win_mask[i] * signal[s + i];
        fftw_execute(p);

        // append new feature
        auto feature = ret.col(s / step);
        
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

AudioFeatureList Features::filter_bank(
    const AudioSamples &signal, int samplerate,
    double winlen, double winstep, int nfft,
    int nfilt, double lowfreq, double highfreq,
    double preemph, AudioMask(*winfunc)(int length))
{
    auto zero_to_eps = [](double x) -> double {
        return (x == 0) ? std::numeric_limits<double>::epsilon() : x;
    };

    int samples = winlen * samplerate;
    int step = winstep * samplerate;
    int N = (nfft > samples) ? nfft : samples;

    if (highfreq == DBL_MAX)
        highfreq = samplerate / 2.0;
    auto pow_spec = power_spectrum(
        pre_emphasis(signal, preemph), samplerate,
        winlen, winstep, N, winfunc);

    auto filters = mel_filters(nfilt, N, samplerate, lowfreq, highfreq);
    AudioFeatureList feat = filters * pow_spec;
    feat.unaryExpr<double(*)(double)>(zero_to_eps);

    return feat;
}

std::pair<AudioFeatureList, AudioFeatureList> Features::filter_bank_energy(
    const AudioSamples &signal, int samplerate,
    double winlen, double winstep, int nfft,
    int nfilt, double lowfreq, double highfreq,
    double preemph, AudioMask(*winfunc)(int length))
{
    auto zero_to_eps = [](double x) -> double {
        return (x == 0) ? std::numeric_limits<double>::epsilon() : x;
    };

    int samples = winlen * samplerate;
    int step = winstep * samplerate;
    int N = (nfft > samples) ? nfft : samples;

    if (highfreq == DBL_MAX)
        highfreq = samplerate / 2.0;
    auto pow_spec = power_spectrum(
        pre_emphasis(signal, preemph), samplerate,
        winlen, winstep, N, winfunc);
    AudioFeatureList energy = pow_spec.colwise().sum();
    energy.unaryExpr<double(*)(double)>(zero_to_eps);

    auto filters = mel_filters(nfilt, N, samplerate, lowfreq, highfreq);
    AudioFeatureList feat = filters * pow_spec;
    feat.unaryExpr<double(*)(double)>(zero_to_eps);

    return {feat, energy};
}

AudioFeatureList Features::mfcc(
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


    return AudioFeatureList();
}


NAMESPACE_END(audio)
NAMESPACE_END(yuki)