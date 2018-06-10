#include "features.h"
#include "fftw3.h"
#include "math/mathutils.h"
#include <iostream>

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(audio)


double  Features::default_winlen            = 0.025;
double  Features::default_winstep           = 0.010;
int     Features::default_nfft              = 0;
int     Features::default_numcep_dct        = 13;
int     Features::default_nfilt_mel         = 40;
int     Features::default_ceplifter_mfcc    = 0;
double  Features::default_preemph           = 0;
double  Features::dafault_lowfreq           = 0;
double  Features::dafault_highfreq          = DBL_MAX;
AudioMask(*Features::default_winfunc)(int)  = Features::Window::ones;

AudioMask Features::Window::ones(int length)
{
    return AudioMask(length, 1.0);
}

Eigen::MatrixXd Features::dct(const Eigen::MatrixXd &in, bool normalization)
{
    const int N = in.rows();
    Eigen::MatrixXd ret(in.rows(), in.cols());
    
    double *fftw_in  = fftw_alloc_real(N);
    double *fftw_out = fftw_alloc_real(N);

    // DCT-II (the well-known DCT)
    auto plan = fftw_plan_r2r_1d(N, fftw_in, fftw_out, FFTW_REDFT10, FFTW_ESTIMATE);

    for (int c = 0; c < in.cols(); ++c)
    {
        for (int r = 0; r < N; ++r)
            fftw_in[r] = in(r, c);
        fftw_execute(plan);
        for (int r = 0; r < N; ++r)
            ret(r, c) = fftw_out[r];
        if (normalization)
        {
            /* f = sqrt(1/(4*N)) if r = 0,
               f = sqrt(1/(2*N)) otherwise. */
            ret(0, c) *= std::sqrt(1.0 / (4.0 * (double)N));
            for (int r = 1; r < N; ++r)
                ret(r, c) *= std::sqrt(1.0 / (2.0 * (double)N));
        }
    }

    fftw_destroy_plan(plan);

    fftw_free(fftw_in);
    fftw_free(fftw_out);

    return ret;
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
    double *        fftw_in  = fftw_alloc_real(N);
    fftw_complex *  fftw_out = fftw_alloc_complex(N / 2 + 1);
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
    fftw_free(fftw_in);
    fftw_free(fftw_out);
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
    feat = feat.unaryExpr<double(*)(double)>(zero_to_eps);

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

    std::cout << "feature " << N << std::endl;

    if (highfreq == DBL_MAX)
        highfreq = samplerate / 2.0;
    auto pow_spec = power_spectrum(
        pre_emphasis(signal, preemph), samplerate,
        winlen, winstep, N, winfunc);
    AudioFeatureList energy = pow_spec.colwise().sum();
    energy = energy.unaryExpr<double(*)(double)>(zero_to_eps);

    auto filters = mel_filters(nfilt, N, samplerate, lowfreq, highfreq);
    AudioFeatureList feat = filters * pow_spec;
    feat = feat.unaryExpr<double(*)(double)>(zero_to_eps);

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
    AudioFeatureList fbank, energy;
    if (!append_energy)
        fbank.noalias() = filter_bank(signal, samplerate, winlen, winstep, nfft, nfilt, lowfreq, highfreq, preemph, winfunc);
    else
    {
        auto pair = filter_bank_energy(signal, samplerate, winlen, winstep, nfft, nfilt, lowfreq, highfreq, preemph, winfunc);
        fbank.noalias() = pair.first;
        energy.noalias() = pair.second;
    }

    int rows = std::min(numcep, (int)fbank.rows());
    fbank = fbank.unaryExpr<double(*)(double)>(std::log);

    AudioFeatureList after_dct = dct(fbank, true).topRows(rows);
    if (ceplifter > 0)
    {
        /* Apply a cepstral lifter the the matrix of cepstra. This has the effect of increasing the
           magnitude of the high frequency DCT coeffs */
        Eigen::MatrixXd lift = Eigen::MatrixXd::Zero(rows, rows);
        for (int i = 0; i < rows; ++i)
            lift(i, i) = 1.0 + (ceplifter / 2.0) * std::sin(M_PI * (double)i / (double)ceplifter);
        after_dct = lift * after_dct;
    }
    if (append_energy)
    {
        energy = energy.unaryExpr<double((*)(double))>(std::log);
        after_dct.row(0) << energy.row(0);
        return after_dct;
    }
    else
    {
        return after_dct.bottomRows(rows - 1);
    }
}

AudioFeatureList Features::autocorrelation(
    const std::vector<AudioSamples> &signal_list,
    bool biased_acorre_estimator)
{
    if (signal_list.size() == 0) return AudioFeatureList();
    int maxlag = (int)signal_list[0].size();
    int N = math::nextpow2(maxlag * 2 - 1);
    
    double *        fftw_in  = fftw_alloc_real(N);
    fftw_complex *  fftw_cp  = fftw_alloc_complex(N);
    double *        fftw_out = fftw_alloc_real(N);

    auto plan_fft  = fftw_plan_dft_r2c_1d(N, fftw_in, fftw_cp,  FFTW_ESTIMATE);
    auto plan_ifft = fftw_plan_dft_c2r_1d(N, fftw_cp, fftw_out, FFTW_ESTIMATE);
    memset(fftw_in,  0, sizeof(double) * N);
    memset(fftw_out, 0, sizeof(double) * N);

    AudioFeatureList ret(maxlag + 1, signal_list.size());
    for (size_t i = 0; i < signal_list.size(); ++i)
    {
        for (int j = 0; j < maxlag; ++j)
            fftw_in[j] = signal_list[i][j];
        fftw_execute(plan_fft);
        // sqaure
        for (int j = 0; j < N; ++j)
        {
            fftw_cp[j][0] = fftw_cp[j][0] * fftw_cp[j][0] + fftw_cp[j][1] * fftw_cp[j][1];
            fftw_cp[j][1] = 0;
        }
        fftw_execute(plan_ifft);
        // copy result
        for (int j = 0; j < maxlag + 1; ++j)
        {
            fftw_out[j] /= (double)N;
            ret(j, i) = (biased_acorre_estimator)
                ? fftw_out[j] / (double)maxlag
                : fftw_out[j];
        }
    }

    fftw_destroy_plan(plan_fft);
    fftw_destroy_plan(plan_ifft);
    fftw_free(fftw_in);
    fftw_free(fftw_cp);
    fftw_free(fftw_out);

    return ret;
}


NAMESPACE_END(audio)
NAMESPACE_END(yuki)