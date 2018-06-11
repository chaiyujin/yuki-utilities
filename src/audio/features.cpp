#include "features.h"
#include "fftw3.h"
#include "math/math.h"
#include "signal_process.h"
#include <iostream>
#include <complex>

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

    for (int s = 0; s + samples < (int)signal.size(); s += step)
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
        DSP::pre_emphasis(signal, preemph), samplerate,
        winlen, winstep, N, winfunc);

    auto filters = DSP::mel_filters(nfilt, N, samplerate, lowfreq, highfreq);
    AudioFeatureList feat = filters * pow_spec;
    feat = feat.unaryExpr<double(*)(double)>(zero_to_eps);

    return feat;
}

AudioFeatureList Features::log_filter_bank(
    const AudioSamples &signal, int samplerate,
    double winlen, double winstep, int nfft,
    int nfilt, double lowfreq, double highfreq,
    double preemph, AudioMask(*winfunc)(int length))
{
    AudioFeatureList feat = filter_bank(signal, samplerate, winlen, winstep, nfft, nfilt, lowfreq, highfreq, preemph, winfunc);
    return feat.unaryExpr<double(*)(double)>(DSP::energy2dB);
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
        DSP::pre_emphasis(signal, preemph), samplerate,
        winlen, winstep, N, winfunc);
    AudioFeatureList energy = pow_spec.colwise().sum();
    energy = energy.unaryExpr<double(*)(double)>(zero_to_eps);

    auto filters = DSP::mel_filters(nfilt, N, samplerate, lowfreq, highfreq);
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

    AudioFeatureList after_dct = DSP::dct(fbank, true).topRows(rows);
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


AudioFeatureList Features::lpc(
    const AudioSamples &signal, int samplerate, int order,
    double  winlen, double  winstep, double  preemph,
    AudioMask(*winfunc)(int length))
{
    auto signal_emph = DSP::pre_emphasis(signal, preemph);

    int win_samples = winlen * samplerate;
    int win_step    = winstep * samplerate;
    // check if order is valid
    if (order > signal.size() || order > win_samples)
        YUKI_ERROR_EXIT("[LPC]: Input signal (or windowed signal) must have length >= order.");

    int frames = 0;
    if ((int)signal.size() >= win_samples)
        frames = ((int)signal.size() - win_samples) / win_step + 1;

    auto win_mask = winfunc(win_samples);
    std::vector<AudioSamples> slices(frames, AudioSamples(win_samples));
    for (int s = 0, fi = 0; s + win_samples < (int)signal.size(); s += win_step, fi++)
    {
        for (int i = 0; i < win_samples; ++i)
            slices[fi][i] = win_mask[i] * signal[s + i];
        // avoid div 0 error.
        if (slices[fi][0] == 0) slices[fi][0] = 1e-10;
    }

    AudioFeatureList acorre = DSP::autocorrelation(slices, false);
    AudioFeatureList acoeff(order + 1, frames);
    AudioFeatureList err(1, frames);
    AudioFeatureList kcoeff(order, frames);
#pragma omp parallel for
    for (int i = 0; i < frames; ++i)
    {
        math::levinson(
            acorre.col(i).data(), order,
            acoeff.col(i).data(),
            err.col(i).data(),
            kcoeff.col(i).data());
    }
    return acoeff;
}


std::vector<std::vector<double>> Features::formants(
    const AudioFeatureList &lpc_A, int samplerate)
{
    std::vector<std::vector<double>> ret(lpc_A.cols());
    double F = samplerate / 2.0;

    const int rows = lpc_A.rows();
// #pragma omp parallel for
    for (int i = 0; i < lpc_A.cols(); ++i)
    {
        auto rts = math::roots(lpc_A.col(i).data(), rows);
        for (int j = 0; j < rts.size(); ++j)
        {
            if (rts[j].imag() < 0) continue;
            double ang = std::arg(rts[j]);
            double frq = ang * (F / (2.0 * M_PI));
            double bwd = -0.5 * (F / (2.0 * M_PI)) * std::log(std::abs(rts[j]));
            if (frq > 90 && bwd < 400)
                ret[i].push_back(frq);
        }
        std::sort(ret[i].begin(), ret[i].end());
    }
    return ret;
}


NAMESPACE_END(audio)
NAMESPACE_END(yuki)