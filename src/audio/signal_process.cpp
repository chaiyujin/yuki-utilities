#include "signal_process.h"
#include "math/math.h"
#include <fftw3.h>

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(audio)


Eigen::MatrixXd DSP::dct(const Eigen::MatrixXd &in, bool normalization)
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


AudioSamples DSP::pre_emphasis(const AudioSamples &signal, double preemph)
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

Eigen::MatrixXd DSP::mel_filters(int nfilt, int nfft, int samplerate, double lowfreq, double highfreq)
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

AudioFeatureList DSP::autocorrelation(
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

std::pair<Eigen::VectorXd, Eigen::VectorXcd> DSP::freqz(
    const Eigen::VectorXd &b,
    const Eigen::VectorXd &a,
    int wor_n,
    bool whole)
{
    if (b.rows() != a.rows())
        YUKI_ERROR_EXIT("[DSP]: freqz asks for `a` and `b` in same size.\n");
    double lastp = (whole) ? 2.0 * M_PI : M_PI;
    std::pair<Eigen::VectorXd, Eigen::VectorXcd> ret;
    Eigen::VectorXd &W = ret.first;
    Eigen::VectorXcd &H = ret.second;
    W.resize(wor_n);
    H.resize(wor_n);
    for (int i = 0; i < wor_n; ++i)
    {
        W(i, 0) = lastp * i / (double)wor_n;
        std::complex<double> zm1 = std::exp(std::complex<double>(0, -1) * W(i, 0));
        auto h = math::polyval(b.data(), b.rows(), zm1, true) / 
                 math::polyval(a.data(), a.rows(), zm1, true);
        H(i, 0) = h;
    }
    return ret;
}

std::pair<Eigen::VectorXd, Eigen::VectorXcd> DSP::freqz(
    double b,
    const Eigen::VectorXd &a,
    int wor_n,
    bool whole)
{
    Eigen::VectorXd b_(a.rows());
    b_.setZero();
    b_[0] = b;
    return freqz(b_, a, wor_n, whole);
}


Eigen::MatrixX2d DSP::WH_to_FreqdB(
    const std::pair<Eigen::VectorXd, Eigen::VectorXcd> &WH,
    double frequency)
{
    /*
        w = frequency * w / (2 * np.pi)
        h = 20 * np.log10(np.abs(h))
    */
   Eigen::MatrixX2d ret(WH.first.rows(), 2);

    for (int i = 0; i < ret.rows(); ++i)
    {
        ret(i, 0) = WH.first(i, 0) * frequency / (2.0 * M_PI);
        ret(i, 1) = 20.0 * std::log10(std::abs(WH.second(i, 0)));
    }

   return ret;
}

NAMESPACE_END(audio)
NAMESPACE_END(yuki)