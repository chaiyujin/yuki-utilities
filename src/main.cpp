#include "core.h"
#include "math/mathutils.h"
#include "audio/io_wav.h"
#include "audio/features.h"
#include "audio/signal_process.h"
#include "image/image.h"

using namespace std;
using namespace yuki;
void test_vec();
void test_audio();
void test_image();

int main()
{
    using namespace yuki::math;
    test_vec();
    test_audio();
    test_image();

    return 0;
}

void test_vec()
{
    using namespace yuki::math;
    cout << float3() << endl;
    cout << "cross (1,2,3) (3,2,1) = " << cross(double3(1,2,3), double3(3,2,1)) << endl;
    cout << "dot   (1,2,3) (3,2,1) = " << dot(double3(1,2,3), double3(3,2,1)) << endl;
    cout << "to_eigen:\n" << to_eigen(double3(1,2,3)).transpose() << endl;
    cout << "double3 size: " << sizeof(double3) << ", 3 double size:" << sizeof(double) * 3 << endl;
}

void test_audio()
{
    using namespace yuki::audio;
    using namespace yuki::image;
    WAVPCM wav_file;

    wav_file.read("../asset/test1.wav");
    /* spectrum related */
    {
        auto save_feat_img = [](const Eigen::MatrixXd &feature, const std::string &filename) -> void
        {
            cout << "save feature into " << filename << endl;
            cout << "  feature size:  " << feature.rows() << " x " << feature.cols() << endl;
            cout << "  feature range: " << feature.minCoeff() << " ~ " << feature.maxCoeff() << endl;
            auto img = ColorMap::colorize_eigen(feature, {0, 0}, ColorMap::Jet, true);
            cv::resize(img, img, cv::Size(640, 320), 0, 0, cv::INTER_NEAREST);
            cv::imwrite(filename, img);
        };
        auto dB = [](double x) -> double
        {
            return std::log10(x) * 20;
        };
        Eigen::MatrixXd spec = Features::power_spectrum(wav_file.track(0), wav_file.samplerate()).leftCols(64);
        Eigen::MatrixXd fbank = Features::filter_bank(wav_file.track(0), wav_file.samplerate()).leftCols(64);
        Eigen::MatrixXd logfbank = fbank.unaryExpr<double(*)(double)>(dB);
        Eigen::MatrixXd mfcc = Features::mfcc(wav_file.track(0), wav_file.samplerate()).leftCols(64);

        save_feat_img(spec, "../asset/spectrum.jpg");
        save_feat_img(fbank, "../asset/fbank.jpg");
        save_feat_img(logfbank, "../asset/log_fbank.jpg");
        save_feat_img(mfcc, "../asset/mfcc.jpg");
    }
    /* linear prediction  */
    {
        AudioSamples sample({2, 3, -1});
        std::vector<AudioSamples> samples({sample});
        auto acorre = DSP::autocorrelation(samples);
        cout << "auto correlation of (2, 3, -1): " << acorre.transpose() << endl;
        Eigen::MatrixXd lpc = Features::lpc(wav_file.track(0), wav_file.samplerate(), 12).leftCols(64);
        cout << "LPC feature:\n";
        cout << "  feature size:  " << lpc.rows() << " x " << lpc.cols() << endl;
        // cout << "  col 0: " << lpc.col(0).transpose() << endl;
        // cout << "  col 10: " << lpc.col(10).transpose() << endl;
        std::vector<std::vector<double>> formants = Features::formants(lpc.leftCols(64), wav_file.samplerate());
        cout << "Formants:\n";
        cout << "  " << formants[0] << endl;
        
        auto WH_ = DSP::freqz(1, lpc.col(0));
        auto WH = DSP::WH_to_FreqdB(WH_, wav_file.samplerate() / 2.0);
        auto H = WH.col(1);
        cout << "Freqz:\n";
        cout << "  " << H.topRows(10).transpose() << endl;
    }
}

void test_image()
{
    // test color map
    using namespace yuki::image;

    auto img = Draw::create_canvas(320, 240, {255, 255, 255});
    cv::imwrite("../asset/white.jpg", img);
}
