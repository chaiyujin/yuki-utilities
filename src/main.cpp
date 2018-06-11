#include "core.h"
#include "math/math.h"
#include "audio/audio.h"
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
    {
        /* draw audio */
        using namespace yuki::image;

        int S = 0;
        int E = Features::num_samples(wav_file.samplerate(), 64);
        AudioSamples wins = slice(wav_file.track(0), S, E);

        auto img = Draw::create_canvas(640, 480, {255, 255, 255});
        Eigen::MatrixXd feat = Features::log_filter_bank(wins, wav_file.samplerate());
        Draw::audio_feature(img, feat, Rect2d({0, 160}, {640, 320}));
        Draw::audio_wav(
            img, wins, Rect2d({0, 0}, {640, 160}),
            0, Features::samples_aligned_to_feature(wav_file.samplerate(), feat.cols()));

        auto haxis = Draw::horizon_axis({0.0, 0.64}, 5, 640, "SEC");
        auto vaxis = 
            concat_row(
                Draw::vertical_axis({-1, 1}, 3, 640, 160, "wav"),
                Draw::vertical_axis({0, 4000}, 5, 640, 320, "freq"),
                "right");
        auto full = concat_col(vaxis, concat_row(img, haxis));

        cv::imwrite("../asset/white.jpg", full);
    }
}

void test_image()
{
    // test color map
}
