#include "math/mathutils.h"
#include "audio/io_wav.h"
#include "audio/features.h"
#include "image/image.h"

using namespace std;
void test_vec();
void test_audio();
void test_image();

int main()
{
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
    cout << "to_eigen:\n" << to_eigen(double3(1,2,3)) << endl;
}

void test_audio()
{
    using namespace yuki::audio;
    using namespace yuki::image;
    WAVPCM wav_file;

    auto calc = [](double x) -> double
    {
        return std::log10(std::abs(x) + 1e-8) * 20;
    };

    try
    {
        wav_file.read("../asset/test1.wav");
        Eigen::MatrixXd feature = Features::filter_bank(wav_file.track(0), wav_file.samplerate());
        cout << feature.rows() << " " << feature.cols() << endl;
        feature = feature.unaryExpr<double(*)(double)>(calc);
        auto min_ = feature.minCoeff();
        auto max_ = feature.maxCoeff();
        cout << min_ << " " << max_ << endl;
        
        auto img = ColorMap::colorize_eigen(feature.leftCols(5000), {-150, 100}, ColorMap::Jet);
        // cv::resize(img, img, cv::Size(900, 240), 0, 0, cv::INTER_NEAREST);
        cout << img.rows << " " << img.cols << endl;
        cv::imwrite("../asset/spectrum.jpg", img);
    }
    catch (std::runtime_error &e)
    {
        cerr << e.what() << endl;
    }
}

void test_image()
{
    // test color map
    using namespace yuki::image;
    std::cout << ColorMap::Jet.get(0.0) << std::endl;
    std::cout << ColorMap::Jet.get(1.0) << std::endl;
}
