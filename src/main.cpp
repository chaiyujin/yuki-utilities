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
    WAVPCM wav_file;

    try
    {
        wav_file.read("../asset/test.wav");
        auto feature = Features::power_spectrum(wav_file.track(0), wav_file.samplerate());
        std::cout << feature[0].size() << std::endl;
        for (int i = 0; i < 10; ++i)
        {
            std::cout << feature[0][i] << " ";
            if (i % 5 == 4) std::cout << endl;
        }
        std::cout << endl;
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