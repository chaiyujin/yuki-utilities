#include "math/mathutils.h"
#include "audio/io_wav.h"

using namespace std;
void test_vec();
void test_audio();

int main()
{
    test_vec();
    test_audio();

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
    }
    catch (std::runtime_error &e)
    {
        cerr << e.what() << endl;
    }
}
