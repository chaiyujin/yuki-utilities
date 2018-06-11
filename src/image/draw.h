#pragma once
#include "common.h"
#include "math/mathutils.h"
#include <vector>

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(image)

using namespace yuki::math;

class Draw
{
public:
    static cv::Mat create_canvas(int w, int h, rgba32 clear_color={255, 255, 255, 255});
    static void audio_wav(cv::Mat &img, const std::vector<double> &signal);

};

NAMESPACE_END(image)
NAMESPACE_END(yuki)