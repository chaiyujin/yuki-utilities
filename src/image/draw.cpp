#include "draw.h"

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(image)


cv::Mat Draw::create_canvas(int w, int h, rgba32 clear_color)
{
    cv::Mat ret = cv::Mat(h, w, CV_8UC4);
    ret.setTo(cv::Scalar(clear_color.r(), clear_color.g(), clear_color.b(), clear_color.a()));
    return ret;
}

NAMESPACE_END(image)
NAMESPACE_END(yuki)