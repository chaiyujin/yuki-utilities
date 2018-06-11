#pragma once
#include "common.h"
#include "math/math.h"
#include <vector>

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(image)

using namespace yuki::math;

class Draw
{
public:
    static cv::Mat create_canvas(int w, int h, rgba32 clear_color={255, 255, 255, 255});
    static void audio_wav(
        cv::Mat &img,                     const std::vector<double> &signal,
        const Rect2d &draw_rect=Rect2d(), int s=0, int e=INT32_MAX);
    static void audio_feature(
        cv::Mat &img,                     Eigen::MatrixXd &feature,
        const Rect2d &draw_rect=Rect2d(), const range_float &feat_minmax = {0, 0});
    static cv::Mat horizon_axis(
        range_double minmax, int num_label,
        int main_w, const std::string &unit="");
    static cv::Mat vertical_axis(
        range_double minmax, int num_label,
        int main_w, int main_h, const std::string &unit="");
};

NAMESPACE_END(image)
NAMESPACE_END(yuki)