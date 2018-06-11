#pragma once

#include "core.h"
#include "opencv2/opencv.hpp"
#include "Eigen/Eigen"
#include "math/mathutils.h"

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(image)
using namespace yuki::math;

inline cv::Mat concat_row(const cv::Mat &a, const cv::Mat &b, const std::string &align="left", rgba32 clear_color={255, 255, 255, 255})
{
    int cols = std::max(a.cols, b.cols);
    int sc_a = 0, sc_b = 0;
    if (align == "right") {
        sc_a = cols - a.cols;
        sc_b = cols - b.cols;
    }
    else if (align == "center") {
        sc_a = (cols - a.cols) / 2;
        sc_b = (cols - b.cols) / 2;
    }
    cv::Mat ret = cv::Mat(a.rows + b.rows, cols, CV_8UC4);
    ret.setTo(cv::Scalar(clear_color.b(), clear_color.g(), clear_color.r(), clear_color.a()));
    a.copyTo(ret(cv::Rect(sc_a, 0, a.cols, a.rows)));
    b.copyTo(ret(cv::Rect(sc_b, a.rows, b.cols, b.rows)));
    return ret;
}

inline cv::Mat concat_col(const cv::Mat &a, const cv::Mat &b, const std::string &align="top", rgba32 clear_color={255, 255, 255, 255})
{
    int rows = std::max(a.rows, b.rows);
    int sr_a = 0, sr_b = 0;
    if (align == "bottom") {
        sr_a = rows - a.rows;
        sr_b = rows - b.rows;
    }
    else if (align == "center") {
        sr_a = (rows - a.rows) / 2;
        sr_b = (rows - b.rows) / 2;
    }
    cv::Mat ret = cv::Mat(rows, a.cols + b.cols, CV_8UC4);
    ret.setTo(cv::Scalar(clear_color.b(), clear_color.g(), clear_color.r(), clear_color.a()));
    a.copyTo(ret(cv::Rect(0,      sr_a, a.cols, a.rows)));
    b.copyTo(ret(cv::Rect(a.cols, sr_b, b.cols, b.rows)));
    return ret;
}

NAMESPACE_END(image)
NAMESPACE_END(yuki)