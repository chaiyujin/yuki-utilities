#include "draw.h"
#include "color_map.h"

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(image)


cv::Mat Draw::create_canvas(int w, int h, rgba32 clear_color)
{
    cv::Mat ret = cv::Mat(h, w, CV_8UC4);
    ret.setTo(cv::Scalar(clear_color.b(), clear_color.g(), clear_color.r(), clear_color.a()));
    return ret;
}

void Draw::audio_wav(
        cv::Mat &img,
        const std::vector<double> &signal,
        const Rect2d &draw_rect,
        int s, int e)
{
    BBox2d rect = (!draw_rect.is_valid())
        ? BBox2d({0, 0}, {(double)img.cols, (double)img.rows})
        : BBox2d(draw_rect) & BBox2d({0, 0}, {(double)img.cols, (double)img.rows});
    if ((int)signal.size() < e) e = (int)signal.size();
    // draw signal[s: e] into img[rect]
    const double MAX = 1.0;
    const double MIN = -1.0;
    auto map_h = [&](double v) -> int {
        return (int)round( (v - MIN) / (MAX - MIN) * (-(rect.height() - 1)) + rect.bottom() - 1 );
    };
    auto map_w = [&](int i) -> int {
        return (int)round( (double)(i - s) / (double)(e - 1 - s) * (rect.width() - 1) + rect.left() );
    };
    for (int i = s; i + 1 < e; ++i)
    {
        int x0 = map_w(i);
        int x1 = map_w(i + 1);
        int y0 = map_h(signal[i]);
        int y1 = map_h(signal[i + 1]);
        cv::line(img, cv::Point(x0, y0), cv::Point(x1, y1), CV_RGB(0, 100, 200), 1, cv::LINE_AA);
    }
}

void Draw::audio_feature(
    cv::Mat &img, Eigen::MatrixXd &feature,
    const Rect2d &draw_rect,
    const range_float &feat_minmax)
{
    BBox2d rect = (!draw_rect.is_valid())
        ? BBox2d({0, 0}, {(double)img.cols, (double)img.rows})
        : BBox2d(draw_rect) & BBox2d({0, 0}, {(double)img.cols, (double)img.rows});
    // draw feature into img[rect]
    
    auto small = ColorMap::colorize_eigen(feature, feat_minmax, ColorMap::Jet, true);
    cv::resize(small, small, cv::Size(rect.width(), rect.height()), 0, 0, cv::INTER_NEAREST);
    small.copyTo(img(cv::Rect(rect.p0().x, rect.p0().y, rect.wh().x, rect.wh().y)));
}

cv::Mat Draw::horizon_axis(
    range_double minmax, int num_label,
    int main_w, const std::string &unit)
{
    const int len = main_w;
    const double H = 12 * len / 320;
    const double W = 5 * len / 320;
    const double Scale = 0.3 * len / 320.0;

    int h = 0;
    cv::Mat img = create_canvas(main_w, H * 5 / 3 + 5);
    cv::line(img, cv::Point(0,          h), cv::Point(main_w,       h),     CV_RGB(50, 50, 50), 2, cv::LINE_AA);
    cv::line(img, cv::Point(0,          h), cv::Point(0,            h + H), CV_RGB(50, 50, 50), 2, cv::LINE_AA);
    cv::line(img, cv::Point(main_w,     h), cv::Point(main_w,       h + H), CV_RGB(50, 50, 50), 2, cv::LINE_AA);
    
    auto num_to_str = [](double x) -> std::string {
        std::string str = std::to_string(x);
        while (str.back() == '0' && str.length() > 1) str = str.substr(0, str.length() - 1);
        if (str.back() == '.') str = str.substr(0, str.length() - 1);
        return str;
    };

    /* labels */
    if (num_label < 2) num_label = 2;
    for (int i = 1; i + 1 < num_label; ++i)
    {
        int label_x = (double)main_w * i / (num_label - 1);
        double label_v = (double)(minmax.end - minmax.start) * i / (num_label - 1) + minmax.start;
        auto str = num_to_str(label_v);
        cv::line(img, cv::Point(label_x, h), cv::Point(label_x, h + H / 3), CV_RGB(50, 50, 50), 1, cv::LINE_AA);
        cv::putText(
            img, str, cv::Point(label_x - str.length() * W / 2, h + H),
            cv::FONT_HERSHEY_DUPLEX, Scale, CV_RGB(50, 50, 50), 1, cv::LINE_AA);
    }
    {
        auto str = num_to_str(minmax.start);
        cv::putText(
            img, str, cv::Point(0, h + H),
            cv::FONT_HERSHEY_DUPLEX, Scale, CV_RGB(50, 50, 50), 1, cv::LINE_AA);
        str = num_to_str(minmax.end);
        cv::putText(
            img, str, cv::Point(main_w - str.length() * W - 3, h + H),
            cv::FONT_HERSHEY_DUPLEX, Scale, CV_RGB(50, 50, 50), 1, cv::LINE_AA);
    }
    /* axis unit */
    if (unit.length()) {
        int x = (main_w) / 2 - unit.length() * W / 2 ;
        int y = h + H * 5 / 3;
        cv::putText(
            img, unit, cv::Point(x, y),
            cv::FONT_HERSHEY_DUPLEX, Scale, CV_RGB(50, 50, 50), 1, cv::LINE_AA);
    }
    return img;
}

cv::Mat Draw::vertical_axis(
    range_double minmax, int num_label,
    int main_w, int main_h, const std::string &unit)
{
    const double H = 6 * main_w / 320.0;
    const double W = 6 * main_w / 320.0;
    const int gap = 3 + W;
    const double Scale = 0.3 * main_w / 320.0;
    auto num_to_str = [](double x) -> std::string {
        std::string str = std::to_string(x);
        while (str.back() == '0' && str.length() > 1) str = str.substr(0, str.length() - 1);
        if (str.back() == '.') str = str.substr(0, str.length() - 1);
        return str;
    };

    std::vector<int> label_y;
    std::vector<std::string> labels;
    label_y.push_back(main_h);
    labels.push_back(num_to_str(minmax.start));
    if (num_label < 2) num_label = 2;
    for (int i = 1; i + 1 < num_label; ++i)
    {
        int y_ = main_h - (double)main_h * i / (num_label - 1);
        double v_ = (double)(minmax.end - minmax.start) * i / (num_label - 1) + minmax.start;
        label_y.push_back(y_);
        labels.push_back(num_to_str(v_));
    }
    label_y.push_back(0);
    labels.push_back(num_to_str(minmax.end));

    int max_len = 0;
    for (auto &s : labels)
        if (s.length() > max_len) max_len = s.length();

    max_len += unit.length() + 1;
    int w = max_len * W + gap;
    cv::Mat img = create_canvas(w, main_h);

    cv::line(img, cv::Point(w - 1, main_h - 1), cv::Point(w - 1,               0),          CV_RGB(50, 50, 50), 1, cv::LINE_AA);
    cv::line(img, cv::Point(w - 1, main_h - 1), cv::Point(w - W * 2 - gap,     main_h - 1), CV_RGB(50, 50, 50), 1, cv::LINE_AA);
    cv::line(img, cv::Point(w - 1, 0),          cv::Point(w - W * 2 - gap,     0),          CV_RGB(50, 50, 50), 1, cv::LINE_AA);
    
    {
        auto str = num_to_str(minmax.start);
        cv::putText(
            img, str, cv::Point(w - str.length() * W - gap, main_h - 3),
            cv::FONT_HERSHEY_DUPLEX, Scale, CV_RGB(50, 50, 50), 1, cv::LINE_AA);
        str = num_to_str(minmax.end);
        cv::putText(
            img, str, cv::Point(w - str.length() * W - gap, H + 3),
            cv::FONT_HERSHEY_DUPLEX, Scale, CV_RGB(50, 50, 50), 1, cv::LINE_AA);
    }
    for (int i = 1; i + 1 < num_label; ++i)
    {
        auto str = labels[i];
        int x = w - gap - W * str.length();
        int y = label_y[i];
        cv::line(img, cv::Point(w, y), cv::Point(w - gap, y), CV_RGB(50, 50, 50), 1, cv::LINE_AA);
        cv::putText(
            img, str, cv::Point(x, y + H / 2),
            cv::FONT_HERSHEY_DUPLEX, Scale, CV_RGB(50, 50, 50), 1, cv::LINE_AA);
    }
    if (unit.length()) {
        int x = 1 * W;
        int y = main_h / 2;
        cv::putText(
            img, unit, cv::Point(x, y + H / 2),
            cv::FONT_HERSHEY_DUPLEX, Scale, CV_RGB(50, 50, 50), 1, cv::LINE_AA);
    }

    return img;
}

NAMESPACE_END(image)
NAMESPACE_END(yuki)