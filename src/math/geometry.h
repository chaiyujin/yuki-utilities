#pragma once
#include "mathutils.h"
namespace yuki {
namespace math {

class BBox2d;
class Point2d
{
    double2 xy_;
public:
    friend class BBox2d;
    Point2d(double2 xy = {0, 0}) : xy_(xy) {}
    double x() const { return xy_.x; }
    double y() const { return xy_.y; }
    double2 xy() const { return xy_; }
};

class BBox2d
{
    double2 xy_;
    double2 wh_;
public:
    BBox2d(double2 left_top = {NAN, NAN}, double2 wh = {0, 0}) : xy_(left_top), wh_(wh) {}
    
    double left()   const { return xy_.x; }
    double top()    const { return xy_.y; }
    double right()  const { return xy_.x + wh_.x; }
    double bottom() const { return xy_.y + wh_.y; }
    double width()  const { return wh_.x; }
    double height() const { return wh_.y; }

    BBox2d operator+(const Point2d &p) const
    {
        if (xy_.isnan())
            return BBox2d(p.xy(), {0, 0});
        double x0 = std::min(xy_.x, p.x());
        double y0 = std::min(xy_.y, p.y());
        double x1 = std::max(xy_.x + wh_.x, p.x());
        double y1 = std::max(xy_.y + wh_.y, p.y());
        return BBox2d({x0, y0}, {x1 - x0, y1 - y0});
    }
    BBox2d &operator+=(const Point2d &p)
    {
        if (xy_.isnan())
        {
            xy_ = p.xy();
            wh_ = {0, 0};
        }
        else
        {
            double x0 = std::min(xy_.x, p.x());
            double y0 = std::min(xy_.y, p.y());
            double x1 = std::max(xy_.x + wh_.x, p.x());
            double y1 = std::max(xy_.y + wh_.y, p.y());
            xy_ = {x0, y0};
            wh_ = {x1 - x0, y1 - y0};
        }
        return (*this);
    }
};

}}