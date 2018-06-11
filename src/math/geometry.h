#pragma once
#include "mathutils.h"
#include <iostream>

namespace yuki {
namespace math {

class Rect2d;
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

class Rect2d
{
protected:
    double2 xy_;
    double2 wh_;
public:
    Rect2d(double2 left_top = {NAN, NAN}, double2 wh = {0, 0}) : xy_(left_top), wh_(wh) {}
    double left()   const { return xy_.x; }
    double top()    const { return xy_.y; }
    double right()  const { return xy_.x + wh_.x; }
    double bottom() const { return xy_.y + wh_.y; }
    double width()  const { return wh_.x; }
    double height() const { return wh_.y; }
    // left top corner
    double2 p0() const { return xy_; }
    // right bottom corner
    double2 p1() const { return xy_ + wh_;}
    double2 wh() const { return wh_; }
    // is a valid rectangle
    bool is_valid() const { return !xy_.isnan() && wh_.x >= 0 && wh_.y >= 0; }
};

class BBox2d : public Rect2d
{
public:
    BBox2d(double2 left_top = {NAN, NAN}, double2 wh = {0, 0}) : Rect2d(left_top, wh) {}
    BBox2d(const Rect2d &rect) : Rect2d(rect) {}

    BBox2d operator|(const Point2d &p) const
    {
        if (!is_valid())
            return BBox2d(p.xy(), {0, 0});
        double x0 = std::min(xy_.x, p.x());
        double y0 = std::min(xy_.y, p.y());
        double x1 = std::max(xy_.x + wh_.x, p.x());
        double y1 = std::max(xy_.y + wh_.y, p.y());
        return BBox2d({x0, y0}, {x1 - x0, y1 - y0});
    }
    BBox2d &operator|=(const Point2d &p)
    {
        if (!is_valid()) {
            xy_ = p.xy();
            wh_ = {0, 0};
        }
        else {
            double x0 = std::min(xy_.x, p.x());
            double y0 = std::min(xy_.y, p.y());
            double x1 = std::max(xy_.x + wh_.x, p.x());
            double y1 = std::max(xy_.y + wh_.y, p.y());
            xy_ = {x0, y0};
            wh_ = {x1 - x0, y1 - y0};
        }
        return (*this);
    }
    BBox2d operator|(const BBox2d &b) const
    {
        // return the valid one
        if (!is_valid() && b.is_valid()) return b;
        if (is_valid() && !b.is_valid()) return (*this);
        // return invalid bbox if both are invalid
        if (!is_valid() && !b.is_valid()) return BBox2d();
        // both are valid, merge them
        double x0 = std::min(xy_.x, b.xy_.x);
        double y0 = std::min(xy_.y, b.xy_.y);
        double x1 = std::max(xy_.x + wh_.x, b.xy_.x + b.wh_.x);
        double y1 = std::max(xy_.y + wh_.y, b.xy_.y + b.wh_.y);
        return BBox2d({x0, y0}, {x1 - x0, y1 - y0});
    }
    BBox2d &operator|=(const BBox2d &b)
    {
        // return the valid one
        if (!is_valid() && b.is_valid())
            *this = b;
        // return invalid bbox if both are invalid
        else if (!is_valid() && !b.is_valid())
            *this = BBox2d();
        // both are valid, merge them
        else {
            double x0 = std::min(xy_.x, b.xy_.x);
            double y0 = std::min(xy_.y, b.xy_.y);
            double x1 = std::max(xy_.x + wh_.x, b.xy_.x + b.wh_.x);
            double y1 = std::max(xy_.y + wh_.y, b.xy_.y + b.wh_.y);
            xy_ = {x0, y0};
            wh_ = {x1 - x0, y1 - y0};
        }
        return (*this);
    }
    BBox2d operator&(const BBox2d &b) const
    {
        // return a invalid bbox, if a or b is invalid.
        if (!is_valid() || !b.is_valid())
            return BBox2d();
        double x0 = std::max(xy_.x, b.xy_.x);
        double y0 = std::max(xy_.y, b.xy_.y);
        double x1 = std::min(xy_.x + wh_.x, b.xy_.x + b.wh_.x);
        double y1 = std::min(xy_.y + wh_.y, b.xy_.y + b.wh_.y);
        return BBox2d({x0, y0}, {x1 - x0, y1 - y0});
    }
    BBox2d &operator&=(const BBox2d &b)
    {
        // return a invalid bbox, if a or b is invalid.
        if (!is_valid() || !b.is_valid())
            *this = BBox2d();
        else {
            double x0 = std::max(xy_.x, b.xy_.x);
            double y0 = std::max(xy_.y, b.xy_.y);
            double x1 = std::min(xy_.x + wh_.x, b.xy_.x + b.wh_.x);
            double y1 = std::min(xy_.y + wh_.y, b.xy_.y + b.wh_.y);
            xy_ = {x0, y0};
            wh_ = {x1 - x0, y1 - y0};
        }
        return (*this);
    }
    
};

inline std::ostream &operator<<(std::ostream &out, const Rect2d &rect)
{
    if (rect.is_valid())
        out << "xy: " << rect.p0() << ", wh: " << rect.wh() << std::endl;
    else
        out << "invalid rect\n";
    return out;
}

}}