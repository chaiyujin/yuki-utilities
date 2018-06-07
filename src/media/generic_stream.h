#pragma once
#include "common.h"
#include <map>
#include <vector>

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(media)

template <class T>
struct GenericFrame
{
    double ts;
    T      data;
};

template <class T>
class GenericStream
{
public:
    enum Type
    {
        None = 0,
        AudioFeature = 1,
        VideoFeature = 2,
        PhonemeLabel = 3,
        Audio = 4,
        Video = 5,
        Undefined
    };

    GenericStream(Type type=Type::None, int size=0) : frames(size), type_(type) {}
    
    size_t size() const { return frames_.size(); }
    void resize(int size) { frames_.resize(size); }
    void push_back(const T &d, double ts) { frames_.push_back(GenericFrame<T>{ts, d}); }
    void push_back(const GenericFrame& frame) { frames_.push_back(frame); }
    GenericFrame &operator[](int idx) { return frames_[idx]; }
    const GenericFrame &operator[](int idx) const { return frames_[idx]; }

private:
    std::vector<GenericFrame<T>> frames_;
    Type                         type_;
};

class GenericStreamGroup
{
public:
    
};

NAMESPACE_END(media)
NAMESPACE_END(yuki)