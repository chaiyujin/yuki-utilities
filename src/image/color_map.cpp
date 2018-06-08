#include "color_map.h"


NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(image)

const ColorMap ColorMap::Jet
{{
    { 0, 0, 100 },
    { 0, 0, 255 },
    { 0, 255, 255 },
    { 255, 255, 0 },
    { 255, 0, 0 },
    { 100, 0, 0 },
}};
const ColorMap ColorMap::Classic
{{
    { 30, 77, 203 },
    { 25, 60, 192 },
    { 45, 117, 220 },
    { 204, 108, 191 },
    { 196, 57, 178 },
    { 198, 33, 24 },
}};
const ColorMap ColorMap::GrayScale
{{
    { 255, 255, 255 },
    { 0, 0, 0 },
}};
const ColorMap ColorMap::InvGrayScale
{{
    { 0, 0, 0 },
    { 255, 255, 255 },
}};
const ColorMap ColorMap::Biomes
{{
    { 0, 0, 204 },
    { 204, 230, 255 },
    { 255, 255, 153 },
    { 170, 255, 128 },
    { 0, 153, 0 },
    { 230, 242, 255 },
}};
const ColorMap ColorMap::Cold
{{
    { 230, 247, 255 },
    { 0, 92, 230 },
    { 0, 179, 179 },
    { 0, 51, 153 },
    { 0, 5, 15 }
}};
const ColorMap ColorMap::Warm
{{
    { 255, 255, 230 },
    { 255, 204, 0 },
    { 255, 136, 77 },
    { 255, 51, 0 },
    { 128, 0, 0 },
    { 10, 0, 0}
}};


NAMESPACE_END(image)
NAMESPACE_END(yuki)