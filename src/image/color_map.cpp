#include "color_map.h"


NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(image)

const ColorMap ColorMap::Jet {{
    { 0, 0, 255 },
    { 0, 255, 255 },
    { 255, 255, 0 },
    { 255, 0, 0 },
    { 50, 0, 0 },
}};

    // const static ColorMap Classic;
    // const static ColorMap GrayScale;
    // const static ColorMap InvGrayScale;
    // const static ColorMap Biomes;
    // const static ColorMap Cold;
    // const static ColorMap Warm;


NAMESPACE_END(image)
NAMESPACE_END(yuki)