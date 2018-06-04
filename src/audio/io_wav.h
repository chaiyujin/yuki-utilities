#pragma once
#include "core.h"
#include <string>
#include <fstream>
#include <vector>

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(audio)

class WAVPCM
{
    std::vector<std::vector<uint8_t>> data_;
public:
    struct Header
    {
        char            chunk_id[4];
        uint32_t        chunk_size;
        char            format[4];
        char            sub_chunk1_id[4];
        uint32_t        sub_chunk1_size;
        uint16_t        audio_format;
        uint16_t        num_channels;
        uint32_t        sample_rate;
        uint32_t        byte_rate;
        uint16_t        block_align;
        uint16_t        bits_per_sample;
    };

    struct Chunk
    {
        char            id[4];
        int32_t         size;
    };

    bool read(const std::string &path);
    bool write(const std::string &path);
};

NAMESPACE_END(audio)
NAMESPACE_END(yuki)