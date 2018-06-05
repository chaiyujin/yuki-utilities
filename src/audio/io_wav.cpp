#include "io_wav.h"
#include <iostream>

NAMESPACE_BEGIN(yuki)
NAMESPACE_BEGIN(audio)

bool WAVPCM::read(const std::string &path)
{
    if (!std::ifstream(path).good())
    {
        YUKI_ERROR("[WAV] No such file: %s\n", path.c_str());
        return false;
    }
    
    std::ifstream fin(path, std::ios::binary);

    fin.read((char *)(&header_), sizeof(header_));
    YUKI_DEBUG("WAV PCM File header_: \n");
    YUKI_DEBUG("File Type: %c%c%c%c\n", header_.chunk_id[0], header_.chunk_id[1], header_.chunk_id[2], header_.chunk_id[3]);
    YUKI_DEBUG("File Size: %d\n", header_.chunk_size);
    YUKI_DEBUG("Format Name: %c%c%c%c\n", header_.format[0], header_.format[1], header_.format[2], header_.format[3]);
    YUKI_DEBUG("Format chunk: %c%c%c%c\n", header_.sub_chunk1_id[0], header_.sub_chunk1_id[1], header_.sub_chunk1_id[2], header_.sub_chunk1_id[3]);
    YUKI_DEBUG("Format Length: %d\n", header_.sub_chunk1_size);
    YUKI_DEBUG("Format Type: %d\n", header_.audio_format);
    YUKI_DEBUG("Number of channels: %d\n", header_.num_channels);
    YUKI_DEBUG("Sample Rate: %d\n", header_.sample_rate);
    YUKI_DEBUG("Sample Rate * Bits/Sample * Channels / 8: %d\n", header_.byte_rate);
    YUKI_DEBUG("Bits per Sample * Channels / 8: %d\n", header_.block_align);
    YUKI_DEBUG("Bits per Sample: %d\n", header_.bits_per_sample);

    /* check if the header_ is valid */ {    
        if (header_.sub_chunk1_size < 16)
        {
            YUKI_ERROR("[WAV]: Wrong header_ %s.\n", path.c_str());
            return false;
        }
        if (header_.audio_format != 1)
        {
            YUKI_ERROR("[WAV]: Not PCM %s\n", path.c_str());
            return false;
        }
    }

    // skip extra infomations
    if (header_.sub_chunk1_size > 16)
        fin.seekg(header_.sub_chunk1_size - 16, std::ios_base::cur); 

    Chunk chunk;
    for (;;)
    {
        fin.read((char *)&(chunk), sizeof(Chunk));
        YUKI_DEBUG("%c%c%c%c\t%d\n", chunk.id[0], chunk.id[1], chunk.id[2], chunk.id[3], chunk.size);
        if (*(uint32_t *)&chunk.id == 0x61746164)  // == "data"
            break;
        fin.seekg(chunk.size, std::ios_base::cur);  // skip not important subchunk
    }

    int samples_count = chunk.size * 8 / header_.bits_per_sample / header_.num_channels;

    YUKI_DEBUG("Samples count = %i\n", samples_count);
    data_.clear();
    for (uint16_t i = 0; i < header_.num_channels; ++i)
    {
        data_.emplace_back(samples_count);
    }

    for (int i = 0; i < samples_count; ++i)
    {
        for (int ch = 0; ch < header_.num_channels; ++ch)
        {
            if (fin.eof())
            {
                YUKI_ERROR("[WAV] Reach end of file.\n");
                return false;
            }

            switch (header_.bits_per_sample)
            {
            case 8:
            {
                uint8_t val;
                fin.read((char *)&val, 1);
                data_[ch][i] = (((double)val - 128.0) / 128.0);
                break;
            }
            case 16:
            {
                int16_t val;
                fin.read((char *)&val, 2);
                data_[ch][i] = ((double)val / 32768.);
                break;
            }
            default:
                YUKI_ERROR("[WAV] unsupport bit_depth %d\n", header_.bits_per_sample);
                return false;
            }
        }
    }

    YUKI_DEBUG("%lu channels\n", data_.size());
    for (size_t i = 0; i < data_.size(); ++i)
        YUKI_DEBUG("%lu samples\n", data_[i].size());

    fin.close();

    return true;
}

bool WAVPCM::write(const std::string &path)
{

    return true;
}    

NAMESPACE_END(audio)
NAMESPACE_END(yuki)