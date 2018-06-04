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
    
    data_.clear();
    std::ifstream fin(path, std::ios::binary);

    Header header;
    fin.read((char *)(&header), sizeof(header));
    YUKI_DEBUG("WAV PCM File Header: \n");
    YUKI_DEBUG("File Type: %c%c%c%c\n", header.chunk_id[0], header.chunk_id[1], header.chunk_id[2], header.chunk_id[3]);
    YUKI_DEBUG("File Size: %ld\n", header.chunk_size);
    YUKI_DEBUG("Format Name: %c%c%c%c\n", header.format[0], header.format[1], header.format[2], header.format[3]);
    YUKI_DEBUG("Format chunk: %c%c%c%c\n", header.sub_chunk1_id[0], header.sub_chunk1_id[1], header.sub_chunk1_id[2], header.sub_chunk1_id[3]);
    YUKI_DEBUG("Format Length: %ld\n", header.sub_chunk1_size);
    YUKI_DEBUG("Format Type: %hd\n", header.audio_format);
    YUKI_DEBUG("Number of channels: %d\n", header.num_channels);
    YUKI_DEBUG("Sample Rate: %ld\n", header.sample_rate);
    YUKI_DEBUG("Sample Rate * Bits/Sample * Channels / 8: %ld\n", header.byte_rate);
    YUKI_DEBUG("Bits per Sample * Channels / 8: %hd\n", header.block_align);
    YUKI_DEBUG("Bits per Sample: %hd\n", header.bits_per_sample);

    /* check if the header is valid */ {    
        if (header.sub_chunk1_size < 16)
        {
            YUKI_ERROR("[WAV]: Wrong header %s.\n", path.c_str());
            return false;
        }
        if (header.audio_format != 1)
        {
            YUKI_ERROR("[WAV]: Not PCM %s\n", path.c_str());
            return false;
        }
    }

    // skip extra infomations
    if (header.sub_chunk1_size > 16)
        fin.seekg(header.sub_chunk1_size - 16, std::ios_base::cur);
        

    Chunk chunk;
    for (;;)
    {
        fin.read((char *)&(chunk), sizeof(Chunk));
        YUKI_DEBUG("%c%c%c%c\t%d\n", chunk.id[0], chunk.id[1], chunk.id[2], chunk.id[3], chunk.size);
        if (*(uint32_t *)&chunk.id == 0x61746164)  // == "data"
            break;
        fin.seekg(chunk.size, std::ios_base::cur);  // skip not important subchunk
    }



    fin.close();

    return true;
}

bool WAVPCM::write(const std::string &path)
{

    return true;
}    

NAMESPACE_END(audio)
NAMESPACE_END(yuki)