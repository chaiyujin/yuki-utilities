#pragma once

extern "C" {
    #ifndef __STDC_CONSTANT_MACROS
    #define __STDC_CONSTANT_MACROS
    #endif

    #include <libavutil/opt.h>
    #include <libavutil/dict.h>
    #include <libavutil/eval.h>
    #include <libavutil/fifo.h>
    #include <libavutil/time.h>
    #include <libavutil/avutil.h>
    #include <libavutil/pixfmt.h>
    #include <libavutil/pixdesc.h>
    #include <libavutil/rational.h>
    #include <libavutil/timestamp.h>
    #include <libavutil/samplefmt.h>
    #include <libavutil/hwcontext.h>
    #include <libavutil/audio_fifo.h>
    #include <libavutil/threadmessage.h>
    #include <libavutil/channel_layout.h>
    #include <libavformat/avio.h>
    #include <libavformat/avformat.h>
    #include <libavcodec/avcodec.h>
    #include <libavfilter/avfilter.h>
    #include <libavfilter/buffersink.h>
    #include <libavfilter/buffersrc.h>
    #include <libswscale/swscale.h>
    #include <libswresample/swresample.h>
    #include <libavdevice/avdevice.h>
};

#include <memory>

namespace yuki {
namespace media {

    enum Type
    {
        Audio = 0,
        Video = 1
    };

    class Frame
    {
    protected:
        std::shared_ptr<uint8_t>    data_;
        double                      ts_;
    public:
        Frame() : data_(NULL), ts_(0) {}
        virtual ~Frame() {}
        bool is_null() const { return !data_; }
        uint8_t *data() { return data_.get(); }
        const uint8_t *data() const { return data_.get(); }
    };

    class VideoFrame : public Frame
    {
        int                            w_, h_, pixels_;
        bool                        is_depth_;
    public:
        VideoFrame() : Frame(), w_(0), h_(0), pixels_(0), is_depth_(false) {}
        VideoFrame(AVFrame *frame, double pts, bool depth)
            : Frame() { from_avframe(frame, pts, depth); }

        int width() const { return w_; }
        int height() const { return h_; }
        int pixels() const { return pixels_; }
        bool is_depth() const { return is_depth_; }

        // only RGBA or DepthRVL
        void from_avframe(AVFrame *frame, double pts, bool depth)
        {
            is_depth_ = depth;
            int size = frame->width * frame->height * 2;    // uint16_t
            if (!depth) size *= 2;                            // RGBA
            data_.reset(new uint8_t[size]);
            memcpy(data_.get(), frame->data[0], size);
            w_ = frame->width;
            h_ = frame->height;
            pixels_ = w_ * h_;
            ts_ = pts;
        }
    };

    struct AudioFrame : public Frame
    {
        int                            nb_samples_;
        int                            byte_per_sample_;
    public:
        AudioFrame() : Frame(), nb_samples_(0), byte_per_sample_(0) {}
        AudioFrame(AVFrame *frame, double pts, int nb_samples, int byte_per_sample)
            : Frame() { from_avframe(frame, pts, nb_samples, byte_per_sample); }

        int nb_samples() const { return nb_samples_; }
        int byte_per_sample() const { return byte_per_sample_; }

        void from_avframe(AVFrame *frame, double pts, int nb_samples, int byte_per_sample)
        {
            nb_samples_ = nb_samples;
            byte_per_sample_ = byte_per_sample;
            data_.reset(new uint8_t[nb_samples_ * byte_per_sample_]);
            memcpy(data_.get(), frame->data[0], nb_samples_ * byte_per_sample_);
            ts_ = pts;
        }
    };

    template <class T> class Stream;
    template <>
    class Stream<AudioFrame>
    {
    public:
        const media::Type Type = media::Type::Audio;

    };
}}
