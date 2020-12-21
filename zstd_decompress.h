#ifndef _ZSTD_DECOMPRESS_
#define _ZSTD_DECOMPRESS_

#include <stdio.h>
#include <zstd.h>

typedef struct zstd_stream zstdStrm;

struct zstd_stream {
        int           back;
        size_t        o_len;
        size_t        o_pos;
        size_t        i_len;
        size_t        i_pos;
        void         *i_buf;
        void         *o_buf;
        ZSTD_DStream *strm;
        FILE         *fp;
};

zstdStrm *initDStream(zstdStrm *s);
zstdStrm *zstdOpen(const char *fn);
zstdStrm *zstdFdOpen(int fd);
int zstdDecompress(zstdStrm *s);
int zstdRead(zstdStrm *s, void *buf, size_t len);
int zstdUngetc(int c, zstdStrm *s);
int zstdGetc(zstdStrm *s);
int zstdClose(zstdStrm *s);
int zstdRewind(zstdStrm *s);

#endif // end _ZSTD_DECOMPRESS_
