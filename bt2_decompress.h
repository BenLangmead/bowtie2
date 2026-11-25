#ifndef _BT2_DECOMPRESS_
#define _BT2_DECOMPRESS_

#include <cstdint>
#include <stdio.h>
#include <zlib.h>
#ifdef WITH_ZSTD
#include <zstd.h>

typedef struct zstd_stream zstdStrm;

struct zstd_stream {
        int           back;
        size_t        err;      // last ZSTD error code, 0 if none
        size_t        o_len;
        size_t        o_pos;
        size_t        i_len;
        size_t        i_pos;
        uint8_t      *i_buf;
        uint8_t      *o_buf;
        uint32_t      data;
        ZSTD_DCtx    *strm;
        FILE         *fp;
};

zstdStrm *initDStream(zstdStrm *s);
zstdStrm *zstdStrmInit(uint32_t *header);
zstdStrm *zstdFileOpen(FILE *fp, uint32_t *header);
size_t zstdRead(zstdStrm *s, uint8_t *buf, size_t len);
int zstdClose(zstdStrm *s);
int zstdRewind(zstdStrm *s);
#endif

struct GzipReader {
        z_stream strm;
        FILE *fp;
        unsigned char *inbuf;   // Input buffer for compressed data
        size_t inbuf_size;
        size_t inbuf_pos;
        size_t inbuf_avail;
        uint32_t data;
        int eof;              // Set when we hit EOF on the file
        int z_err;            // Last zlib error
};

GzipReader *gzipReaderInit(FILE *fp, uint32_t *header);
void gzClose(GzipReader *gr);
int refillInput(GzipReader *gr);
size_t gzRead(GzipReader *gr, uint8_t *buf, size_t buf_size);
int gzGetc(GzipReader *gr);
int gzUnGetc(int c, GzipReader *gr);
void gzRewind(GzipReader *gr);
const char *gzErrorString(GzipReader *gr);

#endif // end _BT2_DECOMPRESS_
