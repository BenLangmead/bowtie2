#include <cstdio>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zstd.h>

#include "zstd_decompress.h"


zstdStrm *zstdStrmInit()
{
	zstdStrm *s;

	if ((s = (zstdStrm *)malloc(sizeof(zstdStrm))) == NULL)
                return NULL;
        s->i_len = 0;
        s->o_len = 0;
        s->i_pos = 0;
        s->o_pos = 0;
        s->i_buf = NULL;
        if ((s->i_buf = malloc(ZSTD_DStreamInSize())) == NULL) {
		zstdClose(s);
                return NULL;
        }
	if ((s->o_buf = malloc(ZSTD_DStreamOutSize())) == NULL) {
		zstdClose(s);
                return NULL;
        }

        if ((s->strm = ZSTD_createDCtx()) == NULL) {
                zstdClose(s);
                return NULL;
        }

        return s;
}

zstdStrm *zstdOpen(const char *fn) {
        FILE *fp;
        zstdStrm *s;

	s = zstdStrmInit();
	if (s == NULL)
		return NULL;
        fp = fopen(fn, "rb");
        if (fp == NULL) {
		zstdClose(s);
		return NULL;
        }
        s->fp = fp;

        return s;
}

zstdStrm *zstdFdOpen(int fd) {
	FILE *fp;
	zstdStrm *s;

	if (fd == -1)
		return NULL;
	s = zstdStrmInit();
	if (s == NULL)
		return NULL;
	fp = fdopen(fd, "rb");
        if (fp == NULL) {
		zstdClose(s);
		return NULL;
        }
        s->fp = fp;

	return s;
}

int zstdDecompress(zstdStrm *s)
{
        int ret;
        ZSTD_inBuffer in = { s->i_buf, s->i_len, s->i_pos };
	ZSTD_outBuffer out = { s->o_buf, s->o_len, s->o_pos };

	ret = ZSTD_decompressStream(s->strm, &out, &in);
	s->i_pos = in.pos;
	s->o_pos = out.pos;

        return ret;
}


int zstdGetUncompressedSize();

int zstdDecompressBuffer(zstdStrm *s, unsigned char *in_buf, size_t in_len, unsigned char *out_buf, size_t out_len) {
  ZSTD_inBuffer in = {in_buf, in_len, 0};
  ZSTD_outBuffer out = {out_buf, out_len, 0};

  return ZSTD_decompressStream(s->strm, &out, &in);
}

int zstdRead(zstdStrm *s, void *buf, size_t len)
{
	ZSTD_outBuffer out = { buf, len, 0 };

	while (out.pos < len) {
		if (s->i_pos == s->i_len) {
			s->i_len = fread(s->i_buf, 1, ZSTD_DStreamInSize(), s->fp);
			s->i_pos = 0;
		}
		size_t before = out.pos;
		ZSTD_inBuffer in = { s->i_buf, s->i_len, s->i_pos };
		size_t ret = ZSTD_decompressStream(s->strm, &out, &in);
		s->i_pos = in.pos;
		if (ZSTD_isError(ret))
			return -1; // EOF: any other negative value makes getc_wrapper retry forever
		// A step may consume input (a frame header, a skippable frame) without
		// output; the stream only ends once the file is read and nothing is left to flush.
		if (s->i_len == 0 && out.pos == before)
			break;
	}

	return out.pos > 0 ? (int)out.pos : -1;
}

int zstdUngetc(int c, zstdStrm *s)
{
	if (s == NULL || c == EOF || s->o_pos == 0)
		return EOF;
	((unsigned char *)s->o_buf)[--s->o_pos] = c;
	s->o_len++;

	return c;
}

int zstdGetc(zstdStrm *s)
{
	if (s->o_len == 0) {
		int ret = zstdRead(s, s->o_buf, ZSTD_DStreamOutSize());
		if (ret < 0)
			return ret;
		s->o_len = ret;
		s->o_pos = 0;
	}
	s->o_len--;
	return ((unsigned char *)s->o_buf)[s->o_pos++];
}

int zstdClose(zstdStrm *s)
{
        if (s == NULL)
                return -1;
        if (s->fp != NULL)
                fclose(s->fp);
        if (s->i_buf != NULL)
                free(s->i_buf);
        if (s->o_buf != NULL)
                free(s->o_buf);
        if (s->strm != NULL)
		ZSTD_freeDStream(s->strm);
#ifndef _WIN32
	bzero(s, sizeof(zstdStrm));
#endif
	free(s);

        return 0;
}

int zstdRewind(zstdStrm *s)
{
        if (s == NULL)
                return -1;
        s->i_len = 0;
        s->i_pos = 0;
        ZSTD_DCtx_reset(s->strm, ZSTD_reset_session_only);
        return fseek(s->fp, 0, SEEK_SET);
}
