#include <cstdio>
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
        s->back  = EOF;
        s->i_len = ZSTD_DStreamInSize();
        s->o_len = 0;
        s->i_pos = 0;
        s->o_pos = 0;
        s->i_buf = NULL;
        if ((s->i_buf = malloc(s->i_len)) == NULL) {
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

int zstdRead(zstdStrm *s, void *buf, size_t len)
{
	if (s->i_pos == s->i_len && feof(s->fp))
		return -1;

        s->o_pos = 0;
        s->o_buf = buf;
        s->o_len = len;

        if (s->back != EOF) {
                *((char *)(s->o_buf) + s->o_pos++) = (char)s->back;
                s->back = EOF;
        }
        while (s->o_pos < len) {
                int ret;

                if ((s->i_pos == 0 || s->i_pos == s->i_len) && !feof(s->fp)) {
                        int nread;
			size_t nleft = ZSTD_DStreamInSize();

			s->i_len = s->i_pos = 0;
			do {
				nread = fread((char *)s->i_buf + s->i_len, 1, nleft, s->fp);
				s->i_len += nread;
				if (feof(s->fp))
					break;
				nleft -= nread;
			} while (nleft > 0);
                }
                ret = zstdDecompress(s);
                if (ret == 0)
                        break;
                if (ZSTD_isError(ret))
                        return ret;
		else if (s->o_pos == 0)
                        return -1; // we have reached EOF

        }

        return s->o_pos;
}

int zstdUngetc(int c, zstdStrm *s)
{
        if (s == NULL || c == EOF || s->back != EOF)
                return EOF;
        s->o_pos--;
        s->back = c;

        return c;
}

int zstdGetc(zstdStrm *s)
{
        int ret;


	if (s->o_len == 0) {
		ret = zstdRead(s, s->o_buf, ZSTD_DStreamOutSize());
		if (ret < 0)
			return ret;
		s->o_len = s->o_pos;
		s->o_pos = 0;
	}
	s->o_len--;
	return ((char *)s->o_buf)[s->o_pos++];
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
	bzero(s, sizeof(zstdStrm));
	free(s);

        return 0;
}

int zstdRewind(zstdStrm *s)
{
        if (s == NULL)
                return -1;
        s->back = EOF;
        s->i_len = 0;
        s->i_pos = 0;
        ZSTD_DCtx_reset(s->strm, ZSTD_reset_session_only);
        return fseek(s->fp, 0, SEEK_SET);
}
