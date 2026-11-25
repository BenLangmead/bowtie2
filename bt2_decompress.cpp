#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>
#include "bt2_decompress.h"

#ifdef WITH_ZSTD

zstdStrm *zstdStrmInit(uint32_t *header)
{
	zstdStrm *s;

	if ((s = (zstdStrm *)malloc(sizeof(zstdStrm))) == NULL)
                return NULL;
        s->back  = EOF;
        s->err   = 0;
        s->i_len = ZSTD_DStreamInSize();
        s->o_len = 0;
        s->i_pos = 0;
        s->o_pos = 0;
        s->i_buf = NULL;
        if (header) {
                s->data = *header;
        }
        if ((s->i_buf = (uint8_t *)malloc(s->i_len)) == NULL) {
		zstdClose(s);
                return NULL;
        }

        if ((s->strm = ZSTD_createDCtx()) == NULL) {
                zstdClose(s);
                return NULL;
        }

        return s;
}

zstdStrm *zstdFileOpen(FILE *fp, uint32_t *header) {
	zstdStrm *s;

	if (fp == NULL)
		return NULL;
	s = zstdStrmInit(header);
	if (s == NULL)
		return NULL;
        if (fp == NULL) {
		zstdClose(s);
		return NULL;
        }
        s->fp = fp;

	return s;
}

size_t zstdRead(zstdStrm *s, uint8_t *buf, size_t len) {
        if (s == NULL || buf == NULL || len == 0) {
                return 0;
        }

        s->o_buf = buf;
        s->o_len = len;
        s->o_pos = 0;

        // If we've already hit EOF on the file and have no more input
        if (s->i_pos >= s->i_len && feof(s->fp)) {
                return 0;
        }

        while (s->o_pos < len) {
                bool need_more = s->i_pos == 0 || s->i_pos == s->i_len || (s->err != 0 && !ZSTD_isError(s->err));
                size_t prev_o_pos = s->o_pos;

                if (need_more) {
                        // Shift any remaining data to the start of the buffer
                        if (s->i_pos > 0 && s->i_pos < s->i_len) {
                                size_t remaining = s->i_len - s->i_pos;
                                memmove(s->i_buf, s->i_buf + s->i_pos, remaining);
                                s->i_pos = 0;
                                s->i_len = remaining;
                        } else {
                                s->i_pos = 0;
                                s->i_len = 0;
                        }

                        // Read more data from file
                        if (!feof(s->fp)) {
                                // Check if we got the header data that was passed in
                                if (s->data != 0 && s->i_len == 0) {
                                        uint8_t *h = (uint8_t *)&s->data;
                                        for (size_t i = 0; i < sizeof(s->data); i++) {
                                                s->i_buf[i] = h[i];
                                        }
                                        s->i_len = sizeof(s->data);
                                        s->data = 0;
                                }

                                size_t bytes_to_read = ZSTD_DStreamInSize() - s->i_len;
                                size_t nread = fread(s->i_buf + s->i_len, 1, bytes_to_read, s->fp);
                                s->i_len += nread;
                        }
                }

                // If we have no input data and are at EOF, stop
                if (s->i_len == 0 && feof(s->fp)) {
                        break;
                }

                ZSTD_inBuffer input = { s->i_buf + s->i_pos, s->i_len - s->i_pos, 0 };
                ZSTD_outBuffer output = { buf + s->o_pos, len - s->o_pos, 0 };

                s->err = ZSTD_decompressStream(s->strm, &output, &input);

                s->i_pos += input.pos;
                s->o_pos += output.pos;

                // Check for errors
                if (ZSTD_isError(s->err)) {
                        return 0;
                }

                // If no progress was made and we didn't get new input, we're done
                if (prev_o_pos == s->o_pos && s->i_pos >= s->i_len && feof(s->fp)) {
                        break;
                }
        }

        return s->o_pos;
}

int zstdClose(zstdStrm *s)
{
        if (s == NULL)
                return -1;
        if (s->fp != NULL && s->fp != stdin)
                fclose(s->fp);
        if (s->i_buf != NULL)
                free(s->i_buf);
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
        s->back = EOF;
        s->i_len = 0;
        s->i_pos = 0;
        ZSTD_DCtx_reset(s->strm, ZSTD_reset_session_only);
        return fseek(s->fp, 0, SEEK_SET);
}

#endif

GzipReader *gzipReaderInit(FILE *fp, uint32_t *header) {
        GzipReader *gr = NULL;
        if ((gr = (GzipReader *)malloc(sizeof(GzipReader))) == NULL)
                return NULL;

        memset(gr, 0, sizeof(*gr));
        gr->fp = fp;
        gr->inbuf_size = 64 * 1024;
        gr->inbuf = (unsigned char *)malloc(gr->inbuf_size);
        gr->inbuf_pos = 0;
        gr->inbuf_avail = 0;
        gr->eof = 0;
        gr->z_err = Z_OK;

        // Initialize inflate
        gr->strm.zalloc = Z_NULL;
        gr->strm.zfree = Z_NULL;
        gr->strm.opaque = Z_NULL;
        gr->strm.avail_in = 0;
        gr->strm.next_in = Z_NULL;
        if (header) {
                gr->data = *header;
        }

        int ret = inflateInit2(&gr->strm, 16 + MAX_WBITS); // +16 for gzip decoding
        if (ret != Z_OK) {
                free(gr->inbuf);
                gr->inbuf = NULL;
                gr->z_err = ret;
                free(gr);
                gr = NULL;
        }

        return gr;
}

void gzClose(GzipReader *gr) {
        if (gr != NULL) {
                if (gr->inbuf) free(gr->inbuf);
                inflateEnd(&gr->strm);
                memset(gr, 0, sizeof(*gr));
                if (gr->fp != NULL && gr->fp != stdin) {
                        fclose(gr->fp);
                }
                free(gr);
                gr = NULL;
        }
}

// Fill the input buffer from the file
int refillInput(GzipReader *gr) {
        if (gr->eof) return Z_STREAM_END;

        uint8_t *buf = gr->inbuf;
        size_t len = 0;

        if (gr->data != 0) {
                uint8_t *h = (uint8_t *)&gr->data;
                for (size_t i = 0; i < sizeof(gr->data); i++) {
                        gr->inbuf[i] = h[i];
                        len += 1;
                }

                buf += len;
                gr->data = 0;
        }
        size_t bytes_read = fread((void *)buf, 1, gr->inbuf_size - len, gr->fp) + len;
        if (bytes_read == 0) {
                gr->eof = 1;
                return Z_BUF_ERROR; // No more data available
        }

        gr->inbuf_pos = 0;
        gr->inbuf_avail = bytes_read;
        gr->strm.avail_in = gr->inbuf_avail;
        gr->strm.next_in = gr->inbuf;

        return Z_OK;
}

size_t gzRead(GzipReader *gr, uint8_t *outbuf, size_t outbuf_size) {
        if (gr->strm.avail_in == 0 && !gr->eof) {
                if (refillInput(gr) != Z_OK) {
                        // If we couldn't fill and have no data, check if we're at EOF
                        if (gr->eof && gr->strm.avail_in == 0) {
                                return Z_STREAM_END;
                        }
                }
        }

        gr->strm.avail_out = outbuf_size;
        gr->strm.next_out = outbuf;
        size_t outbuf_avail = 0;

        int ret = Z_OK;
        while (gr->strm.avail_out > 0) {
                ret = inflate(&gr->strm, Z_NO_FLUSH);

                outbuf_avail = outbuf_size - gr->strm.avail_out;

                if (ret == Z_STREAM_END) {
                        ret = Z_OK;
                        break;
                } else if (ret == Z_BUF_ERROR) {
                        if (gr->eof) {
                                // File is truncated
                                break;
                        }
                        if (refillInput(gr) != Z_OK) {
                                break;
                        }
                } else if (ret != Z_OK) {
                        break;
                }
        }

        gr->z_err = ret;
        return outbuf_avail;
}

// int gzGetc(GzipReader *gr) {
//         if (gr->outbuf_avail == 0) {
//                 int err = gzRead(gr);
//                 if (err != Z_OK && err != Z_STREAM_END) {
//                         return EOF;
//                 }
//         }

//         if (gr->outbuf_pos < gr->outbuf_avail) {
//                 return (unsigned char)gr->outbuf[gr->outbuf_pos++];
//         }

//         return EOF;
// }

// int gzUnGetc(int c, GzipReader *gr) {
//         if (gr->outbuf_avail == 0) {
//                 return EOF;
//         }
//         gr->outbuf[--gr->outbuf_pos] = c;

//         return 0;
// }

void gzRewind(GzipReader *gr) {
        FILE *fp = gr->fp;
        if (fp != NULL && fp != stdin) {
                // if (gr->inbuf) free(gr->inbuf);
                // inflateEnd(&gr->strm);
                gr->inbuf_pos = 0;
                gr->inbuf_avail = 0;

                inflateReset(&gr->strm);
                fseek(gr->fp, 0, SEEK_SET);
        }
}

const char *gzErrorString(GzipReader *gr) {
        if (gr->z_err == Z_STREAM_END) {
                return "End of file";
        } else if (gr->z_err == Z_BUF_ERROR) {
                return "Truncated gzip stream";
        } else if (gr->z_err == Z_DATA_ERROR) {
                return "Corrupt gzip data";
        } else if (gr->z_err == Z_MEM_ERROR) {
                return "Memory error";
        } else if (gr->z_err == Z_VERSION_ERROR) {
                return "zlib version mismatch";
        }
        return "Unknown error";
}
