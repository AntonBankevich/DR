/* The MIT License

   Copyright (c) 2008, 2009, 2011 Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Last Modified: 05MAR2012 */

#ifndef AC_KSEQ_H
#define AC_KSEQ_H

#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <zlib.h>
#include "kstring.h"

#define KS_SEP_SPACE 0 // isspace(): \t, \n, \v, \f, \r
#define KS_SEP_TAB   1 // isspace() && !' '
#define KS_SEP_LINE  2 // line separator: "\n" (Unix) or "\r\n" (Windows)
#define KS_SEP_MAX   2

#ifndef klib_unused
#if (defined __clang__ && __clang_major__ >= 3) || (defined __GNUC__ && __GNUC__ >= 3)
#define klib_unused __attribute__ ((__unused__))
#else
#define klib_unused
#endif
#endif /* klib_unused */

typedef struct {
    char** s;
    char** names;
    int cur;
    int size;
} runtime_sequences;

typedef struct __kstream_t {
	int begin, end;
	int is_eof:2, bufsize:30;
	gzFile f;
	unsigned char *buf;
} kstream_t;

typedef struct {
	kstring_t name, comment, seq, qual;
	int last_char;
    kstream_t *f;
	runtime_sequences *s;
} kseq_t;

#define ks_eof(ks) (((ks)->s != 0 && (ks)->s->cur == (ks)->s->size) || ((ks)->s == 0 && (ks)->f->is_eof && (ks)->f->begin >= (ks)->f->end))
#define ks_rewind(ks) ((ks)->is_eof = (ks)->begin = (ks)->end = 0)

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

extern int ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append);

static inline klib_unused int ks_getc(kstream_t *ks)
{
    if (ks->is_eof && ks->begin >= ks->end) return -1;
    if (ks->begin >= ks->end) {
        ks->begin = 0;
        ks->end = gzread(ks->f, ks->buf, ks->bufsize);
        if (ks->end < ks->bufsize) ks->is_eof = 1;
        if (ks->end == 0) return -1;
    }
    return (int)ks->buf[ks->begin++];
}
static inline int ks_getuntil(kstream_t *ks, int delimiter, kstring_t *str, int *dret) \
{ return ks_getuntil2(ks, delimiter, str, dret, 0); }



/******************
 * FASTA/Q parser *
 ******************/

#define kseq_rewind(ks) ((ks)->last_char = (ks)->f->is_eof = (ks)->f->begin = (ks)->f->end = 0)

kstream_t *ks_init(gzFile f);
void ks_destroy(kstream_t *ks);
kseq_t *kseq_init(gzFile fd, runtime_sequences *seq);
void kseq_destroy(kseq_t *ks);

/* Return value:
   >=0  length of the sequence (normal)
   -1   end-of-file
   -2   truncated quality string
 */
int kseq_read(kseq_t *seq);


#endif
