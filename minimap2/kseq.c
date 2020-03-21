#include "kseq.h"
inline kstream_t *ks_init(gzFile f)
{
    kstream_t *ks = (kstream_t*)calloc(1, sizeof(kstream_t));
    ks->f = f; ks->bufsize = 16384;
    ks->buf = (unsigned char*)malloc(16384);
    return ks;
}

inline void ks_destroy(kstream_t *ks)
{
    if (!ks) return;
    free(ks->buf);
    free(ks);
}

inline kseq_t *kseq_init(gzFile fd, runtime_sequences *seq)
{
    kseq_t *s = (kseq_t*)calloc(1, sizeof(kseq_t));
    s->f = ks_init(fd);
    s->s = seq;
    return s;
}

inline void kseq_destroy(kseq_t *ks)
{
    if (!ks) return;
    free(ks->name.s); free(ks->comment.s); free(ks->seq.s);	free(ks->qual.s);
    ks_destroy(ks->f);
    free(ks);
}


int ks_getuntil2(kstream_t *ks, int delimiter, kstring_t *str, int *dret, int append)
{
    if (dret) *dret = 0;
    str->l = append? str->l : 0;
    if (ks->begin >= ks->end && ks->is_eof) return -1;
    for (;;) {
        int i;
        if (ks->begin >= ks->end) {
            if (!ks->is_eof) {
                ks->begin = 0;
                ks->end = gzread(ks->f, ks->buf, ks->bufsize);
                if (ks->end < ks->bufsize) ks->is_eof = 1;
                if (ks->end == 0) break;
            } else break;
        }
        if (delimiter == KS_SEP_LINE) {
            for (i = ks->begin; i < ks->end; ++i)
                if (ks->buf[i] == '\n') break;
        } else if (delimiter > KS_SEP_MAX) {
            for (i = ks->begin; i < ks->end; ++i)
                if (ks->buf[i] == delimiter) break;
        } else if (delimiter == KS_SEP_SPACE) {
            for (i = ks->begin; i < ks->end; ++i)
                if (isspace(ks->buf[i])) break;
        } else if (delimiter == KS_SEP_TAB) {
            for (i = ks->begin; i < ks->end; ++i)
                if (isspace(ks->buf[i]) && ks->buf[i] != ' ') break;
        } else i = 0; /* never come to here! */
        if (str->m - str->l < (size_t)(i - ks->begin + 1)) {
            str->m = str->l + (i - ks->begin) + 1;
            kroundup32(str->m);
            str->s = (char*)realloc(str->s, str->m);
        }
        memcpy(str->s + str->l, ks->buf + ks->begin, i - ks->begin);
        str->l = str->l + (i - ks->begin);
        ks->begin = i + 1;
        if (i < ks->end) {
            if (dret) *dret = ks->buf[i];
            break;
        }
    }
    if (str->s == 0) {
        str->m = 1;
        str->s = (char*)calloc(1, 1);
    } else if (delimiter == KS_SEP_LINE && str->l > 1 && str->s[str->l-1] == '\r') --str->l;
    str->s[str->l] = '\0';
    return str->l;
}

int kseq_read(kseq_t *seq) {
    if (seq->s != 0) {
        if (seq->s->cur == seq->s->size) {
            return -1;
        }
        seq->name.s = seq->s->names[seq->s->cur];
        seq->name.l = strlen(seq->s->names[seq->s->cur]);
        seq->name.m = (strlen(seq->s->names[seq->s->cur]) + 255) & ~255u;
        seq->seq.s = seq->s->s[seq->s->cur];
        seq->seq.l = strlen(seq->s->s[seq->s->cur]);
        seq->seq.m = (strlen(seq->s->s[seq->s->cur]) + 255) & ~255u;
        seq->s->cur += 1;
        return seq->seq.l;
    }
    int c;
    kstream_t *ks = seq->f;
    if (seq->last_char == 0) { /* then jump to the next header line */
        while ((c = ks_getc(ks)) != -1 && c != '>' && c != '@');
        if (c == -1) return -1; /* end of file */
        seq->last_char = c;
    } /* else: the first header char has been read in the previous call */
    seq->comment.l = seq->seq.l = seq->qual.l = 0; /* reset all members */
    if (ks_getuntil(ks, 0, &seq->name, &c) < 0) return -1; /* normal exit: EOF */
    if (c != '\n') ks_getuntil(ks, KS_SEP_LINE, &seq->comment, 0); /* read FASTA/Q comment */
    if (seq->seq.s == 0) { /* we can do this in the loop below, but that is slower */
        seq->seq.m = 256;
        seq->seq.s = (char*)malloc(seq->seq.m);
    }
    while ((c = ks_getc(ks)) != -1 && c != '>' && c != '+' && c != '@') {
        if (c == '\n') continue; /* skip empty lines */
        seq->seq.s[seq->seq.l++] = c; /* this is safe: we always have enough space for 1 char */
        ks_getuntil2(ks, KS_SEP_LINE, &seq->seq, 0, 1); /* read the rest of the line */
    }
    if (c == '>' || c == '@') seq->last_char = c; /* the first header char has been read */
    if (seq->seq.l + 1 >= seq->seq.m) { /* seq->seq.s[seq->seq.l] below may be out of boundary */
        seq->seq.m = seq->seq.l + 2;
        kroundup32(seq->seq.m); /* rounded to the next closest 2^k */
        seq->seq.s = (char*)realloc(seq->seq.s, seq->seq.m);
    }
    seq->seq.s[seq->seq.l] = 0;	/* null terminated string */
    if (c != '+') return seq->seq.l; /* FASTA */
    if (seq->qual.m < seq->seq.m) {	/* allocate memory for qual in case insufficient */
        seq->qual.m = seq->seq.m;
        seq->qual.s = (char*)realloc(seq->qual.s, seq->qual.m);
    }
    while ((c = ks_getc(ks)) != -1 && c != '\n'); /* skip the rest of '+' line */
    if (c == -1) return -2; /* error: no quality string */
    while (ks_getuntil2(ks, KS_SEP_LINE, &seq->qual, 0, 1) >= 0 && seq->qual.l < seq->seq.l);
    seq->last_char = 0;	/* we have not come to the next header line */
    if (seq->seq.l != seq->qual.l) return -2; /* error: qual string is of a different length */
    return seq->seq.l;
}
