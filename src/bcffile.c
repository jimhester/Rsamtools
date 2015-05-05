#include <string.h>
#include <errno.h>
#include "kstring.h"
#include "bgzf.h"
#include "samtools_patch.h"
#include "vcf.h"
#include "bcffile.h"
#include "utilities.h"

#define smkChar(x) ((x) ? mkChar(x) : R_NaString)

struct typemap {
    uint32_t *id;
    int *type, *mult;
    int len;
};

enum {
    BCF_HDR_REF = 0, BCF_HDR_SAMPLE, BCF_HDR_HEADER, BCF_HDR_LAST
};

static const char *BCF_HDR_NM[] = { "Reference", "Sample", "Header" };

enum {
    BCF_TID = 0, BCF_POS, BCF_ID, BCF_REF, BCF_ALT, BCF_QUAL,
    BCF_FLT, BCF_INFO, BCF_FMT, BCF_GENO, BCF_RECS_PER_RANGE,
    BCF_LAST
};

enum {
    BCF_TYPE_Integer = 1, BCF_TYPE_Float, BCF_TYPE_Character,
    BCF_TYPE_String, BCF_TYPE_Flag, BCF_TYPE_Last
};

static SEXP BCFFILE_TAG = NULL;
static const int BCF_BUFSIZE_GROW = 100000;	/* initial # records */

static vcfFile *_bcf_tryopen(const char *fname, const char *mode)
{
    return vcf_open(fname, mode);
}

static hts_idx_t *_bcf_idx_load(const char *fname)
{
    return bcf_index_load(fname);
}

static void _bcf_close(vcfFile *bcf, int errmsg)
{
    int err = vcf_close(bcf);
    if ((0 != err) && errmsg) {
        if (Z_ERRNO == err) {
            err = errno;
            Rf_error("_bcf_close file system error (%d): %s",
                     err, strerror(err));
        }
        Rf_error("_bcf_close error (%d)", err);
    }
}

static void _bcffile_close(SEXP ext)
{
    _BCF_FILE *bfile = BCFFILE(ext);
    if (NULL != bfile->file)
        vcf_close(bfile->file);
    if (NULL != bfile->index)
        hts_idx_destroy(bfile->index);
    bfile->file = NULL;
    bfile->index = NULL;
}

static void _bcffile_finalizer(SEXP ext)
{
    if (NULL == R_ExternalPtrAddr(ext))
        return;
    _bcffile_close(ext);
    _BCF_FILE *bfile = BCFFILE(ext);
    Free(bfile);
    R_SetExternalPtrAddr(ext, NULL);
}

SEXP bcffile_init()
{
    BCFFILE_TAG = install("BcfFile");
    return R_NilValue;
}

SEXP bcffile_open(SEXP filename, SEXP indexname, SEXP filemode)
{
    _checknames(filename, indexname, filemode);

    _BCF_FILE *bfile = Calloc(1, _BCF_FILE);

    bfile->file = NULL;
    if (0 != Rf_length(filename)) {
        const char *cfile = translateChar(STRING_ELT(filename, 0));
        bfile->file = _bcf_tryopen(cfile, CHAR(STRING_ELT(filemode, 0)));
        if (NULL == bfile->file) {
            Free(bfile);
            Rf_error("'open' BCF failed\n  filename: %s", cfile);
        }
    }

    bfile->index = NULL;
    if (0 != Rf_length(indexname) && bfile->file->format.format != vcf) {
        const char *cindex = translateChar(STRING_ELT(indexname, 0));
        bfile->index = _bcf_idx_load(cindex);
        if (NULL == bfile->index) {
            _bcf_close(bfile->file, 0);
            Free(bfile);
            Rf_error("'open' BCF index failed\n  indexname: %s\n", cindex);
        }
    }

    SEXP ext = PROTECT(R_MakeExternalPtr(bfile, BCFFILE_TAG, filename));
    R_RegisterCFinalizerEx(ext, _bcffile_finalizer, TRUE);
    UNPROTECT(1);

    return ext;
}

SEXP bcffile_close(SEXP ext)
{
    _checkext(ext, BCFFILE_TAG, "close");
    _bcffile_close(ext);
    return ext;
}

SEXP bcffile_isopen(SEXP ext)
{
    SEXP ans = ScalarLogical(FALSE);
    if (NULL != BCFFILE(ext)) {
        _checkext(ext, BCFFILE_TAG, "isOpen");
        if (BCFFILE(ext)->file)
            ans = ScalarLogical(TRUE);
    }
    return ans;
}

SEXP bcffile_isvcf(SEXP ext)
{
    SEXP ans = ScalarLogical(FALSE);
    if (NULL != BCFFILE(ext)) {
        _checkext(ext, BCFFILE_TAG, "isVcf");
        if (BCFFILE(ext)->file && BCFFILE(ext)->file->format.format == vcf)
            ans = ScalarLogical(TRUE);
    }
    return ans;
}

/* implementation */
static int _bcf_ans_grow(SEXP ans, R_len_t sz, int n_smpl)
{
    R_len_t n = sz;
    if (0 <= sz)
        n += Rf_length(VECTOR_ELT(ans, BCF_TID));
    else
        n *= -1;
    for (int i = 0; i < BCF_LAST; ++i) {
        SEXP elt = VECTOR_ELT(ans, i);
        switch (i) {
        case BCF_GENO:
            for (int i = 0; i < Rf_length(elt); ++i) {
                SEXP g = VECTOR_ELT(elt, i);
                SEXP dim = GET_DIM(g);
                if (R_NilValue == dim) {
                    g = Rf_lengthgets(g, n);
                    SET_VECTOR_ELT(elt, i, g);	/* protect */
                } else {
                    PROTECT(dim);
                    g = Rf_lengthgets(g, n_smpl * n);
                    SET_VECTOR_ELT(elt, i, g);	/* protect */
                    INTEGER(dim)[0] = n_smpl;
                    INTEGER(dim)[1] = n;
                    SET_DIM(g, dim);
                    UNPROTECT(1);
                }
            }
            break;
        case BCF_RECS_PER_RANGE:
            break;
        default:
            elt = Rf_lengthgets(elt, n);
            SET_VECTOR_ELT(ans, i, elt);	/* protect */
            break;
        }
    }
    return n;
}

/* static int _bcf_sync1(bcf1_t * b) */
/* { */
/*     /\* called when no FORMAT / GENO fields present; */
/*        from bcftools/bcf.c:bcf_sync *\/ */
/*     char *p, *tmp[5]; */
/*     int n; */
/*     for (p = b->str, n = 0; p < b->str + b->l_str; ++p) { */
/*         if (*p == 0 && p + 1 != b->str + b->l_str) { */
/*             if (n == 5) { */
/*                 ++n; */
/*                 break; */
/*             } else */
/*                 tmp[n++] = p + 1; */
/*         } */
/*     } */
/*     if (n != 4) */
/*         return -1; */
/*     b->ref = tmp[0]; */
/*     b->alt = tmp[1]; */
/*     b->flt = tmp[2]; */
/*     b->info = tmp[3]; */
/*     b->fmt = 0; */
/*     if (*b->alt == 0) */
/*         b->n_alleles = 1; */
/*     else { */
/*         for (p = b->alt, n = 1; *p; ++p) */
/*             if (*p == ',') */
/*                 ++n; */
/*         b->n_alleles = n + 1; */
/*     } */
/*     b->n_gi = 0; */
/*     return 0; */
/* } */

/* NOTE: in bcf2 format / geno fields are encoded as vectors of values
 * by field across samples rather than by sample (e.g., GT:0/0, GQ:48)
 * so that, for example, GT field across all samples in the bcf1_t
 * record are encoded as GT [valforsample1, valforsample2,
 * valforsample3] */

/* convert GENO fields (typically called fmt / format fields in htslib
 * inline documentation) to SEXP; */
static void _bcf_gi2sxp(SEXP geno, const int i_rec, const bcf_hdr_t * h,
                        bcf1_t * b)
{
    SEXP nm = GET_NAMES(geno);
    /* FIXME: more flexible geno not supported by bcftools */

    if (b->n_fmt == 0)
        return;

    int nvals = 0, fieldmemsz = 0;
    char **strdest = NULL;
    int32_t *intdest = NULL;
    float *floatdest = NULL;

    /* from bcftools/bcf.c */
    for (int i = 0; i < b->n_fmt; ++i) {
        const int off = i_rec * bcf_hdr_nsamples(h);
        SEXP g;
        int t;

        for (t = 0; t < Rf_length(nm); ++t)
            if (bcf_hdr_id2int(h, BCF_DT_ID, CHAR(STRING_ELT(nm, t))) == b->d.fmt[i].id)
                break;
        if (Rf_length(nm) <= t) {
            Rf_error("fmt field '%s' not among those allowed by scanBcf",
                     bcf_hdr_int2id(h, BCF_DT_ID, b->d.fmt[i].id));
        }
        g = VECTOR_ELT(geno, t);

        if (b->d.fmt[i].id == bcf_hdr_id2int(h, BCF_DT_ID, "PL")) {
            const int x = b->n_allele * (b->n_allele + 1) / 2;
            SEXP pl = Rf_allocMatrix(INTSXP, x, bcf_hdr_nsamples(h));
            SET_VECTOR_ELT(g, off, pl);	/* protect */
            if((nvals = bcf_get_format_int32(h, b, "PL", &intdest, &fieldmemsz)) < 0)
               Rf_error("internal: problem reading PL fmt field");
            for (int smplnum = 0; smplnum < bcf_hdr_nsamples(h); ++smplnum) {
                for (int k = 0; k < x; ++k)
                    INTEGER(pl)[smplnum * x + k] = intdest[smplnum * x + k];
            }
        } else if (b->d.fmt[i].id == bcf_hdr_id2int(h, BCF_DT_ID, "DP")) {
            int *dp = INTEGER(g) + off;
            if((nvals = bcf_get_format_int32(h, b, "DP", &intdest, &fieldmemsz)) < 0)
                Rf_error("internal: problem reading DP fmt field");
            for (int smplnum = 0; smplnum < bcf_hdr_nsamples(h); ++smplnum)
                *dp++ = intdest[smplnum];
        } else if (b->d.fmt[i].id == bcf_hdr_id2int(h, BCF_DT_ID, "GQ") ||
                   b->d.fmt[i].id == bcf_hdr_id2int(h, BCF_DT_ID, "SP")) {
            int *gq = INTEGER(g) + off;
            if((nvals = bcf_get_format_int32(h, b, "DP", &intdest, &fieldmemsz)) < 0)
                Rf_error("internal: problem reading GQ or SP fmt field");
            for (int smplnum = 0; smplnum < bcf_hdr_nsamples(h); ++smplnum)
                *gq++ = intdest[smplnum];
        } else if (b->d.fmt[i].id == bcf_hdr_id2int(h, BCF_DT_ID, "GT")) {
            if((nvals = bcf_get_genotypes(h, b, &strdest, &fieldmemsz)) < 0)
                Rf_error("internal: problem reading GT field");
            int idx = off;
            for(int smplnum = 0; smplnum < bcf_hdr_nsamples(h); ++smplnum) {
                SET_STRING_ELT(g, idx++, mkChar(strdest[smplnum]));
            }
        } else if (b->d.fmt[i].id == bcf_hdr_id2int(h, BCF_DT_ID, "GL")) {
            const int x = b->n_allele * (b->n_allele + 1) / 2;
            SEXP gl = Rf_allocMatrix(REALSXP, x, bcf_hdr_nsamples(h));
            SET_VECTOR_ELT(g, off, gl);	/* protect */
            if((nvals = bcf_get_format_int32(h, b, "GL", &floatdest, &fieldmemsz)) < 0)
                Rf_error("intenal: problem reading GL field");
            for (int smplnum = 0; smplnum < bcf_hdr_nsamples(h); ++smplnum) {
                for (int k = 0; k < x; ++k)
                    REAL(gl)[smplnum * x + k] = floatdest[smplnum * x + k];
            }
        }
    }
    free(strdest);
    free(intdest);
    free(floatdest);
}

SEXP scan_bcf_header(SEXP ext)
{
    _checkext(ext, BCFFILE_TAG, "scanBcfHeader");
    vcfFile *bcf = BCFFILE(ext)->file;
    if (!bcf->format.format != vcf && 0 != bgzf_seek(bcf->fp.bgzf, 0, SEEK_SET))
        Rf_error("internal: failed to 'seek'");
    bcf_hdr_t *hdr = bcf_hdr_read(bcf);
    if (NULL == hdr)
        Rf_error("no 'header' line \"#CHROM POS ID...\"?");

    int nseqs = 0;
    const char **seqnames = bcf_hdr_seqnames(hdr, &nseqs);
    if(nseqs < 1)
        Rf_error("header contains no seqname info");

    SEXP ans = PROTECT(NEW_LIST(BCF_HDR_LAST));
    SET_VECTOR_ELT(ans, BCF_HDR_REF, NEW_STRING(nseqs));
    SET_VECTOR_ELT(ans, BCF_HDR_SAMPLE, NEW_STRING(bcf_hdr_nsamples(hdr)));
    SET_VECTOR_ELT(ans, BCF_HDR_HEADER, NEW_STRING(hdr->nhrec));

    int i;
    SEXP x = VECTOR_ELT(ans, BCF_HDR_REF);
    for (i = 0; i < nseqs; ++i)
        SET_STRING_ELT(x, i, mkChar(seqnames[i]));
    x = VECTOR_ELT(ans, BCF_HDR_SAMPLE);
    for (i = 0; i < bcf_hdr_nsamples(hdr); ++i)
        SET_STRING_ELT(x, i, mkChar(hdr->samples[i]));
    x = VECTOR_ELT(ans, BCF_HDR_HEADER);

    int hdrlen = 0;
    char *hdr_text = bcf_hdr_fmt_text(hdr, 1, &hdrlen);
    const char *s;
    if (hdrlen > 0) {
        char *txt = (char *) R_alloc(hdrlen+1, sizeof(char));
        strncpy(txt, hdr_text, hdrlen+1);
        s = strtok(txt, "\n");
        for (i = 0; i < hdrlen; ++i) {
            SET_STRING_ELT(x, i, mkChar(s));
            s = strtok(NULL, "\n");
        }
    }

    SEXP nm = NEW_CHARACTER(3);
    SET_NAMES(ans, nm);         /* protect */
    for (i = 0; i < BCF_HDR_LAST; ++i)
        SET_STRING_ELT(nm, i, mkChar(BCF_HDR_NM[i]));

    free(hdr_text);
    bcf_hdr_destroy(hdr);
    UNPROTECT(1);
    return ans;
}

int scan_bcf_range(vcfFile * bcf, bcf_hdr_t * hdr, SEXP ans, int tid, int start,
                   int end, int n)
{
    const int TID_BUFSZ = 8;
    static char *buf = NULL;
    bcf1_t *bcf1 = bcf_init1();
    if (NULL == bcf1)
        Rf_error("scan_bcf_region: failed to allocate memory");
    int sz = Rf_length(VECTOR_ELT(ans, BCF_TID));
    int res;
    if (NULL == buf)
        buf = Calloc(TID_BUFSZ, char);	/* leaks, but oh well */
    int nseqs = 0;
    const char **seqnames = bcf_hdr_seqnames(hdr, &nseqs);
    if(nseqs < 1)
        Rf_error("header contains no seqname info");
    while (0 <= (res = bcf_read(bcf, hdr, bcf1))) {
        bcf_unpack(bcf1, BCF_UN_ALL);
        if (tid >= 0) {
            /* Finding relative position based on length of ref allele */
            int pos = bcf1->rlen;
            pos = bcf1->pos + (pos > 0 ? pos : 1);
            if (bcf1->rid != tid || bcf1->pos > end)
                break;
            if (!(pos >= start && end > bcf1->pos))
                continue;
        }
        if (n >= sz)
            sz = _bcf_ans_grow(ans, BCF_BUFSIZE_GROW, bcf_hdr_nsamples(hdr));
        if (n >= sz) {
            bcf_destroy(bcf1);
            Rf_error("bcf_scan: failed to increase size; out of memory?");
        }
        if (nseqs > 0)
            SET_STRING_ELT(VECTOR_ELT(ans, BCF_TID), n,
                           smkChar(seqnames[bcf1->rid]));
        else {
            snprintf(buf, TID_BUFSZ, "%d", bcf1->rid);
            SET_STRING_ELT(VECTOR_ELT(ans, BCF_TID), n, smkChar(buf));
        }

        INTEGER(VECTOR_ELT(ans, BCF_POS))[n] = bcf1->pos + 1;
        REAL(VECTOR_ELT(ans, BCF_QUAL))[n] = bcf1->qual;
        SET_STRING_ELT(VECTOR_ELT(ans, BCF_ID), n, smkChar(bcf1->d.id));
        SET_STRING_ELT(VECTOR_ELT(ans, BCF_REF), n, smkChar(bcf1->d.allele[0]));
        /* FIX ME: get all alleles, not just the first alt */
        SET_STRING_ELT(VECTOR_ELT(ans, BCF_ALT), n, smkChar(bcf1->d.allele[1]));
        /* FIX ME: placeholder to get all flt values */
        SET_STRING_ELT(VECTOR_ELT(ans, BCF_FLT), n, smkChar(bcf_hdr_int2id(hdr, BCF_DT_ID, bcf1->d.flt[0])));
        /* FIX ME: placeholder to get all info values */
        SET_STRING_ELT(VECTOR_ELT(ans, BCF_INFO), n, smkChar(bcf_hdr_int2id(hdr, BCF_DT_ID, bcf1->d.info[0].key)));
        SET_STRING_ELT(VECTOR_ELT(ans, BCF_FMT), n, smkChar(bcf_hdr_int2id(hdr, BCF_DT_ID, bcf1->d.fmt[0].id)));
        _bcf_gi2sxp(VECTOR_ELT(ans, BCF_GENO), n, hdr, bcf1);
        ++n;
    }
    free(seqnames);
    bcf_destroy(bcf1);
    return n;
}

SEXP scan_bcf(SEXP ext, SEXP space, SEXP tmpl)
{
    _checkparams(space, R_NilValue, R_NilValue);
    _checkext(ext, BCFFILE_TAG, "scanBcf");
    vcfFile *bcf = BCFFILE(ext)->file;
    hts_idx_t *idx = BCFFILE(ext)->index;
    if (bcf->format.format != vcf && 0 != bgzf_seek(bcf->fp.bgzf, 0, SEEK_SET))
        Rf_error("internal: failed to 'seek' on bcf file");
    bcf_hdr_t *hdr = vcf_hdr_read(bcf);
    if (NULL == hdr)
        Rf_error("no 'header' line \"#CHROM POS ID...\"?");

    int n = 0;
    tmpl = PROTECT(Rf_duplicate(tmpl));

    if (R_NilValue == space) {
        SET_VECTOR_ELT(tmpl, BCF_RECS_PER_RANGE, NEW_INTEGER(1));
        n = scan_bcf_range(bcf, hdr, tmpl, -1, -1, -1, n);
        INTEGER(VECTOR_ELT(tmpl, BCF_RECS_PER_RANGE))[0] = n;
    } else {
        SEXP spc = VECTOR_ELT(space, 0);
        const int
        *start = INTEGER(VECTOR_ELT(space, 1)),
            *end = INTEGER(VECTOR_ELT(space, 2)), nspc = Rf_length(spc);
        SEXP nrec = NEW_INTEGER(nspc);
        SET_VECTOR_ELT(tmpl, BCF_RECS_PER_RANGE, nrec);

        for (int i = 0; i < nspc; ++i) {
            int tid = bcf_hdr_name2id(hdr, CHAR(STRING_ELT(spc, i)));
            if (tid < 0) {
                Rf_error("'space' not in file: %s", CHAR(STRING_ELT(spc, i)));
            }
            hts_itr_t *iter = bcf_itr_queryi(idx, tid, start[i], end[i]);
            uint64_t off = (iter != NULL ? iter->curr_off : 0);
            if (off == 0) {
                INTEGER(nrec)[i] = 0;
                continue;
            }
            bgzf_seek(bcf->fp.bgzf, off, SEEK_SET);
            n = scan_bcf_range(bcf, hdr, tmpl, tid, start[i], end[i], n);
            if (i == 0)
                INTEGER(nrec)[i] = n;
            else
                INTEGER(nrec)[i] = n - INTEGER(nrec)[i - 1];
        }
    }
    _bcf_ans_grow(tmpl, -1 * n, bcf_hdr_nsamples(hdr));

    UNPROTECT(1);
    return tmpl;
}

int _as_bcf(vcfFile * fin, vcfFile * fout)
{
    bcf1_t *b = bcf_init1(); /* free'd in bcf_destroy */
    if (NULL == b)
        Rf_error("_as_bcf: failed to allocate memory");
    bcf_hdr_t *hdr = bcf_hdr_read(fin);
    int r, count = 0;

    bcf_hdr_write(fout, hdr);
    while (0 <= (r = bcf_read(fin, hdr, b))) {
        /* FIX ME: what was this trying to protect against? */
        /* if (NULL == b->ref) */
        /*     Rf_error("cannot (yet) coerce VCF files without FORMAT"); */
        bcf_write(fout, hdr, b);
        count++;
    }

    bcf_hdr_destroy(hdr);
    bcf_destroy(b);

    return r >= -1 ? count : -1 * count;
}

SEXP as_bcf(SEXP file, SEXP dictionary, SEXP destination)
{
    if (!IS_CHARACTER(file) || 1 != LENGTH(file))
        Rf_error("'file' must be character(1)");
    if (!IS_CHARACTER(dictionary) || 1 != LENGTH(dictionary))
        Rf_error("'dictionary' must be character(1)");
    if (!IS_CHARACTER(destination) || 1 != LENGTH(destination))
        Rf_error("'destination' must be character(1)");

    vcfFile *fin = _bcf_tryopen(translateChar(STRING_ELT(file, 0)), "r");
    if (NULL == fin)
        Rf_error("failed to open VCF 'file'");

    vcfFile *fout = _bcf_tryopen(translateChar(STRING_ELT(destination, 0)), "wb");
    if (NULL == fout)
        Rf_error("failed to open BCF 'destination'");

    int count = _as_bcf(fin, fout);

    _bcf_close(fin, 0);
    _bcf_close(fout, 0);
    if (count < 0)
        Rf_error("truncated input file at record %d", -1 * count);

    return destination;
}

SEXP index_bcf(SEXP file)
{
    if (!IS_CHARACTER(file) || 1 != LENGTH(file))
        Rf_error("'file' must be character(1)");
    const char *fbcf = translateChar(STRING_ELT(file, 0));
    int status = bcf_index_build(fbcf, 0);
    if (0 != status)
        Rf_error("failed to build index");
    char *fidx = (char *) R_alloc(strlen(fbcf) + 5, sizeof(char));
    sprintf(fidx, "%s.bci", fbcf);
    return mkString(fidx);
}
