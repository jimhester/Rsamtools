#include "sam.h"
#include "bamfile.h"

int _as_bam(samFile * fin, bam_hdr_t *inhdr, samFile * fout)
{
    /* header already written to fout by `as_bam`, so just records */
    bam1_t *b = bam_init1();
    int r, count = 0;

    while (0 <= (r = sam_read1(fin, inhdr, b))) {
        sam_write1(fout, inhdr, b);
        count++;
    }
    bam_destroy1(b);

    return r >= -1 ? count : -1 * count;
}

SEXP as_bam(SEXP file, SEXP destination, SEXP binary)
{
    if (!IS_CHARACTER(file) || 1 != LENGTH(file))
        Rf_error("'file' must be character(1)");
    if (!IS_CHARACTER(destination) || 1 != LENGTH(destination))
        Rf_error("'destination' must be character(1)");
    if (!IS_LOGICAL(binary) || 1 != LENGTH(binary))
        Rf_error("'binary' must be logical(1)");

    samFile *fin, *fout;
    const char *finname = translateChar(STRING_ELT(file, 0)),
        *foutname = translateChar(STRING_ELT(destination, 0));
    bam_hdr_t *hdr;
    if (LOGICAL(binary)[0]) {
        /* SAM --> BAM */
        fin = _bam_tryopen(finname, "r");
        if(fin == NULL) {
            Rf_error("failed to open SAM file '%s'", finname);
        }
        hdr = sam_hdr_read(fin);
        if(hdr == NULL)
            Rf_error("failed to read header of SAM file '%s'", finname);

        fout = _bam_tryopen(foutname, "wb");
    } else {
        /* BAM --> SAM */
        fin = _bam_tryopen(finname, "rb");
        if(fin == NULL) {
            Rf_error("failed to open BAM file '%s'", finname);
        }
        hdr = sam_hdr_read(fin);
        if(hdr == NULL)
            Rf_error("failed to read header of BAM file '%s'", finname);

        fout = _bam_tryopen(foutname, "wh");
    }
    sam_hdr_write(fout, hdr);
    int count = _as_bam(fin, hdr, fout);

    bam_hdr_destroy(hdr);
    sam_close(fin);
    sam_close(fout);
    if (count < 0)
        Rf_error("truncated input file at record %d", -1 * count);

    return destination;
}
