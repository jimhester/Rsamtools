#include "sam.h"
#include "bamfile.h"

int _as_bam(samFile * fin, samFile * fout)
{
    bam1_t *b = bam_init1();
    int r, count = 0;
    bam_hdr_t *hdr = sam_hdr_read(fin);

    while (0 <= (r = sam_read1(fin, hdr, b))) {
        /* less general than sam_write1--maybe bam_write1 faster */
        bam_write1(fout->fp.bgzf, b);
        count++;
    }
    bam_hdr_destroy(hdr);
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
    bam_hdr_t *hdr;
    if (LOGICAL(binary)[0]) {
        /* SAM --> BAM */
        fin = _bam_tryopen(translateChar(STRING_ELT(file, 0)), "r");
        if(fin == NULL) {
            Rf_error("failed to open samfile '%s'",
                     translateChar(STRING_ELT(file, 0)));
        }
        /* header no longer bundled with file (no `header` member) */
        hdr = sam_hdr_read(fin);
        /* check that header opened correctly? */

        fout = _bam_tryopen(translateChar(STRING_ELT(destination, 0)), "wb");
    } else {
        /* BAM --> SAM */
        fin = _bam_tryopen(translateChar(STRING_ELT(file, 0)), "rb");
        if(fin == NULL) {
            Rf_error("failed to open bamfile '%s'",
                     translateChar(STRING_ELT(file, 0)));
        }
        /* FIX ME: Always start by writing header ? */
        hdr = sam_hdr_read(fin);

        fout = _bam_tryopen(translateChar(STRING_ELT(destination, 0)), "wh");
    }
    sam_hdr_write(fout, hdr);
    int count = _as_bam(fin, fout);

    sam_close(fin);
    sam_close(fout);
    if (count < 0)
        Rf_error("truncated input file at record %d", -1 * count);

    return destination;
}
