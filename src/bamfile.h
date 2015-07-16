#ifndef BAMFILE_H
#define BAMFILE_H

#include <Rdefines.h>
#include "sam.h"
#include "bgzf.h"
#include "bambuffer.h"
#include "bam_mate_iter.h"
#include "pbuffer_wrapper.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    samFile *file;
    bam_hdr_t *header;
    hts_idx_t *index;
    uint64_t pos0;
    int irange0;
    bam_mate_iter_t iter;
    void *pbuffer; /* for buffered pileup */
} _BAM_FILE, *BAM_FILE;

#define BAMFILE(b) ((BAM_FILE) R_ExternalPtrAddr(b))

SEXP bamfile_init();
SEXP bamfile_open(SEXP file0, SEXP file1, SEXP mode);
SEXP bamfile_close(SEXP ext);
SEXP bamfile_isopen(SEXP ext);
SEXP bamfile_isincomplete(SEXP ext);

SEXP read_bamfile_header(SEXP ext, SEXP what);
SEXP scan_bamfile(SEXP ext, SEXP space, SEXP keepFlags,
                  SEXP simpleCigar, SEXP tagFilter, 
                  SEXP reverseComplement, SEXP yieldSize,
                  SEXP tmpl, SEXP obeyQname, 
                  SEXP asMates, SEXP qnamePrefix, SEXP qnameSuffix);
SEXP count_bamfile(SEXP ext, SEXP space, SEXP keepFlags, SEXP isSimpleCigar,
                   SEXP tagFilter);
SEXP prefilter_bamfile(SEXP ext, SEXP space, SEXP keepFlags,
		       SEXP isSimpleCigar, SEXP tagFilter, SEXP yieldSize,
                       SEXP obeyQname, SEXP asMates, SEXP qnamePrefix,
                       SEXP qnameSuffix);
SEXP filter_bamfile(SEXP ext, SEXP space, SEXP keepFlags,
                    SEXP isSimpleCigar, SEXP tagFilter, SEXP fout_name,
                    SEXP fout_mode);

void _check_isbamfile(SEXP ext, const char *lbl);
samFile *_bam_tryopen(const char *filename, const char *mode);

#ifdef __cplusplus
}
#endif

#endif
