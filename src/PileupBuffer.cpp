#include "PileupBuffer.h"

int Pileup::insert(uint32_t, uint32_t pos, int n,
                        const bam_pileup1_t *pl, void *data)
{
    //Rprintf("pos: %d\n", pos);
    Pileup *pileup = (Pileup *) data;
    pos = pos + 1;
    int bamBufOffset = 0;
    if(pos >= pileup->start && pos <= pileup->end) {
        pileup->resultMgr->signalPosStart(pos);
        for(bamBufOffset = 0; bamBufOffset != n; ++bamBufOffset) {
            const bam_pileup1_t *curBam = pl + bamBufOffset;
            // Rprintf("distance from end: %u\n", curBam->b->core.l_qseq - curBam->qpos);

            // whole-alignment disqualifiers/filters
            if(curBam->is_refskip) // never pileup a refskip ('N' op in cigar)
                continue;
            const uint8_t mapqual = curBam->b->core.qual;
            if(mapqual < pileup->min_mapq()) continue;
            
            // individual nucleotide disqualifiers
            char nucleotide = 'X', strand = 'X';
            int bin = 0;
            const uint8_t basequal = bam1_qual(curBam->b)[curBam->qpos];
            if(basequal < pileup->min_baseq()) continue;
            bool isDeletion = curBam->is_del;
            if(isDeletion && !pileup->include_deletions())
                continue;
            if(isDeletion && pileup->include_deletions())
                nucleotide = '-';
            else {
                nucleotide =
                    char(bam_nt16_rev_table[bam1_seqi(bam1_seq(curBam->b),
                                                      curBam->qpos)]);
            }

            bool dropNucleotide = (nucleotide == 'N' && pileup->ignoreNs());
            if(dropNucleotide) continue;

            // invariant: alignments that fail strand criterion not included
            if(pileup->hasStrands())
                strand = bam1_strand(curBam->b) ? '-' : '+';

            if(pileup->hasBins()) { // all bin work
                int minBinPoint = pileup->minBinPoint(),
                    maxBinPoint = pileup->maxBinPoint();
                int qpos = curBam->qpos + 1;
                // discard if outside outer range
                if(qpos > minBinPoint && qpos <= maxBinPoint) {
                    bin = pileup->calcBin(qpos);
                    //Rprintf("bin number %d\n", bin);
                } else
                    continue;
            }

            pileup->resultMgr->forwardTuple(BamTuple(nucleotide, strand, bin));
        }
        // extract tuples for pos
        pileup->resultMgr->signalPosEnd();
    }
    return 0;
}

void extract(const ResultMgrInterface * const from, SEXP to, bool hasStrands,
             bool hasNucleotides, bool hasBins) {
    #ifdef PILEUP_DEBUG
    assert(IS_LIST(to));
    for(int i = 0; i != Rf_length(to); ++i) {
        SEXP elt = VECTOR_ELT(to, i);
        assert(IS_LIST(elt));
        assert((unsigned int)Rf_length(elt) == from->size());
    }
    #endif
    int curDim = 1; // 0 is seqnames; must start with pos
    std::copy(from->posBeg(), from->posEnd(),
              INTEGER(VECTOR_ELT(to, curDim++)));

    // Might still be worthwhile to make original vectors ints
    SEXP strand = R_NilValue, nucleotide = R_NilValue;
    if(hasStrands) {
        strand = VECTOR_ELT(to, curDim++);
        std::transform(from->strandBeg(), from->strandEnd(),
                       INTEGER(strand), Pileup::strand_to_lvl);
    }
    if(hasNucleotides) {
        nucleotide = VECTOR_ELT(to, curDim++);
        std::transform(from->nucBeg(), from->nucEnd(),
                       INTEGER(nucleotide), Pileup::nuc_to_lvl);
    }
    if(hasBins) {
        std::copy(from->binBeg(), from->binEnd(),
                  INTEGER(VECTOR_ELT(to, curDim++)));
    }

    std::copy(from->countBeg(), from->countEnd(),
              INTEGER(VECTOR_ELT(to, curDim++)));
    if(hasStrands)
        _as_strand(strand);
    if(hasNucleotides)
        _as_nucleotide(nucleotide);
}

SEXP Pileup::yield() {
    int numDims = 3;
    numDims += hasStrands() ? 1 : 0;
    numDims += hasNucleotides() ? 1 : 0;
    numDims += hasBins() ? 1 : 0;
    uint32_t numResults = resultMgr->size();
    SEXP result = PROTECT(Rf_allocVector(VECSXP, numDims));
    int curDim = 0;
    SET_VECTOR_ELT(result, curDim, Rf_allocVector(INTSXP, numResults));//seqns
    SEXP seqnames = VECTOR_ELT(result, curDim++);
    std::fill_n(INTEGER(seqnames), numResults, getSeqlevelValue(rname));
    SET_VECTOR_ELT(result, curDim++, Rf_allocVector(INTSXP, numResults)); // pos
    if(hasStrands())
        SET_VECTOR_ELT(result, curDim++, Rf_allocVector(INTSXP, numResults));
    if(hasNucleotides())
        SET_VECTOR_ELT(result, curDim++, Rf_allocVector(INTSXP, numResults));
    if(hasBins())
        SET_VECTOR_ELT(result, curDim++, Rf_allocVector(INTSXP, numResults));
    SET_VECTOR_ELT(result, curDim++, Rf_allocVector(INTSXP, numResults));//count

    SEXP nms = PROTECT(Rf_allocVector(STRSXP, numDims));
    curDim = 0;
    SET_STRING_ELT(nms, curDim++, mkChar("seqnames"));
    SET_STRING_ELT(nms, curDim++, mkChar("pos"));
    if(hasStrands())
        SET_STRING_ELT(nms, curDim++, mkChar("strand"));
    if(hasNucleotides())
        SET_STRING_ELT(nms, curDim++, mkChar("nucleotide"));
    if(hasBins())
        SET_STRING_ELT(nms, curDim++, mkChar("cycle_bin"));
    SET_STRING_ELT(nms, curDim++, mkChar("count"));
    SET_ATTR(result, R_NamesSymbol, nms);

    extract(resultMgr, result, hasStrands(), hasNucleotides(), hasBins());
    resultMgr->signalYieldEnd();
    _as_seqlevels(VECTOR_ELT(result, 0), seqnamesLevels);

    UNPROTECT(2);
    return result;
}