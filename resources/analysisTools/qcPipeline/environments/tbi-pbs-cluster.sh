#!/usr/bin/env bash

# The sambamba version used for sorting, viewing. Note that v0.5.9 is segfaulting on convey during view or sort.
export SAMBAMBA_BINARY="sambamba_v${SAMBAMBA_VERSION:?SAMBAMBA_VERSION is not set}"

# The sambamba version used only for making flagstats.
export SAMBAMBA_FLAGSTATS_BINARY="sambamba_v${SAMBAMBA_VERSION:?SAMBAMBA_VERSION is not set}"
# Warning: Currently bwaMemSortSlim uses sambamba flagstats, while mergeAndMarkOrRemoveSlim uses samtools flagstats.

# The sambamba version used only for duplication marking and merging.
export SAMBAMBA_MARKDUP_BINARY="sambamba_v${SAMBAMBA_MARKDUP_VERSION:?SAMBAMBA_MARKDUP_VERSION is not set}"

export PICARD_BINARY="picard-${PICARD_VERSION:?PICARD_VERSION is not set}.sh"
export FASTQC_BINARY="fastqc-${FASTQC_VERSION:?FASTQC_VERSION is not set}"
export SAMTOOLS_BINARY="samtools-${SAMTOOLS_VERSION:?SAMTOOLS_VERSION is not set}"

# biobambam
export BAMSORT_BINARY="bamsort-${BIOBAMBAM_VERSION:?BIOBAMBAM_VERSION is not set}"
export BAMMARKDUPLICATES_BINARY="bammarkduplicates-${BIOBAMBAM_VERSION:?BIOBAMBAM_VERSION is not set}"

# htslib/0.2.5
export BGZIP_BINARY="bgzip-tabix-${HTSLIB_VERSION:?HTSLIB_VERSION is not set}"
export TABIX_BINARY="tabix-${HTSLIB_VERSION}:?HTSLIB_VERSION is not set}"

# WES only
export INTERSECTBED_BINARY="intersectBed-${BEDTOOLS_VERSION:?BEDTOOLS_VERSION is not set}"
export COVERAGEBED_BINARY="coverageBed-${BEDTOOLS_VERSION:?BEDTOOLS_VERSION is not set}"

export PERL_BINARY="perl-${PERL_VERSION:?PERL_VERSION is not set}"
export PYTHON_BINARY="python-${PYTHON_VERSION:?PYTHON_VERSION is not set}"
export RSCRIPT_BINARY="Rscript-${R_VERSION:?R_VERSION is not set}"

export JAVA_BINARY="java8"

if [[ "$WORKFLOW_ID" == "bisulfiteCoreAnalysis" ]]; then
    export BWA_BINARY="bwa-${BWA_VERSION:?BWA_VERSION is not set}-bisulfite"
else
    export BWA_BINARY="bwa-${BWA_VERSION:?BWA_VERSION is not set}"
    export BWA_ACCELERATED_BINARY="/opt/bb/bwa-${BWA_VERSION:?BWA_VERSION is not set}/bin/bwa-bb"
fi

export MBUFFER_BINARY=mbuffer

