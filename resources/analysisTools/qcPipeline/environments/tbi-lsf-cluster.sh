#!/usr/bin/env bash

set -xv

# Given a name of a tool version variable, get the name of the tool version variable. The
# idea is that if the tool name variable is not defined, the name is interpreted as plain tool
# name and an upper-cased _VERSION is appended.
#
# sambamba => SAMBAMBA_VERSION
# SAMBAMBA_VERSION => SAMBAMBA_VERSION    // if SAMBAMBA_VERSION is defined as variable
#
getVersionVariableName () {
    local name="${1:?No binary name given}"
    local versionVariable="${2:-}"
    if [[ -z "$versionVariable" ]]; then
        local ucVersionVariable
        ucVersionVariable=$(echo "$name" | tr '[a-z]' '[A-Z]')
        echo "${ucVersionVariable}_VERSION"
    else
        echo "$versionVariable"
    fi
}
export -f getVersionVariableName

# Load a module given the module name. The module version is implicitly taken from the
# matching tool version variable, optionally given as second parameter. If the second
# parameter is not given, the tool name will be upper-cased and _VERSION appended, to
# get the name of the tool version variable.
#
# R_VERSION=3.3.1, then $(moduleLoad R) will do $(module load R/3.3.1).
# RSCRIPT_VERSION=3.4.0, then $(moduleLoad R RSCRIPT_VERSION) will do $(module load R/3.4.0)
moduleLoad() {
    local name="${1:?No module name given}"
    local versionVariable="${2:-}"
    versionVariable=$(getVersionVariableName "$name" "$versionVariable")
    if [[ -z "${!versionVariable}" ]]; then
        throw 200 "$versionVariable is not set"
    fi
    module load "${name}/${!versionVariable}" || throw 200 "Could not load '${name}/${!versionVariable}'"
}
export -f moduleLoad

# Like moduleLoad, but for unloading the module/version.
moduleUnload() {
    local name="${1:?No module name given}"
    local versionVariable="${2:-}"
    local versionVariable=$(getVersionVariableName "$name" "$versionVariable")
    if [[ -z "${!versionVariable}" ]]; then
        throw 200 "versionVariable for $name is not set"
    fi
    module unload "${name}/${!versionVariable}" || throw 200 "Could not unload '${name}/${!versionVariable}'"
}
export -f moduleUnload

moduleLoad htslib
export BGZIP_BINARY=bgzip
export TABIX_BINARY=tabix

moduleLoad R
export RSCRIPT_BINARY=Rscript

moduleLoad java
export JAVA_BINARY=java

moduleLoad fastqc
export FASTQC_BINARY=fastqc

moduleLoad perl
export PERL_BINARY=perl

moduleLoad python
export PYTHON_BINARY=python

moduleLoad samtools
export SAMTOOLS_BINARY=samtools
export BCFTOOLS_BINARY=bcftools

moduleLoad bedtools
export INTERSECTBED_BINARY=intersectBed
export COVERAGEBED_BINARY=coverageBed
export FASTAFROMBED_BINARY=fastaFromBed

moduleLoad libmaus
moduleLoad biobambam
export BAMSORT_BINARY=bamsort
export BAMMARKDUPLICATES_BINARY=bammarkduplicates

moduleLoad picard
export PICARD_BINARY=picard.sh

moduleLoad vcftools
export VCFTOOLS_SORT_BINARY=vcf-sort

# There are different sambamba versions used for different tasks. The reason has something to do with performance and differences in the versions.
# We define functions here that take the same parameters as the original sambamba, but apply them to the appropriate version by first loading
# the correct version and after the call unloading the version. The _BINARY variables are set to the functions.

# The sambamba version used for sorting, viewing. Note that v0.5.9 is segfaulting on Convey during view or sort.
sambamba_sort_view() {
    moduleLoad sambamba SAMBAMBA_VERSION
    sambamba "$@"
    moduleUnload sambamba
}
export -f sambamba_sort_view
export SAMBAMBA_BINARY=sambamba_sort_view

# The sambamba version used only for flagstats. For the flagstats sambamba 0.4.6 used is equivalent to samtools 0.1.19 flagstats. Newer versions
# use the new way of counting in samtools (accounting for supplementary reads).
# Warning: Currently bwaMemSortSlim uses sambamba flagstats, while mergeAndMarkOrRemoveSlim uses samtools flagstats.
sambamba_flagstat() {
    moduleLoad sambamba SAMBAMBA_FLAGSTATS_VERSION
    sambamba "$@"
    moduleUnload sambamba SAMBAMBA_FLAGSTATS_VERSION
}
export -f sambamba_flagstat
export SAMBAMBA_FLAGSTATS_BINARY=sambamba_flagstat

# The sambamba version used only for duplication marking and merging. Use the bash function here!
# Should be changeable independently also for performance and stability reasons.
sambamba_markdup() {
    moduleLoad sambamba SAMBAMBA_MARKDUP_VERSION
    sambamba "$@"
    moduleUnload sambamba SAMBAMBA_MARKDUP_VERSION
}
export -f sambamba_markdup
export SAMBAMBA_MARKDUP_BINARY=sambamba_markdup


if [[ "$WORKFLOW_ID" == "bisulfiteCoreAnalysis" ]]; then
    ## For bisulfite alignment, we suffix the the value of BINARY_VERSION by '-bisulfite', because that's the name in LSF cluster.
    export BWA_VERSION="${BWA_VERSION:?BWA_VERSION is not set}-bisulfite"
    moduleLoad bwa
    export BWA_BINARY=bwa

elif [[ "$WORKFLOW_ID" == "qcAnalysis" || "$WORKFLOW_ID" == "exomeAnalysis" ]]; then
    if [[ "${useAcceleratedHardware:-false}" == false ]]; then
        moduleLoad bwa
        export BWA_BINARY=bwa
    elif [[ "${useAcceleratedHardware:-true}" == true ]]; then
        moduleLoad bwa-bb BWA_VERSION
        export BWA_ACCELERATED_BINARY=bwa-bb
    else
        throw 200 "Uninterpretable value for boolean 'useAcceleratedHardware': '$useAcceleratedHardware'"
    fi
else
    throw 200 "Unknown workflow ID '$WORKFLOW_ID'"
fi

# Unversioned binaries.
export MBUFFER_BINARY=mbuffer
export CHECKSUM_BINARY=md5sum
