#!/bin/bash

source "$TOOL_WORKFLOW_LIB"

printInfo 

set -o pipefail 

ID=${RUN}_${LANE}
SM=sample_${SAMPLE}_${PID}

# use scratch dir for temp files: samtools sort uses the current working directory for them unless provided with a different prefix
# so far, it's not the local scratch (${RODDY_SCRATCH}) but inside roddyExecutionStore of the PID
# here, using the job ID is a bad idea because if the job fails, the temp files remain
#WORKDIR=${DIR_TEMP}/${RODDY_JOBID}
WORKDIR=${FILENAME_SORTED_BAM}_TEMP_SORTING
mkdir -p ${WORKDIR}

# pipes via local scratch dir
FNPIPE1=${RODDY_SCRATCH}/NAMED_PIPE1
FNPIPE2=${RODDY_SCRATCH}/NAMED_PIPE2
NP_READBINS_IN=${RODDY_SCRATCH}/np_readbins_in
NP_COVERAGEQC_IN=${RODDY_SCRATCH}/np_coverageqc_in
NP_COMBINEDANALYSIS_IN=${RODDY_SCRATCH}/np_combinedanalysis_in
NP_FLAGSTATS=${RODDY_SCRATCH}/np_flagstats_in
NP_SAMTOOLS_INDEX_IN=${RODDY_SCRATCH}/np_samtools_index_in
NP_BWA_OUT=${RODDY_SCRATCH}/np_bwa_out
NP_BCONV_OUT=${RODDY_SCRATCH}/np_bconv_out

MBUF_SMALL="${MBUFFER_BINARY} -m ${MBUFFER_SIZE_SMALL} -q -l /dev/null"
MBUF_LARGE="${MBUFFER_BINARY} -m ${MBUFFER_SIZE_LARGE} -q -l /dev/null"

mkfifo ${NP_READBINS_IN} ${NP_COVERAGEQC_IN} ${NP_COMBINEDANALYSIS_IN} ${NP_FLAGSTATS} ${NP_BWA_OUT} ${NP_BCONV_OUT}

bamname=`basename ${FILENAME_SORTED_BAM}`
INDEX_FILE=${FILENAME_SORTED_BAM}.bai
tempSortedBamFile=${FILENAME_SORTED_BAM}.tmp
tempFileForSort=${WORKDIR}/${bamname}_forsorting
tempBamIndexFile=${FILENAME_SORTED_BAM}.tmp.bai
tempFlagstatsFile=${FILENAME_FLAGSTATS}.tmp

# samtools sort may complain about truncated temp files and for each line outputs
# the error message. This happens when the same files are written at the same time,
# see http://sourceforge.net/p/samtools/mailman/samtools-help/thread/BAA90EF6FE3B4D45A7B2F6E0EC5A8366DA3AB5@USTLMLLYC102.rf.lilly.com/
NP_SORT_ERRLOG="$RODDY_SCRATCH/NP_SORT_ERRLOG"
FILENAME_SORT_LOG="$DIR_TEMP/${bamname}_errlog_sort"
FILENAME_BWA_EC="$DIR_TEMP/${bamname}_ec"

RAW_SEQ=${RAW_SEQ_1}
source ${TOOL_COMMON_ALIGNMENT_SETTINGS_SCRIPT}

NP_BAMSORT=${RODDY_SCRATCH}/NAMED_PIPE_BAMSORT
mkfifo ${NP_BAMSORT}

# Create the following variable for error checking issues
TMP_FILE=${tempSortedBamFile}
# error tracking
FILENAME_BWA_LOG=${DIR_TEMP}/${bamname}_errlog_bwamem

bamFileExists=false
# in case the BAM already exists, but QC files are missing, create these only
if [[ -f ${FILENAME_SORTED_BAM} ]] && [[ -s ${FILENAME_SORTED_BAM} ]]
then
    checkBamIsComplete "$FILENAME_SORTED_BAM"
	bamFileExists=true
fi

# test if one of the fastq files is a fake fastq file to simulate paired end sequencing in PIDs with mixed sequencing (single and paired end)
LENGTH_SEQ_1=`${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_1} 2>/dev/null | head | wc -l`
LENGTH_SEQ_2=`${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_2} 2>/dev/null | head | wc -l`
[[ ${LENGTH_SEQ_1} == 0 || ${LENGTH_SEQ_2} == 0 ]] && useSingleEndProcessing=true

# make biobambam sort default
useBioBamBamSort=${useBioBamBamSort-true}

# Default: Dummy process IDs to simplify downstream logic.
true & procTrim=$!
true & procUnpack1=$!
true & procUnpack2=$!
true & procID_fqconv_r1=$!
true & procID_fqconv_r2=$!

if [[ ${bamFileExists} == false ]]; then	# we have to make the BAM 
	mkfifo ${FNPIPE1} ${FNPIPE2}
	if [[ $useAdaptorTrimming == true ]]; then # [[ "$ADAPTOR_TRIMMING_TOOL" == *.jar ]]
		if [ "${qualityScore}" = "illumina" ]; then
			eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_1} | ${PERL_BINARY} ${TOOL_CONVERT_ILLUMINA_SCORES} - > $i1" & procUnpack1=$!
			eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_2} | ${PERL_BINARY} ${TOOL_CONVERT_ILLUMINA_SCORES} - > $i2" & procUnpack2=$!
		else
			eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_1} > $i1" & procUnpack1=$!
			eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_2} > $i2" & procUnpack2=$!
		fi
		# done in the sourced script designed for bwa aln:
		#     i1=$DIR_SCRATCH/at_i1
		#     i2=$DIR_SCRATCH/at_i2
		#     o1=$DIR_SCRATCH/at_o1
		#     o2=/dev/null
		#     u1=/dev/null
		#     u2=/dev/null
		#     mkfifo $i1 $i2 $o1

		# But o2 now has to contain read2 here:
		o2=${RODDY_SCRATCH}/at_o2
		mkfifo $o2
		eval "$JAVA_BINARY -jar ${TOOL_ADAPTOR_TRIMMING} $ADAPTOR_TRIMMING_OPTIONS_0 $i1 $i2 $o1 $u1 $o2 $u2 $ADAPTOR_TRIMMING_OPTIONS_1" & procTrim=$!

		# trimming with fastx does not work in combination with Trimmomatic!
		# besides, bwa mem automagically reverts mate pair data
		cat $o1 | ${PYTHON_BINARY} ${TOOL_METHYL_C_TOOLS} fqconv -1 - - > $FNPIPE1 & procID_fqconv_r1=$!
		cat $o2 | ${PYTHON_BINARY} ${TOOL_METHYL_C_TOOLS} fqconv -2 - - > $FNPIPE2 & procID_fqconv_r2=$!
	elif [ "${qualityScore}" = "illumina" ]; then	# bwa mem has no possibility to convert Illumina 1.3 scores
	    true & procTrim=$!     # dummy process id
		eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_1} | ${PERL_BINARY} ${TOOL_CONVERT_ILLUMINA_SCORES} - | \
		      ${PYTHON_BINARY} ${TOOL_METHYL_C_TOOLS} fqconv -1 - - > $FNPIPE1" & procID_fqconv_r1=$!
		eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_2} | ${PERL_BINARY} ${TOOL_CONVERT_ILLUMINA_SCORES} - | \
		      ${PYTHON_BINARY} ${TOOL_METHYL_C_TOOLS} fqconf -2 - - > $FNPIPE2" & procID_fqconv_r2=$!
	else
	    true & procTrim=$!     # dummy process id
		eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_1} | ${PYTHON_BINARY} ${TOOL_METHYL_C_TOOLS} fqconv -1 - - > $FNPIPE1" & procID_fqconv_r1=$!
		eval "${UNZIPTOOL} ${UNZIPTOOL_OPTIONS} ${RAW_SEQ_2} | ${PYTHON_BINARY} ${TOOL_METHYL_C_TOOLS} fqconv -2 - - > $FNPIPE2" & procID_fqconv_r2=$!
	fi

	INPUT_PIPES=""
	[[ ${LENGTH_SEQ_1} != 0 ]] && INPUT_PIPES="$FNPIPE1"
	[[ ${LENGTH_SEQ_2} != 0 ]] && INPUT_PIPES="${INPUT_PIPES} $FNPIPE2"
	[[ ${LENGTH_SEQ_1} == 0 ]] && cat $FNPIPE1 >/dev/null
	[[ ${LENGTH_SEQ_2} == 0 ]] && cat $FNPIPE2 >/dev/null
fi

# Try to read from pipes BEFORE they are filled.
# In all cases:
# SAM output is piped to perl script that calculates various QC measures
(${PERL_BINARY} ${TOOL_COMBINED_BAM_ANALYSIS} -i ${NP_COMBINEDANALYSIS_IN} -a ${FILENAME_DIFFCHROM_MATRIX}.tmp -c ${CHROM_SIZES_FILE} -b ${FILENAME_ISIZES_MATRIX}.tmp  -f ${FILENAME_EXTENDED_FLAGSTATS}.tmp  -m ${FILENAME_ISIZES_STATISTICS}.tmp -o ${FILENAME_DIFFCHROM_STATISTICS}.tmp -p ${INSERT_SIZE_LIMIT} ) & procIDCBA=$!

# genome coverage (depth of coverage and other QC measures in one file)
(${TOOL_COVERAGE_QC_D_IMPL} --alignmentFile=${NP_COVERAGEQC_IN} --outputFile=${FILENAME_GENOME_COVERAGE}.tmp --processors=1 --basequalCutoff=${BASE_QUALITY_CUTOFF} --ungappedSizes=${CHROM_SIZES_FILE}) & procIDGenomeCoverage=$!

# read bins
(set -o pipefail; ${TOOL_GENOME_COVERAGE_D_IMPL} --alignmentFile=${NP_READBINS_IN} --outputFile=/dev/stdout --processors=4 --mode=countReads --windowSize=${WINDOW_SIZE} | $MBUF_SMALL | ${PERL_BINARY} ${TOOL_FILTER_READ_BINS} - ${CHROM_SIZES_FILE} > ${FILENAME_READBINS_COVERAGE}.tmp) & procIDReadbinsCoverage=$!

# use sambamba for flagstats
(set -o pipefail; cat ${NP_FLAGSTATS} | ${SAMBAMBA_FLAGSTATS_BINARY} flagstat -t 1 /dev/stdin > $tempFlagstatsFile) & procIDFlagstat=$!

if [[ ${bamFileExists} == true ]]; then
	echo "the BAM file already exists, re-creating other output files."
	# make all the pipes
	(cat ${FILENAME_SORTED_BAM} | ${MBUF_LARGE} | tee ${NP_COVERAGEQC_IN} ${NP_READBINS_IN} ${NP_FLAGSTATS} \
	    | ${SAMBAMBA_BINARY} view /dev/stdin | ${MBUF_LARGE} > $NP_COMBINEDANALYSIS_IN) & procIDOutPipe=$!
else
	if [[ ${useBioBamBamSort} == false ]]; then
		# we use samtools for making the index
		# filter secondary and supplementary alignments (flag 2304) after alignment
		mkfifo $NP_SORT_ERRLOG ${NP_SAMTOOLS_INDEX_IN}

		# Index bam file
		${SAMTOOLS_BINARY} index ${NP_SAMTOOLS_INDEX_IN} ${tempBamIndexFile} & procID_IDX=$!

		# Align converted fastq files
		${BWA_BINARY} mem \
			-t ${BWA_MEM_THREADS} \
			-R "@RG\tID:${ID}\tSM:${SM}\tLB:${LB}\tPL:ILLUMINA" \
			$BWA_MEM_OPTIONS ${INDEX_PREFIX} ${INPUT_PIPES} 2> $FILENAME_BWA_LOG \
		> ${NP_BWA_OUT} & procID_BWA=$!

		# Convert aligned reads back to original state
		${SAMTOOLS_BINARY} view -uSbh -F 2304 ${NP_BWA_OUT} | \
		${PYTHON_BINARY} ${TOOL_METHYL_C_TOOLS} bconv - - \
		> ${NP_BCONV_OUT} & procID_BCONV=$!

		# Sort bam file
		(set -o pipefail; \
		${SAMTOOLS_BINARY} view -h ${NP_BCONV_OUT} | \
		tee ${NP_COMBINEDANALYSIS_IN} | \
		${SAMTOOLS_BINARY} view -uSbh -F 2304 - | \
		$MBUF_LARGE | \
		${SAMTOOLS_BINARY} sort -@ 8 \
			-m ${SAMPESORT_MEMSIZE} \
			-o - ${tempFileForSort} 2>$NP_SORT_ERRLOG | \
		tee ${NP_COVERAGEQC_IN} \
			${NP_READBINS_IN} \
			${NP_FLAGSTATS} \
			${NP_SAMTOOLS_INDEX_IN} > ${tempSortedBamFile}; \
		    echo "$? (Pipe Exit Codes: ${PIPESTATUS[@]})" > "$FILENAME_BWA_EC") & procID_MEMSORT=$!

		# filter samtools error log
		(cat $NP_SORT_ERRLOG | uniq > $FILENAME_SORT_LOG) & procID_logwrite=$!

		# Check for errors
		wait $procID_logwrite	# do we need a check for it?
		wait $procID_BWA || throw 20 "Error in BWA alignment"
		wait $procID_BCONV || throw 21 "Error in BCONV"
		wait $procID_MEMSORT || throw 24 "Error in samtools sort"
		wait $procID_IDX; [[ ! $? -eq 0 ]] && echo "Error from samtools index" && exit 10
	else
		throw 11 "biobambam sort not implemented yet";
	fi;
fi

if [[ ${bamFileExists} == true ]]; then
	wait $procIDOutPipe || throw 13 "Error from sambamba view pipe"
else	# make sure to rename BAM file when it has been produced correctly
	[[ -p $i1 ]] && rm $i1 $i2 $o1 $o2 2> /dev/null
	rm -rf ${WORKDIR} 2> /dev/null # Remove temporary files sorting
	rm $FNPIPE1
	rm $FNPIPE2
	errorString="There was a non-zero exit code in the bwa mem - sort pipeline; exiting..."
	source ${TOOL_BWA_ERROR_CHECKING_SCRIPT}
	checkBamIsComplete "$tempSortedBamFile"
	mv ${tempSortedBamFile} ${FILENAME_SORTED_BAM} || throw 36 "Could not move file"
	# index is only created by samtools or biobambam when producing the BAM, it may be older than the BAM, so update time stamp
	if [[ -f ${tempBamIndexFile} ]]; then
	 	mv ${tempBamIndexFile} ${INDEX_FILE} && touch ${INDEX_FILE} || throw 37 "Could not move file"
	fi
	# clean up adapter trimming pipes if they exist
	[[ -p $i1 ]] && rm $i1 $i2 $o1 $o2 2> /dev/null
	
fi

wait $procTrim || throw 38 "Error from trimming"
wait $procUnpack1 || throw 39 "Error from reading FASTQ 1"
wait $procUnpack2 || throw 40 "Error from reading FASTQ 2"
wait $procID_fqconv_r1 || throw 18 "Error in FQCONV (read 1)"
wait $procID_fqconv_r2 || throw 19 "Error in FQCONV (read 2)"

wait $procIDFlagstat || throw 14 "Error from sambamba flagstats"
wait $procIDReadbinsCoverage || throw 15 "Error from genomeCoverage read bins"
wait $procIDGenomeCoverage || throw 16 "Error from coverageQCD"
wait $procIDCBA || throw 17 "Error from combined QC perl script"

# rename QC files
mv ${FILENAME_DIFFCHROM_MATRIX}.tmp ${FILENAME_DIFFCHROM_MATRIX} || throw 28 "Could not move file"
mv ${FILENAME_ISIZES_MATRIX}.tmp ${FILENAME_ISIZES_MATRIX} || throw 29 "Could not move file"
mv ${FILENAME_EXTENDED_FLAGSTATS}.tmp ${FILENAME_EXTENDED_FLAGSTATS} || throw 30 "Could not move file"
mv ${FILENAME_ISIZES_STATISTICS}.tmp ${FILENAME_ISIZES_STATISTICS} || throw 31 "Could not move file"
mv ${FILENAME_DIFFCHROM_STATISTICS}.tmp ${FILENAME_DIFFCHROM_STATISTICS} || throw 32 "Could not move file"
mv ${FILENAME_READBINS_COVERAGE}.tmp ${FILENAME_READBINS_COVERAGE} || throw 34 "Could not move file"
mv ${FILENAME_GENOME_COVERAGE}.tmp ${FILENAME_GENOME_COVERAGE} || throw 35 "Could not move file"
mv ${tempFlagstatsFile} ${FILENAME_FLAGSTATS} || throw 33 "Could not move file"

# QC summary
# remove old warnings file if it exists (due to errors in run such as wrong chromsizes file)
[[ -f ${FILENAME_QCSUMMARY}_WARNINGS.txt ]] && rm ${FILENAME_QCSUMMARY}_WARNINGS.txt
${PERL_BINARY} $TOOL_WRITE_QC_SUMMARY \
    -p $PID \
    -s $SAMPLE \
    -r $RUN \
    -l $LANE \
    -w ${FILENAME_QCSUMMARY}_WARNINGS.txt \
    -f $FILENAME_FLAGSTATS \
    -d $FILENAME_DIFFCHROM_STATISTICS \
    -i $FILENAME_ISIZES_STATISTICS \
    -c $FILENAME_GENOME_COVERAGE \
    > ${FILENAME_QCSUMMARY}_temp \
    && mv ${FILENAME_QCSUMMARY}_temp $FILENAME_QCSUMMARY \
    || throw 12 "Error from writeQCsummary.pl"

# Produced qualitycontrol.json for OTP.
${PERL_BINARY} ${TOOL_QC_JSON} \
	${FILENAME_GENOME_COVERAGE} \
	${FILENAME_ISIZES_STATISTICS} \
	${FILENAME_FLAGSTATS} \
	${FILENAME_DIFFCHROM_STATISTICS} \
	> ${FILENAME_QCJSON}.tmp \
	|| throw 26 "Error when compiling qualitycontrol.json for ${FILENAME_QCJSON}, stopping here"
mv ${FILENAME_QCJSON}.tmp ${FILENAME_QCJSON} || throw 27 "Could not move file"

# plots are only made for paired end and not on convey
[[ ${useSingleEndProcessing-false} == true ]] && exit 0

${RSCRIPT_BINARY} ${TOOL_INSERT_SIZE_PLOT_SCRIPT} ${FILENAME_ISIZES_MATRIX} ${FILENAME_ISIZES_STATISTICS} ${FILENAME_ISIZES_PLOT}_temp "PE insertsize of ${bamname}" && mv ${FILENAME_ISIZES_PLOT}_temp ${FILENAME_ISIZES_PLOT} || ( echo "Error from insert sizes plotter" && exit 22)

${RSCRIPT_BINARY} ${TOOL_PLOT_DIFFCHROM} -i "$FILENAME_DIFFCHROM_MATRIX" -s "$FILENAME_DIFFCHROM_STATISTICS" -o "${FILENAME_DIFFCHROM_PLOT}_temp" && mv  ${FILENAME_DIFFCHROM_PLOT}_temp ${FILENAME_DIFFCHROM_PLOT} || ( echo "Error from chrom_diff.r" && exit 23)


