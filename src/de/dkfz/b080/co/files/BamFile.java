package de.dkfz.b080.co.files;

import de.dkfz.b080.co.common.ParallelizationHelper;
import de.dkfz.b080.co.methods.Common;
import de.dkfz.roddy.config.Configuration;
import de.dkfz.roddy.core.ExecutionContext;
import de.dkfz.roddy.execution.jobs.*;
import de.dkfz.roddy.knowledge.files.*;
import de.dkfz.roddy.knowledge.methods.GenericMethod;
//import sun.net.www.content.text.Generic;

import java.io.File;
import java.util.*;
import java.util.stream.Stream;

/**
 * A bam file is the binary version of a sam file and contains sequence data.
 *
 * @author heinold
 * @url http://www.broadinstitute.org/igv/bam
 */
public class BamFile extends BasicBamFile implements ITestdataSource {

    private BamIndexFile indexFile;
    private FlagstatsFile flagstatsFile;
    private TextFile extendedFlagstatsFile;
    private ChromosomeDiffFileGroup chromosomeDiffFileGroup;
    private ChromosomeDiffValueFile chromosomeDiffStatisticsFile;
    private ChromosomeDiffTextFile chromosomeDiffMatrixFile;
    private ChromosomeDiffPlotFile chromosomeDiffPlotFile;
    private InsertSizesFileGroup insertSizesFileGroup;
    private InsertSizesValueFile insertSizesStatisticsFile;
    private InsertSizesTextFile insertSizesMatrixFile;
    private InsertSizesPlotFile insertSizesPlotFile;
    private CoverageTextFile genomeCoverageTextFile;
    private CoverageTextFile targetCoverageTextFile;
    private CoverageTextFile rawBamCoverageTextFile;
    private CoverageTextFile readBinsCoverageTextFile;
    private CoverageTextFile targetsWithCoverageTextFile;
    private BamMetricsFile bamMetricsFile;
    private BamFile targetExtractedBamFile;
    private boolean _isTargetExractedFile;
    private OnTargetCoverageTextFile onTargetCoverageTextFile;
    private QCSummaryFile qcSummaryFile;
    private TextFile qcJsonFile;
    private TextFile dipStatisticsFile;
    private InsertSizesPlotFile dipStatisticsPlotFile;
    private MethylationMetaCheckpointFile methylationMetaCheckpointFile;
    private MethylationMetaMetricsCheckpointFile methylationMetaMetricsCheckpointFile;

    public BamFile(ConstructionHelperForBaseFiles helper) {
        super(helper);
    }
//
//    private static ConstructionHelperForBaseFiles adaptHelperBasedOnInput(ConstructionHelperForBaseFiles helper) {
//        if (helper instanceof ConstructionHelperForSourceFiles) {
//            return helper;
//        } else if (helper instanceof ConstructionHelperForGenericCreation) {
//            ConstructionHelperForGenericCreation h = (ConstructionHelperForGenericCreation) (helper);
//            /* This is a bit complicated and I did not know how to solve it in a more convenient way, but:
//             * BamFile fomerly had a lot of different constructors for different input. Most of them
//             * decreased the file stage level of the parent object. Now we have only one constructor, which just checks
//             * if we have a BamFile as the parent file(in this case, we do not decrease levels)
//             */
//            FileStageSettings newFS = h.fileStageSettings;
//            if (!(h.parentObject instanceof BamFile)) {
//                if(h.parentObject instanceof FileGroup) {
//                    newFS = ((BaseFile)((FileGroup)h.parentObject).getFilesInGroup().get(0)).getFileStage().decreaseLevel();
//                } else {
//                    newFS = ((BaseFile)h.parentObject).getFileStage().decreaseLevel();
//                }
//            } else {
//                newFS = ((BamFile)h.parentObject).getFileStage();
//            }
//
//            //Reset the previous construction helper with the new filestage settings. Do a copy?
//
//        } else {
//            throw new RuntimeException("Oh oh, there is an unresolved case for the BamFile constructor.");
//        }
//    }

    /**
     * "Copy" constructor for Bam conversion.
     *
     * @param minorBam
     */
    public BamFile(BasicBamFile minorBam) {
        super(minorBam);
    }

    @Override
    public void setAsTemporaryFile() {
        super.setAsTemporaryFile();
        if (hasIndex()) index().setAsTemporaryFile(); //Also indexFile files created by a or for a temporary bam file are temporary
    }

    public String getType() {
        return _isTargetExractedFile ? "targetExtractedBam" : "commonBam";
    }

    public Sample getSample() {
        return ((COFileStageSettings) getFileStage()).getSample();
    }

    public boolean hasIndex() {
        return indexFile != null;
    }

    public boolean hasFlagstats() {
        return flagstatsFile != null;
    }

    public boolean hasInsertSizes() {
        return insertSizesFileGroup != null;
    }

    public boolean hasChromosomeDiff() {
        return chromosomeDiffFileGroup != null;
    }

    public boolean hasCoverage() {
        return genomeCoverageTextFile != null;
    }

    public boolean hasMetricsFile() {
        return bamMetricsFile != null;
    }

    public void setMetricsFile(BamMetricsFile metricsFile) {
        this.bamMetricsFile = metricsFile;
    }

    public boolean hasOnTargetCoverageTextFile() {
        return onTargetCoverageTextFile != null;
    }

    public boolean isRawBamFile() {
        return !_isTargetExractedFile;
    }

    public void setIsTargetExtracted() {
        _isTargetExractedFile = true;
    }

    public boolean isTargetExtractedBamFile() {
        return _isTargetExractedFile;
    }

    public boolean hasTargetExtractedBamFile() {
        return targetExtractedBamFile != null;
    }

    public boolean hasRawBamCoverageTextFile() {
        return rawBamCoverageTextFile != null;
    }

    public boolean hasTargetCoverageTextFile() {
        return targetCoverageTextFile != null;
    }

    public void setIndexFile(BamIndexFile indexFile) {
        this.indexFile = indexFile;
    }

    public void setFlagstatsFile(FlagstatsFile flagstatsFile) {
        this.flagstatsFile = flagstatsFile;
    }

    public BamIndexFile getIndexFile() {
        return indexFile;
    }

    public ChromosomeDiffFileGroup getChromosomeDiffFileGroup() {
        return chromosomeDiffFileGroup;
    }

    public CoverageTextFile getGenomeCoverageTextFile() {
        return genomeCoverageTextFile;
    }

    public FlagstatsFile getFlagstatsFile() {
        return flagstatsFile;
    }

    public void setExtendedFlagstatsFile(TextFile extendedFlagstatsFile) {
        this.extendedFlagstatsFile = extendedFlagstatsFile;
    }

    public TextFile getExtendedFlagstatsFile() {
        return extendedFlagstatsFile;
    }

    public BamMetricsFile getMetricsFile() {
        return bamMetricsFile;
    }

    public InsertSizesFileGroup getInsertSizesFileGroup() {
        return insertSizesFileGroup;
    }

    public CoverageTextFile getTargetCoverageTextFile() {
        return targetCoverageTextFile;
    }

    public CoverageTextFile getRawBamCoverageTextFile() {
        return rawBamCoverageTextFile;
    }

    public MethylationMetaCheckpointFile getMethylationMetaCheckpointFile() {
        return methylationMetaCheckpointFile;
    }

    public MethylationMetaMetricsCheckpointFile getMethylationMetaMetricsCheckpointFile() { return methylationMetaMetricsCheckpointFile; }

    public QCSummaryFile getQcSummaryFile() {
        return qcSummaryFile;
    }

    public TextFile getQcJsonFile() {
        return qcJsonFile;
    }

    public void setChromosomeDiffStatisticsFile(ChromosomeDiffValueFile chromosomeDiffStatisticsFile) {
        this.chromosomeDiffStatisticsFile = chromosomeDiffStatisticsFile;
    }

    public void setChromosomeDiffMatrixFile(ChromosomeDiffTextFile chromosomeDiffMatrixFile) {
        this.chromosomeDiffMatrixFile = chromosomeDiffMatrixFile;
    }

    public void setChromosomeDiffPlotFile(ChromosomeDiffPlotFile chromosomeDiffPlotFile) {
        this.chromosomeDiffPlotFile = chromosomeDiffPlotFile;
    }

    public void setInsertSizesStatisticsFile(InsertSizesValueFile insertSizesStatisticsFile) {
        this.insertSizesStatisticsFile = insertSizesStatisticsFile;
    }

    public void setInsertSizesMatrixFile(InsertSizesTextFile insertSizesMatrixFile) {
        this.insertSizesMatrixFile = insertSizesMatrixFile;
    }

    public void setInsertSizesPlotFile(InsertSizesPlotFile insertSizesPlotFile) {
        this.insertSizesPlotFile = insertSizesPlotFile;
    }

    public void setGenomeCoverageTextFile(CoverageTextFile genomeCoverageTextFile) {
        this.genomeCoverageTextFile = genomeCoverageTextFile;
    }

    public void setReadBinsCoverageTextFile(CoverageTextFile readBinsCoverageTextFile) {
        this.readBinsCoverageTextFile = readBinsCoverageTextFile;
    }

    public void setBamMetricsFile(BamMetricsFile bamMetricsFile) {
        this.bamMetricsFile = bamMetricsFile;
    }

    public void setQcSummaryFile(QCSummaryFile qcSummaryFile) {
        this.qcSummaryFile = qcSummaryFile;
    }

    public void setQcJsonFile(TextFile qcJsonFile) {
        this.qcJsonFile = qcJsonFile;
    }

    public ChromosomeDiffValueFile getChromosomeDiffStatisticsFile() {
        return chromosomeDiffStatisticsFile;
    }

    public ChromosomeDiffTextFile getChromosomeDiffMatrixFile() {
        return chromosomeDiffMatrixFile;
    }

    public ChromosomeDiffPlotFile getChromosomeDiffPlotFile() {
        return chromosomeDiffPlotFile;
    }

    public InsertSizesValueFile getInsertSizesStatisticsFile() {
        return insertSizesStatisticsFile;
    }

    public InsertSizesTextFile getInsertSizesMatrixFile() {
        return insertSizesMatrixFile;
    }

    public InsertSizesPlotFile getInsertSizesPlotFile() {
        return insertSizesPlotFile;
    }

    public BamMetricsFile getBamMetricsFile() {
        return bamMetricsFile;
    }

    public BamFile getTargetExtractedBamFile() {
        return targetExtractedBamFile;
    }

    public void setTargetCoverageTextFile(CoverageTextFile targetCoverageTextFile) {
        this.targetCoverageTextFile = targetCoverageTextFile;
    }

    public CoverageTextFile getTargetsWithCoverageTextFile() {
        return targetsWithCoverageTextFile;
    }

    public void setTargetsWithCoverageTextFile(CoverageTextFile targetsWithCoverageTextFile) {
        this.targetsWithCoverageTextFile = targetsWithCoverageTextFile;
    }

    public void setOnTargetCoverageTextFile(OnTargetCoverageTextFile tf) {
        this.onTargetCoverageTextFile = tf;
    }

    public OnTargetCoverageTextFile getOnTargetCoverageTextFile() {
        return onTargetCoverageTextFile;
    }

    @Override
    /**
     * Performs the following operations in the necessary order:<br> <ul>
     * <il>indexFile</il> <il>flagstat</il> <il>insertsizes</il> <il>chromosome
     * diff</il> <il>coverage</il> <il>qcsummary</il> </ul>
     */
    public void runDefaultOperations() {
        index();
        flagstat();
        if (!getExecutionContext().getConfiguration().getConfigurationValues().getBoolean("useSingleEndProcessing", false))
            determineInsertSizes();
        diffChroms();
//        if (getExecutionContext().getConfiguration().getConfigurationValues().getBoolean("runGenomeCoverage", true))
//            calcCoverage();
    }

    @ScriptCallingMethod
    public BamIndexFile index() {
        if (indexFile == null)
            indexFile = GenericMethod.callGenericTool(COConstants.TOOL_SAMTOOLS_INDEX, this);
        return indexFile;
    }

    @ScriptCallingMethod
    public FlagstatsFile flagstat() {
        if (flagstatsFile == null)
            flagstatsFile = GenericMethod.callGenericTool(COConstants.TOOL_SAMTOOLS_FLAGSTAT, this);
        return flagstatsFile;
    }

    @ScriptCallingMethod
    public Tuple2<MethylationMetaMetricsCheckpointFile,MethylationMetaCheckpointFile> methylationCallingMeta() {
        if (methylationMetaCheckpointFile == null || methylationMetaMetricsCheckpointFile == null) {
            Tuple2<MethylationMetaMetricsCheckpointFile, MethylationMetaCheckpointFile> res = GenericMethod.callGenericTool("methylationCallingMeta", this);
            methylationMetaCheckpointFile = res.value1;
            methylationMetaMetricsCheckpointFile = res.value0;
        }
        return new Tuple2(methylationMetaCheckpointFile, methylationMetaMetricsCheckpointFile);
    }

    @ScriptCallingMethod
    public Tuple2<MethylationMetaMetricsCheckpointFile,MethylationMetaCheckpointFile> libraryMethylationCallingMeta() {
        return methylationCallingMeta();
    }

    @ScriptCallingMethod
    public Tuple2<MethylationMetaMetricsCheckpointFile,MethylationMetaCheckpointFile>  mergedMethylationCallingMeta() {
        return methylationCallingMeta();
    }


    @ScriptCallingMethod
    public TextFile performPurityAnalysis(BamFile bamControlMerged) {
        return GenericMethod.callGenericTool(COConstants.TOOL_PURITY_ESTIMATION, this, bamControlMerged);
    }

    @ScriptCallingMethod
    public BamMetricsCollectionFileGroup collectMetrics() {
        final ExecutionContext run = getExecutionContext();
        final Configuration cfg = run.getConfiguration();
        BamMetricsCollectionFileGroup resultGroup = new BamMetricsCollectionFileGroup(new LinkedList());

        String metricsOutputDirectory = cfg.getConfigurationValues().getString("metricsOutputDirectory");

        List<BaseFile> pFiles = Arrays.asList((BaseFile) this);

        final String TOOL = COConstants.TOOL_COLLECT_BAM_METRICS;
        Map<String, Object> parameters = run.getDefaultJobParameters(TOOL);
        parameters.put("FILENAME", this.getPath().getAbsolutePath());
        parameters.put("DIR_METRICS", run.getRuntimeService().getDirectory(metricsOutputDirectory, run).getAbsolutePath());
        parameters.put("sample", ((COFileStageSettings) this.getFileStage()).getSample().getName());
        parameters.put("SAMPLE", ((COFileStageSettings) this.getFileStage()).getSample().getName());

        BamMetricsAlignmentSummaryFile alnSummaryFile = new BamMetricsAlignmentSummaryFile(this);
        resultGroup.addFile(alnSummaryFile);

        Job job = new Job(run, run.createJobName(this, TOOL), TOOL, parameters, pFiles, Arrays.asList((BaseFile) alnSummaryFile));
        alnSummaryFile.setCreatingJobsResult(job.run());
        return resultGroup;
    }

    @ScriptCallingMethod
    public InsertSizesFileGroup determineInsertSizes() {
        if (insertSizesFileGroup == null)
            insertSizesFileGroup = GenericMethod.callGenericTool("insertSizes", this);
        return insertSizesFileGroup;
    }

    @ScriptCallingMethod
    public ChromosomeDiffFileGroup diffChroms() {
        if (chromosomeDiffFileGroup == null)
            chromosomeDiffFileGroup = GenericMethod.callGenericTool("chromosomeDiff", this);
        return chromosomeDiffFileGroup;
    }

    @ScriptCallingMethod
    public CoverageTextFile calcCoverage() {
        if (genomeCoverageTextFile == null)
            genomeCoverageTextFile = GenericMethod.callGenericTool(Common.GENOMECOVERAGE, this, this.index(), "COVERAGE_TYPE=" + CoverageTextFile.CoverageType.Default.toString());
        return genomeCoverageTextFile;
    }

    @ScriptCallingMethod
    public CoverageTextFile targetCoverage() {
        if (targetCoverageTextFile == null)
            targetCoverageTextFile = GenericMethod.callGenericTool(Common.GENOMECOVERAGE, this, this.index(), "COVERAGE_TYPE=" + CoverageTextFile.CoverageType.TargetEnrichment.toString());
        return targetCoverageTextFile;
    }

    @ScriptCallingMethod
    public CoverageTextFile rawBamCoverage() {
        if (rawBamCoverageTextFile == null)
            rawBamCoverageTextFile = GenericMethod.callGenericTool(Common.GENOMECOVERAGE, this, this.index(), "COVERAGE_TYPE=" + CoverageTextFile.CoverageType.RawBamTargetEnrichment.toString());
        return rawBamCoverageTextFile;
    }

    public Tuple3<BamFile, BamFile, TextFile> extractTelomeres() {
        return GenericMethod.callGenericTool("telomereExtraction", this);
    }

    @ScriptCallingMethod
    public BamFile extractTargetsCalculateCoverage() {
        if (targetExtractedBamFile == null)
            targetExtractedBamFile = GenericMethod.callGenericTool(COConstants.TARGET_EXTRACTION_AND_COVERAGE_SLIM, this, this.getGenomeCoverageTextFile(), "SAMPLE=" + getSample().getName());
        return targetExtractedBamFile;
    }


    public boolean hasReadBinsCoverageTextFile() {
        return readBinsCoverageTextFile != null;
    }

    @ScriptCallingMethod
    public CoverageTextFile getReadBinsCoverageTextFile() {
        return calcReadBinsCoverage();
    }

    /**
     * Counts reads in Window-Size kb
     */
    @ScriptCallingMethod
    public CoverageTextFile calcReadBinsCoverage() {
        if (readBinsCoverageTextFile == null)
            readBinsCoverageTextFile = GenericMethod.callGenericTool("readBinsCoverage", this);
        return readBinsCoverageTextFile;
    }

    public void performPostMergeQCAnalysis() {
        if (insertSizesPlotFile == null && chromosomeDiffPlotFile == null) {
            Tuple13<BamIndexFile, FlagstatsFile, TextFile, BamMetricsFile, ChromosomeDiffValueFile, ChromosomeDiffTextFile, ChromosomeDiffPlotFile, InsertSizesValueFile, InsertSizesTextFile, InsertSizesPlotFile, CoverageTextFile, CoverageTextFile, QCSummaryFile> results = GenericMethod.callGenericTool("postMergeQCAnalysis", this, "SAMPLE=" + this.getSample().getName());
            setIndexFile(results._a);
            setFlagstatsFile(results._b);
            setExtendedFlagstatsFile(results._c);
            setBamMetricsFile(results._d);
            setChromosomeDiffStatisticsFile(results._e);
            setChromosomeDiffMatrixFile(results._f);
            setChromosomeDiffPlotFile(results._g);
            setInsertSizesStatisticsFile(results._h);
            setInsertSizesMatrixFile(results._i);
            setInsertSizesPlotFile(results._j);
            setGenomeCoverageTextFile(results._k);
            setReadBinsCoverageTextFile(results._l);
            setQcSummaryFile(results._m);
        }
    }

    public List<COBaseFile> getCreatedFiles() {
        List<COBaseFile> files = new LinkedList<>();
        if (hasIndex()) {
            files.add(getIndexFile());
        }
        if (hasFlagstats()) {
            files.add(flagstatsFile);
        }
        if (hasChromosomeDiff()) {
            files.add(getChromosomeDiffFileGroup().getTextFile());
            files.add(getChromosomeDiffFileGroup().getPlotFile());
            files.add(getChromosomeDiffFileGroup().getValueFile());
        }
        if (hasCoverage()) {
            files.add(getGenomeCoverageTextFile());
        }
        if (hasInsertSizes()) {
            files.add(getInsertSizesFileGroup().getTextFile());
            files.add(getInsertSizesFileGroup().getPlotFile());
            files.add(getInsertSizesFileGroup().getValueFile());
        }
        if (hasMetricsFile()) {
            files.add(getMetricsFile());
        }
        //TODO Find out for which files / operations an indexFile is used
        if (hasTargetCoverageTextFile()) {
            files.add(getTargetCoverageTextFile());
        }
        if (hasRawBamCoverageTextFile()) {
            files.add(getRawBamCoverageTextFile());
        }
        if (hasTargetExtractedBamFile()) {
            files.add(getTargetExtractedBamFile());
        }
        if (hasOnTargetCoverageTextFile()) {
            files.add(getOnTargetCoverageTextFile());
        }

        return files;
    }

    @ScriptCallingMethod
    public QCSummaryFile createQCSummaryFile() {
        if (qcSummaryFile == null)
            qcSummaryFile = QCSummaryFile.createFromFileList(getExecutionContext(), this, getCreatedFiles());
        return qcSummaryFile;
    }

    @Override
    public String toString() {
        return String.format("BamFile %s hasIndex=%b hasFlagstats=%b hasChromosomeDiff=%b hasCoverage=%b",
                super.toString(), hasIndex(), hasFlagstats(), hasChromosomeDiff(), hasCoverage());
    }

    @Override
    public boolean createTestData() {
        throw new UnsupportedOperationException("Creating test data for bamfiles not supported yet.");
    }
}
