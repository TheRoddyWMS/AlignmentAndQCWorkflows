/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.methods
import de.dkfz.b080.co.files.*
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.execution.jobs.Job
import de.dkfz.roddy.execution.jobs.BEJobResult
import de.dkfz.roddy.execution.jobs.ScriptCallingMethod
import de.dkfz.roddy.execution.jobs.StaticScriptProviderClass
import de.dkfz.roddy.knowledge.files.BaseFile
import de.dkfz.roddy.knowledge.files.FileObject

/**
 *
 * @author michael
 */
@groovy.transform.CompileStatic
@StaticScriptProviderClass
class Common {

    public static final String QCSUMMARY = "qcSummary";

    public static final String PID = "DataSet";

    private static LaneFile recursivelySearchLaneFile(List<BaseFile> files) {
        LaneFile lf = null;
        for (BaseFile bf : files) {
            if (bf instanceof LaneFile) {
                lf = (LaneFile) bf;
            } else {
                lf = recursivelySearchLaneFile(bf.getParentFiles());
                if (lf != null) {
                    break;
                }
            }
        }
        return lf;
    }

    @ScriptCallingMethod
    public static QCSummaryFile createQCSummaryFileFromList(ExecutionContext run, BamFile bamFile, List<COBaseFile> files) {
        QCSummaryFile qcSummaryFile = BaseFile.constructManual(QCSummaryFile, bamFile, files as List<FileObject>, null, null,null,null,null,null) as QCSummaryFile;

        LaneFile laneFile = recursivelySearchLaneFile([(BaseFile) bamFile]);
        COFileStageSettings bamFileFileStage = (COFileStageSettings) bamFile.getFileStage()

        String sample;
        if(laneFile != null) {
            sample = laneFile.getSample().getName();
        } else {
            sample = bamFileFileStage.getSample().getName();
        }
        String runId = "all_merged";
        String lane = "all_merged";

        if (bamFileFileStage.stage.isMoreDetailedOrEqualTo(COFileStage.RUN)) {
            runId = bamFileFileStage.runID;
            if (bamFileFileStage.stage.isMoreDetailedOrEqualTo(COFileStage.LANE))
                lane = bamFileFileStage.laneId;

        }
        def temp = run.getDefaultJobParameters(QCSUMMARY);
        Map<String, Object> parameters = (Map<String, Object>)temp;
        parameters.putAll([
                "SAMPLE": sample,
                "RUN": runId,
                "LANE": lane,
                "FILENAME_QCSUM": qcSummaryFile.absolutePath
        ]);
        files.each {
            it ->
                BaseFile bf = (BaseFile) it;
                if (bf instanceof FlagstatsFile) {
                    parameters["FILENAME_FLAGSTAT"] = bf.path.absolutePath;
                } else if (bf instanceof ChromosomeDiffValueFile) {
                    parameters["FILENAME_DIFFCHROM"] = bf.path.absolutePath;
                } else if (bf instanceof InsertSizesValueFile) {
                    parameters["FILENAME_ISIZE"] = bf.path.absolutePath;
                } else if (bf instanceof BamMetricsFile) {
                    parameters["FILENAME_METRICS"] = bf.path.absolutePath;
//                } else if (bf instanceof BamFile && ((BamFile) bf).isTargetExtractedBamFile()) {
                } else if (bf instanceof CoverageTextFile && bamFile.getTargetCoverageTextFile() == bf) {
                    parameters["FILENAME_TARGETCAPTURE"] = bf.path.absolutePath;
                } else if (bf instanceof CoverageTextFile) {
                    parameters["FILENAME_COVERAGE"] = bf.path.absolutePath;
                }
        }

        BEJobResult jobResult = new Job(run, run.createJobName(files[0], QCSUMMARY), QCSUMMARY, parameters, new LinkedList<BaseFile>(files), [(COBaseFile) qcSummaryFile] as List<BaseFile>).run();
        qcSummaryFile.setCreatingJobsResult(jobResult);
        return qcSummaryFile;
    }
}
