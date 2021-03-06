/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.common

import de.dkfz.b080.co.files.*
import de.dkfz.roddy.core.ExecutionContext
import de.dkfz.roddy.execution.jobs.BEFakeJobID
import de.dkfz.roddy.execution.jobs.BEJobResult
import de.dkfz.roddy.execution.jobs.FakeBEJob
import de.dkfz.roddy.knowledge.files.BaseFile
import de.dkfz.roddy.tools.LoggerWrapper

@groovy.transform.CompileStatic
class QCPipelineScriptFileServiceHelper {

    private static final LoggerWrapper logger = LoggerWrapper.getLogger(QCPipelineScriptFileServiceHelper.class.name);

    /**
     *
     * @param context
     * @param sample
     * @param runName
     * @param files
     * @return
     */
    static List<LaneFileGroup>  sortAndPairLaneFilesToGroupsForSampleAndRun(ExecutionContext context, Sample sample, String libraryID, String runName, List<File> files) {
        AlignmentAndQCConfig config = new AlignmentAndQCConfig(context)
        LibraryID libId = libraryID == null ? new LibraryID(config.defaultLibraryName) : new LibraryID(libraryID)

        List<File> sortedFiles = [];
        LinkedHashMap<String, LinkedList<File>> sortedFileGroups = new LinkedHashMap<String, LinkedList<File>>();
        List<LaneFileGroup> fileGroups = new LinkedList<LaneFileGroup>();
        if (files.isEmpty()) {
            return fileGroups;
        }
        sortedFiles += files;
        sortedFiles.sort(new Comparator<File>() {
            @Override
            int compare(File o1, File o2) {
                return o1.getAbsolutePath().compareTo(o2.getAbsolutePath());
            }
        })

        boolean[] paired = new boolean[sortedFiles.size()];
        boolean singleEndProcessing = context.getConfiguration().getConfigurationValues().getBoolean("useSingleEndProcessing", false);
        if (singleEndProcessing) {
            for (int i = 0; i < sortedFiles.size(); i++) {
                File _f0 = sortedFiles[i];
                File _f1 = new File(_f0.getAbsolutePath() + "_dummySecondary");
                IndexID index1 = new IndexID("R1");
                IndexID index2 = new IndexID("R2");
                String lane = String.format("L%03d", i);
                LaneID laneId = new LaneID(String.format("%s_%s_%s_%s_%s", context.getDataSet().getId(), sample.getName(), libId, runName, lane, index1));


                BEJobResult result = getFileExistedFakeJobResult()
                LinkedList<LaneFile> filesInGroup = new LinkedList<LaneFile>(Arrays.asList(
                        (LaneFile) BaseFile.constructSourceFile(LaneFile, _f0, context,
                                new COFileStageSettings(laneId, index1, 0, new RunID(runName), libId, sample, context.getDataSet(), COFileStage.INDEXEDLANE),
                                result),
                        (LaneFile) BaseFile.constructSourceFile(LaneFile, _f1, context,
                                new COFileStageSettings(laneId, index2, 1, new RunID(runName), libId, sample, context.getDataSet(), COFileStage.INDEXEDLANE),
                                result)
                ));
                filesInGroup[1].setFileIsValid();
                fileGroups << new LaneFileGroup(context, laneId.toString(), runName, sample, filesInGroup)
            }
        } else {
            for (int i = 0; i < sortedFiles.size() - 1; i++) {
                File _f0 = sortedFiles[i];
                File _f1 = sortedFiles[i + 1];
                String f0 = _f0.name;
                String f1 = _f1.name;
                int diffCount = 0;
                if (f0.size() == f1.size()) {
                    for (int c = 0; c < f0.size(); c++) {
                        if (f0[c] != f1[c])
                            diffCount++;
                    }
                }
                if (diffCount > 1) {
                    continue;   //Files are not equal enough so skip to the next pair
                } else {
                    int j = i + 1; // Idea/IDE failure? Shows an error...
                    paired[i] = true;    //Detect single files with this
                    paired[j] = true;
                    i++;
                    String[] blocks0 = f0.split("_").reverse();
                    String[] blocks1 = f1.split("_").reverse(); //Rightmost non unique block is the indexFile
                    IndexID index0 = new IndexID("");
                    IndexID index1 = new IndexID("");
                    int indexOfIndex = blocks0.size() - 1;
                    for (int b = 0; b < blocks0.size(); b++) {
                        indexOfIndex--;
                        if (blocks0[b] == blocks1[b]) continue;
                        index0 = new IndexID(blocks0[b]);
                        index1 = new IndexID(blocks1[b]);
                        break;
                    }
                    blocks0 = blocks0.reverse().toList()[0..indexOfIndex];
                    LaneID laneId = new LaneID(blocks0.join("_"));

                    LinkedList<LaneFile> filesInGroup = new LinkedList<LaneFile>();

                    BEJobResult result = getFileExistedFakeJobResult()

                    filesInGroup << (LaneFile) BaseFile.constructSourceFile(LaneFile, _f0, context,
                            new COFileStageSettings(laneId, index0, 0, new RunID(runName), libId, sample, context.getDataSet(), COFileStage.INDEXEDLANE),
                            result);
                    filesInGroup << (LaneFile) BaseFile.constructSourceFile(LaneFile, _f1, context,
                            new COFileStageSettings(laneId, index1, 1, new RunID(runName), libId, sample, context.getDataSet(), COFileStage.INDEXEDLANE),
                            result);

                    fileGroups << new LaneFileGroup(context, laneId.toString(), runName, sample, filesInGroup)
                }
            }
        }
        if(fileGroups.size() == 0) {
            logger.postAlwaysInfo("There were no files for sample ${sample.getName()} and run ${runName}" )
        }
        return fileGroups;
    }

    static BEJobResult getFileExistedFakeJobResult() {
        return new BEJobResult(null, new FakeBEJob(new BEFakeJobID(BEFakeJobID.FakeJobReason.FILE_EXISTED)), null, null, null, null)
    }
}