/*
 * Copyright (c) 2018 German Cancer Research Center (DKFZ).
 *
 * Distributed under the MIT License (license terms are at https://github.com/DKFZ-ODCF/AlignmentAndQCWorkflows).
 */

package de.dkfz.b080.co.methods

import de.dkfz.b080.co.files.*;
import de.dkfz.roddy.execution.jobs.*;
import de.dkfz.roddy.knowledge.files.Tuple2;
import de.dkfz.roddy.knowledge.files.Tuple3
import de.dkfz.roddy.knowledge.files.Tuple4;
import de.dkfz.roddy.knowledge.methods.GenericMethod


@StaticScriptProviderClass
public final class ACEseq {

    @ScriptCallingMethod
    public static Tuple2<TextFile, TextFile> annotateCovWindows(CoverageTextFile windowCoverageTextFile, Sample sample) {
        return (Tuple2<TextFile, TextFile>) GenericMethod.callGenericTool("annotateCovWindows", windowCoverageTextFile, "SAMPLE=" + sample.name);
    }

    @ScriptCallingMethod
    public static TextFile mergeAndFilterCovWindows(TextFile annotatedCoverageFile) {
        return (TextFile) GenericMethod.callGenericTool("mergeAndFilterCovWindows", annotatedCoverageFile);
    }

    @ScriptCallingMethod
    public static Tuple4<TextFile, TextFile, TextFile, TextFile> correctGc(TextFile mergedAndFilteredCovWinFile, Sample sample) {
        return (Tuple4<TextFile, TextFile, TextFile, TextFile>) GenericMethod.callGenericTool("correctGc", mergedAndFilteredCovWinFile, "SAMPLE=" + sample.name);
    }

    public static Tuple3<TextFile, TextFile, TextFile> aceSeqQc(CoverageTextFile windowedCoverageTextFile, Sample sample) {
        Tuple2<TextFile, TextFile> annotationResult = ACEseq.annotateCovWindows(windowedCoverageTextFile, sample);
        TextFile mergedAndFilteredCoverageWindowFile = ACEseq.mergeAndFilterCovWindows(annotationResult.value0);
        Tuple3<TextFile, TextFile, TextFile, TextFile> correctedWindowFile = ACEseq.correctGc(mergedAndFilteredCoverageWindowFile, sample);
        return new Tuple3(annotationResult, mergedAndFilteredCoverageWindowFile, correctedWindowFile)
    }

}

