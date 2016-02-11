package de.dkfz.b080.co.files;

/**
 *
 * @author michael
 */
public class GenomeCoveragePlotFile extends COBaseFile {

    public GenomeCoveragePlotFile(CoverageTextFile parentFile) {
        super(new ConstructionHelperForManualCreation(parentFile, null, null,null,null,null,null,null));
    }

    public GenomeCoveragePlotFile(BamFileGroup group) {
        super(new ConstructionHelperForManualCreation(group, null, null,null,null,null,null,null));
    }

}
