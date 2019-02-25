# Alignment and Quality Control Plugin for Roddy

This plugins contains alignment and quality control related [Roddy](https://github.com/TheRoddyWMS/Roddy) workflows.

This is the branch for 1.2.51 maintenance branch with only minimal documentation. Please see the master branch for details.


## Change Logs

* 1.2.51-2 (branch-specific changes)
  - Improved error checking and reporting for BWA and surrounding pipe

* 1.2.51-1 (branch-specific changes)
  - Update to Roddy 3.0 release (official LSF-capable release)
  - Bugfix in tbi-lsf-cluster.sh

* 1.2.51
  - Lifted 1.1.51 to Roddy 2.4 (development-only release)
  - FASTQ quality classification (Xavier)
  - BAM termination sequence check
  - Bugfixes in WGBS (off-by-one, meth-call splitting CG/CH)
  - Further bugfixes

* 1.1.51
  - Improved and error checking
  - Pre-submission executability checks
  - Progress on WGBS workflow
  - Tuning (sambamba flagstat -t 1)

* 1.1.2
  - Rename from QualityControlWorkflows to AlignmentAndQCWorkflows
  - First development version of WGBS workflow
  - Adaptation to Roddy API changes (FileSystemAccessProvider)
