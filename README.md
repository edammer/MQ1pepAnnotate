# MQ1pepAnnotate
Find and Annotate MS/MS Spectra for One-peptide Identifications in Thermo .RAW files searched using MaxQuant

This script leverages the protViz and rawrr R packages to produce annotated MS/MS peptide spectral matches.  Suitable for generating required supplemental data and table for publications in proteomic journals that require quality control and visualization of a dataset's full compendium of such peptides identified, quantified, and contributing to protein isoform-level quantitative rollup.

The parallel processing script version uses foreach to parallelize RAW file reading and output after generating the .csv table and data frame containing scan numbers of all "one-hit wonder" peptides' best-scored (by MaxQuant) peptide spectral match (PSM) across all raw files in the MaxQuant search, with each thread handling a separate raw file and PDF output. The PDF outputs can be merged later using a PDF file editing app like Adobe Acrobat Pro, or online free apps to do this.  Parallelized processing time for 31 .RAW files (~1 GB each) on a 2.2 GHz CPU with as many threads as RAW files completed in 60 seconds, annotating ~500 PSMs.

The string denoting oxidized methionine, deamidated asparagine or glutamine, and acetylated protein N-termini changed in 2024 (MaxQuat v2.4.x+), and an updated script handling these strings, assigning the appropriate mass shifts to annotations, is a second .R script version of the parallelized version, <A href="https://github.com/edammer/MQ1pepAnnotate/blob/main/Annotate%20One-Hit%20Wonders%20RAW%20extracted%20PSMs%20with%20MaxQuant%20v2.4%2B%20Identifications%5BparallelProcessing%5D.R">here</a>.

Please contact <a href="mailto:edammer@emory.edu">Eric Dammer</a> if you have questions about using or further adapting this script for different post-translational modifications of specific known mass shift or for further development to automate this script, or use it with other proteomic search engine outputs than MaxQuant /txt/ folder contents and Thermo .RAW files, and I will be happy to work with you.

-Eric
