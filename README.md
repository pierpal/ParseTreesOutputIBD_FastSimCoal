# Parse Newick trees from FastSimCoal and output IBD segments

This program parses an ancestral recombination graph (ARG) stored in Newick trees produced by FastSimCoal and outputs IBD segments. The user must specify the recombination rate (cM/Mb), the length of the simulated region (in bp), the minimum length of the IBD segments to be output, and the total number of trees in the file. The script below provides an example. IBD segments are defined as regions for which a pair of individuals share the same most recent common ancestor (MRCA). The MRCA is determined based on coalescent time (note that there is a very small possibility that this is not unambiguous).

Dependencies: uses the Forester library (https://code.google.com/p/forester/) to parse Newick trees.

Contact: ppalama AT hsph DOT harvard DOTAGAIN edu

References: this tool was developed as part of the analysis for
- P. F. Palamara, T. Lencz, A. Darvasi, I. Pe'er. "Length distributions of identity by descent reveal fine-scale demographic history". The American Journal of Human Genetics, 2012.
- P. F. Palamara, I. Pe'er. "Inference of historical migration rates via haplotype sharing". Bioinformatics, 2013.

### Example
    # this assumes the output of FastSimCoal is in the folder "fastSimCoalFileRoot"

    OUT=fastSimCoalFileRoot
    recRate=2
    chrLen=100000000
    minLen=0.5

    cat $OUT/${OUT}_1_true_trees.trees | awk '$1=="tree"' | tr '_' ' ' | awk '{printf $7"\t"; for (i=10;i<=NF;i++) printf $i; print ""}' | gzip --best -c -v - > $OUT/$OUT.trees.gz
    nTrees=`gunzip -c $OUT/$OUT.trees.gz | wc -l`
    gunzip -c $OUT/$OUT.trees.gz | java -jar ParseTreesOutputIBD_FastSimCoal.jar $recRate $chrLen $minLen $nTrees | gzip --best -c -v - > $OUT/$OUT.match.gz
