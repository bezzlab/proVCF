# Protein Variant Calling Format (proVCF)

A proVCF file (Figure 1) contains a header and a data section similar to a VCF file. Header lines describe an arbitrary number of items of meta information used in the data section in a standardized format. The meta-information lines start with ## and the header of the tab delimited data section begins with #. Unlike a VCF file, there are six mandatory meta-information lines describing file format, software used to create the proVCF file (source), path to the reference file (reference), source of reference (referenceSource, e.g. Uniprot13), date of reference file collection (referenceDate), taxonomy id of the species of the sample (species, e.g. for human it is 9606 ). The data section of the proVCF file contains eight mandatory columns where CHROM represents reference protein id, POS is 1-based position of the variation in the reference. I store protein ID in the CHROM column as proteins are the reference for the protein variations, hence POS is the position of the variation within the reference protein. ID represents unique identifier of the variation, REF is reference amino acid, ALT contains alternative amino acid sequence, QUAL holds PSM or Peptide level q-value when available, and variation site filtering information is in FILTER, e.g. in the figure TRANS suggests the variation has transcriptomic level evidence and PASS means peptide evidence was found for the variation site. Lastly, INFO is a semicolon-separated, user defined key-value pair of additional information. Post-translated modifications are reported in the INFO column with the MOD keyword and a combination of three components. The first component is the PTM associated amino acid symbol, the second is the position of the amino acid within the ALT sequence and last is the mass difference or UniMod ID or name of the modification (e.g. phosphorylation). Chromosome information about the protein can be added to the INFO column using the Chromosome key.

## Conventions and reserved keywords

VCF files contain a number of standardized keywords with specific meaning; they are not allowed to be re-defined for proVCF file. I have added a number of keywords for the proVCF format related to protein polymorphisms. The following section describes reserved tags of current proVCF version:
MOD: This keyword is reserved to store PTMs. Each MOD key may have multiple values separated by comma. Each value has three components, amino acid symbol, position of the amino acid within ALT sequence and mass shift or Unimod ID of the modification or common name.
APOS: This keyword is reserved to store position of the variation in the alternative sequence. Data type is integer and should be defined in the header when used.
Chromosome: As the name suggest, this key is used for storing chromosome information of the reference when available.
AID: Identification number or sequence header from fasta file of the alternative sequence.
Type: This keyword is reserved to describe type of polymorphisms.


![Alt text](proVCF.png)
![Alt text](proVCFEg.png)

# proVCF Tools

![Alt text](proVCFTools.png)
