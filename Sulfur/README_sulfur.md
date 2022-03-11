---
title: "Sulfur workflow"
output:
  html_document: rmdformats::downcute
  toc: yes
  toc_float: yes
  pdf_document: default
---

# Generating table of sulfur genes that are present from IMG data

1. Subset all of the IMG tables for sulfur cycling genes
```
for i in *mapping_scaffolds.tsv; do for gene in `cat sulfur_genes/sulfur_genes.txt`; do grep -w $gene $i >> sulfur_genes/$i.sulfur; done; done
```

2. Sort for unique rows to remove duplicates grabbed with grep
```
for i in *.sulfur; do sort $i | uniq >> $i.uniq; done
```

3. Subset for bins only
```
for i in *.uniq; do grep -v "NoBin" $i >> $i.binsOnly; done
```

4. Cat all tables
```
cat *.binsOnly > all_SFsamples_sulfur_binsOnly.tsv
```

5. Add header info
```
echo -e "\tOriginal_Contig_Name\tIMG_Contig_Name\tLocus_Tag\tIMG_Gene_ID\tGene_Type\tGene_Start\tGene_Stop\tGene_Length\tHomolog_Gene_ID\tHomolog_Taxon_ID\tLineage_%ID\tLineage\tProduct\tSource\tCOG_ID\tCog_%ID\tPFAM_ID\tKO_Term\tKO_Term.1\tKO_%ID\tEC_Number\tBin\tGC\tLength\tDepth\tBin_scaffold_name\tAbun_bin" | cat - all_SFsamples_sulfur_binsOnly.tsv > all_SFsamples_sulfur_binsOnly.header.tsv
```

# Search for sulfur genes using custom HMMs from Karthik's lab

Justification: IMG sulfur output is limited e.g., no dsrC hits? 

1. Searching using Karthiks 79 sulfur HMMs:
```
for file in /storage1/data12/Plume_Viruses/Marinobacter_MAGs/sulfur_hmmsearch/sulfur_hmms/*.hmm; do hmmsearch --cut_tc --cpu 10 --tblout $file.txt $file 639_SFB_mags.faa; echo "next hmm"; done
```

2. Python script to parse hmm output for the best hits *per protein*
```
python3 /slowdata/scripts/python_scripts/hmm_parser.py -i dsrC.hmm.txt

for i in *.txt; do python3 /slowdata/scripts/python_scripts/hmm_parser.py -i $i; done
```

3. Cat the sorted output
```
cat *parsed-score.tsv > all_sulfur_hmmOutput_SFB.tsv
```



