# build annotation from Ensembl with version control


library(biomaRt)
library(GenomicFeatures)
#hg.gtf.db <- makeTxDbFromGFF("/path/to/Homo_sapiens.GRCh38.111.chr.gtf.gz", format="gtf" )
hg.gtf.db <- makeTxDbFromGFF("/path/to/Homo_sapiens.GRCh38.110.chr.gtf", format="gtf" )
ensembl.genes = genes(hg.gtf.db)
print(length(ensembl.genes))
#Homo_sapiens.GRCh38.111.chr.gtf.gz
#[1] 63187
#Homo_sapiens.GRCh38.110.chr.gtf
#[1] 62700

print(length(unique(ensembl.genes)))
#[1] 63183
#human = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL",host="asia.ensembl.org", dataset="hsapiens_gene_ensembl", version="79")
human = useEnsembl(biomart="ENSEMBL_MART_ENSEMBL",host="https://jul2023.archive.ensembl.org", dataset="hsapiens_gene_ensembl")#host = "mar2015.archive.ensembl.org"

gene_annot = getBM(attributes=c("ensembl_gene_id", "entrezgene_id", "gene_biotype","external_gene_name", "hgnc_symbol", "description"), mart=human, filters="ensembl_gene_id", values=ensembl.genes$gene_id, uniqueRows=TRUE)

ensembl.genes$external_gene_name = gene_annot$external_gene_name[ match(ensembl.genes$gene_id, gene_annot$ensembl_gene_id) ]  

ensembl.genes$hgnc_symbol = gene_annot$hgnc_symbol[ match(ensembl.genes$gene_id, gene_annot$ensembl_gene_id) ]
ensembl.genes$gene_biotype = gene_annot$gene_biotype[ match(ensembl.genes$gene_id, gene_annot$ensembl_gene_id) ]
ensembl.genes$status = gene_annot$status[ match(ensembl.genes$gene_id, gene_annot$ensembl_gene_id) ]
ensembl.genes$description = gene_annot$description[ match(ensembl.genes$gene_id, gene_annot$ensembl_gene_id) ]
ensembl.genes$entrezgene = gene_annot$entrezgene[ match(ensembl.genes$gene_id, gene_annot$ensembl_gene_id) ]

# saveRDS(ensembl.genes, file = "/Users/jing/nBox/Projects/bin/CommonFunctionR/ensembl_genes_gtf111.RDS")
saveRDS(ensembl.genes, file = "/Users/jing/nBox/Projects/bin/CommonFunctionR/ensembl_genes_gtf110.RDS")