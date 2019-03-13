## ============================================== ##
## 6. SPir Germline :: CGP Script for Gene Panel
##
##		6.1 Find germline mutations in top 50 genes using CGP script
##		6.2 Filter noncoding germline mutations. Promoter regions for these genes are annotated separetely 
##		6.3 Annotate these germline mutations using annovar
##		6.4 Filter based on Depth + known polymorphic ones using EXac
##		6.5 Summarise the germline mutations
##
## ============================================== ##

	## ------- ##
	## 6.1 Find germline mutations in top 50 genes using CGP script	
	## ------- ##

		rec_genes_50 = c("CYLD","COL11A1","DNMT3A","SPTA1","DLG2","TDRD9","TPTE","UGGT2","IFT80","TTN","TP53","CD1D","CSMD1","AKT1","DNAH17","ALPK1","KMT2E","ARHGAP21","CUX1","MAGEC3","BTBD19","EIF2B3","BRDT","F5","TNN","TNR","HMCN1","NLRP3","USP54","CCSER2","EML3","EFCAB4B","SLCO1C1","BAZ2A","NXPH4","PTPRR","STARD13","FAM81A","NEO1","E4F1","ZNF267","GLG1","POLR2A","MYH8","ITGA2B","SMG8","TANC2","LOXHD1","ZNF701","PPP1R12C","BLM")

		library(biomaRt)
		ens37mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host = "grch37.ensembl.org", path = "/biomart/martservice", dataset = "hsapiens_gene_ensembl") 
		filters <- listFilters(ens37mart)
		attributes <- listAttributes(ens37mart)
		rec_genes.coords37 <- getBM(attributes=c( "hgnc_symbol",  "start_position", "end_position", "chromosome_name"), filters  = 'hgnc_symbol',  values = rec_genes_50, mart = ens37mart) 
		#rec_genes.coords37_BLM <- getBM(attributes=c( "hgnc_symbol",  "start_position", "end_position", "chromosome_name"), filters  = 'hgnc_symbol',  values = "BLM", mart = ens37mart) 
		rec_genes.coords37 = rec_genes.coords37[ -grep("HG",rec_genes.coords37$chromosome_name), ]
		gene_locs = paste(rec_genes.coords37[,4], ":", rec_genes.coords37[,2], "-", rec_genes.coords37[,3], sep = "" )
		gene_locs_parameter = paste( paste("-r ", gene_locs, sep = ""), collapse = " " )
		

		command = paste("perl /software/CGP/projects/vcfutil/perl/scripts/GermlineSummary.pl -c /nfs/cancer_ref02/human/GRCh37d5/vagrent/e75/vagrent.cache.gz -p 1524 -a caveman -a pindel ", gene_locs_parameter ," -o /lustre/scratch116/casm/team113/projects/4205_Spiradenocarcinoma/analysis/germline/CGP_Target_Gene/CGP_Target_Gene.csv", sep = "" )
		system( command )

		bsub -e CGP_Germline_top_50_RecGenes.e -o CGP_Germline_top_50_RecGenes.o -R "select[type==X86_64 && mem > 4000] rusage[mem=4000]" -M4000 'perl /software/CGP/projects/vcfutil/perl/scripts/GermlineSummary.pl -c /nfs/cancer_ref02/human/GRCh37d5/vagrent/e75/vagrent.cache.gz -p 1524 -a caveman -r 14:105235686-105262088 -r 4:113206665-113363776 -r 10:24872538-25012597 -r 12:56989380-57030600 -r 1:92414928-92479983 -r 1:45274154-45281257 -r 10:86088342-86278273 -r 1:158149737-158154686 -r 1:103342023-103574052 -r 8:2792875-4852494 -r 7:101458959-101927249 -r 16:50775961-50835846 -r 11:83166055-85338966 -r 17:76419778-76573476 -r 2:25455845-25565459 -r 16:2273567-2285743 -r 12:3715799-3873985 -r 1:45316450-45452282 -r 11:62369690-62380237 -r 1:169483404-169555826 -r 15:59664892-59815748 -r 16:74485856-74641012 -r 1:185703683-186160085 -r 3:159974774-160117668 -r 17:42449548-42466873 -r 7:104654626-104754808 -r 18:44056935-44236996 -r X:140926102-140985618 -r 17:10293639-10325267 -r 15:73344051-73597547 -r 1:247579458-247612410 -r 12:57610578-57620232 -r 17:7387685-7417933 -r 19:55602281-55628927 -r 12:71031853-71314623 -r 12:20848289-20906320 -r 17:57286761-57292608 -r 1:158580496-158656488 -r 13:33677272-33924767 -r 17:61086917-61505060 -r 14:104394799-104519004 -r 1:175036994-175117202 -r 1:175284330-175712906 -r 17:7565097-7590856 -r 21:10906201-11029719 -r 2:179390716-179695529 -r 13:96453834-96705736 -r 10:75257296-75385711 -r 16:31885079-31928668 -r 19:53059075-53090427 -r 15:91260558-91358859 -o /lustre/scratch116/casm/team113/projects/4205_Spiradenocarcinoma/analysis/germline/CGP_Target_Gene/CGP_Target_50_Rec_Genes_and_BLM.csv'
		
		
		/lustre/scratch116/casm/team113/projects/4205_Spiradenocarcinoma/analysis/germline/CGP_Target_Gene/CGP_Target_50_Rec_Genes_and_BLM.csv

		command = paste("perl /software/CGP/projects/vcfutil/perl/scripts/GermlineSummary.pl -c /nfs/cancer_ref02/human/GRCh37d5/vagrent/e75/vagrent.cache.gz -p 1524 -a caveman -a pindel ", gene_locs_parameter ," -o /lustre/scratch116/casm/team113/projects/4205_Spiradenocarcinoma/analysis/germline/CGP_Target_Gene/CGP_Target_Gene.csv", sep = "" )



		## Read the CGP target germline file 
		require(data.table)
		setwd("/lustre/scratch116/casm/team113/projects/4205_Spiradenocarcinoma/analysis/germline/CGP_Target_Gene")
		target_gene_germline = as.data.frame( fread("CGP_Target_50_Rec_Genes_and_BLM.csv", sep = "\t", header = T, stringsAsFactors = F, check.names = F , skip = 48 ) );


		interesting_effects = c( "missense","silent","ess_splice","nonsense","start_lost" );
		target_gene_germline_exonic = target_gene_germline[ (target_gene_germline$Effect %in% interesting_effects), ]
