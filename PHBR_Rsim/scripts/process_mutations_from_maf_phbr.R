#!/usr/bin/env Rscript

# === Libraries ===
library(Biostrings)
library(dplyr)
library(tidyr)
library(stringr)

# === Input arguments from Nextflow ===
args <- commandArgs(trailingOnly = TRUE)
input_file    <- args[1]
output_file   <- args[2]
sample_id     <- args[3]
mapping_file  <- args[4]
mane_file     <- args[5]
fasta_file    <- args[6]

message("🔧 Input file: ", input_file)
message("🔧 Output file: ", output_file)
message("🔧 Sample ID: ", sample_id)
message("🔧 Mapping file: ", mapping_file)
message("🔧 MANE file: ", mane_file)
message("🔧 FASTA file: ", fasta_file)

# === Load RefSeq ↔ Ensembl Mapping ===
message("📥 Loading RefSeq-Ensembl mapping...")
mapping <- read.delim(mapping_file, stringsAsFactors = FALSE)
mapping$refseq_mrna <- sub("\\.\\d+$", "", mapping$refseq_mrna)

# === Load MANE Canonical Transcripts ===
message("📥 Loading MANE Select transcripts...")
header_line <- readLines(mane_file, n = 1)
mane_df <- read.delim(mane_file, header = FALSE, sep = "\t", comment.char = "")
colnames(mane_df) <- strsplit(sub("^#", "", header_line), "\t")[[1]]
mane_df <- subset(mane_df, MANE_status == "MANE Select")
mane_df <- mane_df[, c("symbol", "RefSeq_nuc", "Ensembl_nuc")]
mane_df$RefSeq_nuc <- sub("\\.\\d+$", "", mane_df$RefSeq_nuc)
mane_df$Ensembl_nuc <- sub("\\.\\d+$", "", mane_df$Ensembl_nuc)
colnames(mane_df) <- c("gene", "refseq", "ensembl")
canonical_transcripts <- setNames(mane_df$ensembl, mane_df$gene)

# === Load Ensembl protein FASTA ===
message("📥 Reading Ensembl protein FASTA...")
ensembl_proteins <- readAAStringSet(fasta_file)
transcript_ids <- sub(".*transcript:([A-Z0-9]+)(\\.\\d+)?\\b.*", "\\1", names(ensembl_proteins))
protein_lookup <- setNames(as.character(ensembl_proteins), transcript_ids)
missing_transcripts <- c()

get_protein_seq <- function(ensembl_transcript_id) {
  ensembl_transcript_id <- sub("\\.\\d+$", "", ensembl_transcript_id)
  if (is.na(ensembl_transcript_id) || ensembl_transcript_id == "") return(NULL)
  if (!(ensembl_transcript_id %in% names(protein_lookup))) {
    missing_transcripts <<- unique(c(missing_transcripts, ensembl_transcript_id))
    return(NULL)
  }
  return(protein_lookup[[ensembl_transcript_id]])
}

get_longest_isoform <- function(gene) {
  candidates <- ensembl_proteins[grepl(paste0("gene_symbol:", gene), names(ensembl_proteins))]
  if (length(candidates) == 0) return(NULL)
  longest <- which.max(width(candidates))
  return(as.character(candidates[longest]))
}

extract_driver_mutations <- function(maf_df) {
  if (!"all_effects" %in% colnames(maf_df)) return(data.frame())
  if (all(is.na(maf_df$all_effects))) return(data.frame())
  if (!is.character(maf_df$all_effects)) return(data.frame())

  all_entries <- unlist(strsplit(maf_df$all_effects[!is.na(maf_df$all_effects)], ";"))

  parsed <- lapply(all_entries, function(entry) {
    fields <- strsplit(entry, ",")[[1]]
    if (length(fields) >= 4 && grepl("^p\\.[A-Z]\\d+[A-Z]$", fields[3])) {
      data.frame(
        Hugo_Symbol = fields[1],
        Variant_Classification = fields[2],
        HGVSp_Short = fields[3],
        Transcript_ID = ifelse(fields[1] %in% names(canonical_transcripts),
                       canonical_transcripts[[fields[1]]],
                       fields[4]),
        stringsAsFactors = FALSE
      )
    } else {
      NULL
    }
  })

  do.call(rbind, parsed)
}


apply_aa_change <- function(seq, pos, new_aa) {
  if (pos > nchar(seq)) return(NULL)
  substr(seq, pos, pos) <- new_aa
  return(seq)
}

generate_peptides <- function(mutated_seq, position, to_aa, from_aa, gene, mut_id) {
  n <- nchar(mutated_seq)

  pep_mhc1 <- NA
  pep_mhc2 <- NA

  # MHC-I: 21mer, center at 11th (position)
  if (position > 10 && position + 10 <= n) {
    pep <- substr(mutated_seq, position - 10, position + 10)
    if (substr(pep, 11, 11) != to_aa) {
      message("⚠️ MHC-I peptide center != mutated AA for ", mut_id,
              " | expected: ", to_aa, ", got: ", substr(pep, 11, 11))
    }
    pep_mhc1 <- pep
  }

  # MHC-II: 29mer, center at 15th
  if (position > 14 && position + 14 <= n) {
    pep <- substr(mutated_seq, position - 14, position + 14)
    if (substr(pep, 15, 15) != to_aa) {
      message("⚠️ MHC-II peptide center != mutated AA for ", mut_id,
              " | expected: ", to_aa, ", got: ", substr(pep, 15, 15))
    }
    pep_mhc2 <- pep
  }

  return(list(Peptide_MHC1 = pep_mhc1, Peptide_MHC2 = pep_mhc2))
}



hotspot_enst <- list(
  "BRAF"   = "ENST00000288602",     # use .11 if your FASTA includes version
  "IDH1"   = "ENST00000345146",
  "PIK3CA" = "ENST00000263967",
  "KRAS"   = "ENST00000256078",
  "TP53"   = "ENST00000269305",
  "FBXW7"  = "ENST00000357628"
)

hotspot_aachanges <- c(
  "BRAF_p.V600E", "IDH1_p.R132H", "PIK3CA_p.E545K", "PIK3CA_p.E542K",
  "KRAS_p.G12D", "KRAS_p.G12C", "KRAS_p.G12V", "KRAS_p.G13D",
  "TP53_p.R175H", "TP53_p.R282W", "TP53_p.Y220C", "TP53_p.H179R",
  "FBXW7_p.R465C"
)


# === Load MAF ===
message("📥 Reading MAF file...")
maf <- read.delim(gzfile(input_file), comment.char = "#", stringsAsFactors = FALSE)

missense <- maf %>%
  filter(Variant_Classification == "Missense_Mutation",
         !is.na(HGVSp_Short), !is.na(Transcript_ID)) %>%
  distinct(Hugo_Symbol, HGVSp_Short, Transcript_ID)

message("🔍 Found ", nrow(missense), " initial missense mutations")

# === Filter for hotspot mutations only (based on mut_matrix column names) ===
hotspot_set <- c("BRAF_V600E", "IDH1_R132H", "PIK3CA_E545K", "PIK3CA_H1047R", 
"KRAS_G12D", "KRAS_G12V", "PIK3CA_E542K", "TP53_R175H", "TP53_R273C", 
"TP53_R248Q", "NRAS_Q61R", "KRAS_G12C", "TP53_R273H", "TP53_R248W", 
"TP53_R282W", "KRAS_G13D", "PIK3CA_R88Q", "TP53_Y220C", "NRAS_Q61K", 
"AKT1_E17K", "PTEN_R130Q", "IDH1_R132C", "TP53_G245S", "FGFR3_S249C", 
"ERBB2_S310F", "PTEN_R130G", "KRAS_G12A", "TP53_H179R", "BRAF_V600M", 
"FBXW7_R465C", "PIK3CA_N345K", "PIK3CA_H1047L", "GNA11_Q209L", 
"KRAS_G12R", "TP53_H193R", "TP53_V157F", "BCOR_N1459S", "HRAS_Q61R", 
"PPP2R1A_P179R", "FBXW7_R465H", "FBXW7_R505G", "GNAQ_Q209P", 
"PIK3CA_E726K", "PIK3CA_G118D", "KRAS_G12S", "TP53_R249S", "TP53_Y163C", 
"POLE_P286R", "TP53_C176F", "TP53_R273L", "U2AF1_S34F", "ERBB3_V104M", 
"FGFR2_S252W", "PIK3CA_C420R", "TP53_M237I", "CTNNB1_S37F", "TP53_R158L", 
"PIK3CA_Q546R", "TP53_E285K", "CTNNB1_T41A", "NRAS_Q61L", "SMAD4_R361H", 
"TP53_S241F", "EGFR_A289V", "EGFR_L858R", "KRAS_Q61H", "MB21D2_Q311E", 
"KRAS_A146T", "PIK3CA_E545A", "TP53_I195T", "TP53_Y205C", "TP53_Y234C", 
"EGFR_G598V", "NFE2L2_R34G", "PIK3CA_E453K", "RAC1_P29S", "TP53_R280T", 
"TP53_V272M", "CTNNB1_S45P", "FBXW7_R479Q", "IDH1_R132G", "POLE_V411L", 
"TP53_C176Y", "TP53_K132N", "TP53_R158H", "NFE2L2_D29H", "PIK3CA_M1043I", 
"PIK3CA_R108H", "PPP2R1A_R183W", "TP53_G245D", "TP53_G245V", 
"TP53_H214R", "ATM_R337C", "CTNNB1_S37C", "GNAS_R844C", "MAPK1_E322K", 
"TP53_C275Y", "FBXW7_R505C", "GNAS_R844H", "TP53_C238Y", "TP53_E286K", 
"TP53_R280K", "TP53_Y236C", "CTNNB1_S33C", "ERBB2_R678Q", "ERBB2_V842I", 
"FGFR2_N550K", "IDH2_R172K", "PCBP1_L100Q", "PIK3CA_Q546K", "SF3B1_K700E", 
"TP53_A159V", "TP53_G266V", "TP53_H179Y", "TP53_L194R", "TP53_P151S", 
"TP53_R337C", "CDKN2A_H83Y", "CTNNB1_S33F", "FBXW7_R689W", "KRAS_G13C", 
"NFE2L2_E79Q", "NFE2L2_E82D", "TP53_C238F", "TP53_P278S", "TP53_V173M", 
"TP53_V216M", "CTNNB1_S45F", "EP300_D1399N", "MAP2K1_P124S", 
"TP53_C141Y", "TP53_E271K", "TP53_H193Y", "TP53_R110L", "TP53_V173L", 
"CDKN2A_P114L", "CNTNAP2_R483Q", "CTNNB1_G34R", "CTNNB1_T41I", 
"ERBB2_L755S", "HRAS_G13V", "MAX_R60Q", "NRAS_G12D", "PIK3CA_H1047Y", 
"PIK3CA_R93Q", "PIK3R1_G376R", "PTEN_R173C", "SF3B1_R625H", "SMAD4_R361C", 
"BRAF_K601E", "CHD4_R975H", "CTNNB1_D32G", "CTNNB1_G34E", "MTOR_E1799K", 
"MUC16_P5119S", "PIK3CA_E81K", "PIK3CA_M1043V", "PIK3CA_Y1021C", 
"PPP6C_R301C", "PTEN_R142W", "RAF1_S257L", "RGS7_R44C", "STAG1_A116T", 
"TP53_A159P", "TP53_C135Y", "TP53_C242F", "TP53_C275F", "TP53_G245C", 
"TP53_G266R", "TP53_H193L", "TP53_P250L", "CIC_R215W", "CTNNB1_S33Y", 
"HRAS_Q61K", "NBEA_E1710K", "PIK3CA_K111E", "PIK3CA_Q546P", "PIK3CA_R38C", 
"PIK3CA_R38H", "PIK3R1_N564D", "TP53_G266E", "TP53_P152L", "TP53_R248L", 
"TP53_R337L", "TP53_S241C", "XPO1_E571K", "BRAF_G469A", "CASP8_R127Q", 
"CDKN2A_D108Y", "CNTNAP2_G362E", "CTCF_R377C", "CTNNB1_D32N", 
"ERBB2_V777L", "FGFR3_Y375C", "GNAQ_Q209L", "HRAS_G13R", "IDH1_R132S", 
"ISX_R86C", "KRAS_Q61L", "MTOR_S2215Y", "NFE2L2_E79K", "PIK3CA_E545G", 
"PIK3CA_K111N", "SPOP_M117V", "TP53_P278A", "TP63_R379C", "TRRAP_S722F", 
"ARID2_S297F", "BCLAF1_A306T", "BMP5_R321Q", "CDKN2A_D84N", "CHD4_R1105W", 
"CHD4_R1338I", "CNBD1_L396R", "CTNNB1_D32Y", "CTNNB1_S33P", "DCAF12L2_P334L", 
"ERCC2_N238S", "FUBP1_R430C", "HRAS_G12S", "HRAS_Q61L", "KDR_R1032Q", 
"KNSTRN_S24F", "MYC_S161L", "MYCN_P44L", "NRAS_Q61H", "PIK3CA_G106V", 
"PIK3CA_R93W", "RUNX1T1_E249K", "RUNX1T1_R254Q", "SF3B1_R625C", 
"SMARCA4_R1192C", "SPOP_F133L", "TP53_G244S", "TP53_I195F", "TP53_N239D", 
"TP53_R175G", "TP53_R213Q", "TP53_R249G", "TP53_R283P", "TP53_S127F", 
"BCLAF1_R283Q", "BRAF_D594N", "BRAF_G466V", "BRAF_N581S", "CASP8_R292W", 
"CREBBP_R1446C", "CTNNB1_G34V", "CUX1_E47K", "EGFR_A289T", "FBXW7_G423V", 
"FGFR2_C383R", "KRAS_Q61K", "KRAS_Q61R", "NFE2L2_R34P", "NRAS_G13R", 
"PDE4DIP_R634Q", "PIK3CA_D350N", "PIK3CA_E365K", "PIK3CA_E545Q", 
"PIK3CA_G1007R", "PIK3CA_T1025A", "POLQ_R860Q", "PPP2R1A_S256F", 
"PREX2_E1295K", "PTEN_C136R", "PTEN_D92E", "PTEN_G132D", "PTEN_R130L", 
"PTPRD_S431L", "RGPD3_E755K", "RSPO2_R64Q", "RUNX1T1_D246N", 
"SMARCA4_T910M", "SPOP_W131G", "STIL_S76L", "TNC_S1725L", "TP53_C135F", 
"TP53_D281Y", "TP53_G105C", "TP53_G154V", "TP53_G244D", "TP53_G262V", 
"TP53_P151H", "TP53_P190L", "TP53_Q136E", "TP53_S241Y", "TPR_S2155L", 
"XPO1_R749Q", "ZNF479_R485I", "ARID2_R314C", "ATM_R337H", "ATRX_R781Q", 
"BCL6_R594Q", "BCLAF1_E163K", "CDC73_R91Q", "CNTNAP2_E149K", 
"CRNKL1_S128F", "CSMD3_R3127Q", "CTNNB1_D32V", "DCC_R443Q", "DICER1_R944Q", 
"EGFR_L62R", "EGFR_R222C", "EP300_R580Q", "ERBB4_R711C", "ETV1_R187C", 
"FBXW7_R658Q", "FBXW7_Y545C", "GATA3_A396T", "GRIN2A_G1322E", 
"GRIN2A_R1022C", "GRIN2A_R1067W", "HRAS_G12D", "IDH2_R140Q", 
"KRAS_K117N", "KRAS_L19F", "LRP1B_R2443C", "LRP1B_R295Q", "MAP3K1_S1330L", 
"NFE2L2_G31A", "NFE2L2_R34Q", "NT5C2_E390K", "P2RY8_R234H", "PIK3CA_G1049R", 
"PIK3CA_V344M", "PREX2_E414K", "PTEN_F341V", "PTEN_R173H", "PTPN11_G503V", 
"RABEP1_R75Q", "RAC1_A178V", "RHOA_E47K", "RPL5_E82K", "SMARCA4_R1192H", 
"SPOP_F133V", "TP53_A161T", "TP53_G199V", "TP53_G279E", "TP53_H193P", 
"TP53_K132E", "TP53_K132R", "TP53_P278L", "TP53_P278R", "TP53_R175C", 
"TP53_R273S", "TP53_T155P", "TP53_V274F", "TP53_Y163H", "TP63_E609K", 
"ZNF521_F513L", "A1CF_E42K", "ABI1_K445N", "ARHGEF10_A459V", 
"ARID1A_G2087R", "BAZ1A_R1480C", "BMP5_R449C", "BRAF_G469V", 
"CACNA1D_E2107K", "CACNA1D_R116Q", "CACNA1D_R462Q", "CCND1_T286I", 
"CDKN2A_P48L", "CHD4_R1162W", "CIC_R201W", "CIITA_E727K", "CNBD1_L396P", 
"CNOT3_E20K", "CNTNAP2_R389W", "CREBBP_R1446H", "CRLF2_S128L", 
"CSMD3_R100Q", "CSMD3_R334Q", "CSMD3_S1090Y", "CTCF_S354F", "CTNNA2_S891L", 
"CTNNB1_H36P", "CTNNB1_K335I", "CTNNB1_R587Q", "CTNNB1_S37A", 
"CTNNB1_S45Y", "CTNND2_V428I", "DICER1_D1709N", "ECT2L_R685Q", 
"EGFR_L861Q", "EGFR_R108K", "EGFR_R252C", "EGFR_V774M", "ERBB3_D297Y", 
"ERBB3_V104L", "ESR1_E247K", "EZH2_E745K", "FAM135B_S11L", "FAT1_E4454K", 
"FAT3_E299K", "FAT3_S2675L", "FGFR1_N577K", "FOXA1_I176M", "GRIN2A_E962K", 
"HIST1H3B_E74K", "HRAS_K117N", "JAK3_R445Q", "KEAP1_R470C", "KEAP1_V271L", 
"KIT_D816V", "KIT_K642E", "LIFR_E319K", "LIFR_E391K", "LRP1B_L3714F", 
"LZTR1_G248R", "MAX_H28R", "MET_R1166Q", "MUC16_S10659L", "MYH9_E530K", 
"NFE2L2_D29N", "NFE2L2_L30F", "NTRK1_G169R", "NUMA1_A1102T", 
"PAX3_T424M", "PAX7_S155L", "PDE4DIP_E1147K", "PIK3CA_D350G", 
"PIK3CA_E453Q", "PIK3CA_E542A", "PIK3CA_E600K", "PIK3CA_V344G", 
"PIK3R1_K567E", "POLG_R386C", "PREX2_H1404Y", "PTEN_P246L", "PTPRB_D1778N", 
"PTPRB_R570Q", "PTPRT_D1266N", "PTPRT_E324K", "PTPRT_M881I", 
"RHOA_E40Q", "ROBO2_R261C", "SF3B1_E902K", "SFRP4_R232Q", "SMAD2_R120Q", 
"SMARCA4_G1232S", "SMARCA4_R885C", "TCF7L2_R471C", "TGFBR2_D471N", 
"TGFBR2_R553C", "TP53_D259V", "TP53_D281E", "TP53_E286G", "TP53_F270S", 
"TP53_G244C", "TP53_G334V", "TP53_L130V", "TP53_M246I", "TP53_Q331H", 
"TP53_R156P", "TP53_R249M", "TP53_R267W", "TP53_R280G", "TP53_S127Y", 
"TP53_S215I", "ZMYM2_E64K", "A1CF_A203V", "ACVR1_R206H", "ACVR1_R375C", 
"AFF3_R902W", "AKAP9_R1947C", "AKT3_R367Q", "ASXL1_A611T", "BCL11B_D885N", 
"BCL11B_E123K", "BCL2L12_R18W", "BCLAF1_G390E", "BCLAF1_R544H", 
"BRAF_G466E", "BRCA2_A1393V", "BRD3_A373T", "CACNA1D_R1594W", 
"CACNA1D_R533H", "CACNA1D_S1755L", "CAMTA1_S287L", "CBL_E276K", 
"CBL_R206Q", "CCNC_R25H", "CDC73_R139Q", "CDH1_D254Y", "CDH10_R638Q", 
"CDK12_R890H", "CHD4_R1068C", "CHD4_R877W", "CIC_R1512H", "CIC_R1515C", 
"CIC_R202W", "CNTNAP2_P137S", "CSMD3_P258S", "CSMD3_R3025C", 
"CSMD3_R683C", "CTNNA2_R449W", "CTNNB1_D32H", "CTNND2_R1063C", 
"DCAF12L2_A430V", "DCAF12L2_R337H", "DCAF12L2_R378H", "DCC_T1255M", 
"DICER1_S1344L", "EGFR_A289D", "EGFR_E114K", "EGFR_S768I", "ELF3_E262Q", 
"EP300_A1629V", "ERBB2_S310Y", "ERBB3_A232V", "ERBB3_E332K", 
"ERBB3_E928G", "ERBB3_R475W", "ERG_V162I", "ETV5_A192V", "FAM135B_R1211Q", 
"FAM135B_R719Q", "FAM135B_S645R", "FAT3_A3348V", "FAT3_E4502K", 
"FAT3_R258C", "FAT3_R567Q", "FAT4_D1790N", "FAT4_H2514Y", "FAT4_R1671C", 
"FAT4_R1815C", "FAT4_R2685Q", "FBLN2_A1156T", "FBXW7_R14Q", "FBXW7_R505H", 
"FBXW7_S582L", "FGFR3_G372C", "FGFR3_G382R", "FGFR3_R248C", "FNBP1_R293C", 
"FOXP1_S450L", "FUS_D46H", "GATA3_M294K", "GLI1_R380Q", "GOLGA5_R674C", 
"GRIN2A_R1159H", "GRIN2A_S1425L", "GRM3_E49K", "GRM3_E724K", 
"GRM3_M518I", "HIST1H3B_E134Q", "ITK_E196K", "KAT6A_E1408K", 
"KAT7_R127Q", "KDM5A_R719C", "KDR_S1021L", "KIAA1549_R1275Q", 
"KIT_E257D", "KIT_N822K", "KMT2C_F357L", "KMT2C_R4693Q", "LCP1_R488H", 
"LMO2_A193T", "LRIG3_R224Q", "LSM14A_R272C", "LZTR1_F243L", "MAP2K4_R134W", 
"MAP2K4_S184L", "MAP3K1_R306H", "MECOM_S237L", "MECOM_S594L", 
"MECOM_S888L", "MED12_D23Y", "MET_M1268T", "MLLT11_E86K", "MTOR_R2152C", 
"MUC16_A2118T", "MUC16_E5116K", "MUC16_R4603Q", "MUC16_R8606H", 
"MUC16_S2253F", "MUC16_S5646F", "MUC16_S6608L", "MYC_P74L", "MYD88_L273P", 
"NBEA_R2083C", "NBEA_V2417I", "NBN_P198Q", "NCOA4_R578Q", "NCOR2_P20Q", 
"NF1_R1870Q", "NF2_R418C", "NFATC2_S808L", "NFE2L2_D29Y", "NFE2L2_G81S", 
"NFE2L2_W24C", "NOTCH1_A465T", "NRAS_G12C", "NUTM1_S116L", "PABPC1_R518C", 
"PBRM1_R876C", "PDE4DIP_R430Q", "PIK3CA_C378F", "PIK3CA_C901F", 
"PIK3CA_E418K", "PIK3CA_E545D", "PIK3CA_E970K", "PIK3CA_M1004I", 
"PIK3CA_N1044K", "PIK3CA_R357Q", "PIK3CB_D1067V", "PIK3CB_E1051K", 
"PIK3CB_E470K", "PIK3CB_R321Q", "POLE_R1826W", "POLE_S459F", 
"PPP2R1A_R183Q", "PPP2R1A_S256Y", "PREX2_D269N", "PREX2_E553K", 
"PRKCB_S660F", "PTEN_G129E", "PTEN_G36R", "PTEN_R14M", "PTPN11_F71L", 
"PTPN11_T468M", "PTPN13_S887L", "PTPRB_R747Q", "PTPRC_R1115Q", 
"PTPRD_L1053I", "PTPRD_R561Q", "PTPRK_E358K", "PTPRT_E548K", 
"PTPRT_L695F", "PTPRT_R328C", "RAC1_P29L", "REL_S440L", "RGPD3_R1394C", 
"RGS7_R273W", "RNF43_S607L", "ROBO2_D1030N", "ROS1_R1311Q", "RUNX1T1_A211T", 
"RUNX1T1_H177Y", "RUNX1T1_R405W", "RUNX1T1_R484Q", "SETD2_R2109Q", 
"SF3B1_G740E", "SF3B1_K666T", "SF3B1_R957Q", "SFPQ_A672V", "SMAD2_R321Q", 
"SMAD2_S276L", "SPEN_R637Q", "SSX4_R88Q", "STAG2_V465F", "TCEA1_R153Q", 
"TCF7L2_R420W", "TET1_R1158W", "TFRC_R174C", "TP53_C141W", "TP53_C242S", 
"TP53_C242Y", "TP53_C277F", "TP53_D259Y", "TP53_E224D", "TP53_F109C", 
"TP53_F134L", "TP53_H179D", "TP53_L130F", "TP53_L265P", "TP53_N239S", 
"TP53_P151T", "TP53_P278H", "TP53_P278T", "TP53_R213L", "TP53_R267P", 
"TP53_R273P", "TP53_R280I", "TP53_S106R", "TP53_S215G", "TP53_T125M", 
"TP53_T211I", "TP53_Y126C", "TPR_R1698C", "VAV1_A263T", "ZFHX3_R1439Q", 
"ZMYM2_R435C", "ZNF429_R221I", "ZNF479_R261I", "ZNF479_R429I", 
"ZNF521_D889N")

# Add Simple_Mutation column for matching
# === Flag hotspot mutations ===
missense <- missense %>%
  mutate(Simple_Mutation = paste0(Hugo_Symbol, "_", sub("p\\.", "", HGVSp_Short)),
         is_hotspot = Simple_Mutation %in% hotspot_set)



# === Add driver mutations from all_effects ===
driver_muts <- extract_driver_mutations(maf)


if (!is.null(driver_muts) && nrow(driver_muts) > 0) {
  message("➕ Parsed ", nrow(driver_muts), " mutations from all_effects column")

  driver_muts <- driver_muts %>%
    mutate(Simple_Mutation = paste0(Hugo_Symbol, "_", sub("p\\.", "", HGVSp_Short)),
          is_hotspot = Simple_Mutation %in% hotspot_set)


  # Keep only hotspot mutations from the annotation field
  hotspot_from_all_effects <- driver_muts %>%
    filter(Simple_Mutation %in% hotspot_set) %>%
    rowwise() %>%
    mutate(Transcript_ID = ifelse(
      Hugo_Symbol %in% names(canonical_transcripts),
      canonical_transcripts[[Hugo_Symbol]],
      Transcript_ID
    )) %>%
    ungroup()

  if (nrow(hotspot_from_all_effects) > 0) {
    message("🔥 Found ", nrow(hotspot_from_all_effects), " hotspot mutations in all_effects")

    missense <- bind_rows(missense, hotspot_from_all_effects) %>%
      distinct(Hugo_Symbol, HGVSp_Short, Transcript_ID)
  }
}


message("🔎 Total unique missense/driver mutations: ", nrow(missense))

# === Process mutations ===
# === Process mutations ===
all_peptides <- list()

for (i in seq_len(nrow(missense))) {
  gene         <- missense$Hugo_Symbol[i]
  aa_change    <- missense$HGVSp_Short[i]
  ensembl_id   <- missense$Transcript_ID[i]
  ensembl_id_clean <- sub("\\.\\d+$", "", ensembl_id)
  mut_id       <- paste0(gene, "_", sub("p\\.", "", aa_change))

  message("➡️ Processing ", gene, " | ", aa_change, " | ", ensembl_id_clean)

  # Parse mutation string (e.g., p.V600E)
  match <- regexec("p\\.([A-Z*]+)(\\d+)([A-Z*]+)", aa_change)
  parsed <- regmatches(aa_change, match)[[1]]

  if (length(parsed) != 4) {
    message("❌ Could not parse mutation: ", aa_change)
    next
  }

  from_aa <- parsed[2]
  pos     <- as.integer(parsed[3])
  to_aa   <- parsed[4]

  # --- Get protein sequence ---
  prot_seq <- get_protein_seq(ensembl_id_clean)
  if (is.null(prot_seq)) {
    message("⚠️ Transcript not found: ", ensembl_id_clean, " — trying fallback")
    prot_seq <- get_longest_isoform(gene)
    ensembl_id_clean <- "LONGEST_ISOFORM"
  }

  if (is.null(prot_seq)) {
    message("❌ No sequence available for ", gene, " — skipping")
    next
  }

  # --- Confirm reference AA matches ---
  actual_aa <- substr(prot_seq, pos, pos)
  if (actual_aa != from_aa) {
    message("⚠️ Ref AA mismatch for ", ensembl_id_clean, ": expected ", from_aa, ", found ", actual_aa,
            " — trying longest isoform for ", gene)
    prot_seq <- get_longest_isoform(gene)
    ensembl_id_clean <- "LONGEST_ISOFORM"
    if (is.null(prot_seq) || substr(prot_seq, pos, pos) != from_aa) {
      message("❌ Still no match — skipping ", mut_id)
      next
    }
  }


  # --- Apply mutation ---
  mut_seq <- apply_aa_change(prot_seq, pos, to_aa)
  if (substr(mut_seq, pos, pos) != to_aa) {
    message("❌ Mutation failed to apply at pos ", pos, " for ", mut_id)
    next
  }

  message("✅ Mutation applied correctly: ", from_aa, " → ", to_aa, " at pos ", pos)
  message("📏 Protein length: ", nchar(mut_seq))

  # --- Peptide generation ---
  peptides <- generate_peptides(mut_seq, pos, to_aa, from_aa, gene, mut_id)

  message("📍 MHC-I window: ", pos - 10, "-", pos + 10, " | Peptide: ", peptides$Peptide_MHC1)
  message("📍 MHC-II window: ", pos - 14, "-", pos + 14, " | Peptide: ", peptides$Peptide_MHC2)

  # --- Save output ---
# --- Save output ---
  all_peptides <- append(all_peptides, list(data.frame(
    Gene = gene,
    Ensembl_Transcript = ensembl_id_clean,
    AAChange = aa_change,
    Simple_Mutation = mut_id,
    Peptide_MHC1 = peptides$Peptide_MHC1,
    Peptide_MHC2 = peptides$Peptide_MHC2,
    AA_Pos = pos,
    Ref_AA = from_aa,
    Alt_AA = to_aa,
    is_hotspot = mut_id %in% hotspot_set  
  )))

}


# === Output ===
if (length(all_peptides) > 0) {
  peptides_df <- bind_rows(all_peptides) %>%
    filter(!is.na(Peptide_MHC1) | !is.na(Peptide_MHC2))

  if (nrow(peptides_df) > 0) {
    peptides_df$sample_id <- sample_id
    write.table(peptides_df, file = output_file, sep = "\t", row.names = FALSE, quote = FALSE)
    message("✅ Peptides written for sample: ", sample_id, " (", nrow(peptides_df), " entries)")
  } else {
    message("⚠️ No valid peptides found (after NA filtering) — writing empty file.")
    file.create(output_file)
  }

} else {
  message("⚠️ No peptides found for sample: ", sample_id, " — writing empty file.")
  file.create(output_file)
}


# === Log missing transcripts ===
if (length(missing_transcripts) > 0) {
  writeLines(unique(missing_transcripts), "missing_transcripts.txt")
  message("📝 Logged missing transcripts to missing_transcripts.txt")
}
