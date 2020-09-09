#written by Noah Friedman
import pandas as pd
import numpy as np
import sys
import re
import os
from collections import Counter

#RETURN DIFFERENT GENE SETS
#GIVES different sets of genes
def get_im3_genes():
	im3Genes = set(['ABL1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXP1', 'FUBP1', 'GATA1', 'GATA2', 'GATA3', 'GNA11', 'GNAQ', 'GNAS', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3B', 'HNF1A', 'HRAS', 'ICOSLG', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAPK1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOR1', 'NF1', 'NF2', 'NFE2L2', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLE', 'PPP2R1A', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAC1', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'STAG2', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'YAP1', 'YES1'])
	return im3Genes

def get_im5_genes():
	return set(['ABL1', 'ACVR1', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'EP300', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'MLL', 'MLL2', 'MLL3', 'MPL', 'MRE11A', 'MSH2', 'MSH6', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPM1D', 'PPP2R1A', 'PPP6C', 'PRDM1', 'PRKAR1A', 'PTCH1', 'PTEN', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAD21', 'RAD50', 'RAD51', 'RAD51B', 'RAD51C', 'RAD51D', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RUNX1', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHQ1', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SOCS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'VEGFA', 'VHL', 'VTCN1', 'WT1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])

def get_im6_genes():
	im6Genes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'KMT2A', 'KMT2B', 'KMT2C', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
	return im6Genes

def get_tsgs():
	tsgs = set(['ERRFI1', 'ASXL2', 'PMAIP1', 'ACTG1', 'SUFU', 'FBXO11', 'MEN1', 'FAM58A', 'B2M', 'RB1', 'DUSP22', 'SESN1', 'GPS2', 'RAD51D', 'SMG1', 'CDC73', 'MAP3K1', 'SMARCB1', 'INPP4B', 'PARK2', 'SMAD4', 'CBFB', 'CDH1', 'PPP6C', 'SETDB1', 'SETDB2', 'NF2', 'CDKN2B', 'CDKN2C', 'CDKN2A', 'DDX3X', 'PIK3R1', 'BARD1', 'PDS5B', 'KLF4', 'SPRED1', 'VHL', 'SMAD2', 'PMS1', 'PMS2', 'SETD2', 'GATA3', 'TBL1XR1', 'MUTYH', 'SOCS1', 'FAM175A', 'ROBO1', 'ARID1B', 'ARID1A', 'TCF7L2', 'STK11', 'FOXA1', 'PTEN', 'FAT1', 'FAS', 'CYLD', 'MAX', 'SH2D1A', 'APC', 'NTHL1', 'CTCF', 'KDM5C', 'KMT2C', 'ZFHX3', 'FOXP1', 'PIGA', 'CDKN1B', 'CDKN1A', 'FUBP1', 'MSH2', 'ID3', 'TNFRSF14', 'TRAF3', 'EP400', 'BRIP1', 'ARID4A', 'ARID4B', 'XRCC2', 'DAXX', 'SDHAF2', 'ASXL1', 'AMER1', 'RASA1', 'EGR1', 'MST1', 'SOX17', 'RUNX1', 'PIK3R3', 'NCOR1', 'NF1', 'JAK1', 'PTPRD', 'CHEK2', 'CHEK1', 'SMC1A', 'TMEM127', 'STAG1', 'RAD51', 'TCF3', 'STAG2', 'ARID2', 'RAD50', 'RNF43', 'PARP1', 'BLM', 'CUX1', 'RECQL', 'RAD21', 'PTPN2', 'PTPN1', 'SLX4', 'INHA', 'PAX5', 'IRF1', 'TP53', 'HLA-A', 'IRF8', 'CBL', 'TOP1', 'SHQ1', 'PRDM1', 'NSD1', 'ATXN2', 'CREBBP', 'HDAC4', 'SESN2', 'PPP2R1A', 'EPHA7', 'ATM', 'EPHA3', 'POT1', 'SMAD3', 'MOB3B', 'TBX3', 'POLE', 'ATR', 'FANCD2', 'FH', 'BCORL1', 'SOX9', 'IKZF3', 'TSC1', 'TP63', 'MRE11A', 'SDHC', 'BTG1', 'POLD1', 'CIITA', 'SMC3', 'SAMHD1', 'RTEL1', 'ECT2L', 'PIK3R2', 'CRBN', 'FANCC', 'NBN', 'FANCA', 'HLA-B', 'RECQL4', 'DUSP4', 'ERCC2', 'FBXW7', 'TGFBR2', 'TGFBR1', 'MSH3', 'RBM15', 'TET1', 'TET3', 'SESN3', 'MGA', 'LTB', 'FOXL2', 'SH2B3', 'BCOR', 'HIST1H1D', 'ATRX', 'EP300', 'RAD51C', 'RAD51B', 'HIST1H1B', 'TNFAIP3', 'DICER1', 'ARID5B', 'LATS2', 'FOXO1', 'KEAP1', 'EZH2', 'SP140', 'NKX3-1', 'PBRM1', 'PALB2', 'CIC', 'BRCA1', 'DTX1', 'FLCN', 'SPEN', 'CD58', 'ERCC3', 'ERCC4', 'MSH6', 'BCL11B', 'BMPR1A', 'ERF', 'BRCA2', 'NOTCH2', 'EED', 'MITF', 'ELF3', 'SMARCA4', 'BBC3', 'ANKRD11', 'CEBPA', 'BCL2L11', 'AXIN2', 'AXIN1', 'CDK12', 'ESCO2', 'MLH1', 'SDHB', 'MED12', 'HNF1A', 'RYBP', 'ATP6V1B2', 'DNMT3B', 'KMT2B', 'KMT2A', 'DNMT3A', 'NFKBIA', 'TRAF5', 'KMT2D', 'SPOP', 'RBM10', 'P2RY8', 'TP53BP1', 'TSC2', 'KDM6A', 'EPCAM', 'PHOX2B', 'NPM1', 'BCL10', 'LATS1', 'HOXB13', 'ARID3A', 'PTPRT', 'PTPRS', 'INPPL1', 'NOTCH4', 'TET2', 'NOTCH1', 'CASP8', 'NOTCH3', 'GRIN2A', 'MAP2K4', 'WT1', 'BACH2', 'SDHA', 'BAP1', 'PTCH1', 'SDHD'])
	return tsgs

def get_oncogenes():
	impactGenes = get_im6_genes()
	tsgs = get_tsgs()
	return impactGenes - tsgs

def get_gene_size_info():
	geneLengthInfo = pd.read_table('~/Desktop/hypermutationProjectFinal/files/infoFiles/all_gene_coding_sequences.tsv')
	geneLengthDict = dict(zip(geneLengthInfo['hgnc_symbol'], geneLengthInfo['nt.length']))
	return geneLengthDict


#pathways taken from https://www.cell.com/cell/fulltext/S0092-8674(18)30359-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867418303593%3Fshowall%3Dtrue
#francisco cell paper
def get_pathway_genes(pathway):
	pathwayFileMap = {'WNT': '~/Desktop/hypermutationProjectFinal/files/infoFiles/wntBetaCateninPathwayData.csv',
					  'PI3K': '~/Desktop/hypermutationProjectFinal/files/infoFiles/PI3KPathwayData.csv'}
	df = pd.read_csv(pathwayFileMap[pathway])
	return set(df['Gene'])

def set_up_dep_map_data(depMapPath = 'files/infoFiles/projectAchillesDepMapData.csv', cancerTypeLines = [], loadedDf = pd.DataFrame([{}])):
	
	#if we are iterating over this function we have the option to only load the dataframe once
	depMapData = None
	if loadedDf.shape[0] > 0:
		depMapData = loadedDf
	else:
		depMapData = pd.read_csv(depMapPath)

	if len(cancerTypeLines) > 0:
		depMapData = depMapData[depMapData['Unnamed: 0'].isin(cancerTypeLines)]

	listOfDicts = []
	for col in depMapData.columns.values:
	    gene = col.split(' ')[0]
	    if col != 'Unnamed: 0':
	        listOfDicts.append({'Hugo_Symbol': gene, 'score': np.nanmean(depMapData[col])})
	depMapDf = pd.DataFrame(listOfDicts)
	return depMapDf

#two modes: getEssentialGenes, and getAllGenes which just returns the df of scores
def get_essential_genes(depMapPath = 'files/infoFiles/projectAchillesDepMapData.csv', mode='getEssentialGenes', cancerTypeLines = []):
	depMapDf = load_dep_map_data(depMapPath = depMapPath, cancerTypeLines = cancerTypeLines)

	if mode == 'getEssentialGenes':
		essentialGenes = set(depMapDf[depMapDf['score'] <= -1]['Hugo_Symbol'])
		return essentialGenes
	elif mode == 'getAllGenes':
		return depMapDf
	elif mode == 'getModerateEssentialGenes':
		moderateEssentialGenes = set(depMapDf[(depMapDf['score'] > -1)  & (depMapDf['score'] <= -.5)]['Hugo_Symbol'])
		return moderateEssentialGenes
	else:
		print 'error invalid mode specified'

#returns all genes that have dep map scores >.5 and are not in impact 468
def get_cancer_neutral_genes(depMapPath = 'files/infoFiles/projectAchillesDepMapData.csv', cancerTypeLines = []):
	
	depMapDf = load_dep_map_data(depMapPath = depMapPath, cancerTypeLines = cancerTypeLines)

	impactGenes = set(['ABL1', 'ACVR1', 'AGO2', 'AKT1', 'AKT2', 'AKT3', 'ALK', 'ALOX12B', 'ANKRD11', 'APC', 'AR', 'ARAF', 'ARID1A', 'ARID1B', 'ARID2', 'ARID5B', 'ASXL1', 'ASXL2', 'ATM', 'ATR', 'ATRX', 'AURKA', 'AURKB', 'AXIN1', 'AXIN2', 'AXL', 'B2M', 'BABAM1', 'BAP1', 'BARD1', 'BBC3', 'BCL10', 'BCL2', 'BCL2L1', 'BCL2L11', 'BCL6', 'BCOR', 'BIRC3', 'BLM', 'BMPR1A', 'BRAF', 'BRCA1', 'BRCA2', 'BRD4', 'BRIP1', 'BTK', 'CALR', 'CARD11', 'CARM1', 'CASP8', 'CBFB', 'CBL', 'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CD274', 'CD276', 'CD79A', 'CD79B', 'CDC42', 'CDC73', 'CDH1', 'CDK12', 'CDK4', 'CDK6', 'CDK8', 'CDKN1A', 'CDKN1B', 'CDKN2A', 'CDKN2B', 'CDKN2C', 'CEBPA', 'CENPA', 'CHEK1', 'CHEK2', 'CIC', 'CREBBP', 'CRKL', 'CRLF2', 'CSDE1', 'CSF1R', 'CSF3R', 'CTCF', 'CTLA4', 'CTNNB1', 'CUL3', 'CXCR4', 'CYLD', 'CYSLTR2', 'DAXX', 'DCUN1D1', 'DDR2', 'DICER1', 'DIS3', 'DNAJB1', 'DNMT1', 'DNMT3A', 'DNMT3B', 'DOT1L', 'DROSHA', 'DUSP4', 'E2F3', 'EED', 'EGFL7', 'EGFR', 'EIF1AX', 'EIF4A2', 'EIF4E', 'ELF3', 'EP300', 'EPAS1', 'EPCAM', 'EPHA3', 'EPHA5', 'EPHA7', 'EPHB1', 'ERBB2', 'ERBB3', 'ERBB4', 'ERCC2', 'ERCC3', 'ERCC4', 'ERCC5', 'ERF', 'ERG', 'ERRFI1', 'ESR1', 'ETV1', 'ETV6', 'EZH1', 'EZH2', 'FAM123B', 'FAM175A', 'FAM46C', 'FAM58A', 'FANCA', 'FANCC', 'FAT1', 'FBXW7', 'FGF19', 'FGF3', 'FGF4', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FH', 'FLCN', 'FLT1', 'FLT3', 'FLT4', 'FOXA1', 'FOXL2', 'FOXO1', 'FOXP1', 'FUBP1', 'FYN', 'GATA1', 'GATA2', 'GATA3', 'GLI1', 'GNA11', 'GNAQ', 'GNAS', 'GPS2', 'GREM1', 'GRIN2A', 'GSK3B', 'H3F3A', 'H3F3B', 'H3F3C', 'HGF', 'HIST1H1C', 'HIST1H2BD', 'HIST1H3A', 'HIST1H3B', 'HIST1H3C', 'HIST1H3D', 'HIST1H3E', 'HIST1H3F', 'HIST1H3G', 'HIST1H3H', 'HIST1H3I', 'HIST1H3J', 'HIST2H3C', 'HIST2H3D', 'HIST3H3', 'HLA-A', 'HLA-B', 'HNF1A', 'HOXB13', 'HRAS', 'ICOSLG', 'ID3', 'IDH1', 'IDH2', 'IFNGR1', 'IGF1', 'IGF1R', 'IGF2', 'IKBKE', 'IKZF1', 'IL10', 'IL7R', 'INHA', 'INHBA', 'INPP4A', 'INPP4B', 'INPPL1', 'INSR', 'IRF4', 'IRS1', 'IRS2', 'JAK1', 'JAK2', 'JAK3', 'JUN', 'KDM5A', 'KDM5C', 'KDM6A', 'KDR', 'KEAP1', 'KIT', 'KLF4', 'KMT2B', 'KMT5A', 'KNSTRN', 'KRAS', 'LATS1', 'LATS2', 'LMO1', 'LYN', 'MALT1', 'MAP2K1', 'MAP2K2', 'MAP2K4', 'MAP3K1', 'MAP3K13', 'MAP3K14', 'MAPK1', 'MAPK3', 'MAPKAP1', 'MAX', 'MCL1', 'MDC1', 'MDM2', 'MDM4', 'MED12', 'MEF2B', 'MEN1', 'MET', 'MGA', 'MITF', 'MLH1', 'KMT2A', 'KMT2B', 'KMT2C', 'MPL', 'MRE11A', 'MSH2', 'MSH3', 'MSH6', 'MSI1', 'MSI2', 'MST1', 'MST1R', 'MTOR', 'MUTYH', 'MYC', 'MYCL1', 'MYCN', 'MYD88', 'MYOD1', 'NBN', 'NCOA3', 'NCOR1', 'NEGR1', 'NF1', 'NF2', 'NFE2L2', 'NFKBIA', 'NKX2-1', 'NKX3-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NPM1', 'NRAS', 'NSD1', 'NTHL1', 'NTRK1', 'NTRK2', 'NTRK3', 'NUF2', 'NUP93', 'PAK1', 'PAK7', 'PALB2', 'PARK2', 'PARP1', 'PAX5', 'PBRM1', 'PDCD1', 'PDCD1LG2', 'PDGFRA', 'PDGFRB', 'PDPK1', 'PGR', 'PHOX2B', 'PIK3C2G', 'PIK3C3', 'PIK3CA', 'PIK3CB', 'PIK3CD', 'PIK3CG', 'PIK3R1', 'PIK3R2', 'PIK3R3', 'PIM1', 'PLCG2', 'PLK2', 'PMAIP1', 'PMS1', 'PMS2', 'PNRC1', 'POLD1', 'POLE', 'PPARG', 'PPM1D', 'PPP2R1A', 'PPP4R2', 'PPP6C', 'PRDM1', 'PRDM14', 'PREX2', 'PRKAR1A', 'PRKCI', 'PRKD1', 'PTCH1', 'PTEN', 'PTP4A1', 'PTPN11', 'PTPRD', 'PTPRS', 'PTPRT', 'RAB35', 'RAC1', 'RAC2', 'RAD21', 'RAD50', 'RAD51', 'RAD51C', 'RAD51L1', 'RAD51L3', 'RAD52', 'RAD54L', 'RAF1', 'RARA', 'RASA1', 'RB1', 'RBM10', 'RECQL', 'RECQL4', 'REL', 'RET', 'RFWD2', 'RHEB', 'RHOA', 'RICTOR', 'RIT1', 'RNF43', 'ROS1', 'RPS6KA4', 'RPS6KB2', 'RPTOR', 'RRAGC', 'RRAS', 'RRAS2', 'RTEL1', 'RUNX1', 'RXRA', 'RYBP', 'SDHA', 'SDHAF2', 'SDHB', 'SDHC', 'SDHD', 'SESN1', 'SESN2', 'SESN3', 'SETD2', 'SF3B1', 'SH2B3', 'SH2D1A', 'SHOC2', 'SHQ1', 'SLX4', 'SMAD2', 'SMAD3', 'SMAD4', 'SMARCA4', 'SMARCB1', 'SMARCD1', 'SMO', 'SMYD3', 'SOCS1', 'SOS1', 'SOX17', 'SOX2', 'SOX9', 'SPEN', 'SPOP', 'SPRED1', 'SRC', 'SRSF2', 'STAG2', 'STAT3', 'STAT5A', 'STAT5B', 'STK11', 'STK19', 'STK40', 'SUFU', 'SUZ12', 'SYK', 'TAP1', 'TAP2', 'TBX3', 'TCEB1', 'TCF3', 'TCF7L2', 'TEK', 'TERT', 'TET1', 'TET2', 'TGFBR1', 'TGFBR2', 'TMEM127', 'TMPRSS2', 'TNFAIP3', 'TNFRSF14', 'TOP1', 'TP53', 'TP53BP1', 'TP63', 'TRAF2', 'TRAF7', 'TSC1', 'TSC2', 'TSHR', 'U2AF1', 'UPF1', 'VEGFA', 'VHL', 'VTCN1', 'WHSC1', 'WHSC1L1', 'WT1', 'WWTR1', 'XIAP', 'XPO1', 'XRCC2', 'YAP1', 'YES1', 'ZFHX3', 'ZRSR2'])
	neutralGenes = set(depMapDf[depMapDf['score'] >= -.5]['Hugo_Symbol'])
	neutralGenes = neutralGenes - impactGenes
	return neutralGenes

#note cancer types must be passed in proper format
def get_dict_of_dep_map_data_by_cancer_type(cancerTypes):
	depMapInfo = pd.read_csv('~/Desktop/hypermutationProjectFinal/files/infoFiles/dep_map_sample_info.csv')
	d = {}
	for cancerType in cancerTypes:
		relevantCellLines = depMapInfo[depMapInfo['primary_disease'] == cancerType]['DepMap_ID']
		d[cancerType] = relevantCellLines
	return d
		
def get_cancer_type_specific_dep_map_data(cancerTypes):
	dictOfCellLines = get_dict_of_dep_map_data_by_cancer_type(cancerTypes)
	dictOfDepMapData = {}

	depMapData = pd.read_csv('~/Desktop/hypermutationProjectFinal/files/infoFiles/projectAchillesDepMapData.csv')

	for cancerType in cancerTypes:
		print 'analyzing ', cancerType
		depMapDf = set_up_dep_map_data(cancerTypeLines = dictOfCellLines[cancerType], loadedDf=depMapData)
		dictOfDepMapData[cancerType] = dict(zip(depMapDf['Hugo_Symbol'], depMapDf['score']))
	return dictOfDepMapData


#######
###########
#HYPERMUTATOR COHORTS

#note before running this function I added the n non-synonymous mutations to the signatures file.  I did it using the following

"""nonsynonMutTypes = ["Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", 
"Frame_Shift_Ins", "Frame_Shift_Del","In_Frame_Del",
"In_Frame_Ins","Translation_Start_Site","Splice_Site"]
mc3NonSynom = mc3maf[mc3maf['Variant_Classification'].isin(nonsynonMutTypes)]
nmutDict = dict(mc3NonSynom['SAMPLE_ID'].value_counts())
tcgaSigs = pd.read_table(pathPrefix + '/juno/work/taylorlab/friedman/myAdjustedDataFiles/tcgaSigsCombined.txt')
tcgaSigs['Sample_Name_Short'] = tcgaSigs['Sample Name'].apply(lambda x: x[:15])
tcgaSigs['nNonSynonymous'] = tcgaSigs['Sample_Name_Short'].apply(lambda x: nmutDict[x] if x in nmutDict else None)
tcgaSigs.to_csv(pathPrefix + '/juno/work/taylorlab/friedman/myAdjustedDataFiles/tcgaSigsCombined.txt',
index = False, sep='\t')
"""

#TODO: write a function to get exome signature cohorts
def get_exome_signature_cohorts(tcgaSigPath = '', exomeSigPath = ''):
	def dominant_signature_lambda_function(row, signatureColumns):
	    rowAsDict = row.to_dict()
	    listOfSigs = [(key, value) for key, value in rowAsDict.items() if key in signatureColumns]
	    sortedList = sorted(listOfSigs, key = lambda x: -1*x[1])
	    return sortedList[0][0]

	tcgaSignatures = pd.read_table(tcgaSigPath)
	exomeRecaptureSignatures = pd.read_table(exomeSigPath)

	exomeRecaptureSignatures['Signature.MMR'] = exomeRecaptureSignatures['Signature.6'] + exomeRecaptureSignatures['Signature.15'] + exomeRecaptureSignatures['Signature.20'] +  + exomeRecaptureSignatures['Signature.21'] + exomeRecaptureSignatures['Signature.26']
	exomeRecaptureSignatures['Signature.APOBEC'] = exomeRecaptureSignatures['Signature.2'] + exomeRecaptureSignatures['Signature.13']
	exomeRecaptureSignatures['Signature.SMOKING'] = exomeRecaptureSignatures['Signature.4'] + exomeRecaptureSignatures['Signature.24'] + exomeRecaptureSignatures['Signature.29']
	tcgaSignatures['Signature.MMR'] = tcgaSignatures['Signature.6'] + tcgaSignatures['Signature.15'] + tcgaSignatures['Signature.20'] +  + tcgaSignatures['Signature.21'] + tcgaSignatures['Signature.26']
	tcgaSignatures['Signature.APOBEC'] = tcgaSignatures['Signature.2'] + tcgaSignatures['Signature.13']
	tcgaSignatures['Signature.SMOKING'] = tcgaSignatures['Signature.4'] + tcgaSignatures['Signature.24'] + tcgaSignatures['Signature.29']

	signatureColumns = [col for col in tcgaSignatures.columns.values if 'Signature' in col]
	exomeRecaptureSignatures['dominantSignature'] = exomeRecaptureSignatures.apply(lambda row: dominant_signature_lambda_function(row, signatureColumns), axis=1)
	tcgaSignatures['dominantSignature'] = tcgaSignatures.apply(lambda row: dominant_signature_lambda_function(row, signatureColumns), axis=1)

	return dict(zip(exomeRecaptureSignatures['Sample Name'], exomeRecaptureSignatures['dominantSignature'])), dict(zip(tcgaSignatures['Sample Name'], tcgaSignatures['dominantSignature']))

nmutAttributionThresh = 400
def get_tcga_pole_mmr_hypermutator_ids(tcgaSigsPath = 'files/infoFiles/tcgaSigDecomposition.txt'):
	tcgaSigs = pd.read_table(tcgaSigsPath)
	
	tcgaSigs = tcgaSigs[tcgaSigs['nNonSynonymous'].notnull()]
	tcgaSigs['Signature.MMR'] = tcgaSigs['Signature.6'] + tcgaSigs['Signature.15'] + tcgaSigs['Signature.20']+ tcgaSigs['Signature.21'] + tcgaSigs['Signature.26']
	tcgaSigs['Signature.POLE'] = tcgaSigs['Signature.10'] + tcgaSigs['Signature.14']
	tcgaSigs['mmrAttributed'] = tcgaSigs['Signature.MMR'] * tcgaSigs['nNonSynonymous']
	tcgaSigs['poleAttributed'] = tcgaSigs['Signature.POLE'] * tcgaSigs['nNonSynonymous']

	mmrCases = set(tcgaSigs[tcgaSigs['mmrAttributed'] > nmutAttributionThresh]['Sample_Name_Short'])
	poleCases = set(tcgaSigs[tcgaSigs['poleAttributed'] > nmutAttributionThresh]['Sample_Name_Short'])

	return mmrCases, poleCases

#alert for now we are basing this on total TMB not 
def get_exome_recapture_pole_mmr_hypermutator_ids(exomeRecaptureSigsPath = 'files/infoFiles/exomeRecaptureSigDecomposition.tsv'):
	exomeSigs = pd.read_table(exomeRecaptureSigsPath)
	exomeSigs['Signature.MMR'] = exomeSigs['Signature.6'] + exomeSigs['Signature.15'] + exomeSigs['Signature.20']+ exomeSigs['Signature.21'] + exomeSigs['Signature.26']
	exomeSigs['Signature.POLE'] = exomeSigs['Signature.10'] + exomeSigs['Signature.14']
	exomeSigs['mmrAttributed'] = exomeSigs['Signature.MMR'] * exomeSigs['nNonSynonymous']
	exomeSigs['poleAttributed'] = exomeSigs['Signature.POLE'] * exomeSigs['nNonSynonymous']

	mmrCases = set(exomeSigs[exomeSigs['mmrAttributed'] > nmutAttributionThresh]['Sample Name'])
	poleCases = set(exomeSigs[exomeSigs['poleAttributed'] > nmutAttributionThresh]['Sample Name'])

	return mmrCases, poleCases

#function to get ids for getting ids for different groups of hypermutator cohorts
def get_hypermutator_signature_cohorts(impactSigsPath=''):

	impactSigs = pd.read_table(impactSigsPath)
	return dict(zip(impactSigs['Tumor_Sample_Barcode'], impactSigs['dominantSignature']))


###CANCER TYPE BASED cohorts

#all TCGA cancer types fyi: set(['TCGA-TGCT', 'TCGA-SARC', 'TCGA-KIRP', 'TCGA-UCS', 'TCGA-KICH', 'TCGA-LAML', 'TCGA-KIRC', 'TCGA-BRCA', 'TCGA-PRAD', 'TCGA-LUAD', 'TCGA-OV', 'TCGA-UVM', 'TCGA-PCPG', 'TCGA-PAAD', 'TCGA-READ', 'TCGA-DLBC', 'TCGA-ESCA', 'TCGA-THCA', 'TCGA-UCEC', 'TCGA-STAD', 'TCGA-CHOL', 'TCGA-LGG', 'TCGA-LIHC', 'TCGA-CESC', 'TCGA-COAD', 'TCGA-HNSC', 'TCGA-SKCM', 'TCGA-GBM', 'TCGA-THYM', 'TCGA-LUSC', 'TCGA-ACC', 'TCGA-MESO', 'TCGA-BLCA'])
#GIVES YOU THE IDS OF TCGA BY CANCER TYPE
def get_tcga_cancer_type_info(tcgaInfoPath = 'files/infoFiles/tcgaCancerTypeInfo.txt', cancerTypes = []):
	tcgaDf = pd.read_table(tcgaInfoPath)
	return dict(zip(tcgaDf['TCGA_ID'], tcgaDf['Clinical_type']))
	
def get_impact_cancer_type_info(impactCancerTypeInfoPath = '~/Desktop/hypermutationProjectFinal/files/infoFiles/cancerTypeInfo.txt'):
	impactCancerTypeDf = pd.read_table(impactCancerTypeInfoPath)
	return dict(zip(impactCancerTypeDf['#Sample Identifier'], impactCancerTypeDf['Cancer Type']))


def get_ids_by_hypermutant_status(hypermutantIdDir='/Users/friedman/Desktop/hypermutationProjectFinal/files/hypermutationStatusIds',
	cancerType='', hypermutantStatus = 'Hypermutated'):

	cancerTypeAdj = re.sub(' ', '_', cancerType)
	path = os.path.join(hypermutantIdDir, cancerTypeAdj + '.tsv')
	df = pd.read_table(path)
	if hypermutantStatus == 'all':
		return set(df['Tumor_Sample_Barcode'])
	else:
		return set(df[df['hypermutantClassification'] == hypermutantStatus]['Tumor_Sample_Barcode'])

#ENUMERATES EVERY SINGLE HYPERMUTATED CASES
def get_all_hypermutant_ids(hypermutantIdDir='/Users/friedman/Desktop/hypermutationProjectFinal/files/hypermutationStatusIds'):
	hypermutatedCases = set([])
	#these are the cancer types that we will ignore here
	for f in os.listdir(hypermutantIdDir):
	    cType = re.sub('_', ' ', f)[:-4]
	    hypermutatedCases = hypermutatedCases | get_ids_by_hypermutant_status(hypermutantIdDir=hypermutantIdDir, cancerType=cType, hypermutantStatus = 'Hypermutated')
	return hypermutatedCases

#ENUMERATES EVERY SINGLE NON-HYPERMUTATED CASES
def get_all_normal_ids(hypermutantIdDir='/Users/friedman/Desktop/hypermutationProjectFinal/files/hypermutationStatusIds'):
	hypermutatedCases = set([])
	#these are the cancer types that we will ignore here
	for f in os.listdir(hypermutantIdDir):
	    cType = re.sub('_', ' ', f)[:-4]
	    hypermutatedCases = hypermutatedCases | get_ids_by_hypermutant_status(hypermutantIdDir=hypermutantIdDir, cancerType=cType, hypermutantStatus = 'Normal')
	return hypermutatedCases

#enumerates msi case
def get_msi_cases(msiInfoFilePath = '', msiScoreThresh=10):
	msiInfoDf = pd.read_table(msiInfoFilePath)
	msiIds = set(msiInfoDf[msiInfoDf['MSI_SCORE'] > msiScoreThresh]['Tumor_Sample_Barcode'])
	return msiIds

def get_tcga_msi_cases(tcgaMsiScoresPath = '', msiScoreThresh=10):
	msiInfoDf = pd.read_table(tcgaMsiScoresPath, skiprows=[0,1])
	msiIds = set(msiInfoDf[msiInfoDf['MSIsensor score'] > msiScoreThresh]['Participant Barcode'])
	return msiIds


def get_all_tmb_info(tmbFilePath = '/juno/work/taylorlab/friedman/myAdjustedDataFiles/mutations_TMB_and_MSI_stats.txt'):
	df = pd.read_table(tmbFilePath)
	return dict(zip(df['Tumor_Sample_Barcode'], df['tmb']))

#SIGNATURE BASED COHORTS
#get all cases in impact that fit a signature cohort
def get_impact_signature_cohort(impactSigsPath, signatureName):
	sigsDf = pd.read_table(impactSigsPath)
	return set(sigsDf[sigsDf['dominantSignature'] == signatureName]['Tumor_Sample_Barcode'])

def get_pan_impact_signature_mapping(impactSigsPath):
	sigsDf = pd.read_table(impactSigsPath)
	return dict(zip(sigsDf['Tumor_Sample_Barcode'], sigsDf['dominantSignature']))

#returns a mapping from exome ids to dmp ids
def get_exome_id_mapping(exomeIdMappingPath = '~/Desktop/hypermutationProjectFinal/files/infoFiles/dmpToSampleNameMapping.tsv'):
	df = pd.read_table(exomeIdMappingPath)
	return dict(zip(df['Tumor_Sample_Barcode'], df['DMP'])), dict(zip(df['DMP'], df['Tumor_Sample_Barcode']))


##########N SAMPLE BASED COHORTS

###RETURNS ALL CASES IN MAF WHERE THERE ARE MULTIPLE SAMPLES
def enumerate_cases_with_multiple_mutations(maf):
	maf['pid'] = maf['Tumor_Sample_Barcode'].apply(lambda x: x[:9])
	limitedDf = maf.drop_duplicates(subset='Tumor_Sample_Barcode')
	multiplePatientCases = [patient for patient, count in dict(limitedDf['pid'].value_counts()).items() if count > 1]
	return multiplePatientCases

def get_hypermutator_pairs_samples():
	samples = set(['P-0022520-T02-IM6', 'P-0044415-T02-IM6','P-0044415-T01-IM6', 'P-0036882-T01-IM6',
	'P-0036882-T02-IM6','P0008345-T01-IM5', 'P-0008345-T02-IM6','P-0040250-T02-IM6', 'P-0040250-T01-IM6',
	'P-0002049-T02-IM6','P-0030205-T03-IM6', 
	'P-0030205-T02-IM6', 'P-0002463-T01-IM3','P-0002463-T02-IM5','P-0001685-T02-IM3','P-0001685-T01-IM3','P-0017986-T01-IM6',
	'P-0017986-T02-IM6', 'P-0001237-T01-IM3', 'P-0001237-T03-IM6', 'P-0001237-T02-IM6', 'P-0026297-T02-IM6', 'P-0026297-T01-IM6', 'P-0018106-T01-IM6',
	 'P-0018106-T02-IM6', 'P-0004379-T01-IM5', 'P-0004379-T02-IM6','P-0026221-T02-IM6','P-0026221-T01-IM6','P-0003650-T01-IM5',
	 'P-0003650-T02-IM6', 'P-0043518-T01-IM6','P-0043518-T02-IM6','P-0013991-T02-IM6','P-0013991-T01-IM5','P-0042357-T02-IM6',
	 'P-0042357-T01-IM6','P-0019556-T02-IM6','P-0019556-T01-IM6','P-0001882-T03-IM6','P-0001882-T02-IM5','P-0035180-T01-IM6',
	 'P-0002265-T04-IM5','P-0002265-T02-IM5','P-0001420-T02-IM5','P-0008682-T02-IM6','P-0008682-T01-IM5','P-0019199-T02-IM6',
	 'P-0019199-T01-IM6','P-0001325-T02-IM5', 'P-0004910-T03-IM5','P-0004910-T01-IM5', 'P-0004910-T04-IM5'])
	return samples

def get_tcga_whitelist():
	return 0


#GET RELATED GENES
def get_related_genes_by_cancer_type(
	thresh = 1.0/30.0, impactMafPath = '~/Desktop/hypermutationProjectFinal/files/mafs/impact_mutations_maf_filtered'):

	def enumerate_related_genes(maf, thresh):
	    nCases = len(set(maf['Tumor_Sample_Barcode']))
	    oncMaf = maf[maf['oncogenic'].notnull()]
	    oncMaf['caseGene'] = oncMaf['Tumor_Sample_Barcode'] + '_' + oncMaf['Hugo_Symbol']
	    geneCounter = Counter(oncMaf.drop_duplicates(subset=['caseGene'])['Hugo_Symbol'])
	    
	    return [gene for gene, count in geneCounter.items() if (1.0*count)/nCases > thresh]


	allImpactMutsMaf = pd.read_table(impactMafPath) #todo get this in a better way
	cancerTypeDict = dict(get_impact_cancer_type_info())
	allImpactMutsMaf['cancerType'] = allImpactMutsMaf['Tumor_Sample_Barcode'].apply(lambda x: cancerTypeDict[x] if x in cancerTypeDict else None)
	hypermutantIds = get_all_hypermutant_ids()
	normalIds = get_all_normal_ids()
	allImpactMutsMaf['hypermutationStatus'] = allImpactMutsMaf['Tumor_Sample_Barcode'].apply(lambda x:
	    'hypermutated' if x in hypermutantIds else 'normal' if x in normalIds else 'Intermediate')

	cancerTypes = set(allImpactMutsMaf['cancerType'])
	d = {}
	for ct in cancerTypes:
	    ctMaf = allImpactMutsMaf[allImpactMutsMaf['cancerType'] == ct]
	    if ctMaf.shape[0] > 0:
			normalMaf = ctMaf[ctMaf['hypermutationStatus'] == 'normal']
			relatedGenes = enumerate_related_genes(normalMaf, thresh)
			d[ct] = relatedGenes
	return d






