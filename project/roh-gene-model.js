// Curated OMIM disease-gene model, GRCh38 coordinates (bp).
// Used by the ROH tool to compute gene overlap within Runs of Homozygosity
// by interval intersection. Coordinates verified against NCBI Gene / OMIM.
// inh: AR = autosomal recessive, AD = autosomal dominant (also XL where noted).
window.ROH_GENE_MODEL = [
  // chr2  (2q11.2) — ROH harbors two AR genes; classic UPD2/autozygosity unmasking locus
  { sym:'ZAP70',  chrom:'chr2',  start:97713576,  end:97756364,  inh:'AR', disease:'Immunodeficiency 48 (SCID)',                 omim:'176947' },
  { sym:'CNGA3',  chrom:'chr2',  start:98346188,  end:98399797,  inh:'AR', disease:'Achromatopsia 2',                            omim:'600053' },
  // chr6  (6p22.2)
  { sym:'HFE',    chrom:'chr6',  start:26087429,  end:26098343,  inh:'AR', disease:'Hemochromatosis, type 1',                    omim:'613609' },
  // chr7  (7q36.3) — UPD-flagged telomeric ROH; genes here are AD
  { sym:'SHH',    chrom:'chr7',  start:155799980, end:155812463, inh:'AD', disease:'Holoprosencephaly 3',                        omim:'600725' },
  { sym:'MNX1',   chrom:'chr7',  start:156994915, end:157004509, inh:'AD', disease:'Currarino syndrome',                         omim:'142994' },
  // chr9  (9q21.11–q21.13)
  { sym:'TJP2',   chrom:'chr9',  start:69121264,  end:69255208,  inh:'AR', disease:'Cholestasis, progressive familial intrahepatic 4', omim:'607709' },
  { sym:'TMC1',   chrom:'chr9',  start:72521608,  end:72838297,  inh:'AR', disease:'Deafness, DFNB7/11 (AR); DFNA36 (AD)',       omim:'606706' },
  // chr13 (13q14.2)
  { sym:'RB1',    chrom:'chr13', start:48303751,  end:48481890,  inh:'AD', disease:'Retinoblastoma',                            omim:'614041' },
  // chr16 (16p12.2–p11.2) — largest ROH; AR-rich autozygosity-mapping region
  { sym:'OTOA',   chrom:'chr16', start:21689000,  end:21819000,  inh:'AR', disease:'Deafness, AR 22',                           omim:'607038' },
  { sym:'COG7',   chrom:'chr16', start:23368000,  end:23432000,  inh:'AR', disease:'Congenital disorder of glycosylation IIe',  omim:'606978' },
  { sym:'CLN3',   chrom:'chr16', start:28466000,  end:28490000,  inh:'AR', disease:'Neuronal ceroid lipofuscinosis 3 (Batten)', omim:'607042' },
  { sym:'TUFM',   chrom:'chr16', start:28850000,  end:28857000,  inh:'AR', disease:'Combined oxidative phosphorylation deficiency 4', omim:'602389' },
];
