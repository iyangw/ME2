setwd("/data/wangyy/CUT_Tag/reanalysis/11_chipqc/H3K4me3")

# Annotation of latest genome version generation
library(GenomicFeatures)
library(GenomicRanges)
library(IRanges)

specie = "mmu"
txdb=switch(specie,
            mmu=TxDb.Mmusculus.UCSC.mm39.knownGene::TxDb.Mmusculus.UCSC.mm39.knownGene,
            hsa=TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
)

All5utrs <- reduce(unique(unlist(fiveUTRsByTranscript(txdb))))
All3utrs <- reduce(unique(unlist(threeUTRsByTranscript(txdb))))
Allcds <- reduce(unique(unlist(cdsBy(txdb,"tx"))))
Allintrons <- reduce(unique(unlist(intronsByTranscript(txdb))))
Alltranscripts <- reduce(unique(transcripts(txdb)))

posAllTranscripts <- Alltranscripts[strand(Alltranscripts) == "+"]
posAllTranscripts <- posAllTranscripts[!(start(posAllTranscripts)-20000 < 0)]
negAllTranscripts <- Alltranscripts[strand(Alltranscripts) == "-"]
chrLimits <- seqlengths(negAllTranscripts)[as.character(seqnames(negAllTranscripts))]      
negAllTranscripts <- negAllTranscripts[!(end(negAllTranscripts)+20000 > chrLimits)]      
Alltranscripts <- c(posAllTranscripts,negAllTranscripts)
Promoters500 <-  reduce(flank(Alltranscripts,500))    
Promoters2000to500 <-  reduce(flank(Promoters500,1500))
LongPromoter20000to2000  <- reduce(flank(Promoters2000to500,18000))

annotation = list(
  version="mm39",
  LongPromoter20000to2000=LongPromoter20000to2000,
  Promoters2000to500=Promoters2000to500,
  Promoters500=Promoters500,
  All5utrs=All5utrs,
  Alltranscripts=Alltranscripts,
  Allcds=Allcds,
  Allintrons=Allintrons,
  All3utrs=All3utrs
)


library(ChIPQC)

# patch QCmetrics bug
m <- getMethod("QCmetrics", "ChIPQCsample")
f <- m@.Data
txt <- paste(deparse(f), collapse = "\n")

pat <- "(?s)names\\s*\\(\\s*res\\s*\\)\\s*(<-|=)\\s*c\\(.*?\\)"

rep <- paste(
  'nm <- c("Reads","Map%","Filt%","Dup%","ReadL","FragL","RelCC","SSD","RiP%")',
  'out <- setNames(rep(NA_real_, length(nm)), nm)',
  
  'k <- min(5L, length(res))',
  'if (k > 0) out[seq_len(k)] <- as.numeric(res[seq_len(k)])',
  
  'frag_val <- try(fragmentlength(object), silent = TRUE)',
  'relcc_val <- try(RelativeCrossCoverage(object), silent = TRUE)',
  'ssd_val <- try(ssd(object), silent = TRUE)',
  'rip_val <- try(frip(object) * 100, silent = TRUE)',
  
  'is_missing <- c(',
  '  FragL = inherits(frag_val, "try-error") || length(frag_val) == 0 || all(is.na(frag_val)),',
  '  RelCC = inherits(relcc_val, "try-error") || length(relcc_val) == 0 || all(is.na(relcc_val)),',
  '  SSD   = inherits(ssd_val, "try-error")  || length(ssd_val)  == 0 || all(is.na(ssd_val)),',
  '  `RiP%`= inherits(rip_val, "try-error")  || length(rip_val)  == 0 || all(is.na(rip_val))',
  ')',
  
  'if (!is_missing["FragL"]) out["FragL"] <- as.numeric(frag_val)[1]',
  'if (!is_missing["RelCC"]) out["RelCC"] <- as.numeric(relcc_val)[1]',
  'if (!is_missing["SSD"])   out["SSD"]   <- as.numeric(ssd_val)[1]',
  'if (!is_missing["RiP%"])  out["RiP%"]  <- as.numeric(rip_val)[1]',
  
  'res <- out',
  sep = "\n"
)

txt2 <- sub(pat, rep, txt, perl = TRUE)

f2 <- eval(parse(text = txt2))

setMethod("QCmetrics", "ChIPQCsample", f2)

# Now we can start analysis
samples <- read.csv('H3K4me3.csv')
H3K4me3 <- ChIPQC(samples, annotation=annotation, chromosomes = NULL) 
ChIPQCreport(H3K4me3, reportName="ChIP QC report H3K4me3", reportFolder="report")

