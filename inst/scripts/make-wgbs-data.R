#Example WGBS data downloaded from ENCODE. It is converted to GRanges format,
#so it can be directly used as example data.
#source from https://www.encodeproject.org/reference-epigenomes/ENCSR527XKV/

#Read the BED methylation calls as methylKit object
wgbs <- methylKit::methRead("ENCFF059FZP.bed",
                           sample.id = "sample1", assembly = "hg38",
                           header = FALSE,context = "CpG",
                           resolution = "base", 
                           pipeline = list(fraction = FALSE,
                                           chr.col = 1, start.col = 3, 
                                           end.col = 3,coverage.col = 5,
                                           strand.col = 6,freqC.col = 11))

#Subset the first 25000 rows
WGBS_GRanges <- methylKit::select(wgbs,1:25000)

#Convert the raw calls to GRanges and extract methylation values with convertMe
WGBS_GRanges <- deconvR:::convertMe(WGBS_GRanges)

#Save the example data as rda in extdata folder
save(WGBS_GRanges,file = "extdata/WGBS_GRanges.rda")

#Recompress it 
tools::resaveRdaFiles("extdata/WGBS_GRanges.rda","auto")

