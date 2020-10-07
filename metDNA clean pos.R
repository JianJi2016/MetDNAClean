library(readr)
library(dplyr)
metDNA.file <- read_csv("MRN.annotation.result.csv")
metDNA.file %>% filter(!is.na(ID)) -> metDNA.file
Final.data <- NULL
total <- nrow(metDNA.file)
tpb <- txtProgressBar(min = 0,
                      max = 100,
                      # initial = 0,
                      char = ">",
                      width = 50,
                      style = 3
)
for (i in 1:total) {

demo <- slice(metDNA.file,i)
demo %>% select(name:rt) -> compound.info
demo %>% select(Annotation.type,ID,compound.name,isotope,adduct,score,confidence) -> select.file

cbind(
as.character(unlist(strsplit(pull(select.file[,1]), split = ";"))),
as.character(unlist(strsplit(pull(select.file[,2]), split = ";"))),
as.character(unlist(strsplit(pull(select.file[,3]), split = ";"))),
as.character(unlist(strsplit(pull(select.file[,4]), split = ";"))),
as.character(unlist(strsplit(pull(select.file[,5]), split = ";"))),
as.character(unlist(strsplit(pull(select.file[,6]), split = ";"))),
as.character(unlist(strsplit(pull(select.file[,7]), split = ";")))
) -> tobe.select

colnames(tobe.select) <- colnames(select.file)
tobe.select <- as.data.frame(tobe.select)

tobe.select$score <- as.numeric(as.matrix(tobe.select$score))

tobe.select %>% 
  filter(score > 0.4 &
         confidence %in% c("grade1","grade2")  &
         isotope == "[M]" &
         adduct %in% c("M+H","2M+H","2M+NH4","M+NH4")) %>% 
  arrange(desc(score)) %>% slice(1) -> selected.result

if(nrow(selected.result) == 1){
  final.selected <- cbind(compound.info,selected.result)
  Final.data <- rbind(final.selected,Final.data)
 } 

setTxtProgressBar(pb = tpb,
                  value = i/total*100)
}



write.csv(Final.data,"Final.data.pos.csv")
