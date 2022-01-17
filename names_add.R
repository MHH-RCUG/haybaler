# Script to filter names out of name information
# Lisa Hollstein, Dec 2021 - Jan 2022

args <- commandArgs()

print(args)

setwd(args[6])

input_file <- args[7]

bacteria_data <- read.csv(input_file, sep = "\t", dec = ",", row.names = 1, header = TRUE)

# create dataframe in which rownames are shortened
sn_bacteria_data <- bacteria_data


# check if new name already exists and rename the row
check_existing <- function(new_name, old_name) {
  if(new_name %in% rownames(sn_bacteria_data)){
    # add values of this row to values of existing row with the same short name
    sn_bacteria_data[new_name,] <- as.numeric(sn_bacteria_data[new_name,]) + as.numeric(sn_bacteria_data[old_name,])
    sn_bacteria_data <- sn_bacteria_data[rownames(sn_bacteria_data) != old_name,]
    warning(paste(old_name,"has been added to", new_name))
  } else{
    rownames(sn_bacteria_data)[rownames(sn_bacteria_data) == old_name] <- new_name  # rename row
  }
  return(sn_bacteria_data)
}

# add subspecies to new name
subspecies <- function(new_name, split_name) {
  a <- 3
  add <- split_name[n+a]
  while(nchar(add) == 0){
    a <- a + 1
    add <- split_name[n+a]
  }
  new_name <- paste(new_name, "subsp", add)
  return(new_name)
}


row_names <- rownames(bacteria_data)
final_names <- c()
# extract name of species for every rowname
for (i in row_names){
  split_name <- unlist(strsplit(i, split = '_', fixed = TRUE))
  # find the species name
  for (n in 1:(length(split_name)-1)){ 
    first_letter <- utf8ToInt(substring(split_name[n], 1, 1))
    rest <- utf8ToInt(substring(split_name[n], 2))
    second_letter <- utf8ToInt(substring(split_name[n+1], 1, 1))
    if (split_name[n] == "organism" || split_name[n] == "candidatus"){
      name <- paste(split_name[n], split_name[n+1], split_name[n+2])
      if (split_name[n+3] == "subsp") {  # check for subspecies
        name <- subspecies(name, split_name)
      }
      final_names <- c(final_names, name)
      sn_bacteria_data <- check_existing(name, i)
      break
    }else if(length(first_letter) == 0 || length(second_letter) == 0 || length(rest) == 0){
    }else if(64 < first_letter && first_letter < 91 && 96 < second_letter && second_letter < 123){
      name <- paste(split_name[n], split_name[n+1])
      if (split_name[n+2] == "subsp"){  # check for subspecies
        name <- subspecies(name, split_name)
      }
      
      final_names <- c(final_names,name)
      sn_bacteria_data <- check_existing(name, i)
    
      break
    }
  }
}


setwd(paste0(args[6],"/short_names"))

# extract filename from input_file 
filename <- gsub("\\.[^.]*$", "", tail(unlist(strsplit(input_file, split="/")), n = 1))

name_check <- data.frame(row_names, final_names)
colnames(name_check) <- c("original names", "short names")
write.csv(name_check, paste0(filename, "_nametbl.csv"), row.names = FALSE)

if (length(row_names) != length(final_names)){
  stop("A problem occured while cutting the names")
} else {
  write.csv(sn_bacteria_data,paste0(filename,"_short.csv"))
}

warnings()



