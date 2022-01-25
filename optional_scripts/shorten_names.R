# Script to filter names out of name information
# Lisa Hollstein, Dec 2021 - Jan 2022

args <- commandArgs()

setwd(args[6])

input_file <- args[7]

bacteria_data <- read.csv(input_file, sep = "\t", dec = ",", row.names = 1, header = TRUE)

# create dataframe in which rownames are shortened
sn_bacteria_data <- bacteria_data

duplicates <- c()


# check if new name already exists and rename the row
check_existing <- function(new_name, old_name, duplicates) {
  if(new_name %in% rownames(sn_bacteria_data)){
    # keep long name if short row name wouldn't be unique
    place <- match("Botrytis cinerea", rownames(sn_bacteria_data))
    rownames(sn_bacteria_data)[place] <- rownames(bacteria_data)[place]
    duplicates <- c(duplicates, new_name)
  } else if (!(new_name %in% duplicates)){
    rownames(sn_bacteria_data)[rownames(sn_bacteria_data) == old_name] <- new_name  # rename row
  }
  return_items <- list("duplicates" = duplicates, "data" = sn_bacteria_data)
  return(return_items)
}

# add subspecies to new name
subspecies <- function(new_name, split_name) {
  a <- 3
  add <- split_name[n+a]
  while(nchar(add) == 0){
    a <- a + 1
    add <- split_name[n+a]
  }
  new_name <- paste(new_name, "subsp", add, sep = "_")
  return(new_name)
}


row_names <- rownames(bacteria_data)
final_names <- c()
# extract name of species for every rowname
for (i in row_names){
  split_name <- unlist(strsplit(i, split = '_', fixed = TRUE))
  # find the species name
  for (n in 1:(length(split_name)-1)){ 
    # Using Unicode to find the species 
    # Must be a capital letter:
    first_letter <- utf8ToInt(substring(split_name[n], 1, 1))
    # Must be a lowercase letter:
    rest <- utf8ToInt(substring(split_name[n], 2)) # second letter in first word
    # Must be a lowercase letter:
    second_letter <- utf8ToInt(substring(split_name[n+1], 1, 1))
    # if "organism" or "candidatus" is in the name, 
    # the species name consists of the two following words
    if (split_name[n] == "organism" || split_name[n] == "candidatus"){
      new_name <- paste(split_name[n], split_name[n+1], split_name[n+2], sep = "_")
      if (split_name[n+3] == "subsp") {  # check for subspecies
        new_name <- subspecies(new_name, split_name)
      }
      return_items <- check_existing(new_name, i, duplicates)
      sn_bacteria_data <- return_items$data
      duplicates <- return_items$duplicates
      if (new_name %in% duplicates) {
        final_names <- c(final_names, i)
      } else {
        final_names <- c(final_names, new_name)
      }
      
      break
    }else if(length(first_letter) == 0 || length(second_letter) == 0 || length(rest) == 0){
      # Do nothing if one of the letters doesn't exist
    }else if(64 < first_letter && first_letter < 91 && 96 < second_letter && second_letter < 123){
      # Unicode value of lowercase letters is between 97 and 122, 
      # for capital letters it is between 65 and 90
      new_name <- paste(split_name[n], split_name[n+1], sep = "_")
      if (split_name[n+2] == "subsp"){  # check for subspecies
        new_name <- subspecies(new_name, split_name)
      }
      
      return_items <- check_existing(new_name, i, duplicates)
      sn_bacteria_data <- return_items$data
      duplicates <- return_items$duplicates
      if (new_name %in% duplicates) {
        final_names <- c(final_names, i)
      } else {
        final_names <- c(final_names, new_name)
      }
      
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



