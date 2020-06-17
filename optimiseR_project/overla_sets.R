library(Biostrings)
library(dplyr)

# Define the sets of overlaps to choose from

set1 <- c("TGCC", "ACTA", "TTAC", "CAGA", "TGTG", "GAGC", "AGGA", "CGAA", "ATAG", "AAGG", "AACT", "AAAA", "ACCG")

set2 <- c("AGTG", "CAGG", "ACTC", "AAAA", "AGAC","CGAA","ATAG","AACC","TACA","TAGA","ATGC","GATA","CTCC","GTAA","CTGA","ACAA","AGGA","ATTA", "ACCG","GCGA")

set3 <- c("CCTC", "CTAA", "GACA", "GCAC", "AATC", "GTAA", "TGAA", "ATTA", "CCAG", "AGGA", "ACAA", "TAGA", "CGGA", "CATA", "CAGC", "AACG","AAGT", "CTCC", "AGAT", "ACCA", "AGTG", "GGTA", "GCGA", "AAAA", "ATGA")

# define the gene to be split up
dppr_gene <- "tttcagggcgcCATGGGCAACGATGTGGTTACCTACACCACCCTGATCGATGGTCTGGCGAAAGCGGGCCGTCTGGAGGAAGCGCTGCAGCTGTTCCAAGAGATGAAGGAAAAAGGTGTGAAGCCGGACGTGGTTACCTATAACACCCTGATTGACGGTCTGGCGAAGGCTGGCCGTCTGGAGGAAGCGCTGCAACTGTTTCAGGAAATGAAGGAAAAAGGCGTTAAACCGGATGTTGTGACCTATACCACCCTGATTGATGGTCTGGCGAAGGCCGGCCGTCTGGAGGAAGCGCTGTGCCTGTTCCAGGAGATGAAAGAGAAAGGTGTGAAGCCGAATGTTGTGACCTACAATACCCTGATCGACGGTCTGGCGAAGGCAGGCCGTCTGGAGGAAGCGCTGCAGTTATTTCAGGAGATGAAAGAAAAAGGCGTTAAGCCGGATGTTGTTACCTACAACACTTTAATTGACGGTCTGGCGAAAGCTGGCCGTCTGGAGGAAGCGCTGCAGCTATTTCAGGAGATGAAGGAGAAAGGTGTGAAACCGAGCGTGGTGACCTACAACACTTTAATCGATGGTCTGGCGAAGGCGGGCCGTCTGGAGGAAGCGCTGCAGCTTTTTCAGGAGATGAAGGAAAAAGGCGTTAAGCCGAGCGTGGTTACCTACAACACTCTTATTGATGGTCTGGCGAAAGCTGGCCGTCTGGAGGAAGCGCTGCAGCTCTTTCAGGAGATGAAGGAAAAAGGTGTGAAACCGGATGTTGTTACCTATAATACTTTGATAGATGGTCTGGCGAAGGCGGGCCGTCTGGAGGAAGCGCTGCAGCTATTCCAGGAAATGAAAGAGAAAGGCGTTAAACCGGACGTTGTGACCTACACCACTTTAATTGATGGTCTGGCGAAAGCCGGCCGTCTGGAGGAAGCGCTGCAGCTATTCCAAGAAATGAAAGAAAAAGGTGTGAAGCCTAACGTGGTTACCTATACCACTTTAATCGACGGTCTGGCGAAAGCCGGCCGTCTGGAGGAAGCGCTGCAGCTATTCCAAGAGATGAAAGAGAAAGGCGTTAAGCCGAATGTGGTGACCTACAACACATTAATTGATGGTCTGGCGAAAGCCGGCCGTCTGGAGGAAGCGCTGTGTCTGTTCCAGGAGATGAAGGAAAAAGGTGTGAAGCCGAGCGTGGTGACCTATAATACTTTGATCGATGGTCTGGCGAAAGCAGGCCGTCTGGAGGAAGCGCTGCAGCTATTTCAAGAAATGAAAGAAAAAGGCGTTAAGCCGAGCGTGGTAACCTATACCACTCTTATCGACGGTCTGGCGAAAGCAGGCCGTCTGGAGGAAGCGCTGCAGCTATTTCAAGAGATGAAAGAAAAAGGTGTGAAACCTAATGTTGTGACCTACAACACTCTTATCGATGGTCTGGCGAAAGCAGGCCGTCTGGAGGAAGCGCTGCAGCTATTTCAAGAGATGAAGGAGAAAGGCGTTAAGCCTGATGTGGTTACCTATAATACTCTCATCGATGGTCTGGCGAAAGCAGGCCGTCTGGAGGAAGCGCTGCAGCTATTTCAAGAGATGAAGGAAAAAGGTGTGAAGCCAGATGTTGTGACCTATAATACTCTAATAGATGGTCTGGCGAAGGCGGGCCGTCTGGAGGAAGCGCTGCAGCTATTTCAAGAGATGAAGGAAAAAGGTGTTAAACCGGATGTGGTTACCTACAACACACTTATTGATGGTCTGGCGAAAGCAGGCCGTCTGGAGGAAGCGCTGCAGCTATTTCAAGAGATGAAGGAAAAGGGTGTGAAACCGAGCGTTGTGACCAATAACACCCTGAAGGACGGCGCGAGCAAAGCGGGCTAAGAATTCGAGCTCCGTCGACAAG"


# define functions for finding the positions
find_overlap_pos <- function(overlap, gene){
  matches = matchPattern(overlap, gene)
  match.df = data.frame(
    overlap = overlap,
    start = matches@ranges@start
  )
  return(match.df)
}



find_all_overlaps <- function(overlap_list, some_sequence){
  results_list <- list()
  for(i in 1:length(overlap_list)){
    try(
      temp1 <- find_overlap_pos(overlap = overlap_list[i], 
                       gene = some_sequence)
        )
    try(
      results_list[[i]] <- temp1
    )
      
  }
  return(do.call(rbind, results_list))
}

#find initially, all potential positions in the gene.
found_positions <- find_all_overlaps(set1, dppr_gene)

#set parameters for the selection criteria of the fragments
total_length <- stringr::str_length(dppr_gene)
min_length <- 300
max_length <- 500
overlap_length <- 4


potential_starting_pos <- found_positions %>% 
  filter(start > min_length & start < max_length)
potential_starting_pos



#from a list of potential positions, select one (some_integer) and finds the next potential position

remaining_positions <- found_positions

find_next_positons <- function(current_position, some_integer){
  
  #renames coloumds of the input dataframe
  colnames(current_position) <- c("overlap", "start")

  #select a starting position to work from. Uses the sequence and location to calculate further on.
  position = current_position$overlap[some_integer]
  
  n = current_position$start[some_integer]
  

  #update object that contains the remaining positions available
  remaining_positions <- remaining_positions %>% filter(overlap != position)
  assign("remaining_positions", remaining_positions, envir = .GlobalEnv)

  #filters the set of all potential positions for overlaps that have not been used (position) and
  #that they are a certain length away from the current position (more than min_length, less than max_length)

  df = remaining_positions %>%
    filter(start > (n + min_length) & start < (n + max_length))

  #creates a dataframe with results, using
  temp1 = data.frame(
    current.end.overlap = position,
    current.end = n,
    next.end.overlap = df$overlap,
    next.end = df$start
    )


  return(temp1)
  
}

second_pos <- find_next_positons(potential_starting_pos, 5)

third_pos <- find_next_positons(third_pos[,4-1:0], 10)

third_pos %>% filter((total_length - next.end) > 300)



























































