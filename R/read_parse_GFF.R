#' Reads and Parses GFF File Into Data.Frame Object
#'
#' This function reads in a GFF file, parses it to extract necessary info,
#' and returns a data.frame object for use with downstream functions.
#' See ?associate_scaffolds_GFF()
#'
#' Function To Read a GFF File into a data.frame
#' @param file Input GFF file
#' @noRd
read_gff <- function(file){
  # Read gff file into dataframe. Uses readr.
  readr::read_tsv(
    file,
    col_names = c(
      "seqid",
      "source",
      "type",
      "start",
      "stop",
      "score",
      "strand",
      "phase",
      "attr"
    ),
    na        = ".",
    comment   = "#",
    col_types = "ccciidcic"
  )
}

#' Function to Read and Parse a GFF file
#'
#' This function takes a GFF file path as input. It will read it, parse it,
#' and return a data.frame for use with associate_scaffolds_GFF().
#' @param gff.filepath Path to GFF file
#' @return data.frame object with useful GFF info
#' @export
#' @examples
#' gff <- parseGFF("./genes.gff")
parseGFF <- function(gff.filepath){

  `%>%` <- dplyr::`%>%`

  ##############################################################
  ### Read and parse the gff file
  ##############################################################
  # Read in the gff file.
  g <- read_gff(gff.filepath)

  # Get these attributes from attr column.
  tags <- c("ID", "Parent", "gene", "product", "Name")

  # Parse and Add the tags from above as separate columns to the gff dataframe.
  # Removes e.g. "ID="; "gene=".
  # Source:
  #https://cran.r-project.org/web/packages/rmonad/vignettes/gff-processing.html
  gff <- dplyr::data_frame(
    attr  = stringr::str_split(g$attr, ";"),
    order = 1:nrow(g)
  ) %>%
    dplyr::mutate(ntags = sapply(attr, length)) %>%
    tidyr::unnest(attr) %>%
    dplyr::mutate(attr = ifelse(grepl('=', attr),
                                attr, paste(".U", attr, sep="="))) %>%
    tidyr::separate_(
      col   = "attr",
      into  = c("tag", "value"),
      sep   = "=",
      extra = "merge"
    ) %>%
    dplyr::filter(tag %in% c(tags, ".U")) %>%
    {
      if(nrow(.) > 0){
        tidyr::spread(., key="tag", value="value")
      } else {
        .$tag   = NULL
        .$value = NULL
        .
      }
    } %>%
    {
      if("Parent" %in% names(.)){
        .$Parent <- ifelse(.$Parent == "-", NA, .$Parent)
      }
      .
    } %>% {
      for(tag in c(tags, ".U")){
        if(! tag %in% names(.))
          .[[tag]] = NA_character_
      }
      .
    } %>%
    {
      if("ID" %in% names(.))
        .$ID <- ifelse(is.na(.$ID) & !is.na(.$.U) & .$ntags == 1, .$.U, .$ID)
      .
    } %>%
    merge(dplyr::data_frame(order=1:nrow(g)), all=TRUE) %>%
    dplyr::arrange(order) %>%
    { cbind(g, .) } %>%
    dplyr::select(-.U, -order, -ntags, -attr) %>%
    {
      if(all(c("ID", "Parent", "gene", "product", "Name") %in% names(.))){
        parents <- subset(., type %in% c("CDS", "exon"))$Parent
        parent_types <- subset(., ID %in% parents)$type

        if(any(parent_types == "gene"))
          warning("Found CDS or exon directly inheriting from a gene, this may
                  be fine.")

        if(! all(parent_types %in% c("gene", "mRNA")))
          warning("Found CDS or exon with illegal parent")

        if( any(is.na(parents)) )
          warning("Found CDS or exon with no parent")

        if(! any(duplicated(.$ID, incomparables=NA)))
          warning("IDs are not unique, this is probably bad")
      }
      .
    }

  # Strip GenBank version number from TranscriptId. E.g. "XM3333.1" to "XM3333"
  gff$TranscriptId <- gsub("\\..*", "", gff$Name)

  # Strip type prefix from TranscriptId. E.g. "exon-XM3333" to "XM3333".
  #gff$TranscriptId <- gsub("^[^-]*-", "", gff$TranscriptId)

  # Filter down to have only mRNA entries.
  # This allows me to extract the full gene length.
  gff <- subset(gff, tolower(type) == "mrna")

  return(gff)
}
