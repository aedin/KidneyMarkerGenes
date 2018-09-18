
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

downloadData<- function() {
  url="https://www.biorxiv.org/highwire/filestream/105699/field_highwire_adjunct_files/0/348615-1.xlsx"
  download.file(url, destfile = "./data/Clarke_supplement.xlsx")
}



#' Title
#'
#' @param MGI character vector of mouse gene symbols
#'
#' @return data.frame with columns c( MGI.symbol, HGNC.symbol,  Gene.stable.ID)
#' which are mouse gene symbols, human gene symbols and human EnsEMBL gene ids
#' @export
#'
#' @examples
#' mapMouse2Human(c("Cd4", "Cd8a", "Cd44"))
#'
mapMouse2Human<-function(MGI){
  library(AnnotationDbi)
  library(biomaRt)

  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  hugenes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = MGI ,mart = mouse, attributesL = c("hgnc_symbol","ensembl_gene_id"), martL = human, uniqueRows=T)
  return(hugenes)
}

#' Create RData file from Clarke Genes Supplementary materials
#'
#' @return
#' @export
#'
#' @examples
#' .clarkeGenes()
#' data(clarkeGenes)
.clarkeGenes<-function(){
  clarkeGenes<- readxl::read_xlsx("./data/Clarke_supplement.xlsx", sheet = "Complete Database")
  clarkeGenes <-as.data.frame(clarkeGenes)
  hugenes=mapMouse2Human(clarkeGenes$'Marker Gene Symbol'[1:323])
  clarkeGenes= merge(clarkeGenes, hugenes, by.x="Marker Gene Symbol", by.y="MGI.symbol")
  save(clarkeGenes,file="./data/clarkeGenes.rda")
}


#' Title
#'
#' @param id  Character either with value "Symbol" or "Ensembl". Return signature as gene symbols or EnsEMBL gene identifers
#' @param matrix  logical, default TRUE, return matrix, if FALSE return list
#' @param convertNA  Logical. Default is FALSE Replace NA with 0
#'
#' @return
#' @export
#'
#' @examples
#' clarkeGenesSigs(id= "Symbol",matrix=TRUE)
#' clarkeGenesSigs(id= "Ensembl",matrix=FALSE)
#'
clarkeGenesSigs<-function(id="Symbol", matrix=TRUE,  convertNA=FALSE){
  require(reshape2)
  if(id=="Ensembl")sigs<-reshape2::acast(clarkeGenes, Gene.stable.ID~`Cell Type`,median, na.rm=TRUE, value.var="Whole Kidney TPM")
  if(id=="Symbol")sigs<-reshape2::acast(clarkeGenes, HGNC.symbol~`Cell Type`,median, na.rm=TRUE, value.var="Whole Kidney TPM")

  if (!matrix) sigs<-apply(sigs, 2, function(x) x[!is.na(x)])
  if (convertNA) sigs[is.na(sigs)]<-0

  return(sigs)
}


#' excludeZeroSumRowCol: Remove rows and columns where the rowSum or colSum is zero
#'
#' @param x a matrix  or data.frame
#'
#' @return matrix where rowSum(x) and colSum(x) >0
#' @export
#' mm<-matrix(c(1:30, rep(0, 10)), ncol=4, nrow=10)
#' mm[4,]<-0
#' mm
#' excludeZeroSumRowCol(mm)
#'
#' @examples
excludeZeroSumRowCol<-function(x){
  x<-x[!rowSums(x)==0,]
  x<-x[,!colSums(x)==0]
  return(x)
}
