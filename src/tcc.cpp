#include <RcppEigen.h>
#include "data_manipulation.h"
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]


//[[Rcpp::export]]
StringVector TCCMapC(StringVector from, std::string from_type,
                     std::string to, Environment ec_to_enst,
                     Environment enst_to_ec, Environment enst_to_ensg,
                     Environment ensg_to_enst, Environment ensg_to_gene,
                     Environment gene_to_ensg){
  StringVector ec;
  StringVector enst;
  StringVector ensg;
  StringVector gene;
  if (from_type == "EC") {
    ec = from;
    enst = HashTableLookup(ec_to_enst, from);
    if(to == "ENST") {
      return(enst);
    }
    from = enst;
    from_type = "ENST";
  }
  if (from_type == "ENST") {
    enst = from;
    if (to == "EC") {
      return(HashTableLookup(enst_to_ec, from));
    }
    if (to == "ENSG" || to == "GENE") {
      ensg = HashTableLookup(enst_to_ensg, enst);
      if (to == "ENSG") {
        return(ensg);
      }
    }
    from = ensg;
    from_type = "ENSG";
  }
  if (from_type == "ENSG") {
    ensg = from;
    if (to == "GENE") {
      return(HashTableLookup(ensg_to_gene, from));
    } else {
      enst = HashTableLookup(ensg_to_enst, from);
      if(to == "ENST"){
        return(enst);
      }
      if(to == "EC"){
        return(HashTableLookup(enst_to_ec, enst));
      }
    }
  }
  if(from_type == "GENE") {
    gene = from;
    ensg = HashTableLookup(gene_to_ensg, from);
    if (to == "ENSG") {
      return(ensg);
    }
    enst = HashTableLookup(ensg_to_enst, ensg);
    if (to == "ENST") {
      return(enst);
    }
    return(HashTableLookup(enst_to_ec, enst));
  }
  return(ec);
}

//Return only ECs that map uniquely to one gene
//[[Rcpp::export]]
CharacterVector ECUniqueGene(CharacterVector ecs, bool ensg, bool verbose,
                             Environment ec_to_enst, Environment enst_to_ec, 
                             Environment enst_to_ensg, Environment ensg_to_enst, 
                             Environment ensg_to_gene, Environment gene_to_ensg) {
  CharacterVector unambig_ecs;
  std::string to_type;
  if (ensg){
    to_type = "ENSG";
  } else {
    to_type = "GENE";
  }
  for(int i=0; i<ecs.size(); ++i){
    CharacterVector ec = Rcpp::as<CharacterVector>(ecs[i]);
    CharacterVector genes;
    try {
      genes = TCCMapC(ec, "EC", to_type, ec_to_enst, enst_to_ec,
                      enst_to_ensg, ensg_to_enst, ensg_to_gene,
                      gene_to_ensg);
    } catch(std::exception& e) {
      if (verbose) {
        Rcout << ec << " doesn't map to any" << to_type << std::endl;
      }
    }

    if (genes.size() == 1){
      unambig_ecs.push_back(ecs[i]);
    }
  }
  return(unambig_ecs);
}

//[[Rcpp::export]]
CharacterVector GeneToECMapC(CharacterVector gene, bool ambig, bool ensg, bool verbose,
                             Environment ec_to_enst, Environment enst_to_ec,
                             Environment enst_to_ensg, Environment ensg_to_enst,
                             Environment ensg_to_gene, Environment gene_to_ensg)
{
  std::string from_type = "GENE";
  if (ensg) {
    from_type = "ENSG";
  }
  CharacterVector ecs = TCCMapC(gene, from_type, "EC", ec_to_enst, enst_to_ec,
                                enst_to_ensg, ensg_to_enst, ensg_to_gene,
                                gene_to_ensg);
  if(ambig) {
    return(ecs);
  }
  CharacterVector unambig_ecs = ECUniqueGene(ecs, ensg, verbose, ec_to_enst, 
                                             enst_to_ec, enst_to_ensg, 
                                             ensg_to_enst, ensg_to_gene, 
                                             gene_to_ensg);
  return(unambig_ecs);
}

//[[Rcpp::export]]
List GenesToECMap(CharacterVector genes, bool ambig, bool ensg, bool verbose,
                  Environment ec_to_enst, Environment enst_to_ec,
                  Environment enst_to_ensg, Environment ensg_to_enst,
                  Environment ensg_to_gene, Environment gene_to_ensg)
{
  List mappings(genes.size());
  for(int i=0; i<genes.size(); ++i){
    CharacterVector gene = Rcpp::as<CharacterVector>(genes[i]);
    mappings[i] = GeneToECMapC(gene, ambig, ensg, verbose, ec_to_enst, enst_to_ec,
                  enst_to_ensg, ensg_to_enst, ensg_to_gene, gene_to_ensg);
  }
  return(mappings);
  
}
