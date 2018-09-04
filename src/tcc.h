#ifndef TCC
#define TCC

#include <RcppEigen.h>
#include "data_manipulation.h"
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>

using namespace Rcpp;

//----------------------------------------------------
StringVector TCCMapC(StringVector from, std::string from_type,
                     std::string to, Environment ec_to_enst,
                     Environment enst_to_ec, Environment enst_to_ensg,
                     Environment ensg_to_enst, Environment ensg_to_gene,
                     Environment gene_to_ensg);
CharacterVector ECUniqueGene(CharacterVector ecs, Environment ec_to_enst,
                             Environment enst_to_ec, Environment enst_to_ensg,
                             Environment ensg_to_enst, Environment ensg_to_gene,
                             Environment gene_to_ensg);
CharacterVector GeneToECMapC(CharacterVector gene, bool ambig, bool ensg,
                             Environment ec_to_enst, Environment enst_to_ec,
                             Environment enst_to_ensg, Environment ensg_to_enst,
                             Environment ensg_to_gene, Environment gene_to_ensg);
//----------------------------------------------------

#endif//TCC
