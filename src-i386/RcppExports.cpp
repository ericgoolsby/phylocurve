// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// C_ultrafast
arma::mat C_ultrafast(int nedge, arma::mat L, arma::mat E, arma::vec pic_ace, arma::vec ace_hat, arma::vec var_hat, arma::vec len_vec, arma::vec pic_len_vec);
RcppExport SEXP phylocurve_C_ultrafast(SEXP nedgeSEXP, SEXP LSEXP, SEXP ESEXP, SEXP pic_aceSEXP, SEXP ace_hatSEXP, SEXP var_hatSEXP, SEXP len_vecSEXP, SEXP pic_len_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type nedge(nedgeSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type E(ESEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pic_ace(pic_aceSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type ace_hat(ace_hatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type var_hat(var_hatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type len_vec(len_vecSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type pic_len_vec(pic_len_vecSEXP);
    __result = Rcpp::wrap(C_ultrafast(nedge, L, E, pic_ace, ace_hat, var_hat, len_vec, pic_len_vec));
    return __result;
END_RCPP
}
// multipic
List multipic(int ntip, int nnode, arma::vec edge1, arma::vec edge2, arma::vec edge_len, arma::mat phe, arma::mat contr, arma::vec var_contr, int scaled, int pic_len, int pic_recon);
RcppExport SEXP phylocurve_multipic(SEXP ntipSEXP, SEXP nnodeSEXP, SEXP edge1SEXP, SEXP edge2SEXP, SEXP edge_lenSEXP, SEXP pheSEXP, SEXP contrSEXP, SEXP var_contrSEXP, SEXP scaledSEXP, SEXP pic_lenSEXP, SEXP pic_reconSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type ntip(ntipSEXP);
    Rcpp::traits::input_parameter< int >::type nnode(nnodeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type edge1(edge1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type edge2(edge2SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type edge_len(edge_lenSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type phe(pheSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type contr(contrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type var_contr(var_contrSEXP);
    Rcpp::traits::input_parameter< int >::type scaled(scaledSEXP);
    Rcpp::traits::input_parameter< int >::type pic_len(pic_lenSEXP);
    Rcpp::traits::input_parameter< int >::type pic_recon(pic_reconSEXP);
    __result = Rcpp::wrap(multipic(ntip, nnode, edge1, edge2, edge_len, phe, contr, var_contr, scaled, pic_len, pic_recon));
    return __result;
END_RCPP
}
// multipic2
List multipic2(int ntip, int nnode, arma::vec edge1, arma::vec edge2, arma::mat edge_len, arma::mat phe, arma::mat contr, arma::mat var_contr, int scaled, int pic_len, int pic_recon);
RcppExport SEXP phylocurve_multipic2(SEXP ntipSEXP, SEXP nnodeSEXP, SEXP edge1SEXP, SEXP edge2SEXP, SEXP edge_lenSEXP, SEXP pheSEXP, SEXP contrSEXP, SEXP var_contrSEXP, SEXP scaledSEXP, SEXP pic_lenSEXP, SEXP pic_reconSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type ntip(ntipSEXP);
    Rcpp::traits::input_parameter< int >::type nnode(nnodeSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type edge1(edge1SEXP);
    Rcpp::traits::input_parameter< arma::vec >::type edge2(edge2SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type edge_len(edge_lenSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type phe(pheSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type contr(contrSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type var_contr(var_contrSEXP);
    Rcpp::traits::input_parameter< int >::type scaled(scaledSEXP);
    Rcpp::traits::input_parameter< int >::type pic_len(pic_lenSEXP);
    Rcpp::traits::input_parameter< int >::type pic_recon(pic_reconSEXP);
    __result = Rcpp::wrap(multipic2(ntip, nnode, edge1, edge2, edge_len, phe, contr, var_contr, scaled, pic_len, pic_recon));
    return __result;
END_RCPP
}
// logdet
double logdet(arma::mat A);
RcppExport SEXP phylocurve_logdet(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    __result = Rcpp::wrap(logdet(A));
    return __result;
END_RCPP
}
// Rinv4
double Rinv4(unsigned int nspecies, int nedge, int ngroups, arma::uvec ind, arma::vec len_vec, arma::uvec anc, arma::uvec des, arma::vec R, arma::mat painted_edges, arma::mat phylocovs, int REML);
RcppExport SEXP phylocurve_Rinv4(SEXP nspeciesSEXP, SEXP nedgeSEXP, SEXP ngroupsSEXP, SEXP indSEXP, SEXP len_vecSEXP, SEXP ancSEXP, SEXP desSEXP, SEXP RSEXP, SEXP painted_edgesSEXP, SEXP phylocovsSEXP, SEXP REMLSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< unsigned int >::type nspecies(nspeciesSEXP);
    Rcpp::traits::input_parameter< int >::type nedge(nedgeSEXP);
    Rcpp::traits::input_parameter< int >::type ngroups(ngroupsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type len_vec(len_vecSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type anc(ancSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type des(desSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type painted_edges(painted_edgesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type phylocovs(phylocovsSEXP);
    Rcpp::traits::input_parameter< int >::type REML(REMLSEXP);
    __result = Rcpp::wrap(Rinv4(nspecies, nedge, ngroups, ind, len_vec, anc, des, R, painted_edges, phylocovs, REML));
    return __result;
END_RCPP
}
// Rinv6
double Rinv6(unsigned int nspecies, int nedge, int ngroups, arma::uvec ind, arma::vec len_vec, arma::uvec anc, arma::uvec des, arma::mat L, arma::vec R, arma::mat painted_edges, arma::mat phylocovs, int REML);
RcppExport SEXP phylocurve_Rinv6(SEXP nspeciesSEXP, SEXP nedgeSEXP, SEXP ngroupsSEXP, SEXP indSEXP, SEXP len_vecSEXP, SEXP ancSEXP, SEXP desSEXP, SEXP LSEXP, SEXP RSEXP, SEXP painted_edgesSEXP, SEXP phylocovsSEXP, SEXP REMLSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< unsigned int >::type nspecies(nspeciesSEXP);
    Rcpp::traits::input_parameter< int >::type nedge(nedgeSEXP);
    Rcpp::traits::input_parameter< int >::type ngroups(ngroupsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type len_vec(len_vecSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type anc(ancSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type des(desSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type painted_edges(painted_edgesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type phylocovs(phylocovsSEXP);
    Rcpp::traits::input_parameter< int >::type REML(REMLSEXP);
    __result = Rcpp::wrap(Rinv6(nspecies, nedge, ngroups, ind, len_vec, anc, des, L, R, painted_edges, phylocovs, REML));
    return __result;
END_RCPP
}
// theta_Rinv6
arma::mat theta_Rinv6(unsigned int nspecies, int nedge, int ngroups, arma::uvec ind, arma::vec len_vec, arma::uvec anc, arma::uvec des, arma::mat L, arma::vec R, arma::mat painted_edges, arma::mat phylocovs, int REML);
RcppExport SEXP phylocurve_theta_Rinv6(SEXP nspeciesSEXP, SEXP nedgeSEXP, SEXP ngroupsSEXP, SEXP indSEXP, SEXP len_vecSEXP, SEXP ancSEXP, SEXP desSEXP, SEXP LSEXP, SEXP RSEXP, SEXP painted_edgesSEXP, SEXP phylocovsSEXP, SEXP REMLSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< unsigned int >::type nspecies(nspeciesSEXP);
    Rcpp::traits::input_parameter< int >::type nedge(nedgeSEXP);
    Rcpp::traits::input_parameter< int >::type ngroups(ngroupsSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type ind(indSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type len_vec(len_vecSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type anc(ancSEXP);
    Rcpp::traits::input_parameter< arma::uvec >::type des(desSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type L(LSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type R(RSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type painted_edges(painted_edgesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type phylocovs(phylocovsSEXP);
    Rcpp::traits::input_parameter< int >::type REML(REMLSEXP);
    __result = Rcpp::wrap(theta_Rinv6(nspecies, nedge, ngroups, ind, len_vec, anc, des, L, R, painted_edges, phylocovs, REML));
    return __result;
END_RCPP
}
// fast_transform
arma::mat fast_transform(arma::mat m);
RcppExport SEXP phylocurve_fast_transform(SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< arma::mat >::type m(mSEXP);
    __result = Rcpp::wrap(fast_transform(m));
    return __result;
END_RCPP
}
