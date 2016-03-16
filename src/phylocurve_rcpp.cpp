// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat C_ultrafast(int nedge,arma::mat L,arma::mat E,arma::vec pic_ace,arma::vec ace_hat,arma::vec var_hat,arma::vec len_vec,arma::vec pic_len_vec)
{
  arma::vec tempL = arma::zeros(2);
  int root_node = 0;
  double root_val = 0;
  double var_val = 0;
  double root_edge = 0;
  arma::vec es = arma::zeros(2);
  arma::vec ds = arma::zeros(2);
  arma::mat ret = arma::zeros(ace_hat.size(),2);
  int i = 0;
  for(i=nedge-1;i>=0;i--)
  {
    tempL(0) = L(i,0)-1;
    tempL(1) = L(i,1)-1;
    if(tempL(0)!=tempL(1))
    {
      root_node = E(i,0)-1;
      root_val = ace_hat(root_node);
      var_val = var_hat(root_node);
      root_edge = len_vec(i);
      es(0) = pic_len_vec(tempL(0));
      es(1) = pic_len_vec(tempL(1));
      ds(0) = pic_ace(E(tempL(0),1)-1);
      ds(1) = pic_ace(E(tempL(1),1)-1);
      ace_hat(E(i,1)-1) = (es(0)*es(1)*root_val + root_edge*es(1)*ds(0) + root_edge*es(0)*ds(1)) / (es(0)*es(1) + root_edge*es(1) + root_edge*es(0));
      var_hat(E(i,1)-1) = 1/((-(pic_len_vec(i)*pic_len_vec(i)))/((root_edge-pic_len_vec(i))*(var_val*pic_len_vec(i)+root_edge*pic_len_vec(i)-root_edge*var_val)));
    }
  }
  ret.col(0) = ace_hat;
  ret.col(1) = var_hat;
  return ret;
}

// [[Rcpp::export]]
List multipic(int ntip, int nnode, arma::vec edge1, arma::vec edge2,
              arma::vec edge_len, arma::mat phe, arma::mat contr,
              arma::vec var_contr, int scaled, int pic_len, int pic_recon)
{
  /*  Code adapted from ape package
      Paradis, E. (2012) Analysis of Phylogenetics and
      Evolution with R (Second Edition). New York: Springer.
   
      Paradis, E., Claude, J. and Strimmer, K. (2004) APE:
      analyses of phylogenetics and evolution in R language.
      Bioinformatics, 20, 289–290.
   */
  /* The tree must be in pruningwise order */
  int anc, d1, d2, ic, i, j, k;
  double sumbl;
  
  for (i = 0; i < ntip * 2 - 3; i += 2) {
    j = i + 1;
    anc = edge1(i);
    d1 = edge2(i) - 1;
    d2 = edge2(j) - 1;
    sumbl = edge_len(i) + edge_len(j);
    ic = anc - ntip - 1;
    contr.row(ic) = phe.row(d1) - phe.row(d2);
    if (scaled==1) contr.row(ic) = contr.row(ic)/sqrt(sumbl);
    var_contr(ic) = sumbl;
    phe.row(anc - 1) = (phe.row(d1)*edge_len(j) + phe.row(d2)*edge_len(i))/sumbl;
    /* find the edge where `anc' is a descendant (except if at the root):
    it is obviously below the j'th edge */
    if (j != ntip * 2 - 3) {
      k = j + 1;
      while (edge2(k) != anc) k++;
      edge_len(k) += edge_len(i)*edge_len(j)/sumbl;
    }
  }
  double sum_invV = (edge_len(ntip+nnode-2)+edge_len(ntip+nnode-3))/(edge_len(ntip+nnode-2)*edge_len(ntip+nnode-3));
  double log_detV = log(1/sum_invV) + sum(log(var_contr));
  arma::rowvec root = phe.row(ntip);
  //Rcout << sum_invV << " " << log_detV;
  if((pic_len==1) & (pic_recon==0))
  {
    return List::create(_["contrasts"] = contr,
                        _["sum_invV"] = sum_invV,
                        _["log_detV"] = log_detV,
                        _["root"] = root,
                        _["pic_len"] = edge_len);
  } else if(pic_recon==1)
  {
    return List::create(_["contrasts"] = contr,
                        _["sum_invV"] = sum_invV,
                        _["log_detV"] = log_detV,
                        _["root"] = root,
                        _["pic_len"] = edge_len,
                        _["pic_recon"] = phe);
  }
  return List::create(_["contrasts"] = contr,
                      _["sum_invV"] = sum_invV,
                      _["log_detV"] = log_detV,
                      _["root"] = root);
}

// [[Rcpp::export]]
List multipic2(int ntip, int nnode, arma::vec edge1, arma::vec edge2,
               arma::mat edge_len, arma::mat phe, arma::mat contr,
               arma::mat var_contr, int scaled, int pic_len, int pic_recon)
{
  /*  Code adapted from ape package
      Paradis, E. (2012) Analysis of Phylogenetics and
      Evolution with R (Second Edition). New York: Springer.
  
      Paradis, E., Claude, J. and Strimmer, K. (2004) APE:
      analyses of phylogenetics and evolution in R language.
      Bioinformatics, 20, 289–290.
  */
  
  /* The tree must be in pruningwise order */
  int anc, d1, d2, ic, i, j, k;
  arma::rowvec sumbl = edge_len.row(0);
  
  for (i = 0; i < ntip * 2 - 3; i += 2) {
    j = i + 1;
    anc = edge1(i);
    d1 = edge2(i) - 1;
    d2 = edge2(j) - 1;
    sumbl = edge_len.row(i) + edge_len.row(j);
    ic = anc - ntip - 1;
    contr.row(ic) = phe.row(d1) - phe.row(d2);
    if (scaled==1) contr.row(ic) = contr.row(ic)/sqrt(sumbl);
    var_contr.row(ic) = sumbl;
    phe.row(anc - 1) = (phe.row(d1)%edge_len.row(j) + phe.row(d2)%edge_len.row(i))/sumbl;
    /* find the edge where `anc' is a descendant (except if at the root):
    it is obviously below the j'th edge */
    if (j != ntip * 2 - 3) {
      k = j + 1;
      while (edge2(k) != anc) k++;
      edge_len.row(k) += edge_len.row(i)%edge_len.row(j)/sumbl;
    }
  }
  arma::rowvec sum_invV = (edge_len.row(ntip+nnode-2)+edge_len.row(ntip+nnode-3))/(edge_len.row(ntip+nnode-2)%edge_len.row(ntip+nnode-3));
  arma::rowvec log_detV = log(1/sum_invV) + sum(log(var_contr),0);
  arma::rowvec root = phe.row(ntip);
  //Rcout << sum_invV << " " << log_detV;
  if((pic_len==1) & (pic_recon==0))
  {
    return List::create(_["contrasts"] = contr,
                        _["sum_invV"] = sum_invV,
                        _["log_detV"] = log_detV,
                        _["root"] = root,
                        _["pic_len"] = edge_len);
  } else if(pic_recon==1)
  {
    return List::create(_["contrasts"] = contr,
                        _["sum_invV"] = sum_invV,
                        _["log_detV"] = log_detV,
                        _["root"] = root,
                        _["pic_len"] = edge_len,
                        _["pic_recon"] = phe);
  }
  return List::create(_["contrasts"] = contr,
                      _["sum_invV"] = sum_invV,
                      _["log_detV"] = log_detV,
                      _["root"] = root);
}

// [[Rcpp::export]]
double logdet(arma::mat A)
{
  double val = 0;
  double sign = 1;
  log_det(val,sign,A);
  return val;
}

// [[Rcpp::export]]
double Rinv4(unsigned int nspecies,int nedge,int ngroups,arma::uvec ind,arma::vec len_vec,
             arma::uvec anc,arma::uvec des,arma::vec R,arma::mat painted_edges,arma::mat phylocovs,int REML)
{
  int nvar = 2;
  arma::mat RR = arma::zeros(1,nedge+1);
  arma::mat R1 = arma::zeros(nvar,nedge+1);
  arma::vec logd = arma::zeros(nedge+1);
  arma::mat p = arma::zeros(nvar,nvar*(nedge+1));
  arma::mat pA = arma::zeros(nvar,nvar*(nedge+1));
  arma::mat pA_des_value = arma::zeros(nvar,nvar);
  arma::mat len = arma::zeros(nvar,nvar);
  arma::mat id = arma::eye(nvar,nvar);
  arma::mat itpa = arma::zeros(nvar,nvar);
  arma::mat itpainv = arma::zeros(nvar,nvar);
  int i=0;
  int z=0;
  unsigned int des_node=0;
  unsigned int anc_node=0;
  for(i=0;i<nedge;i++)
  {
    des_node = des[i];
    anc_node = anc[i];
    len = len*0;
    for(z=0;z<ngroups;z++)
    {
      len = len + painted_edges(i,z)*len_vec(i)*phylocovs.rows(nvar*z,nvar*(z+1)-1);
    }
    if(des_node<nspecies)
    {
      logd(des_node) = logdet(len);
      logd(anc_node) += logd(des_node);
      p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1) = inv(len);
      pA.submat(0,anc_node*nvar,nvar-1,(anc_node+1)*nvar-1) += p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1);
      R1.col(des_node) = p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1) * R.elem(des_node+nspecies*ind);
      RR.col(des_node) = trans(R.elem(des_node+nspecies*ind))*R1.col(des_node);
      R1.col(anc_node) += R1.col(des_node);
      RR(anc_node) += RR(des_node);
    } else
    {
      pA_des_value = pA.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1);
      itpa = id + len*pA_des_value;
      itpainv = inv(itpa);
      logd(des_node) += logdet(itpa);
      logd(anc_node) += logd(des_node);
      p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1) = pA_des_value * itpainv;
      pA.submat(0,anc_node*nvar,nvar-1,(anc_node+1)*nvar-1) += p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1);
      RR.col(des_node) += (-trans(R1.col(des_node)) * itpainv * len * R1.col(des_node));
      RR(anc_node) += RR(des_node);
      R1.col(des_node) = trans(trans(R1.col(des_node)) * itpainv);
      R1.col(anc_node) += R1.col(des_node);
    }
  }
  i = nedge;
  len = len*0;
  for(z=0;z<ngroups;z++)
  {
    len = len + painted_edges(i,z)*len_vec(i)*phylocovs.rows(nvar*z,nvar*(z+1)-1);
  }
  pA_des_value = pA.submat(0,anc_node*nvar,nvar-1,(anc_node+1)*nvar-1);
  itpa = id + len*pA_des_value;
  itpainv = inv(itpa);
  logd(nedge) = logd(nspecies) + logdet(itpa);
  p.submat(0,(nedge)*nvar,nvar-1,((nedge)+1)*nvar-1) = pA_des_value * itpainv;
  RR.col(nedge) = RR(anc_node) - (trans(R1.col(anc_node)) * itpainv * len * R1.col(anc_node));
  R1.col(nedge) = trans(trans(R1.col(anc_node)) * itpainv);
  arma::mat pfinal = p.submat(0,(nedge)*nvar,nvar-1,((nedge)+1)*nvar-1);
  arma::mat theta = solve(pfinal,R1.col(nedge));
  double loglik = 0;
  arma::mat ret1 = 2*trans(theta)*R1.col(nedge);
  arma::mat ret2 = trans(theta)*pfinal*theta;
  loglik = (nspecies*nvar-(nvar*REML))*log(2*3.1415926535898) + logd(nedge) + RR(0,nedge) - ret1(0,0) + ret2(0,0);
  if(REML==1) loglik += logdet(pfinal);
  loglik /= -2;
  return loglik;
}

// [[Rcpp::export]]
double Rinv6(unsigned int nspecies,int nedge,int ngroups,arma::uvec ind,arma::vec len_vec,
             arma::uvec anc,arma::uvec des,arma::mat L,arma::vec R,arma::mat painted_edges,arma::mat phylocovs,int REML)
{
  int nvar = phylocovs.n_cols;
  arma::mat RR = arma::zeros(1,nedge+1);
  arma::mat R1 = arma::zeros(nvar,nedge+1);
  
  int mL = L.n_cols;
  arma::mat LL = arma::zeros(mL*nvar,mL*nvar*(nedge+1));
  arma::mat L1 = arma::zeros(nvar*mL,nvar*(nedge+1));
  arma::mat LR = arma::zeros(mL*nvar,nedge+1);
  
  arma::vec logd = arma::zeros(nedge+1);
  arma::mat p = arma::zeros(nvar,nvar*(nedge+1));
  arma::mat pA = arma::zeros(nvar,nvar*(nedge+1));
  arma::mat pA_des_value = arma::zeros(nvar,nvar);
  arma::mat len = arma::zeros(nvar,nvar);
  arma::mat id = arma::eye(nvar,nvar);
  arma::mat itpa = arma::zeros(nvar,nvar);
  arma::mat itpainv = arma::zeros(nvar,nvar);
  int i=0;
  int z=0;
  unsigned int des_node=0;
  unsigned int anc_node=0;
  for(i=0;i<nedge;i++)
  {
    des_node = des[i];
    anc_node = anc[i];
    len = len*0;
    for(z=0;z<ngroups;z++)
    {
      len = len + painted_edges(i,z)*len_vec(i)*phylocovs.rows(nvar*z,nvar*(z+1)-1);
    }
    if(des_node<nspecies)
    {
      logd(des_node) = logdet(len);
      logd(anc_node) += logd(des_node);
      p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1) = inv(len);
      pA.submat(0,anc_node*nvar,nvar-1,(anc_node+1)*nvar-1) += p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1);
      R1.col(des_node) = p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1) * R.elem(des_node+nspecies*ind);
      RR.col(des_node) = trans(R.elem(des_node+nspecies*ind))*R1.col(des_node);
      R1.col(anc_node) += R1.col(des_node);
      RR(anc_node) += RR(des_node);
      L1.cols(des_node*nvar,(des_node+1)*(nvar)-1) = kron(p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1),trans(L.row(des_node)));
      LL.cols(des_node*(nvar*mL),(des_node+1)*(nvar*mL)-1) = kron(L1.cols(des_node*nvar,(des_node+1)*(nvar)-1),L.row(des_node));
      L1.cols(anc_node*nvar,(anc_node+1)*(nvar)-1) += L1.cols(des_node*nvar,(des_node+1)*(nvar)-1);
      LL.cols(anc_node*(nvar*mL),(anc_node+1)*(nvar*mL)-1) += LL.cols(des_node*(nvar*mL),(des_node+1)*(nvar*mL)-1);
      LR.col(des_node) = L1.cols(des_node*nvar,(des_node+1)*(nvar)-1) * R.elem(des_node+nspecies*ind);
      LR.col(anc_node) += LR.col(des_node);
    } else
    {
      pA_des_value = pA.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1);
      itpa = id + len*pA_des_value;
      itpainv = inv(itpa);
      logd(des_node) += logdet(itpa);
      logd(anc_node) += logd(des_node);
      p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1) = pA_des_value * itpainv;
      pA.submat(0,anc_node*nvar,nvar-1,(anc_node+1)*nvar-1) += p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1);
      LR.col(des_node) += (-L1.cols(des_node*nvar,(des_node+1)*(nvar)-1) * itpainv * len * R1.col(des_node));
      LR.col(anc_node) += LR.col(des_node);
      LL.cols(des_node*(nvar*mL),(des_node+1)*(nvar*mL)-1) += (-L1.cols(des_node*nvar,(des_node+1)*(nvar)-1)*itpainv*len*trans(L1.cols(des_node*nvar,(des_node+1)*(nvar)-1))); 
      LL.cols(anc_node*(nvar*mL),(anc_node+1)*(nvar*mL)-1) += LL.cols(des_node*(nvar*mL),(des_node+1)*(nvar*mL)-1);
      L1.cols(des_node*nvar,(des_node+1)*(nvar)-1) = L1.cols(des_node*nvar,(des_node+1)*(nvar)-1) * itpainv;
      L1.cols(anc_node*nvar,(anc_node+1)*(nvar)-1) += L1.cols(des_node*nvar,(des_node+1)*(nvar)-1);
      RR.col(des_node) += (-trans(R1.col(des_node)) * itpainv * len * R1.col(des_node));
      RR(anc_node) += RR(des_node);
      R1.col(des_node) = trans(trans(R1.col(des_node)) * itpainv);
      R1.col(anc_node) += R1.col(des_node);
    }
  }
  i = nedge;
  len = len*0;
  for(z=0;z<ngroups;z++)
  {
    len = len + painted_edges(i,z)*len_vec(i)*phylocovs.rows(nvar*z,nvar*(z+1)-1);
  }
  pA_des_value = pA.submat(0,anc_node*nvar,nvar-1,(anc_node+1)*nvar-1);
  itpa = id + len*pA_des_value;
  itpainv = inv(itpa);
  logd(nedge) = logd(nspecies) + logdet(itpa);
  p.submat(0,(nedge)*nvar,nvar-1,((nedge)+1)*nvar-1) = pA_des_value * itpainv;
  
  /*
  LR.col(nedge) = LR.col(anc_node) - (L1.cols(anc_node*(nvar-1),anc_node*(nvar)-1) * itpainv * len * R1.col(anc_node));
  LL.cols(nedge*((nvar*mL)-1),nedge*(nvar*mL)-1) = LL.cols(anc_node*((nvar*mL)-1),anc_node*(nvar*mL)-1) - 
  (L1.cols(anc_node*(nvar-1),anc_node*(nvar)-1)*itpainv*len*trans(L1.cols(anc_node*(nvar-1),anc_node*(nvar)-1))); 
  L1.cols(nedge*(nvar-1),nedge*(nvar)-1) = L1.cols(anc_node*(nvar-1),anc_node*(nvar)-1) * itpainv;
  */
  LR.col(nedge) = LR.col(anc_node) - (L1.cols(anc_node*nvar,(anc_node+1)*(nvar)-1) * itpainv * len * R1.col(anc_node));
  LL.cols(nedge*(nvar*mL),(nedge+1)*(nvar*mL)-1) = LL.cols(anc_node*(nvar*mL),(anc_node+1)*(nvar*mL)-1) -
    (L1.cols(anc_node*nvar,(anc_node+1)*(nvar)-1)*itpainv*len*trans(L1.cols(anc_node*nvar,(anc_node+1)*(nvar)-1))); 
  L1.cols(nedge*nvar,(nedge+1)*(nvar)-1) = L1.cols(anc_node*nvar,(anc_node+1)*(nvar)-1) * itpainv;
  
  RR.col(nedge) = RR(anc_node) - (trans(R1.col(anc_node)) * itpainv * len * R1.col(anc_node));
  R1.col(nedge) = trans(trans(R1.col(anc_node)) * itpainv);
  arma::mat pfinal = p.submat(0,(nedge)*nvar,nvar-1,((nedge)+1)*nvar-1);
  arma::mat LLfinal = LL.cols(nedge*(nvar*mL),(nedge+1)*(nvar*mL)-1);
  arma::mat theta = solve(LLfinal,LR.col(nedge));
  double loglik = 0;
  arma::mat ret1 = 2*trans(theta)*LR.col(nedge);
  arma::mat ret2 = trans(theta)*LLfinal*theta;
  loglik = (nspecies*nvar-(nvar*REML))*log(2*3.1415926535898) + logd(nedge) + RR(0,nedge) - ret1(0,0) + ret2(0,0);
  if(REML==1) loglik += logdet(LLfinal);
  loglik /= -2;
  return loglik;
}

// [[Rcpp::export]]
arma::mat theta_Rinv6(unsigned int nspecies,int nedge,int ngroups,arma::uvec ind,arma::vec len_vec,
                      arma::uvec anc,arma::uvec des,arma::mat L,arma::vec R,arma::mat painted_edges,arma::mat phylocovs,int REML)
{
  int nvar = phylocovs.n_cols;
  arma::mat RR = arma::zeros(1,nedge+1);
  arma::mat R1 = arma::zeros(nvar,nedge+1);
  
  int mL = L.n_cols;
  arma::mat LL = arma::zeros(mL*nvar,mL*nvar*(nedge+1));
  arma::mat L1 = arma::zeros(nvar*mL,nvar*(nedge+1));
  arma::mat LR = arma::zeros(mL*nvar,nedge+1);
  
  arma::vec logd = arma::zeros(nedge+1);
  arma::mat p = arma::zeros(nvar,nvar*(nedge+1));
  arma::mat pA = arma::zeros(nvar,nvar*(nedge+1));
  arma::mat pA_des_value = arma::zeros(nvar,nvar);
  arma::mat len = arma::zeros(nvar,nvar);
  arma::mat id = arma::eye(nvar,nvar);
  arma::mat itpa = arma::zeros(nvar,nvar);
  arma::mat itpainv = arma::zeros(nvar,nvar);
  int i=0;
  int z=0;
  unsigned int des_node=0;
  unsigned int anc_node=0;
  for(i=0;i<nedge;i++)
  {
    des_node = des[i];
    anc_node = anc[i];
    len = len*0;
    for(z=0;z<ngroups;z++)
    {
      len = len + painted_edges(i,z)*len_vec(i)*phylocovs.rows(nvar*z,nvar*(z+1)-1);
    }
    if(des_node<nspecies)
    {
      logd(des_node) = logdet(len);
      logd(anc_node) += logd(des_node);
      p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1) = inv(len);
      pA.submat(0,anc_node*nvar,nvar-1,(anc_node+1)*nvar-1) += p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1);
      R1.col(des_node) = p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1) * R.elem(des_node+nspecies*ind);
      RR.col(des_node) = trans(R.elem(des_node+nspecies*ind))*R1.col(des_node);
      R1.col(anc_node) += R1.col(des_node);
      RR(anc_node) += RR(des_node);
      L1.cols(des_node*nvar,(des_node+1)*(nvar)-1) = kron(p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1),trans(L.row(des_node)));
      LL.cols(des_node*(nvar*mL),(des_node+1)*(nvar*mL)-1) = kron(L1.cols(des_node*nvar,(des_node+1)*(nvar)-1),L.row(des_node));
      L1.cols(anc_node*nvar,(anc_node+1)*(nvar)-1) += L1.cols(des_node*nvar,(des_node+1)*(nvar)-1);
      LL.cols(anc_node*(nvar*mL),(anc_node+1)*(nvar*mL)-1) += LL.cols(des_node*(nvar*mL),(des_node+1)*(nvar*mL)-1);
      LR.col(des_node) = L1.cols(des_node*nvar,(des_node+1)*(nvar)-1) * R.elem(des_node+nspecies*ind);
      LR.col(anc_node) += LR.col(des_node);
    } else
    {
      pA_des_value = pA.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1);
      itpa = id + len*pA_des_value;
      itpainv = inv(itpa);
      logd(des_node) += logdet(itpa);
      logd(anc_node) += logd(des_node);
      p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1) = pA_des_value * itpainv;
      pA.submat(0,anc_node*nvar,nvar-1,(anc_node+1)*nvar-1) += p.submat(0,des_node*nvar,nvar-1,(des_node+1)*nvar-1);
      LR.col(des_node) += (-L1.cols(des_node*nvar,(des_node+1)*(nvar)-1) * itpainv * len * R1.col(des_node));
      LR.col(anc_node) += LR.col(des_node);
      LL.cols(des_node*(nvar*mL),(des_node+1)*(nvar*mL)-1) += (-L1.cols(des_node*nvar,(des_node+1)*(nvar)-1)*itpainv*len*trans(L1.cols(des_node*nvar,(des_node+1)*(nvar)-1))); 
      LL.cols(anc_node*(nvar*mL),(anc_node+1)*(nvar*mL)-1) += LL.cols(des_node*(nvar*mL),(des_node+1)*(nvar*mL)-1);
      L1.cols(des_node*nvar,(des_node+1)*(nvar)-1) = L1.cols(des_node*nvar,(des_node+1)*(nvar)-1) * itpainv;
      L1.cols(anc_node*nvar,(anc_node+1)*(nvar)-1) += L1.cols(des_node*nvar,(des_node+1)*(nvar)-1);
      RR.col(des_node) += (-trans(R1.col(des_node)) * itpainv * len * R1.col(des_node));
      RR(anc_node) += RR(des_node);
      R1.col(des_node) = trans(trans(R1.col(des_node)) * itpainv);
      R1.col(anc_node) += R1.col(des_node);
    }
  }
  i = nedge;
  len = len*0;
  for(z=0;z<ngroups;z++)
  {
    len = len + painted_edges(i,z)*len_vec(i)*phylocovs.rows(nvar*z,nvar*(z+1)-1);
  }
  pA_des_value = pA.submat(0,anc_node*nvar,nvar-1,(anc_node+1)*nvar-1);
  itpa = id + len*pA_des_value;
  itpainv = inv(itpa);
  logd(nedge) = logd(nspecies) + logdet(itpa);
  p.submat(0,(nedge)*nvar,nvar-1,((nedge)+1)*nvar-1) = pA_des_value * itpainv;
  
  LR.col(nedge) = LR.col(anc_node) - (L1.cols(anc_node*nvar,(anc_node+1)*(nvar)-1) * itpainv * len * R1.col(anc_node));
  LL.cols(nedge*(nvar*mL),(nedge+1)*(nvar*mL)-1) = LL.cols(anc_node*(nvar*mL),(anc_node+1)*(nvar*mL)-1) -
    (L1.cols(anc_node*nvar,(anc_node+1)*(nvar)-1)*itpainv*len*trans(L1.cols(anc_node*nvar,(anc_node+1)*(nvar)-1))); 
  L1.cols(nedge*nvar,(nedge+1)*(nvar)-1) = L1.cols(anc_node*nvar,(anc_node+1)*(nvar)-1) * itpainv;
  
  RR.col(nedge) = RR(anc_node) - (trans(R1.col(anc_node)) * itpainv * len * R1.col(anc_node));
  R1.col(nedge) = trans(trans(R1.col(anc_node)) * itpainv);
  arma::mat pfinal = p.submat(0,(nedge)*nvar,nvar-1,((nedge)+1)*nvar-1);
  arma::mat LLfinal = LL.cols(nedge*(nvar*mL),(nedge+1)*(nvar*mL)-1);
  arma::mat theta = solve(LLfinal,LR.col(nedge));
  return theta;
}

// [[Rcpp::export]]
arma::mat fast_transform(arma::mat m)
{
  unsigned int i = 0;
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, m, "dc");
  
  eigval = (eigval+abs(eigval))/2;
  arma::vec eigval2 = 1/sqrt(eigval);
  arma::uvec subst = arma::find(eigval==0);
  unsigned int s = subst.size();
  for(i=0;i<s;i++)
  {
    eigval2(i) = 0;
  }
  return (eigvec * diagmat(eigval2) * eigvec.t());
}
