/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2018 Robert McCulloch and Rodney Sparapani
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */

#include <BART3.h>
#include "rmvnorm.h"

#ifndef NoRcpp

#define TRDRAW(a, b) trdraw(a, b)
#define TEDRAW(a, b) tedraw(a, b)
#define BDRAW(a, b) bdraw(a, b)

RcppExport SEXP cplbart(
   SEXP _in,            //number of observations in training data
   SEXP _ip,		//dimension of x
   SEXP _il,            //dimension of linear x
   SEXP _inp,		//number of observations in test data
   SEXP _ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   SEXP _ixl,           //x, train, (p+1)xn, added 1 for linear beta_0
   SEXP _iy,		//y, train,  nx1
   SEXP _ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   SEXP _ixlp,          //x, test, (p+1)xnp, added 1 for linear beta_0
   SEXP _im,		//number of trees
   SEXP _inc,		//number of cut points
   SEXP _ind,		//number of kept draws (except for thinnning ..)
   SEXP _iburn,		//number of burn-in draws skipped
   SEXP _ipower,
   SEXP _ibase,
   SEXP _imub,        //prior mean of beta
   SEXP _ig,       //g prior of beta
   SEXP _itau,
   SEXP _idart,
   SEXP _itheta,
   SEXP _iomega,
   SEXP _igrp,
   SEXP _ia,
   SEXP _ib,
   SEXP _irho,
   SEXP _iaug,
   SEXP _inkeeptrain,
   SEXP _inkeeptest,
//   SEXP _inkeeptestme,
   SEXP _inkeeptreedraws,
   SEXP _inprintevery,
//   SEXP _treesaslists,
   SEXP _Xinfo
)
{

   //--------------------------------------------------
   //process args
   size_t n = Rcpp::as<int>(_in);
   size_t p = Rcpp::as<int>(_ip);
   size_t l = Rcpp::as<int>(_il);
   size_t np = Rcpp::as<int>(_inp);
   Rcpp::NumericVector  xv(_ix);
   double *ix = &xv[0];
   Rcpp::NumericVector xlv(_ixl);
   double *ixl = &xlv[0];
   Rcpp::IntegerVector  yv(_iy); // binary
   int *iy = &yv[0];
   Rcpp::NumericVector  xpv(_ixp);
   double *ixp = &xpv[0];
   Rcpp::NumericVector xlpv(_ixlp);
   double *ixlp = &xlpv[0];
   size_t m = Rcpp::as<int>(_im);
   //size_t nc = Rcpp::as<int>(_inc);
   Rcpp::IntegerVector _nc(_inc);
   int *numcut = &_nc[0];
   size_t nd = Rcpp::as<int>(_ind);
   size_t burn = Rcpp::as<int>(_iburn);
   double mybeta = Rcpp::as<double>(_ipower);
   double alpha = Rcpp::as<double>(_ibase);
   Rcpp::NumericVector mub(_imub);
   double *imub = &mub[0];
   const double g = Rcpp::as<double>(_ig);

   double tau = Rcpp::as<double>(_itau);
//   double rootM = sqrt(Rcpp::as<double>(_iM));
   bool dart;
   if(Rcpp::as<int>(_idart)==1) dart=true;
   else dart=false;
   double a = Rcpp::as<double>(_ia);
   double b = Rcpp::as<double>(_ib);
   double rho = Rcpp::as<double>(_irho);
   bool aug;
   if(Rcpp::as<int>(_iaug)==1) aug=true;
   else aug=false;
   double theta = Rcpp::as<double>(_itheta);
   double omega = Rcpp::as<double>(_iomega);
   Rcpp::IntegerVector _grp(_igrp);
//   int *grp = &_grp[0];
   size_t nkeeptrain = Rcpp::as<int>(_inkeeptrain);
   size_t nkeeptest = Rcpp::as<int>(_inkeeptest);
//   size_t nkeeptestme = Rcpp::as<int>(_inkeeptestme);
   size_t nkeeptreedraws = Rcpp::as<int>(_inkeeptreedraws);
   size_t printevery = Rcpp::as<int>(_inprintevery);
//   int treesaslists = Rcpp::as<int>(_treesaslists);
   Rcpp::NumericMatrix Xinfo(_Xinfo);
   //Rcpp::List Xinfo(_Xinfo);
//   Rcpp::IntegerMatrix varcount(nkeeptreedraws, p);

   //return data structures (using Rcpp)
/*
   Rcpp::NumericVector trmean(n); //train
   Rcpp::NumericVector temean(np);
*/
   Rcpp::NumericMatrix trdraw(nkeeptrain,n);
   Rcpp::NumericMatrix tedraw(nkeeptest,np);
   Rcpp::NumericMatrix bdraw(nd+burn,l);
   Rcpp::NumericVector sdraw(nd+burn);
   Rcpp::NumericVector tadraw(nd+burn);
   Rcpp::NumericVector tbdraw(nd+burn);
//   Rcpp::List list_of_lists(nkeeptreedraws*treesaslists);
   Rcpp::NumericMatrix varprb(nkeeptreedraws,p);
   Rcpp::IntegerMatrix varcnt(nkeeptreedraws,p);

   //random number generation
   arn gen;

   bart bm(m);

   if(Xinfo.size()>0) {
     xinfo _xi;
     _xi.resize(p);
     for(size_t i=0;i<p;i++) {
       _xi[i].resize(numcut[i]);
       //Rcpp::IntegerVector cutpts(Xinfo[i]);
       for(size_t j=0;j<numcut[i];j++) _xi[i][j]=Xinfo(i, j);
     }
     bm.setxinfo(_xi);
   }
#else

#define TRDRAW(a, b) trdraw[a][b]
#define TEDRAW(a, b) tedraw[a][b]
#define BDRAW(a, b) bdraw[a][b]

void cplbart(
   size_t n,            //number of observations in training data
   size_t p,		//dimension of x
   size_t l,            //dimension of linear x
   size_t np,		//number of observations in test data
   double* ix,		//x, train,  pxn (transposed so rows are contiguous in memory)
   double* ixl,
   int* iy,		//y, train,  nx1
   double* ixp,		//x, test, pxnp (transposed so rows are contiguous in memory)
   double* ixlp,
   size_t m,		//number of trees
   int *numcut, //size_t nc,		//number of cut points
   size_t nd,		//number of kept draws (except for thinnning ..)
   size_t burn,		//number of burn-in draws skipped
   double mybeta,
   double alpha,
   double mub,
   const double g,
   double tau,
   bool dart,
   double theta,
   double omega,
   int *grp,
   double a,
   double b,
   double rho,
   bool aug,
   size_t nkeeptrain,
   size_t nkeeptest,
   //size_t nkeeptestme,
   size_t nkeeptreedraws,
   size_t printevery,
//   int treesaslists,
   unsigned int n1, // additional parameters needed to call from C++
   unsigned int n2,
/*
   double* trmean,
   double* temean,
*/
   double* _trdraw,
   double* _tedraw,
   double* _bdraw,
   double* _wtsdraw
)
{
   //return data structures (using C++)
   std::vector<double*> trdraw(nkeeptrain);
   std::vector<double*> tedraw(nkeeptest);
   std::vector<double*> bdraw(nkeeptrain+nskip);
   vector<double> sdraw(nd+burn),tadraw(nd+burn),tbdraw(nd+burn);

   for(size_t i=0; i<nkeeptrain; ++i) trdraw[i]=&_trdraw[i*n];
   for(size_t i=0; i<nkeeptest; ++i) tedraw[i]=&_tedraw[i*np];
   for(size_t i=0; i<nkeeptrain+nskip; ++i) bdraw[i]=&_bdraw[i*l];

   std::vector< std::vector<size_t> > varcnt;
   std::vector< std::vector<double> > varprb;

   //random number generation
   arn gen(n1, n2);

   bart bm(m);
#endif

/*
   for(size_t i=0; i<n; ++i) trmean[i]=0.0;
   for(size_t i=0; i<np; ++i) temean[i]=0.0;
*/

   double* iz = new double[n];
   double* ibf = new double[n];

   std::stringstream treess;  //string stream to write trees to
   treess.precision(10);
   treess << nkeeptreedraws << " " << m << " " << p << endl;
   // dart iterations
   std::vector<double> ivarprb (p,0.);
   std::vector<size_t> ivarcnt (p,0);

   printf("*****Into main of plbart\n");

   size_t skiptr,skipte,skiptreedraws;
   //size_t skiptr,skipte,skipteme,skiptreedraws;
   if(nkeeptrain) {skiptr=nd/nkeeptrain;}
   else skiptr = nd+1;
   if(nkeeptest) {skipte=nd/nkeeptest;}
   else skipte=nd+1;
/*
   if(nkeeptestme) {skipteme=nd/nkeeptestme;}
   else skipteme=nd+1;
*/
   if(nkeeptreedraws) {skiptreedraws = nd/nkeeptreedraws;}
   else skiptreedraws=nd+1;

   //--------------------------------------------------
   //print args
   printf("*****Data:\n");
   printf("data:n,p,l,np: %zu, %zu, %zu\n",n,p,l,np);
   printf("y1,yn: %d, %d\n",iy[0],iy[n-1]);
   printf("x1,x[n*p]: %lf, %lf\n",ix[0],ix[n*p-1]);
   if(np) printf("xp1,xp[np*p]: %lf, %lf\n",ixp[0],ixp[np*p-1]);
   printf("*****Number of Trees: %zu\n",m);
   printf("*****Number of Cut Points: %d ... %d\n", numcut[0], numcut[p-1]);
   printf("*****burn and ndpost: %zu, %zu\n",burn,nd);
   printf("*****Prior:mybeta,alpha,tau: %lf,%lf,%lf\n",
                   mybeta,alpha,tau);
   cout << "*****Dirichlet:sparse,theta,omega,a,b,rho,augment: " 
	<< dart << ',' << theta << ',' << omega << ',' << a << ',' 
	<< b << ',' << rho << ',' << aug << endl;
   printf("*****nkeeptrain,nkeeptest,nkeeptreedraws: %zu,%zu,%zu\n",
               nkeeptrain,nkeeptest,nkeeptreedraws);
/*
   printf("*****nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws: %zu,%zu,%zu,%zu\n",
               nkeeptrain,nkeeptest,nkeeptestme,nkeeptreedraws);
*/
   printf("*****printevery: %zu\n",printevery);
   //printf("*****skiptr,skipte,skipteme,skiptreedraws: %zu,%zu,%zu,%zu\n",skiptr,skipte,skipteme,skiptreedraws);
   printf("*****skiptr,skipte,skiptreedraws: %zu,%zu,%zu\n",skiptr,skipte,skiptreedraws);

   //--------------------------------------------------
   //bart bm(m);
   bm.setprior(alpha,mybeta,tau);
   bm.setdata(p,n,ix,iz,numcut);
   bm.setdart(a,b,rho,aug,dart,theta,omega);
   //--------------------------------------------------
//init
   Eigen::Map<const Eigen::MatrixXd> xm(ix,p,n);
   Eigen::Map<const Eigen::MatrixXd> xpm(ixp,p,np);
   Eigen::Map<const Eigen::MatrixXd> xlm(ixl,l,n);
   Eigen::Map<const Eigen::MatrixXd> xlpm(ixlp,l,np);
   Eigen::MatrixXd xmt(xm.transpose());
   Eigen::MatrixXd xpmt(xpm.transpose());
   Eigen::MatrixXd xlmt(xlm.transpose());
   Eigen::MatrixXd xlpmt(xlpm.transpose());

   Eigen::VectorXd bv(l);   // beta for linear
   Eigen::Map<Eigen::VectorXd> mu(imub,l);
   //mu=mub;
   //Rcout << "Vector mu: " << mu << std::endl;
   //Eigen::MatrixXd Sigma(Eigen::MatrixXd::Identity(p,p)/taub);
   Eigen::MatrixXd Sigma(xlm*xlmt);
   Sigma = g*Sigma.inverse();
   //bv=rnormXd(mu, Sigma);
   bv=mu;
   
   Eigen::VectorXd xb(xlmt*bv);
   Eigen::VectorXd xpb(xlpmt*bv);
   
   //update Sigma
   //Eigen::MatrixXd Tau((Eigen::MatrixXd(p,p).setZero().selfadjointView<Eigen::Lower>().rankUpdate(xm)));
   //Tau = Tau + Sigma.inverse();
   //Sigma = Tau.inverse();
   Sigma = Sigma/(g+1);


  for(size_t k=0; k<n; k++) {
    if(iy[k]==0) iz[k]= -rtnorm(0., xb(k), 1., gen);
    else iz[k]=rtnorm(0., -xb(k), 1., gen);
  }

  Eigen::Map<Eigen::VectorXd> zv(iz,n), bf(ibf,n);

   //--------------------------------------------------

   //--------------------------------------------------
   //temporary storage
   //out of sample fit
   double* fhattest=0; 
   if(np) { fhattest = new double[np]; }
   Eigen::VectorXd res(n);    // linear residual y-f(x)

   //--------------------------------------------------
   //mcmc
   printf("\nMCMC\n");
   //size_t index;
   size_t trcnt=0; //count kept train draws
   size_t tecnt=0; //count kept test draws
   bool keeptest,/*keeptestme,*/keeptreedraw;

   time_t tp;
   int time1 = time(&tp);
   xinfo& xi = bm.getxinfo();

   std::vector<double> kdraw(nkeeptrain);

   for(size_t i=0;i<(nd+burn);i++) {
      if(i%printevery==0) printf("done %zu (out of %lu)\n",i,nd+burn);
      if(i==(burn/2)&&dart) bm.startdart();
      //draw bart
      bm.draw(1., gen);
      bm.predict(p,n,ix,ibf);
      //keep track of linear residual
      res=zv+xlmt*bv-bf;
      //update mu for beta
      //mu=Sigma*(xm*res + mu*taub);
      mu = Sigma*xlm*res;
// + mu/(g+1);
      
      //draw beta
      bv = rnormXd(mu, Sigma);
      for(size_t k=0;k<l;k++) BDRAW(i,k)=bv[k];

      xb = xlmt*bv;

      for(size_t k=0; k<n; k++){
	if(iy[k]==0) iz[k]= -rtnorm(-bm.f(k), xb(k), 1., gen);
	else iz[k]=rtnorm(bm.f(k), -xb(k), 1., gen);
      }
      
      //draw tau
      std::vector<tree*> bnodes;
      std::vector<double> vtheta;
      for(size_t j=0;j<m;j++){
	bm.gettree(j).getbots(bnodes);
      }
      for(size_t l=0;l<bnodes.size();l++){
	vtheta.push_back(bnodes[l]->gettheta());
      }
      double nb=static_cast<int>(vtheta.size());
      double musq=std::inner_product( vtheta.begin(), vtheta.end(), vtheta.begin(), 0.0L );
      double ta=12.05+ nb/2, tb=0.0406+musq/2;  //gamma prior for tau
      tau = 1/sqrt(gen.gamma(ta,tb));  //gamma posterior for tau
      tadraw[i]=ta;
      tbdraw[i]=tb;
	
      bm.settau(tau);
      
      sdraw[i] = tau;
  
      if(i>=burn) {
         if(nkeeptrain && (((i-burn+1) % skiptr) ==0)) {
	    kdraw[trcnt] = 3/(tau*sqrt(50));
            for(size_t k=0;k<n;k++) TRDRAW(trcnt,k)=bm.f(k)+xb(k); 
            trcnt+=1;
         }
         keeptest = nkeeptest && (((i-burn+1) % skipte) ==0) && np;
         if(keeptest) {
	   bm.predict(p,np,ixp,fhattest);
	   xpb = xlpmt*bv;
            for(size_t k=0;k<np;k++) TEDRAW(tecnt,k)=fhattest[k]+xpb(k);
            tecnt+=1;
         }
         keeptreedraw = nkeeptreedraws && (((i-burn+1) % skiptreedraws) ==0);
         if(keeptreedraw) {

            for(size_t j=0;j<m;j++) {
	      treess << bm.gettree(j);
	    }
	    #ifndef NoRcpp
	    ivarcnt=bm.getnv();
	    ivarprb=bm.getpv();
	    size_t k=(i-burn)/skiptreedraws;
	    for(size_t j=0;j<p;j++){
	      varcnt(k,j)=ivarcnt[j];
	      varprb(k,j)=ivarprb[j];
	    }
            #else
	    varcnt.push_back(bm.getnv());
	    varprb.push_back(bm.getpv());
	    #endif
         }
      }
      
   }
   int time2 = time(&tp);
   printf("time: %ds\n",time2-time1);
   printf("check counts\n");
   printf("trcnt,tecnt: %zu,%zu\n",trcnt,tecnt);
   //printf("trcnt,tecnt,temecnt,treedrawscnt: %zu,%zu,%zu,%zu\n",trcnt,tecnt,temecnt,treedrawscnt);
   //--------------------------------------------------
   //PutRNGstate();

   if(fhattest) delete[] fhattest;
   delete[] iz;
   delete[] ibf;

#ifndef NoRcpp
   //--------------------------------------------------
   //return
   Rcpp::List ret;
   //ret["yhat.train.mean"]=trmean;
   ret["yhat.train"]=trdraw;
   //ret["yhat.test.mean"]=temean;
   ret["yhat.test"]=tedraw;
   ret["beta.train"]=bdraw;
   ret["sigma.mu"]=sdraw;
   //ret["discrete.p"]=wtsdraw;
   ret["gamma.a"]=tadraw;
   ret["gamma.b"]=tbdraw;
   ret["kdraw"]=kdraw;
   //ret["varcount"]=varcount;
   ret["varcount"]=varcnt;
   ret["varprob"]=varprb;

   Rcpp::List xiret(xi.size());
   for(size_t i=0;i<xi.size();i++) {
      Rcpp::NumericVector vtemp(xi[i].size());
      std::copy(xi[i].begin(),xi[i].end(),vtemp.begin());
      xiret[i] = Rcpp::NumericVector(vtemp);
   }

   Rcpp::List treesL;
   treesL["cutpoints"] = xiret;
   treesL["trees"]=Rcpp::CharacterVector(treess.str());
   //if(treesaslists) treesL["lists"]=list_of_lists;
   ret["treedraws"] = treesL;

   return ret;
#else

#endif

}
