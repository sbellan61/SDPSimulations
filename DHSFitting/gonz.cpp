#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
NumericVector cplC(NumericVector seros, IntegerVector active) {
  int n = seros.size();
  // NumericVector out(n);
  for(IntegerVector::iterator jj = active.begin(); jj != active.end(); ++jj) {
    // Rprintf("active %i \n",*jj);
    seros(*jj-1) = seros(*jj-1)*2;
  }
  return seros;
}


// [[Rcpp::export]]
NumericMatrix prec(NumericMatrix seros, IntegerVector active, NumericVector pmb, NumericVector pfb, NumericVector pmbA, NumericVector pfbA) {
  NumericMatrix serosO = clone(seros); // Old seros (last state from which to iterate)
  if (min(active) < 1) stop("active < 0");
  int numcouples = seros.nrow();
  if (max(active) > numcouples) stop("max(active) > numcouples");

  int s__ = 0;
  int mb_a1 = 1;
  int mb_a2 = 2;
  int mb_ = 3;
  int me_a1 = 4; 
  int me_a2 = 5;
  int me_ = 6;
  int f_ba1 = 7;
  int f_ba2 = 8;
  int f_b = 9;
  int f_ea1 = 10;
  int f_ea2 = 11;
  int f_e = 12;
  int hb1b2 = 13;
  int hb2b1 = 14;
  int hbe = 15;
  int heb = 16;
  int hbpa = 17;
  int hpba = 18;
  int hepa = 19;
  int hpea = 20;
  int hbp = 21;
  int hpb = 22;
  int hep = 23;
  int hpe = 24;
  int he1e2 = 25;
  int he2e1 = 26;
  int mb_a1A = 27;
  int mb_a2A = 28;
  int mb_A = 29;
  int me_a1A = 30;
  int me_a2A = 31;
  int me_A = 32;
  int f_ba1A = 33;
  int f_ba2A = 34;
  int f_bA = 35;
  int f_ea1A = 36;
  int f_ea2A = 37;
  int f_eA = 38;
  int hb1b2A = 39;
  int hb2b1A = 40;
  int hbeA = 41;
  int hebA = 42;
  int hbpaA = 43;
  int hpbaA = 44;
  int hepaA = 45;
  int hpeaA = 46;
  int hbpA = 47;
  int hpbA = 48;		// package as a dll
  int hepA = 49;
  int hpeA = 50;
  int he1e2A = 51;
  int he2e1A = 52;
  NumericVector p_mfirstA = clone(pmb);
  NumericVector p_ffirstA = clone(pmb);
  NumericVector p_mfirst = clone(pmb);
  NumericVector p_ffirst = clone(pmb);

  for(IntegerVector::iterator jj = active.begin()-1; jj != active.end()-1; ++jj) {

    seros(*jj,s__)    = serosO(*jj,s__)* (1-pmb(*jj)) * (1 - pfb(*jj));
    // conditional on survival
    p_mfirstA(*jj) = pmbA(*jj) / (pmbA(*jj)+pfbA(*jj)); // Need to deal with NaN still
    p_ffirstA(*jj) = 1-p_mfirstA(*jj); // Need to deal with NaN still			     
    seros(*jj,mb_a1A)  = serosO(*jj,s__) * pmbA(*jj) * (1 - pfb(*jj));
    seros(*jj,mb_a2A) = serosO(*jj,mb_a1A)*(1-pfb(*jj));
    seros(*jj,mb_A)   = serosO(*jj,mb_a2A)*(1-pfb(*jj)) + serosO(*jj,mb_A)*(1 - pfb(*jj));
    seros(*jj,f_ba1A) = serosO(*jj,s__)* pfbA(*jj) *(1-pmb(*jj));
    seros(*jj,f_ba2A) = serosO(*jj,f_ba1A)*(1 - pmb(*jj));
    seros(*jj,f_bA)   = serosO(*jj,f_ba2A)*(1 - pmb(*jj)) + serosO(*jj,f_bA)*(1 - pmb(*jj));
    seros(*jj,hb1b2A) = serosO(*jj,hb1b2A) + p_mfirstA(*jj) * serosO(*jj,s__)* pmbA(*jj) * pfbA(*jj) + (serosO(*jj,mb_a1A) + serosO(*jj,mb_a2A) + serosO(*jj,mb_A)) * pfbA(*jj);
    seros(*jj,hb2b1A) = serosO(*jj,hb2b1A) + p_ffirstA(*jj) * serosO(*jj,s__)* pmbA(*jj) * pfbA(*jj) + (serosO(*jj,f_ba1A) + serosO(*jj,f_ba2A) + serosO(*jj,f_bA)) * pmbA(*jj);
    // unconditional on survival
    p_mfirst(*jj) = pmb(*jj) / (pmb(*jj)+pfb(*jj)); // Need to deal with NaN still    
    p_ffirst(*jj) = 1-p_mfirst(*jj);	     // Need to deal with NaN still		     
    seros(*jj,mb_a1) = serosO(*jj,s__) * pmb(*jj) * (1-pfb(*jj));
    seros(*jj,mb_a2) = serosO(*jj,mb_a1) * (1 - pfb(*jj));
    seros(*jj,mb_) = serosO(*jj,mb_a2) * (1 - pfb(*jj)) + serosO(*jj,mb_) * (1 - pfb(*jj));
    seros(*jj,f_ba1) = serosO(*jj,s__) * pfb(*jj) * (1-pmb(*jj));
    seros(*jj,f_ba2) = serosO(*jj,f_ba1) * (1 - pmb(*jj));
    seros(*jj,f_b) = serosO(*jj,f_ba2) * (1 - pmb(*jj)) + serosO(*jj,f_b) * (1 - pmb(*jj));
    seros(*jj,hb1b2) = serosO(*jj,hb1b2) + p_mfirst(*jj)  *  serosO(*jj,s__) * pmb(*jj) * pfb(*jj) + (serosO(*jj,mb_a1) + serosO(*jj,mb_a2) + serosO(*jj,mb_))  *  pfb(*jj);
    seros(*jj,hb2b1) = serosO(*jj,hb2b1) + p_ffirst(*jj)  *  serosO(*jj,s__) * pmb(*jj) * pfb(*jj) + (serosO(*jj,f_ba1) + serosO(*jj,f_ba2) + serosO(*jj,f_b))  *  pmb(*jj);
    // Rprintf("\n %d", *jj);    
  }
  return seros;
}
 


// [[Rcpp::export]]
NumericVector examp() {
  NumericMatrix out(6,6);
  for(int i = 0; i < 6; ++i) out(i,0) = 1;
  return out;
}




// [[Rcpp::export]]
NumericMatrix precLoop(int max_bd, NumericMatrix pre_fprev, NumericMatrix pre_mprev, NumericMatrix pre_msurv, NumericMatrix pre_fsurv, List activeList,
		       double bmb, double bfb) {

  Rcpp::List xlist(activeList); 
  int n = xlist.size(); 
  std::vector<double> res(n);
    
  Rprintf("num active %d \n", xlist[0](1));

  // Column names for sero
  int s__ = 0; int mb_a1 = 1; int mb_a2 = 2; int mb_ = 3; int me_a1 = 4;  int me_a2 = 5; int me_ = 6; int f_ba1 = 7; int f_ba2 = 8;
  int f_b = 9; int f_ea1 = 10; int f_ea2 = 11; int f_e = 12; int hb1b2 = 13; int hb2b1 = 14; int hbe = 15; int heb = 16; int hbpa = 17;
  int hpba = 18; int hepa = 19; int hpea = 20; int hbp = 21; int hpb = 22; int hep = 23; int hpe = 24; int he1e2 = 25; int he2e1 = 26;
  int mb_a1A = 27; int mb_a2A = 28; int mb_A = 29; int me_a1A = 30; int me_a2A = 31; int me_A = 32; int f_ba1A = 33; int f_ba2A = 34;
  int f_bA = 35; int f_ea1A = 36; int f_ea2A = 37; int f_eA = 38; int hb1b2A = 39; int hb2b1A = 40; int hbeA = 41; int hebA = 42;
  int hbpaA = 43; int hpbaA = 44; int hepaA = 45; int hpeaA = 46; int hbpA = 47; int hpbA = 48; int hepA = 49; int hpeA = 50;
  int he1e2A = 51; int he2e1A = 52;

  int numcouples = pre_fsurv.nrow(); // # of couples
  NumericMatrix seros(numcouples, 53); // Initialize couples by serostate matrix
  for(int ii = 0; ii < numcouples; ++ii) seros(ii,s__) = 1; // All couples start both susceptible
  // Initialize transmission hazards & probabilities & joint probability of transmission & survival (*a*live)
  NumericVector m_haz(numcouples); // Note, these vectors change within the loop over time below.
  NumericVector f_haz(numcouples);
  NumericVector pmb(numcouples);
  NumericVector pfb(numcouples);
  NumericVector pmbA(numcouples);
  NumericVector pfbA(numcouples);
  NumericVector p_mfirstA(numcouples);
  NumericVector p_ffirstA(numcouples);
  NumericVector p_mfirst(numcouples);
  NumericVector p_ffirst(numcouples);

  // For month of marriage
  for(int tt = 0; tt < max_bd; ++tt) { 
    NumericMatrix serosO = clone(seros); // Copy old seros (last state from which to iterate)
    NumericVector active = activeList[tt]; // Currently active couples
    int numactive = active.size();		 // Number of active couples
    if (min(active) < 1) stop("active < 0"); // Check indices are within arrays
    if (max(active) > numcouples) stop("max(active) > numcouples");

    // For each active couple
    // for(IntegerVector::iterator jj = active.begin()-1; jj != active.end()-1; ++jj) {
      
    //   m_haz(*jj) = bmb*pre_fprev(*jj,tt); // hazards to sexually active men
    //   f_haz(*jj) = bfb*pre_mprev(*jj,tt); // hazards to sexually active women
    //   pmb(*jj) = 1 - exp(-m_haz(*jj));              // transmission probabilities
    //   pfb(*jj) = 1 - exp(-f_haz(*jj));
    //   pmbA(*jj) = pmb(*jj)*pre_msurv(*jj,tt); // joint transmission & survival probabilities 
    //   pfbA(*jj) = pfb(*jj)*pre_fsurv(*jj,tt);
    //   seros(*jj,s__)    = serosO(*jj,s__)* (1-pmb(*jj))*(1 - pfb(*jj));
    //   // conditional on survival
    //   p_mfirstA(*jj) = pmbA(*jj) / (pmbA(*jj)+pfbA(*jj)); // Need to deal with NaN still
    //   p_ffirstA(*jj) = 1-p_mfirstA(*jj); // Need to deal with NaN still			     
    //   seros(*jj,mb_a1A)  = serosO(*jj,s__)*pmbA(*jj)*(1 - pfb(*jj));
    //   seros(*jj,mb_a2A) = serosO(*jj,mb_a1A)*(1-pfb(*jj));
    //   seros(*jj,mb_A)   = serosO(*jj,mb_a2A)*(1-pfb(*jj)) + serosO(*jj,mb_A)*(1 - pfb(*jj));
    //   seros(*jj,f_ba1A) = serosO(*jj,s__)* pfbA(*jj) *(1-pmb(*jj));
    //   seros(*jj,f_ba2A) = serosO(*jj,f_ba1A)*(1 - pmb(*jj));
    //   seros(*jj,f_bA)   = serosO(*jj,f_ba2A)*(1 - pmb(*jj)) + serosO(*jj,f_bA)*(1 - pmb(*jj));
    //   seros(*jj,hb1b2A) = serosO(*jj,hb1b2A) + p_mfirstA(*jj)*serosO(*jj,s__)* pmbA(*jj)*pfbA(*jj) + (serosO(*jj,mb_a1A) + serosO(*jj,mb_a2A) + serosO(*jj,mb_A))*pfbA(*jj);
    //   seros(*jj,hb2b1A) = serosO(*jj,hb2b1A) + p_ffirstA(*jj)*serosO(*jj,s__)* pmbA(*jj)*pfbA(*jj) + (serosO(*jj,f_ba1A) + serosO(*jj,f_ba2A) + serosO(*jj,f_bA))*pmbA(*jj);
    //   // unconditional on survival
    //   p_mfirst(*jj) = pmb(*jj) / (pmb(*jj)+pfb(*jj)); // Need to deal with NaN still    
    //   p_ffirst(*jj) = 1-p_mfirst(*jj);	     // Need to deal with NaN still		     
    //   seros(*jj,mb_a1) = serosO(*jj,s__)*pmb(*jj)*(1-pfb(*jj));
    //   seros(*jj,mb_a2) = serosO(*jj,mb_a1)*(1 - pfb(*jj));
    //   seros(*jj,mb_) = serosO(*jj,mb_a2)*(1 - pfb(*jj)) + serosO(*jj,mb_)*(1 - pfb(*jj));
    //   seros(*jj,f_ba1) = serosO(*jj,s__)*pfb(*jj)*(1-pmb(*jj));
    //   seros(*jj,f_ba2) = serosO(*jj,f_ba1)*(1 - pmb(*jj));
    //   seros(*jj,f_b) = serosO(*jj,f_ba2)*(1 - pmb(*jj)) + serosO(*jj,f_b)*(1 - pmb(*jj));
    //   seros(*jj,hb1b2) = serosO(*jj,hb1b2) + p_mfirst(*jj)*serosO(*jj,s__)*pmb(*jj)*pfb(*jj) + (serosO(*jj,mb_a1) + serosO(*jj,mb_a2) + serosO(*jj,mb_))*pfb(*jj);
    //   seros(*jj,hb2b1) = serosO(*jj,hb2b1) + p_ffirst(*jj)*serosO(*jj,s__)*pmb(*jj)*pfb(*jj) + (serosO(*jj,f_ba1) + serosO(*jj,f_ba2) + serosO(*jj,f_b))*pmb(*jj);
    //   // Rprintf("\n %d", *jj);    
    // }
  }
  return seros;
}

/*** R

*/

				// new = prec(test, active = 2:4, pmb = .05, pfb = .03, pmbA = .05*.5, pfbA = .03*.5)
     // print(new)


// pre.couple <- function(seros.active, pmb, pfb, pmb.a, pfb.a, uncond.mort=T) {
//     if(class(seros.active)=='numeric')  seros.active <- t(as.matrix(seros.active))
//     tp <- seros.active ## old temporary array to update from
//     seros.active[,'s..']   <- tp[,'s..'] * (1-pmb) * (1-pfb)
//     ## transmission and alive, these states are used for fitting
//     seros.active[,'mb.a1A'] <- tp[,'s..'] * pmb.a * (1-pfb)
//     seros.active[,'mb.a2A'] <- tp[,'mb.a1A']*(1-pfb)
//     seros.active[,'mb.A'] <- tp[,'mb.a2A']*(1-pfb) + tp[,'mb.A']*(1 - pfb)
//     seros.active[,'f.ba1A'] <- tp[,'s..']* pfb.a *(1-pmb)
//     seros.active[,'f.ba2A'] <- tp[,'f.ba1A']*(1 - pmb)
//     seros.active[,'f.bA'] <- tp[,'f.ba2A']*(1 - pmb) + tp[,'f.bA']*(1 - pmb)
//     p.mfirst.a <- pmb.a / (pmb.a+pfb.a)
//     p.ffirst.a <- 1-p.mfirst.a
//     p.mfirst.a[is.na(p.mfirst.a)] <- 0
//     p.ffirst.a[is.na(p.ffirst.a)] <- 0                
//     seros.active[,'hb1b2A'] <- tp[,'hb1b2A'] + p.mfirst.a * tp[,'s..']* pmb.a * pfb.a + (tp[,'mb.a1A'] + tp[,'mb.a2A'] + tp[,'mb.A']) * pfb.a
//     seros.active[,'hb2b1A'] <- tp[,'hb2b1A'] + p.ffirst.a * tp[,'s..']* pmb.a * pfb.a + (tp[,'f.ba1A'] + tp[,'f.ba2A'] + tp[,'f.bA']) * pmb.a
//     if(uncond.mort) { ## if calculating probabilities of transmission with or without death (not used for fitting)
//         seros.active[,'mb.a1'] <- tp[,'s..'] * pmb * (1-pfb)
//         seros.active[,'mb.a2'] <- tp[,'mb.a1'] * (1 - pfb)
//         seros.active[,'mb.'] <- tp[,'mb.a2'] * (1 - pfb) + tp[,'mb.'] * (1 - pfb)
//         seros.active[,'f.ba1'] <- tp[,'s..'] * pfb * (1-pmb)
//         seros.active[,'f.ba2'] <- tp[,'f.ba1'] * (1 - pmb)
//         seros.active[,'f.b'] <- tp[,'f.ba2'] * (1 - pmb) + tp[,'f.b'] * (1 - pmb)
//         p.mfirst <- pmb / (pmb+pfb)
//         p.ffirst <- 1-p.mfirst
//         p.mfirst[is.na(p.mfirst)] <- 0
//         p.ffirst[is.na(p.ffirst)] <- 0                
//         seros.active[,'hb1b2'] <- tp[,'hb1b2'] + p.mfirst  *  tp[,'s..'] * pmb * pfb + (tp[,'mb.a1'] + tp[,'mb.a2'] + tp[,'mb.'])  *  pfb
//         seros.active[,'hb2b1'] <- tp[,'hb2b1'] + p.ffirst  *  tp[,'s..'] * pmb * pfb + (tp[,'f.ba1'] + tp[,'f.ba2'] + tp[,'f.b'])  *  pmb
//     }
//     return(seros.active)
// }


     // library(microbenchmark)
     // sroR1 = function(seros, active) {
     // out = seros + seros*5
     // return(out)
     // }
     // cplR = function(seros, active) {
     // seros[active] = seros[active]*5
     // return(seros)
     // }
     // seros = 1:10e5
     // active = seq(1, length(seros), by = 2)
     // microbenchmark(
     // cplC(seros, active)[1:5],
     // cplR(seros, active)[1:5]
     // )
     // cplC(seros, active)[1:5]
     // cplR(seros, active)[1:5]

 
