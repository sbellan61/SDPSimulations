#include <Rcpp.h>
using namespace Rcpp;

// Column names for seros (serostate matrix)
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
int hpbA = 48;
int hepA = 49;
int hpeA = 50;
int he1e2A = 51;
int he2e1A = 52;
static const float inv_sqrt_2pi = 1/sqrt(2*M_PI); // for prnormC below

// [[Rcpp::export]]
double dnormC(const double & x, const double & mean, const double & sd) {
    double a = (x - mean) / sd;
    return inv_sqrt_2pi / sd * std::exp(-0.5f * a * a);
}

// [[Rcpp::export]]
Rcpp::List precLoop(const int & max_bd, const int & max_cd, const IntegerVector & datser,
		       const NumericMatrix & pre_fprev, const NumericMatrix & pre_mprev, const NumericMatrix & within_fprev, const NumericMatrix & within_mprev,
		       const NumericMatrix & pre_msurv, const NumericMatrix & pre_fsurv, const NumericMatrix & within_msurv, const NumericMatrix & within_fsurv,
		       const NumericMatrix & within_art_cov,
		       const List & PreActiveList, const List & WithinActiveList,
		       const double &  bmb, const double &  bfb, const double &  bme, const double &  bfe, const double &  bmp,
		       const double &  lrho, const double & lrho_sd, const double & trans_ratio,
		       const double & acute_sc, int uselessnum,
		       const bool & partner_arv, const double & cov_scalar) {

  int numcouples = pre_fsurv.nrow(); // # of couples
  NumericMatrix seros(numcouples, 53); // Initialize couples by serostate matrix
  for(int jj = 0; jj < numcouples; ++jj) seros(jj,s__) = 1; // All couples start both susceptible
  // Initialize transmission hazards & probabilities & joint probability of transmission & survival (*a*live)
  double bfp = bmp * exp(lrho); // Conver M->F/F->M within-couple ratio (rho) to bfp
  NumericVector m_haz(numcouples); // Note, these vectors change within the loop over time below.
  NumericVector f_haz(numcouples);
  NumericVector pmb(numcouples);
  NumericVector pfb(numcouples);
  NumericVector pmbA(numcouples);
  NumericVector pfbA(numcouples);
  NumericVector pme(numcouples);
  NumericVector pfe(numcouples);
  NumericVector pmeA(numcouples);
  NumericVector pfeA(numcouples);
  NumericVector pmp(numcouples);
  NumericVector pfp(numcouples);
  NumericVector pmpA(numcouples);
  NumericVector pfpA(numcouples);
  NumericVector pmp_ac(numcouples);
  NumericVector pfp_ac(numcouples);
  NumericVector pmpA_ac(numcouples);
  NumericVector pfpA_ac(numcouples);
  NumericVector p_mfirstA(numcouples);
  NumericVector p_ffirstA(numcouples);
  NumericVector denomA(numcouples);
  NumericVector p_mfirst(numcouples);
  NumericVector p_ffirst(numcouples);
  NumericVector denom(numcouples);
  NumericVector within_art_scalar(numcouples);
  int numactive(1);

  // for(int useless = 0; useless < uselessnum; ++useless) {

    NumericMatrix serosO = clone(seros); 
  // For month of pre-couple duration
  for(int tt = 0; tt < max_bd; ++tt) { 
    IntegerVector active = PreActiveList[tt]; // Currently active couples, -1 is to switch from R to C array indicing
    serosO = seros;// Copy old seros (last state from which to iterate)
    numactive = active.size();		 // Number of active couples
    if (min(active) < 0) stop("active < 0"); // Check indices are within arrays
    if (max(active) > numcouples) stop("max(active) > numcouples");
    // Rprintf("\n time step %d", tt);
    // For each active couple
    for(IntegerVector::iterator jj = active.begin(); jj != active.end(); ++jj) { 
      // Rprintf(" %d", *jj+1);
      m_haz(*jj) = bmb*pre_fprev(*jj,tt); // hazards to sexually active men
      f_haz(*jj) = bfb*pre_mprev(*jj,tt); // hazards to sexually active women
      pmb(*jj) = 1 - exp(-m_haz(*jj));              // transmission probabilities
      pfb(*jj) = 1 - exp(-f_haz(*jj));
      pmbA(*jj) = pmb(*jj)*pre_msurv(*jj,tt); // joint transmission & survival probabilities 
      pfbA(*jj) = pfb(*jj)*pre_fsurv(*jj,tt);
      seros(*jj,s__)    = serosO(*jj,s__)* (1-pmb(*jj))*(1 - pfb(*jj));
      // conditional on survival
      denomA(*jj) =pmbA(*jj)+pfbA(*jj); // Competing rates (who is infected first if they happen in same month)
      if(denomA(*jj)==0) {		       // If both rates are 0, then avoid divide by 0 error
	p_mfirstA(*jj) = 0;
	p_ffirstA(*jj) = 0;
      }else{
	p_mfirstA(*jj) = pmbA(*jj) / denomA(*jj); // Otherwise use competing risk framework
	p_ffirstA(*jj) = 1-p_mfirstA(*jj);
      }
      seros(*jj,mb_a1A)  = serosO(*jj,s__)*pmbA(*jj)*(1 - pfb(*jj));
      seros(*jj,mb_a2A) = serosO(*jj,mb_a1A)*(1-pfb(*jj));
      seros(*jj,mb_A)   = serosO(*jj,mb_a2A)*(1-pfb(*jj)) + serosO(*jj,mb_A)*(1 - pfb(*jj));
      seros(*jj,f_ba1A) = serosO(*jj,s__)* pfbA(*jj) *(1-pmb(*jj));
      seros(*jj,f_ba2A) = serosO(*jj,f_ba1A)*(1 - pmb(*jj));
      seros(*jj,f_bA)   = serosO(*jj,f_ba2A)*(1 - pmb(*jj)) + serosO(*jj,f_bA)*(1 - pmb(*jj));
      seros(*jj,hb1b2A) = serosO(*jj,hb1b2A) + p_mfirstA(*jj)*serosO(*jj,s__)* pmbA(*jj)*pfbA(*jj) + (serosO(*jj,mb_a1A) + serosO(*jj,mb_a2A) + serosO(*jj,mb_A))*pfbA(*jj);
      seros(*jj,hb2b1A) = serosO(*jj,hb2b1A) + p_ffirstA(*jj)*serosO(*jj,s__)* pmbA(*jj)*pfbA(*jj) + (serosO(*jj,f_ba1A) + serosO(*jj,f_ba2A) + serosO(*jj,f_bA))*pmbA(*jj);
      // unconditional on survival
      p_mfirst(*jj) = pmb(*jj) / (pmb(*jj)+pfb(*jj)); 
      p_ffirst(*jj) = 1-p_mfirst(*jj);
      // conditional on survival
      denom(*jj)=pmb(*jj)+pfb(*jj);
      if(denom(*jj)==0) {
	p_mfirst(*jj) = 0;
	p_ffirst(*jj) = 0;
      }else{
	p_mfirst(*jj) = pmb(*jj) / denom(*jj); 
	p_ffirst(*jj) = 1-p_mfirst(*jj);
      }
      seros(*jj,mb_a1) = serosO(*jj,s__)*pmb(*jj)*(1-pfb(*jj));
      seros(*jj,mb_a2) = serosO(*jj,mb_a1)*(1 - pfb(*jj));
      seros(*jj,mb_) = serosO(*jj,mb_a2)*(1 - pfb(*jj)) + serosO(*jj,mb_)*(1 - pfb(*jj));
      seros(*jj,f_ba1) = serosO(*jj,s__)*pfb(*jj)*(1-pmb(*jj));
      seros(*jj,f_ba2) = serosO(*jj,f_ba1)*(1 - pmb(*jj));
      seros(*jj,f_b) = serosO(*jj,f_ba2)*(1 - pmb(*jj)) + serosO(*jj,f_b)*(1 - pmb(*jj));
      seros(*jj,hb1b2) = serosO(*jj,hb1b2) + p_mfirst(*jj)*serosO(*jj,s__)*pmb(*jj)*pfb(*jj) + (serosO(*jj,mb_a1) + serosO(*jj,mb_a2) + serosO(*jj,mb_))*pfb(*jj);
      seros(*jj,hb2b1) = serosO(*jj,hb2b1) + p_ffirst(*jj)*serosO(*jj,s__)*pmb(*jj)*pfb(*jj) + (serosO(*jj,f_ba1) + serosO(*jj,f_ba2) + serosO(*jj,f_b))*pmb(*jj);
      // Rprintf("\n %d", *jj);    
    } // end couple iteration
  } // end month or pre-couple duration iteration

  /////////////////////////////////////////////////// 
  int max_withinloop = max_cd - 1;
  // //////////////////////////////////////////////////
  for(int ttw = 0; ttw < max_withinloop; ttw++) { // Now loop through marriage
  IntegerVector within_active = WithinActiveList[ttw]; // Currently active couples, -1 is to switch from R to C array indicing
    serosO = seros;// Copy old seros (last state from which to iterate)
  numactive = within_active.size();		 // Number of active couples
  if (min(within_active) < 0) stop("within_active < 0"); // Check indices are within arrays
  if (max(within_active) > numcouples) stop("max(within-active) > numcouples");
  for(IntegerVector::iterator jj = within_active.begin(); jj != within_active.end(); ++jj) {
    // Within-Couple Hazards  
    pmp(*jj) = 1 - exp(-bmp); // probability of being infected by partner (constant, used inside loop)
    pfp(*jj) = 1 - exp(-bfp);
    pmp_ac(*jj) = 1 - exp(-acute_sc * bmp); //  during acute
    pfp_ac(*jj) = 1 - exp(-acute_sc * bfp);
    // Extra-Couple Hazards
    m_haz(*jj) = bme*within_fprev(*jj,ttw); // hazards to sexually active men
    f_haz(*jj) = bfe*within_mprev(*jj,ttw); // hazards to sexually active women
    pme(*jj) = 1 - exp(-m_haz(*jj));              // transmission probabilities
    pfe(*jj) = 1 - exp(-f_haz(*jj));
    // adjust probability of being infected by partner by probability partner is infectious (i.e. not on ART)
    if(partner_arv) { 		// NEED TO CHECK THIS STATEMENT WORKS STILL
    within_art_scalar(*jj) = (1 - cov_scalar * within_art_cov(*jj,ttw));
    pmp(*jj) = 1 - exp(-bmp * within_art_scalar(*jj));
    pfp(*jj) = 1 - exp(-bfp * within_art_scalar(*jj));
    pmp_ac(*jj) = 1 - exp(-acute_sc * bmp * within_art_scalar(*jj)); // acute
    pfp_ac(*jj) = 1 - exp(-acute_sc * bfp * within_art_scalar(*jj));
  }
    // Joint Transmission and survival for all routes 
    pmpA(*jj) = pmp(*jj) * within_msurv(*jj,ttw);; // Transmission probabilities from partner (jointly with survival)
    pfpA(*jj) = pfp(*jj) * within_fsurv(*jj,ttw);;
    pmpA_ac(*jj) = pmp_ac(*jj) * within_msurv(*jj,ttw);; // acute from partner (jointly with survival)
    pfpA_ac(*jj) = pfp_ac(*jj) * within_fsurv(*jj,ttw);;
    pmeA(*jj) = pme(*jj) * within_msurv(*jj,ttw);;	// extra-couple
    pfeA(*jj) = pfe(*jj) * within_fsurv(*jj,ttw);;
    // Within-couple transmission iterator;
    seros(*jj,s__) = serosO(*jj,s__)*(1-pme(*jj))*(1-pfe(*jj));
    // //////////////////////////////////////////////////;
    // Joint with survival;
    seros(*jj,mb_a1A) = 0;
    seros(*jj,mb_a2A) = serosO(*jj,mb_a1A)*(1-pfe(*jj))*(1-pfp_ac(*jj));
    seros(*jj,mb_A)  = serosO(*jj,mb_a2A)*(1-pfe(*jj))*(1-pfp_ac(*jj)) + serosO(*jj,mb_A)*(1-pfe(*jj))*(1-pfp(*jj)) ;
    seros(*jj,me_a1A) = serosO(*jj,s__)*pmeA(*jj)*(1-pfe(*jj));
    seros(*jj,me_a2A) = serosO(*jj,me_a1A)*(1-pfe(*jj))*(1-pfp_ac(*jj));
    seros(*jj,me_A)  = serosO(*jj,me_a2A)*(1-pfe(*jj))*(1-pfp_ac(*jj)) + serosO(*jj,me_A)*(1-pfe(*jj))*(1-pfp(*jj));
    seros(*jj,f_ba1A) = 0;
    seros(*jj,f_ba2A) = serosO(*jj,f_ba1A)*(1-pme(*jj))*(1-pmp_ac(*jj));
    seros(*jj,f_bA)  = serosO(*jj,f_ba2A)*(1-pme(*jj))*(1-pmp_ac(*jj)) + serosO(*jj,f_bA)*(1-pme(*jj))*(1-pmp(*jj));
    seros(*jj,f_ea1A) = serosO(*jj,s__)*pfeA(*jj)*(1-pme(*jj));
    seros(*jj,f_ea2A) = serosO(*jj,f_ea1A)*(1-pme(*jj))*(1-pmp_ac(*jj));
    seros(*jj,f_eA)  = serosO(*jj,f_ea2A)*(1-pme(*jj))*(1-pmp_ac(*jj)) + serosO(*jj,f_eA)*(1-pme(*jj))*(1-pmp(*jj));
    // hb1b2A) = hb1b2A) // Doesn't change during couple duration;
    // hb2b1A) = hb2b1A) // Doesn't change during couple duration                ;
    seros(*jj,hbeA) = serosO(*jj,hbeA)  + (serosO(*jj,mb_a1A) + serosO(*jj,mb_a2A))*(1-pfp_ac(*jj))*pfeA(*jj) + serosO(*jj,mb_A)*(1-pfp(*jj))*pfeA(*jj);
    seros(*jj,hebA) = serosO(*jj,hebA)  + (serosO(*jj,f_ba1A) + serosO(*jj,f_ba2A))*(1-pmp_ac(*jj))*pmeA(*jj) + serosO(*jj,f_bA)*(1-pmp(*jj))*pmeA(*jj);
    seros(*jj,hbpaA) = serosO(*jj,hbpaA) + (serosO(*jj,mb_a1A) + serosO(*jj,mb_a2A))*pfpA_ac(*jj);
    seros(*jj,hbpA) = serosO(*jj,hbpA)  + serosO(*jj,mb_A)*pfpA(*jj);
    seros(*jj,hpbaA) = serosO(*jj,hpbaA) + (serosO(*jj,f_ba1A) + serosO(*jj,f_ba2A))*pmpA_ac(*jj);
    seros(*jj,hpbA) = serosO(*jj,hpbA)  + serosO(*jj,f_bA)*pmpA(*jj);
    seros(*jj,hepaA) = serosO(*jj,hepaA) + (serosO(*jj,me_a1A) + serosO(*jj,me_a2A))*pfpA_ac(*jj);
    seros(*jj,hepA) = serosO(*jj,hepA)  + serosO(*jj,me_A)*pfpA(*jj);
    seros(*jj,hpeaA) = serosO(*jj,hpeaA) + (serosO(*jj,f_ea1A) + serosO(*jj,f_ea2A))*pmpA_ac(*jj);
    seros(*jj,hpeA) = serosO(*jj,hpeA)  + serosO(*jj,f_eA)*pmpA(*jj);
    denomA(*jj)=pmeA(*jj)+pfeA(*jj); // Competing rates (who is infected first if they happen in same month)
    if(denomA(*jj)==0) {		       // If both rates are 0, then avoid divide by 0 error
      p_mfirstA(*jj) = 0;
      p_ffirstA(*jj) = 0;
    }else{
      p_mfirstA(*jj) = pmeA(*jj) / denomA(*jj); // Otherwise use competing risk framework
      p_ffirstA(*jj) = 1-p_mfirstA(*jj);
    }
    seros(*jj,he1e2A) = serosO(*jj,he1e2A) + p_mfirstA(*jj) * serosO(*jj,s__)*pmeA(*jj)*pfeA(*jj) +
			(serosO(*jj,me_a1A) + serosO(*jj,me_a2A))*(1-pfp_ac(*jj))*pfeA(*jj) + serosO(*jj,me_A)*(1-pfp(*jj))*pfeA(*jj);
    seros(*jj,he2e1A) = serosO(*jj,he2e1A) + p_ffirstA(*jj) * serosO(*jj,s__)*pmeA(*jj)*pfeA(*jj) +
			(serosO(*jj,f_ea1A) + serosO(*jj,f_ea2A))*(1-pmp_ac(*jj))*pmeA(*jj) + serosO(*jj,f_eA)*(1-pmp(*jj))*pmeA(*jj);
    seros(*jj,mb_a1) = 0 ;
    seros(*jj,mb_a2) = serosO(*jj,mb_a1) * (1-pfe(*jj))*(1-pfp_ac(*jj));
    seros(*jj,mb_)   = serosO(*jj,mb_a2) * (1-pfe(*jj))*(1-pfp_ac(*jj)) + serosO(*jj,mb_) * (1-pfe(*jj))*(1-pfp(*jj)) ;
    seros(*jj,me_a1) = serosO(*jj,s__) * pme(*jj)*(1-pfe(*jj));
    seros(*jj,me_a2) = serosO(*jj,me_a1)*(1-pfe(*jj))*(1-pfp_ac(*jj)) ;
    seros(*jj,me_)   = serosO(*jj,me_a2)*(1-pfe(*jj))*(1-pfp_ac(*jj)) + serosO(*jj,me_)*(1-pfe(*jj))*(1-pfp(*jj)) ;
    seros(*jj,f_ba1) = 0;
    seros(*jj,f_ba2) = serosO(*jj,f_ba1)*(1-pme(*jj))*(1-pmp_ac(*jj));
    seros(*jj,f_b)   = serosO(*jj,f_ba2)*(1-pme(*jj))*(1-pmp_ac(*jj)) + serosO(*jj,f_b)*(1-pme(*jj))*(1-pmp(*jj))            ;
    seros(*jj,f_ea1) = serosO(*jj,s__)*pfe(*jj)*(1-pme(*jj));
    seros(*jj,f_ea2) = serosO(*jj,f_ea1)*(1-pme(*jj))*(1-pmp_ac(*jj));
    seros(*jj,f_e)   = serosO(*jj,f_ea2)*(1-pme(*jj))*(1-pmp_ac(*jj)) + serosO(*jj,f_e)*(1-pme(*jj))*(1-pmp(*jj));
    //  hb1b2/b2b1 not here b/c don't change during marriage;
    seros(*jj,hbe)  = serosO(*jj,hbe) + (serosO(*jj,mb_a1)+ serosO(*jj,mb_a2))*(1-pfp_ac(*jj))*pfe(*jj) + serosO(*jj,mb_)*(1-pfp(*jj))*pfe(*jj);
    seros(*jj,heb)  = serosO(*jj,heb) + (serosO(*jj,f_ba1)+ serosO(*jj,f_ba2))*(1-pmp_ac(*jj))*pme(*jj) + serosO(*jj,f_b)*(1-pmp(*jj))*pme(*jj);
    seros(*jj,hbpa) = serosO(*jj,hbpa)+ (serosO(*jj,mb_a1)+ serosO(*jj,mb_a2))*pfp_ac(*jj);
    seros(*jj,hbp)  = serosO(*jj,hbp) + serosO(*jj,mb_)*pfp(*jj);
    seros(*jj,hpba) = serosO(*jj,hpba)+ (serosO(*jj,f_ba1)+ serosO(*jj,f_ba2))*pmp_ac(*jj);
    seros(*jj,hpb)  = serosO(*jj,hpb) + serosO(*jj,f_b)*pmp(*jj);
    seros(*jj,hepa) = serosO(*jj,hepa)+ (serosO(*jj,me_a1)+ serosO(*jj,me_a2))*pfp_ac(*jj);
    seros(*jj,hep)  = serosO(*jj,hep) + serosO(*jj,me_)*pfp(*jj);
    seros(*jj,hpea) = serosO(*jj,hpea)+ (serosO(*jj,f_ea1)+ serosO(*jj,f_ea2))*pmp_ac(*jj);
    seros(*jj,hpe)  = serosO(*jj,hpe) + serosO(*jj,f_e)*pmp(*jj);
    denom(*jj)=pme(*jj)+pfe(*jj);
    if(denom(*jj)==0) {
      p_mfirst(*jj) = 0;
      p_ffirst(*jj) = 0;
    }else{
      p_mfirst(*jj) = pme(*jj) / denom(*jj); 
      p_ffirst(*jj) = 1-p_mfirst(*jj);
    }
    seros(*jj,he1e2) = serosO(*jj,he1e2) + p_mfirst(*jj) * serosO(*jj,s__)*pme(*jj)*pfe(*jj) +
		      (serosO(*jj,me_a1) + serosO(*jj,me_a2))*(1-pfp_ac(*jj))*pfe(*jj) + serosO(*jj,me_)*(1-pfp(*jj))*pfe(*jj);
    seros(*jj,he2e1) = serosO(*jj,he2e1) + p_ffirst(*jj) * serosO(*jj,s__)*pme(*jj)*pfe(*jj) +
		      (serosO(*jj,f_ea1) + serosO(*jj,f_ea2))*(1-pmp_ac(*jj))*pme(*jj) + serosO(*jj,f_e)*(1-pmp(*jj))*pme(*jj);


  }  // end active couple loop
  } // end couple duration loop
  //////////////////////////////////////////////////
  // Get likelihood from serostatuses
  LogicalVector impossible(numcouples);
  NumericMatrix pser(numcouples,4);
  NumericVector probs(numcouples); // probability each couple had their observed serostatus
  for(int jj = 0; jj < numcouples; jj++) { 
    pser(jj,0) = seros(jj,hb1b2A) + seros(jj,hb2b1A) + seros(jj,hbeA) + seros(jj,hebA) + seros(jj,hbpaA) + seros(jj,hpbaA) // ++ (1)
      + seros(jj,hepaA) + seros(jj,hpeaA) + seros(jj,hbpA) + seros(jj,hpbA) + seros(jj,hepA) + seros(jj,hpeA) + seros(jj,he1e2A) + seros(jj,he2e1A); // +- (2)
    pser(jj,1) = seros(jj,mb_a1A) + seros(jj,me_a1A) + seros(jj,mb_a2A) + seros(jj,me_a2A) + seros(jj,mb_A) + seros(jj,me_A); // -+ (3)
    pser(jj,2) = seros(jj,f_ba1A) + seros(jj,f_ea1A) + seros(jj,f_ba2A) + seros(jj,f_ea2A) + seros(jj,f_bA) + seros(jj,f_eA); // -- (4)
    pser(jj,3) = seros(jj,s__);
    // datser gives which of the serostatuses the couple wa observed in, we divide by the total
    // because this may not equal 1 since probability flows into mortality states too
    probs(jj)  = pser(jj, datser(jj)-1) / (pser(jj,1)+pser(jj,2)+pser(jj,3)+pser(jj,4)); // -1 deals with R to C indexing
    impossible(jj) = probs(jj)==0;
  }

  double lprob=0;
  bool do_likelihood_sums;
  do_likelihood_sums=is_true(any(impossible));
  if(do_likelihood_sums) {		// If any have zero probabilitie return -Inf
    lprob = -INFINITY;
  }else{			// Otherwise
    lprob=log(dnormC(lrho, log(trans_ratio), lrho_sd)); // take log prior
    for(int jj = 0; jj < numcouples; jj++) {
      lprob += log(probs(jj));	// Add each couple's log probability
  	}
  }
  
  // } // end useless loop 

   return Rcpp::List::create(_["lprob"] = lprob, _["seros"] = seros);
}


// NumericMatrix mcmcsampler(NumericMatrix sd_props = sd_props, Rcpp::NumericVector inits, 
//                     acute_sc,
//                     multiv = F, covar = NULL, # if multiv, sample from multivariate distribution (calculated during adaptive phase)
//                     verbose = T, tell = 100, seed = 1, lrho_sd = 1/2,
//                     niter = 6*1000, survive = T, uncond_mort = F,
//                     keep_seros = F, ## trace all serostate probabilities
//                     nthin = 1,
//                     nburn = 1000, browse=F)

// 			  const int & max_bd, const int & max_cd,
// 			  const NumericMatrix & pre_fprev, const NumericMatrix & pre_mprev, const NumericMatrix & within_fprev, const NumericMatrix & within_mprev,
// 			  const NumericMatrix & pre_msurv, const NumericMatrix & pre_fsurv, const NumericMatrix & within_msurv, const NumericMatrix & within_fsurv,
// 			  const NumericMatrix & within_art_cov,
// 			  const List & PreActiveList, const List & WithinActiveList,
// 			  const double &  bmb, const double &  bfb, const double &  bme, const double &  bfe, const double &  bmp, const double &  lrho,
// 			  const double & acute_sc, int uselessnum,
// 			  const bool & partner_arv, const double & cov_scalar) {

// }




// sampler <- function(sd.props = sd.props, inits, dat,
//                     acute.sc,
//                     multiv = F, covar = NULL, # if multiv, sample from multivariate distribution (calculated during adaptive phase)
//                     verbose = T, tell = 100, seed = 1, lrho.sd = 1/2,
//                     niter = 6*1000, survive = T, uncond.mort = F,
//                     keep.seros = F, ## trace all serostate probabilities
//                     nthin = 1,
//                     nburn = 1000, browse=F)

//   {
//     if(browse)  browser()
//     set.seed(seed)
//     pars <- inits
//     vv <- 2
//     accept <- 0 ## track each parameters acceptance individually
//     cur <- pcalc(pars, acute.sc = acute.sc, dat = dat, survive = survive, uncond.mort = uncond.mort, lrho.sd = lrho.sd)     #calculate first log probability
//     lprob.cur <- cur$lprob
//     out <- t(as.matrix(c(pars, bfp = as.numeric(pars["bmp"]*exp(pars["lrho"])))))
//     if(keep.seros)      seros.out <- cur$seros
//     last.it <- 0
//     start <- Sys.time()
//     while(vv < niter + 1) {
//         if(verbose & vv%%tell+1==1) print(paste("on iteration",vv,"of",last.it + niter + 1))
//         pars.prop <- pars              #initialize proposal parameterr vector
//         ## propose new parameter vector
//         if(multiv)          {
//             pars.prop <- pars.prop + rmnorm(1, mean = 0, varcov = covar)
//             pars.prop <- as.vector(pars.prop) #otherwise is a matrix
//             names(pars.prop) <- parnames
//           }else{
//             pars.prop <- pars.prop + rnorm(length(pars), mean = 0, sd = sd.props)
//           }
//         ## trace = T if in non-thinned iteration, or the previous one (in case of rejection)
//         ## calculate proposal par log probability
//         prop <- pcalc(pars.prop, acute.sc = acute.sc, dat = dat, survive = survive, uncond.mort = uncond.mort, lrho.sd = lrho.sd)
//         lprob.prop <- prop$lprob
//         lmh <- lprob.prop - lprob.cur       # log Metropolis-Hastings ratio
//         ## if MHR >= 1 or a uniform random # in [0,1] is <= MHR, accept otherwise reject
//         if(lmh >= 0 | runif(1,0,1) <= exp(lmh)) {
//             pars <- pars.prop
//             if(vv>nburn) accept <- accept + 1 #only track acceptance after burn-in
//             lprob.cur <- lprob.prop
//             cur <- prop
//           }
//         if(vv%%nthin + 1 ==1) {
//             out <- rbind(out,t(as.matrix(c(pars, bfp = as.numeric(pars["bmp"]*exp(pars["lrho"]))))))
//             if(keep.seros)      seros.out <- abind(seros.out, cur$seros, along = 3)
//         }
//         vv <- vv+1
//     }
//     if(verbose) print(paste("took", difftime(Sys.time(),start, units = "mins"),"mins"))
//     aratio <- accept/((vv-nburn))
//     give <- 1:nrow(out)>(nburn+1)/nthin
//     if(keep.seros) seros.out <- seros.out[,,give] else seros.out <- NULL
//     return(list(out = out[give,], aratio = aratio, inits = inits, seros.arr=seros.out))
// }



// [[Rcpp::export]]
NumericMatrix cplC(const int & max_bd,
		       const NumericMatrix & pre_fprev, const NumericMatrix & pre_msurv, 
		       const List & PreActiveList,	
		       const double &  bmb) {
  // Column names for sero
  int s__ = 0; int mb_a1A = 1; int mb_a2A = 2; // 53 state variables for each individual
  int numcouples = pre_fprev.nrow(); // get number of couples from one of the input matrices
  NumericMatrix seros(numcouples, 3); // Initialize couples by state matrix
  for(int jj = 0; jj < numcouples; ++jj) seros(jj,s__) = 1; // All couples start in first state
  // Initialize transmission hazards & probabilities & joint probability of transmission & survival (*a*live)
  int numactive(1);	      // # of active couples (changes between time steps)
  NumericVector m_haz(numcouples); // intermediate couple-specific parameters
  NumericVector pmb(numcouples);
  NumericVector pmbA(numcouples);
// For each time step up until max time
  NumericMatrix serosO = clone(seros); 
  for(int tt = 0; tt < max_bd; ++tt) { 
    serosO = seros; // Copy old seros (last state from which to iterate forward)
    IntegerVector active = PreActiveList[tt]; // Currently active couples,
    numactive = active.size();		 // Number of active couples
    if (min(active) < 0) stop("active < 0"); // Check indices are contained within arrays
    if (max(active) > numcouples) stop("max(active) > numcouples");
    // For each active couple
    for(IntegerVector::iterator jj = active.begin(); jj != active.end(); ++jj) {
      // ///////////////////////////////////////////////////
      // Assign temporary parameters: about 40 lines of code like these 3 lines
      m_haz(*jj) = bmb*pre_fprev(*jj,tt); // temporary parameter 1
      pmb(*jj) = 1 - exp(-m_haz(*jj)); // temporary parameter 2
      pmbA(*jj) = pmb(*jj)*pre_msurv(*jj,tt); // temporary parameter 3
      // /////////////////////////////////////////////////
      // Update state variables: about 120 lines of code like these 3 lines
      seros(*jj,s__)    = serosO(*jj,s__)* (1-pmbA(*jj));
      seros(*jj,mb_a1A)  = serosO(*jj,s__)*pmbA(*jj)*(1-pmb(*jj));
      seros(*jj,mb_a2A) = serosO(*jj,mb_a1A)*pmbA(*jj)*pmb(*jj);
    } // end couple iteration
  }   // end month or pre-couple duration iteration
return seros;
}


/*** R

*/

// numcouples <- 10^c(2:5)
// times <- numeric(length(numcouples))
// for(cc in 1:length(numcouples)) {
//     nn <- numcouples[cc]
//     maxT <- 500
//     examplist <- list(NA)
//     for(ii in 1:maxT) examplist[[ii]] <- sample(0:(nn-1), max(101-ii,1))
//     pre_fprev <- matrix(runif(nn*maxT),nn,maxT)
//     pre_msurv <- matrix(runif(nn*maxT),nn,maxT)
//     bmb <-  .05
//     times[cc] <- system.time(cplC(max_bd = maxT, pre_fprev = pre_fprev, pre_msurv = pre_msurv, PreActiveList = examplist, bmb = bmb))[3]
// }
// cbind(numcouples,times, times/numcouples)
