/**
* @mainpage
*
* Space-efficient sparse variant of an RNA (loop-based) free energy
* minimization algorithm (RNA folding equivalent to the Zuker
* algorithm).
*
* The results are equivalent to RNAfold -d0.
*
* Demonstration of space-efficient sparsification with trace-back.
*
* Since many matrix entries can not be efficiently recomputed in
* trace back, we store trace arrows to such entries. To save space,
* trace arrows are gc'ed and trace arrows to candidates are omitted
* and reconstructed in trace back.
*
* ----------------------------------------
* Specific recursions:

W(i,j) = min { W(i,j-1);
				min_i<k<j  W(i,k-1) + V(k,j) <-- W (same i), CLW;
				V(i,j);
		0 if i>=j-m
		}

V(i,j) = min { HairpinE(i,j);
		min_kl V(i,j)+ILoopE(i,j,k,l) <-- TAs;
		WM2(i+1,j-1) + a <-- WM2, no TAs;
		}

WM(i,j) = min { V(i,j)+b      <-- candidate in recomp;
		WM(i,j-1) + c <-- ! not via candidate list;
				min_i<k<j (k-i)*c + V(k,j) + b <-- CLWM   ( trick to save trace arrows );
				min_i<k<j  WM(i,k-1) + V(k,j) + b  <-- WM, CLWM;
		INF if i>=j-m
		}

WM2(i,j) = min{ WM2(i,j-1) + c;
				min_i<k<j  WM(i,k-1) + V(k,j) + b }  <-- WM, CLWM, no TAs;

* ----------------------------------------
* Candidate criteria:
*
(i,j) is a candidate for the split in W if
V(i,j)      < min {
					W(i,j-1);
			min_i<k<j  W(i,k-1) + V(k,j)
					}

(i,j) is a candidate for the split in WM if
V(i,j) + b  < min {
					WM(i,j-1)+c;
			min_i<k<j (k-i)*c + V(k,j) + b;
			min_i<k<j  WM(i,k-1) + V(k,j) + b
					}

*
* For simplicity and space savings, we store all candidates that
* meet either criterion in the same list.
*/

#include <iostream>
#include <iomanip>

#include <sstream>

#include <LocARNA/matrix.hh>

#include <limits>

#include <vector>
#include <iterator>

#include <cstring>
#include <string>
#include <cassert>
#include <numeric>
#include <tuple>

#include "base.hh"
#include "trace_arrow.hh"

extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/params/io.h"
}

#include "cmdline.hh"
#include <omp.h>

#include <fstream>


typedef unsigned short int cand_pos_t;

struct triplet
{
    cand_pos_t first; 
    energy_t second;
    energy_t third;
	triplet(){
		first = 1;
		second = 2;
		third = 3;
	}
	triplet(cand_pos_t x, energy_t y , energy_t z){
		first = x;
		second = y;
		third = z;
	}
};

typedef std::pair<cand_pos_t,energy_t> cand_entry_t;
typedef std::vector< cand_entry_t > cand_list_t;

typedef triplet cand_entry_td1;
typedef std::vector< cand_entry_td1 > cand_list_td1;

class SparseMFEFold;

namespace unrolled {
template < typename Iterator, class Functor >
void for_each( Iterator start, Iterator end, Functor f, size_t num_times ) {
    for(Iterator cur = start; cur != end; cur += num_times ) {
        for( int i = 0; i < num_times; ++i ) {
            f( *(cur + i) );
        }
    }
}
}

typedef struct sparse_features{
	int pair;
	char type; // 'H' 'S' 'M' 'I', and 'O' and 'U' for outer unpaired base and unpaired base inside pair 

	int last_j;
	int in_pair;
	sparse_features(){
		pair=-2;
		type = 'N';
		last_j = -1;
		in_pair = -1;
	}

} sparse_features;



energy_t ILoopE(auto const& S_,auto const& S1_, auto const& params_, int ptype_closing,size_t i, size_t j, size_t k,  size_t l);
energy_t MbLoopE(auto const& S_, auto const& params_, int ptype_closing,size_t i, size_t j);
energy_t Mlstem(auto const& S_, auto const& params_, int ptype_closing,size_t i, size_t j);
void trace_V(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e,sparse_features *fres);
void trace_W(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j,sparse_features *fres);
void trace_WM(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e,sparse_features *fres) ;
void trace_WM2(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,size_t i, size_t j,sparse_features *fres);

bool evaluate_restriction(int i, int j, sparse_features *fres);

/**
* Space efficient sparsification of Zuker-type RNA folding with
* trace-back. Provides methods for the evaluation of dynamic
* programming recursions and the trace-back.
*/
class SparseMFEFold {

public:
	std::string seq_;
	size_t n_;

	short *S_;
	short *S1_;

	paramT *params_;

	std::string structure_;
	std::string restricted_;
	

	bool garbage_collect_;

	LocARNA::Matrix<energy_t> V_; // store V[i..i+MAXLOOP-1][1..n]
	
	std::vector<energy_t> W_;
	std::vector<energy_t> WM_;
	std::vector<energy_t> WM2_;

	std::vector<energy_t> dmli1_; // WM2 from 1 iteration ago
	std::vector<energy_t> dmli2_; // WM2 from 2 iterations ago

	// Pseudoknot portion
	LocARNA::Matrix<energy_t> VP_; // store VP[i..i+MAXLOOP-1][1..n]
	std::vector<energy_t> WMB_;
	std::vector<energy_t> dwmbi_; // WMB from 1 iteration ago
	std::vector<energy_t> WMBP_;
	std::vector<energy_t> WI_;
	std::vector<energy_t> dwib1_; // WI from 1 iteration ago
	std::vector<energy_t> WIP_;
	std::vector<energy_t> dwip1_; // WIP from 1 iteration ago


	bool mark_candidates_;


	TraceArrows ta_;
	
	std::vector< cand_list_td1 > CL_;
	std::vector< cand_list_t > CLVP_j_;
	std::vector< cand_list_t > CLVP_i_;
	std::vector< cand_list_t > CLVPP_;
	std::vector< cand_list_t > CLWMB_;
	std::vector< cand_list_t > CLBE_;

	// Holds restricted info
	sparse_features *fres;
	int *B;
	int *b;
	

	/**
	candidate list for decomposition in W or WM

	@note Avoid separate candidate lists CLW and CLWM for split cases in W and
	WM to save even more space; here, this works after
	reformulating the recursions such that both split-cases recurse to
	V-entries. (compare OCTs)
	*/
	

	// compare candidate list entries by keys (left index i) in descending order
	struct {
	bool operator ()(const cand_entry_td1 &x, size_t y) const {
		return x.first > y;
	}
	}
	cand_comp;

	


	SparseMFEFold(const std::string &seq, bool garbage_collect, std::string restricted)
	: seq_(seq),
	n_(seq.length()),
	params_(scale_parameters()),
	ta_(n_),
		garbage_collect_(garbage_collect)
	{
	make_pair_matrix();

	S_ = encode_sequence(seq.c_str(),0);
	S1_ = encode_sequence(seq.c_str(),1);

	V_.resize(MAXLOOP+1,n_+1);
	W_.resize(n_+1,0);

	WM_.resize(n_+1,INF);

	WM2_.resize(n_+1,INF);

	dmli1_.resize(n_+1,INF);

	dmli2_.resize(n_+1,INF);

	// Pseudoknot portion

	VP_.resize(MAXLOOP+1,n_+1);
	WMB_.resize(n_+1,INF);
	dwmbi_.resize(n_+1,INF);
	WMBP_.resize(n_+1,INF);
	WI_.resize(n_+1,0);
	dwib1_.resize(n_+1,0);
	WIP_.resize(n_+1,INF);
	dwip1_.resize(n_+1,INF);

	// init candidate lists
	CL_.resize(n_+1);
	CLWMB_.resize(n_+1);
	CLVP_j_.resize(n_+1);
	CLVP_i_.resize(n_+1);
	CLVPP_.resize(n_+1);
	CLBE_.resize(n_+1);

	resize(ta_,n_+1);

	fres = new sparse_features[n_+1];
	B = (int*) malloc(sizeof(int)*n_+1);
	b = (int*) malloc(sizeof(int)*n_+1);

	

	restricted_ = restricted;
	
	}

	

	~SparseMFEFold() {
	free(params_);
	free(S_);
	free(S1_);
	delete [] fres;
	free(B);
	free(b);
	}
};


// ! TRANSLATED: -----------------------------------------------------------------------------------

energy_t HairpinE(auto const& seq, auto const& S, auto const& S1, auto const& params, size_t i, size_t j) {

	assert(1<=i);
	assert(i<j);
	
	//assert(j<=len); // don't know len here

	const int ptype_closing = pair[S[i]][S[j]];

	if (ptype_closing==0) return INF;

	return E_Hairpin(j-i-1,ptype_closing,S1[i+1],S1[j-1],&seq.c_str()[i-1], const_cast<paramT *>(params));
	}


/**
* @brief Rotate WM2 arrays to store the previous and previous previous iterations
* 
* @param WM WM array
* @param WM2 WM2 array
* @param dmli1 WM2 from one iteration ago
* @param dmli2 WM2 from two iterations ago
* @param n Length
*/
void rotate_arrays( auto const &WM2, auto &dmli1, auto &dmli2, auto const &WMB, auto &dwmbi, auto const &WI, auto &dwib1, auto const &WIP, auto &dwip1, auto n){
	

	for (int j = 1; j <= n; j++){
		dmli2[j] = dmli1[j];
		dmli1[j] = WM2[j];
		dwmbi[j] = WMB[j];
		dwib1[j] = WI[j];
		dwip1[j] = WIP[j];
	} 
}

/**
* @brief Computes the multiloop V contribution (in essence VM)
* 
* @param dmli1 Row of WM2 from one iteration ago
* @param dmli2 Row of WM2 from two iterations ago
* @param S Sequence Encoding
* @param params Parameters
* @param i Current i
* @param j Current j
* @param p_table Restricted Array
* @return energy_t 
*/
energy_t E_MbLoop(auto const& dmli1, auto const& dmli2, auto const& S, auto const& params, size_t i, size_t j, sparse_features *fres){

	int e = INF;
	int en = INF;
  	unsigned int tt;
	tt  = pair[S[j]][S[i]];

	/* double dangles */
	switch(params->model_details.dangles){
		case 2:
			if ((fres[i].pair <-1 && fres[j].pair <-1) || (fres[i].pair == j and fres[j].pair == i)) {
			e = dmli1[j - 1];

			if (e != INF) {

				int si1 = S[i + 1];
				int sj1 = S[j - 1];

				e += E_MLstem(tt, sj1, si1, params) + params->MLclosing;
				// e += E_MLstem(tt, -1, -1, params) + params->MLclosing;
			}

			}
			break;

		case 1:
			/**
			* ML pair D0
			*  new closing pair (i,j) with mb part [i+1,j-1]  
			*/
			if ((fres[i].pair <-1 && fres[j].pair <-1) || (fres[i].pair == j and fres[j].pair == i)) {
        		e = dmli1[j - 1];

        		if (e != INF) {

          			e += E_MLstem(tt, -1, -1, params) + params->MLclosing;

        		}
      		}
     		tt  = pair[S[j]][S[i]];
      		/** 
			* ML pair 5
			* new closing pair (i,j) with mb part [i+2,j-1] 
			*/
      		if (((fres[i].pair <-1 && fres[j].pair <-1) || (fres[i].pair == j and fres[j].pair == i)) && fres[i+1].pair < -1) {
        		en = dmli2[j - 1];

        		if (en != INF) {

          			int si1 =  S[i + 1];

          			en += E_MLstem(tt, -1, si1, params) + params->MLclosing + params->MLbase;
      
        		}
      		}
      		e   = MIN2(e, en);

			/** 
			* ML pair 3
			* new closing pair (i,j) with mb part [i+1, j-2] 
			*/
			if (((fres[i].pair <-1 && fres[j].pair <-1) || (fres[i].pair == j and fres[j].pair == i)) && fres[j-1].pair < 0) {
				en = dmli1[j - 2];

				if (en != INF) {
					int sj1 = S[j - 1];

					en += E_MLstem(tt, sj1, -1, params) + params->MLclosing + params->MLbase; 
				}
			}
			e   = MIN2(e, en);

			/** 
			* ML pair 53
			* new closing pair (i,j) with mb part [i+2.j-2]
			*/
			if (((fres[i].pair <-1 && fres[j].pair <-1) || (fres[i].pair == j and fres[j].pair == i)) && fres[i+1].pair < -1 && fres[j-1].pair <-1) {
				e = dmli2[j - 2];

				if (e != INF) {

					int si1 = S[i + 1];
					int sj1 = S[j - 1];

					e += E_MLstem(tt, sj1, si1, params) + params->MLclosing + 2 * params->MLbase;
				}
			}
			e   = MIN2(e, en);
      		break;
		case 0:
			if ((p_table[i] <-1 && p_table[j] <-1) || (p_table[i] == j and p_table[j] == i)) {
				e = dmli1[j - 1];

				if (e != INF) {
					e += E_MLstem(tt, -1, -1, params) + params->MLclosing;
				}
			}
			break;  
	}


	return e;
}
/**
* @brief Computes the Multiloop WM contribution 
* 
* @param vkj V at k and j
* @param vk1j V at k+1 and j - INF if not Candidate in traceback
* @param vkj1 V at k and j-1 - INF if not Candidate in traceback
* @param vk1j1 V at k+1 and j-1 - INF if not Candidate in traceback
* @param WM WM array
* @param CL Candidate List
* @param S Sequence Encoding
* @param params Parameters
* @param i Current i
* @param j Current j
* @param n Length
* @param p_table Restricted array
* @return energy_t 
*/
energy_t E_MLStem(auto const& vkj,auto const& vk1j,auto const& vkj1,auto const& vk1j1, auto const& WM, auto const& CL,auto const& S, auto const& params,size_t i, size_t j, auto const& n, sparse_features *fres){

	int e = INF,en=INF;

	int type = pair[S[i]][S[j]];
	


	if ((fres[i].pair < -1 && fres[j].pair < -1) || (fres[i].pair == j && fres[j].pair == i)) {
		en = vkj;
		if (en != INF) {
			if (params->model_details.dangles == 2)
				en += E_MLstem(type, (i == 1) ? S[n] : S[i - 1], S[j + 1], params);
			else
				en += E_MLstem(type, -1, -1, params);

			e = MIN2(e, en);
		}
	}

	if(params->model_details.dangles == 1){
		int mm5 = S[i], mm3 = S[j];
		if (((fres[i].pair < -1 && fres[j].pair < -1) || (fres[i].pair == j && fres[j].pair == i)) && fres[i+1].pair < -1) {
      		en = (j-(i+1) >TURN+1) ? vk1j : INF;
      		if (en != INF) {
        		en += params->MLbase;

            	type = pair[S[i+1]][S[j]];
            	en += E_MLstem(type, mm5, -1, params);

        		e = MIN2(e, en);
      		}
    	}

		if (((fres[i].pair < -1 && fres[j].pair < -1) || (fres[i].pair == j && fres[j].pair == i)) && fres[j-1].pair < -1) {
      		en = (j-1-i>TURN+1) ? vkj1 : INF; 
      		if (en != INF) {
       			en += params->MLbase;

            	type = pair[S[i]][S[j-1]];
            	en += E_MLstem(type, -1, mm3, params);
 
        		e = MIN2(e, en);
      		}
    	}

    	if (((fres[i].pair < -1 && fres[j].pair < -1) || (fres[i].pair == j && fres[j].pair == i)) && fres[i+1].pair < -1 && fres[j-1].pair<-1) {
      		en = (j-1-(i+1)>TURN+1) ? vk1j1 : INF; // i+1 j-1
      		if (en != INF) {
        		en += 2 * params->MLbase;

        		type = pair[S[i+1]][S[j-1]];
        		en += E_MLstem(type, mm5, mm3, params);
        
				e = MIN2(e, en);
      		}
    	} 
		
	}


    return e;
}


auto const recompute_WI(auto &WI, auto const &CL, auto const &CLWMB, auto const& S, auto const &params, auto const& n, size_t i, size_t max_j, sparse_features *fres) {
	

	assert(i>=1);
	assert(max_j<=n);

	
	
	for ( size_t j=i; j<=max_j; j++ ) {
		energy_t wi = 0;
		bool paired;
		#pragma omp parallel for num_threads(6);
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			size_t k = it->first;
			if(k<i) continue;
			paired = (fres[k].pair == j && fres[j].pair == k);
			const energy_t v_kj = it->second + params->bp_penalty;
			wi = std::min( wi, WI[k-1]  + v_kj );
			if(paired) break;
		}

		for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i ; ++it ) {
			if(paired) continue;
			size_t k = it->first;
			if(k<i) continue;
			paired = (fres[k].pair == j && fres[j].pair == k);
			const energy_t wmb_kj = it->second + params->bp_penalty + params->PPS_penalty;
			wi = std::min( wi, WI[k-1]  + wmb_kj );	
		}
		if(wi == 0) WI[j] = (j-i)*params->PUP_penalty;
		WI[j] = wi;
		
	}
	return WI;
}

auto const recompute_WIP(auto &WIP, auto const &CL, auto const &CLWMB, auto const& S, auto const &params, auto const& n, size_t i, size_t max_j, sparse_features *fres) {
	

	assert(i>=1);
	assert(max_j<=n);

	
	
	for ( size_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wip = INF;
		bool paired;
		#pragma omp parallel for num_threads(6);
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			size_t k = it->first;
			if(k<i) continue;
			paired = (fres[k].pair == j && fres[j].pair == k);
			const energy_t v_kj = it->second + params->bp_penalty;
			bool can_pair = true;
			for(int m = i;m<k;++m){
				if(fres[m].pair>-1){
					can_pair = false;
					break;
				} 
		
			}
			if(can_pair) wip = std::min( wip, static_cast<energy_t>(params->MLbase*(k-i)) + v_kj );
			wip = std::min( wip, WIP[k-1]  + v_kj );
			if(paired) break;
		}

		for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i ; ++it ) {
			if(paired) continue;
			size_t k = it->first;
			if(k<i) continue;
			paired = (fres[k].pair == j && fres[j].pair == k);
			const energy_t wmb_kj = it->second + params->bp_penalty + params->PPS_penalty;
			wip = std::min( wip, static_cast<energy_t>(params->MLbase*(k-i)) + wmb_kj );
			wip = std::min( wip, WIP[k-1]  + wmb_kj );	
		}
		if(fres[j].pair<0) wip = std::min(wip, WIP[j-1] + params->MLbase);
		WIP[j] = wip;
	}
	return WIP;
}



/**
* @brief Recompute row of WM 
* 
* @param WM WM array
* @param CL Candidate List
* @param S Sequence Encoding
* @param params Parameters
* @param n length
* @param i Current i
* @param max_j Current j
* @param p_table Restricted array
* @return auto const 
*/
auto const recompute_WM(auto const& WM, auto const &CL, auto const& S, auto const &params, auto const& n, size_t i, size_t max_j, sparse_features *fres) {
	

	assert(i>=1);
	assert(max_j<=n);

	std::vector<energy_t> temp = WM;

	for ( size_t j=i-1; j<=std::min(i+TURN,max_j); j++ ) { temp[j]=INF; }
	
	for ( size_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wm = INF;
		bool paired;
		int mm3 = S[j-1];
		#pragma omp parallel for num_threads(6);
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			size_t k = it->first;
			paired = (fres[k].pair == j && fres[j].pair == k);
			int mm5 = S[k+1];
			const energy_t v_kj = it->third;
			bool can_pair = true;
			for(int m = i;m<k;++m){
				if(fres[m].pair>-1){
					can_pair = false;
					break;
				} 
		
			}
			if(can_pair) wm = std::min( wm, static_cast<energy_t>(params->MLbase*(k-i)) + v_kj );
			wm = std::min( wm, temp[k-1]  + v_kj );
			if(paired) break;
		}
		if(fres[j].pair<0) wm = std::min(wm, temp[j-1] + params->MLbase);
		temp[j] = wm;
	}
	return temp;
}

/**
* @brief Recompute row of WM2 
* 
* @param WM WM array
* @param WM2 WM2 array
* @param CL Candidate List
* @param S Sequence Encoding
* @param params parameters
* @param n length
* @param i current i
* @param max_j current j
* @param p_table restricted array
* @param last_j_array restricted array
* @param in_pair_array restricted array
* @return auto const 
*/
auto const recompute_WM2(auto const& WM, auto const& WM2, auto const CL, auto const& S, auto const &params, auto const& n, size_t i, size_t max_j, sparse_features *fres) {
	

	assert(i>=1);
	//assert(i+2*TURN+3<=max_j);
	assert(max_j<= n);

	std::vector<energy_t> temp = WM2;

	for ( size_t j=i-1; j<=std::min(i+2*TURN+2,max_j); j++ ) { temp[j]=INF; }

	#pragma omp parallel for num_threads(6);
	for ( size_t j=i+2*TURN+3; j<=max_j; j++ ) {
		energy_t wm2 = INF;
		bool paired;
		int mm3 = S[j-1];
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>i+TURN+1 ; ++it ) {
			
			size_t k = it->first;
			paired = (fres[k].pair == j && fres[j].pair == k);
			int mm5 = S[k+1];
			energy_t v_kl = it->third;
			wm2 = std::min( wm2, WM[k-1]  + v_kl );
			if(paired) break;
		}
		if(fres[j].pair<0) wm2 = std::min(wm2, temp[j-1] + params->MLbase);
		// if(evaluate_restriction(i,j,last_j_array,in_pair_array)) wm2=INF;
		temp[j] = wm2;
	}
	return temp;
}

/**
 * @brief Test existence of candidate
 * 
 * @param CL Candidate List
 * @param cand_comp 
 * @param i start
 * @param j end
 * @return true 
 * @return whether (i,j) is candidate for W/WM splits 
 */
bool is_candidate(auto const& CL,auto const& cand_comp,size_t i, size_t j) {
	const cand_list_td1 &list = CL[j];

	auto it = std::lower_bound(list.begin(),list.end(),i,cand_comp);

	return it!=list.end() && it->first==i;
	}

/**
 * @brief Trace from W entry
 * 
 * @param seq Sequence
 * @param CL Candidate List
 * @param cand_comp Candidate Comparator
 * @param structure Final structure
 * @param params Parameters
 * @param S Sequence Encoding
 * @param S1 Sequence Encoding
 * @param ta trace arrows
 * @param W W array
 * @param WM WM array
 * @param WM2 WM2 array
 * @param n Length
 * @param mark_candidates Whether candidates are marked as [ ]
 * @param i row index
 * @param j column index
 * @param p_table Restricted Array
 * @param last_j_array Restricted Array
 * @param in_pair_array Restricted Array
 * pre: W contains values of row i in interval i..j
 */
void trace_W(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j,sparse_features *fres) {
	if (i+TURN+1>=j) return;
	// case j unpaired
	if (W[j] == W[j-1]) {
		trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,i,j-1,fres);
		return;
	}
	
	size_t k=j+1;
	energy_t v=INF;
	int sj1 = (j<n) ? S[j+1] : -1;
	energy_t w;
	for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i;++it ) {
		k = it->first;
		int sk1 = (k>1) ? S[k-1] : -1;
		const energy_t v_kj = it->second + vrna_E_ext_stem(pair[S[k]][S[j]],sk1,sj1,params);
		w = W[k-1] + v_kj;
		
		if (W[j] == w) {
		v = it->second;
		break;
		}
	}

	assert(i<=k && k<j);
	assert(v<INF);

	// don't recompute W, since i is not changed
	trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,i,k-1,fres);
	trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,j,v,fres);
}

/**
* @brief Trace from V entry
* 
* @param seq Sequence
* @param CL Candidate List
* @param cand_comp Candidate Comparator
* @param structure Final Structure
* @param params Parameters
* @param S Sequence Encoding
* @param S1 Sequence Encoding
* @param ta Trace Arrows
* @param WM WM array
* @param WM2 WM2 array
* @param n Length
* @param mark_candidates Whether Candidates should be [ ]
* @param i row index
* @param j column index
* @param e energy in V[i,j]
* @param p_table Restricted Array
* @param last_j_array Restricted Array
* @param in_pair_array Restricted Array
* pre: structure is string of size (n+1)
*/
void trace_V(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e,sparse_features *fres) {
	assert( i+TURN+1<=j );
	assert( j<=n );

	if (mark_candidates && is_candidate(CL,cand_comp,i,j)) {
		structure[i]='{';
		structure[j]='}';
	} else {
		structure[i]='(';
		structure[j]=')';
	}
	const int ptype_closing = pair[S[i]][S[j]];

	if (exists_trace_arrow_from(ta,i,j)) {
		// trace arrows may exist for interior loop case
		const TraceArrow &arrow = trace_arrow_from(ta,i,j);

		const size_t k=arrow.k(i,j);
		const size_t l=arrow.l(i,j);
		assert(i<k);
		assert(l<j);
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l, arrow.target_energy(),fres);
		return;

	} else {

		assert(ptype_closing>0);

		// try to trace back to a candidate: (still) interior loop case
		for ( size_t l=i; l<j; l++) {
		for ( auto it=CL[l].begin(); CL[l].end()!=it && it->first>i; ++it ) {
			const size_t k=it->first;
			if (  e == it->second + ILoopE(S,S1,params,ptype_closing,i,j,k,l) ) {
				trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,it->second,fres);
			return;
			}
		}
		}
	}
	
	// is this a hairpin?
	if ( e == HairpinE(seq,S,S1,params,i,j) ) {
		return;
	}
	
	// if we are still here, trace to wm2 (split case);
	// in this case, we know the 'trace arrow'; the next row has to be recomputed
	auto const temp = recompute_WM(WM,CL,S,params,n,i+1,j-1,fres);
	WM = temp;
	auto const temp2 = recompute_WM2(WM,WM2,CL,S,params,n,i+1,j-1,fres);
	WM2 = temp2;
	
	trace_WM2(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i+1,j-1,fres);
}

/**
* @brief Trace from WM
* 
* @param seq Sequence
* @param CL Candidate List
* @param cand_comp Candidate Comparator
* @param structure Final Structure
* @param params Parameters
* @param S Sequence Encoding
* @param S1 Sequence Encoding
* @param ta Trace Arrows
* @param WM WM array
* @param WM2 Wm2 array
* @param n Length
* @param mark_candidates Whether Candidates should be [ ]
* @param i row index
* @param j column index 
* @param e energy in WM[i,j] 
* @param p_table Restricted array
* @param last_j_array Restricted array
* @param in_pair_array Restricted array
* @param dangles Determines Multiloop Contribution
* pre: vector WM is recomputed for row i
*/
void trace_WM(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,size_t i, size_t j, energy_t e, sparse_features *fres) {
	if (i+TURN+1>j) {return;}

	if ( e == WM[j-1] + params->MLbase ) {
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,j-1,WM[j-1],fres);
		return;
	}
	int mm3 = S[j-1];
	for ( auto it=CL[j].begin();CL[j].end() != it && it->first>=i;++it ) {
		const size_t k = it->first;
		int mm5 = S[k+1];
		const energy_t v_kj = it->third;
		if ( e == WM[k-1] + v_kj ) {
		// no recomp, same i
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,k-1,WM[k-1],fres);
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,j,it->second,fres);
		return;
		} else if ( e == static_cast<energy_t>((k-i)*params->MLbase) + v_kj ) {
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,j,it->second,fres);
		return;
		}
	}
	assert(false);
}

/**
* @brief Trace from WM2
* 
* @param seq Sequence
* @param CL Candidate List
* @param cand_comp Candidate Comparator
* @param structure Final Structure
* @param params Parameters
* @param S Sequence Encoding
* @param S1 Sequence Encoding
* @param ta Trace Arrows
* @param WM WM array
* @param WM2 Wm2 array
* @param n Length
* @param mark_candidates Whether Candidates should be [ ]
* @param i row index
* @param j column index
* @param p_table Restricted array
* @param last_j_array Restricted array
* @param in_pair_array Restricted array
* pre: vectors WM and WM2 are recomputed for row i
 */
void trace_WM2(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,size_t i, size_t j,sparse_features *fres) {
	if (i+2*TURN+3>j) {return;}

	const energy_t e = WM2[j];

	// case j unpaired
	if ( e == WM2[j-1] + params->MLbase ) {
		
		// same i, no recomputation
		trace_WM2(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,j-1,fres);
		return;
	}
	int mm3 = S[j-1];
	for ( auto it=CL[j].begin();CL[j].end() != it  && it->first>=i+TURN+1;++it ) {
		size_t k = it->first;
		int mm5 = S[k+1];
		const energy_t v_kj = it->third;
		if ( e == WM[k-1] + v_kj ) {
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,k-1,WM[k-1],fres);
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,j,it->second,fres);
		return;
		}
	}
	assert(false);
}
/**
* @brief Trace back
* pre: row 1 of matrix W is computed
* @return mfe structure (reference)
*/
const std::string & trace_back(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n,sparse_features *fres,auto const& mark_candidates=false) {

	structure.resize(n+1,'.');

	/* Traceback */
	trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,1,n,fres);
	structure = structure.substr(1,n);

	return structure;
}

/* pre: ptype_closing>0 */
energy_t ILoopE(auto const& S, auto const& S1, auto const& params, int ptype_closing,size_t i, size_t j, size_t k,  size_t l)  {
	assert(ptype_closing>0);
	assert(1<=i);
	assert(i<k);
	assert(k<l);
	assert(l<j);
	//assert(l<=len); // don't know len here

	// note: enclosed bp type 'turned around' for lib call
	const int ptype_enclosed = rtype[pair[S[k]][S[l]]];

	if (ptype_enclosed==0) return INF;

	return E_IntLoop(k-i-1,j-l-1,ptype_closing,ptype_enclosed,S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params));
}


/**
* @brief Register a candidate
* @param i start
* @param j end
* @param e energy of candidate "V(i,j)"
* @param ml energy + ml contribution
*/
void register_candidatetd1(auto &CL, size_t i, size_t j, energy_t e, energy_t ml) {
	// assert(i<=j+TURN+1);
	CL[j].push_back( cand_entry_td1(i, e, ml) );
}
/**
* @brief Register a candidate
* @param i start
* @param j end
* @param e energy of candidate "V(i,j)"
*/
void register_candidate(auto &CL, size_t i, size_t j, energy_t e) {
	// assert(i<=j+TURN+1);
	CL[j].push_back( cand_entry_t(i, e) );
}
/**
 * @brief Get the outer right side of the band
 * 
 * @param B 
 * @param fres 
 * @param l 
 * @param j 
 * @return int 
 */
int getB(int *B, sparse_features *fres, int l, int j){
    if(fres[l].pair>-1 || fres[l].in_pair == 0) return -2;
    if((fres[l].in_pair<fres[j].in_pair && fres[l].last_j > fres[j].last_j) || (fres[l].in_pair==fres[j].in_pair  && fres[l].last_j >= fres[j].last_j && (fres[j].pair < 0 || j< fres[j].pair))) return -1;
    if(j<=l) return 0;
   
    int temp = fres[B[j]].pair;
    while(temp>l){
       temp = fres[B[temp]].pair;
    }
    return fres[temp].pair; 
}
/**
 * @brief Get the outer left side of the band
 * 
 * @param b 
 * @param fres 
 * @param i 
 * @param l 
 * @return int 
 */
int getb(int *b, sparse_features *fres, int i, int l){
    if(fres[l].pair>-1 || fres[l].in_pair == 0) return -2;
    if((fres[i].in_pair>fres[l].in_pair && fres[fres[l].last_j].pair < fres[fres[i].last_j].pair) || (fres[i].in_pair==fres[l].in_pair && fres[fres[l].last_j].pair <= fres[fres[i].last_j].pair && (fres[i].pair < 0 || i> fres[i].pair))) return -1;
    if(i>=l) return 0;
    // cout << "here" << endl;
    int temp = fres[b[i]].pair;
    while(temp<l){
       temp = fres[b[temp]].pair;
    }
    return fres[temp].pair; 
}
/**
 * @brief Get the inner left side of the band
 * 
 * @param fres 
 * @param i 
 * @param l 
 * @return int 
 */
int getbp(sparse_features *fres, int i, int l){
    if(fres[l].pair>-1 || fres[l].in_pair == 0) return -2;
    // if((fres[i].in_pair>fres[l].in_pair && fres[fres[l].last_j].pair < fres[fres[i].last_j].pair) || (fres[i].in_pair==fres[l].in_pair && fres[l].last_j >= fres[i].last_j && (fres[i].pair < 0 || i> fres[i].pair))) return -1;
    if(i>fres[fres[l].last_j].pair) return -1;
	if(i>=l) return 0;

    return fres[fres[l].last_j].pair;
}
/**
 * @brief Get the inner right side of the band
 * 
 * @param fres 
 * @param l 
 * @param j 
 * @return int 
 */
int getBp(sparse_features *fres, int l, int j){
    if(fres[l].pair>-1 || fres[l].in_pair == 0) return -2;
    // if((fres[l].in_pair<fres[j].in_pair && fres[l].last_j > fres[j].last_j) || (fres[l].in_pair==fres[j].in_pair  && fres[l].last_j >= fres[j].last_j && (fres[j].pair < 0 || j< fres[j].pair))) return -1;
    if(j<fres[l].last_j) return -1;
	// if(j<=l) return 0;
    return fres[l].last_j;
}
/**
 * @brief Set the array for B which getB is based off of
 * 
 * @param structure 
 * @param B 
 */
void setB(std::string structure, int *B){
   int n = structure.length();
   int prev_j = -1;
   for(int i = 1;i<=n;++i){
       if(structure[i-1] == ')'){
            prev_j = i;
            B[i] = prev_j;
       }
       else if(prev_j != -1) B[i] = prev_j;
       else B[i] = -1;  
   }  
}
/**
 * @brief Set the b array which the getb is based off of
 * 
 * @param structure 
 * @param b 
 */
void setb(std::string structure, int *b){
   int n = structure.length();
   int prev_j = -1;
   for(int i = n;i>=1;--i){
       if(structure[i-1] == '('){
            prev_j = i;
            b[i] = prev_j;
       }
       else if(prev_j != -1) b[i] = prev_j;
       else b[i] = -1;    
   } 
}
/**
 * @brief Find if [i,j] is empty
 * 
 * @param fres 
 * @param B 
 * @param b 
 * @param i 
 * @param j 
 * @return int 
 */
int is_empty_region(sparse_features* fres, int *B, int *b, int i, int j){
    if(j<i) return 0;
    int B_j = B[j];
    int b_i = b[i];
    if((b_i>j || b_i == -1) && B_j<i) return 1;
    return 0;
}

/**
 * @brief Finds whether [i,j] is a weakly closed region
 * 
 * @param fres 
 * @param B 
 * @param b 
 * @param i 
 * @param j 
 * @return int 
 */
int is_weakly_closed(sparse_features* fres, int *B, int *b, int i, int j){
    if(j<i) return 0;
    int B_j = B[j];
    int b_i = b[i];
    bool b_ij = b_i>j || b_i == -1;
    bool B_ji = B_j<i;
    if(b_ij && B_ji) return 1;
    else if(b_ij) return 0;
    else if(B_ji) return 0;
    else{
        if(fres[B_j].pair < i || fres[B_j].pair > j || fres[b_i].pair > j || fres[b_i].pair < i) return 0;
        int temp = B_j;
        while(temp>=i && temp != -1){
            if(fres[temp].pair<i) return 0;
            temp = B[fres[temp].pair];
        }
        temp = b_i;
        while(temp<=j && temp != -1){
            if(fres[temp].pair>j) return 0;
            temp = b[fres[temp].pair];
        }
        return 1;
    }
}

// std::pair< energy_t, energy_t > split_cases( auto const& CL, auto const& WM, auto const& S, auto const& params, int i, int j, auto &km1, int n, sparse_features *fres) {
// 	energy_t wm_split = INF;
// 	energy_t wm2_split = INF;
// 	int mm3 = S[j-1];

// 	for ( auto const [key,val] : CL) {
// 		size_t k = key;
// 		int mm5 = S[k+1];
// 		bool paired = (fres[k].pair == j && fres[j].pair == k);
// 		energy_t v_kj = E_MLStem(val,INF,INF,INF,WM,CL,S,params,k,j,n,fres);
// 		wm_split = std::min( wm_split, WM[k-1] + v_kj );
// 		bool can_pair = true;
// 		// checks to see if the unpaired bases till k can happen
// 		for(int m = i;m<k;++m){
// 			if(fres[m].pair>-1) {
// 				can_pair = false;
// 				break;
// 			}
// 		}
// 		if(can_pair) wm_split = std::min( wm_split,static_cast<energy_t>((k-i)*params->MLbase) + v_kj );
// 		wm2_split = std::min( wm2_split, WM[k-1] + v_kj );
// 		if(wm2_split==WM[k-1] + v_kj) km1 = k-1;
		
// 		if(paired) return std::make_pair( wm_split, wm2_split );
// 	}	
// 	return std::make_pair( wm_split, wm2_split );

// }
/**
 * @brief Evaluates whether a pairing can occur based on the restriction
 * 
 * @param i Current i
 * @param j Current j
 * @param p_table Restricted array
 * @param last_j_array Restricted array
 * @param in_pair_array Restricted array
 * @param multiloop Boolean to check if we are looking at WM and WM2
 * @return whether i and j can be non INF 
 */
bool evaluate_restriction(int i, int j, sparse_features *fres){
	bool evaluate = 1;
	if(fres[i].in_pair>fres[j].in_pair) evaluate = 0;

	if(fres[i].in_pair<fres[j].in_pair) evaluate = 0;

	if(fres[i].in_pair==fres[j].in_pair){
		if(j>fres[i].last_j) evaluate = 0;
	}
	// Resolves the cases where k-1 is the end of a restricted pair but i is less than the beginning of the k-1 pair
	// And where i is the beginning of the restricted pair but k-1 is past the end of the pair 
	// if(multiloop){
		if((fres[j].pair >0 && i<fres[j].pair && j>fres[j].pair) || (fres[i].pair>0 && j > fres[i].pair && i<fres[i].pair)) evaluate = 1;
	// }
	return evaluate;
}

energy_t fold(auto const& seq, auto &V, auto const& cand_comp, auto &CL, auto &CLWMB,auto &CLVP_j,auto &CLVP_i,auto &CLVPP, auto &CLBE, auto const& S, auto const& S1, auto const& params, auto &ta, auto &W, auto &WM, auto &WM2, auto &dmli1, auto &dmli2, auto &VP, auto &WMB, auto &dwmbi,auto &WMBP,auto &WI,auto &dwibi,auto &WIP, auto &dwip1, auto const& n, auto const& garbage_collect, sparse_features *fres, int *B, int *b) {
	for (size_t i=n; i>0; --i) {
		int si1 = (i>1) ? S[i-1] : -1;
		int mm5 = S[i+1];
		energy_t VP_i_split = INF;
		for(size_t j=i;j<i+TURN+1 && j<=n ;++j){
			WI[j] = (j-i+1)*params->PUP_penalty;
		}
		for ( size_t j=i+TURN+1; j<=n; j++ ) {

			int sj1 = (j<n) ? S[j+1] : -1;
			
			int mm3 = S[j-1];
			bool evaluate = evaluate_restriction(i,j,fres);
			// ------------------------------
			// W: split case
			bool pairedkj = 0;
			energy_t w_split = INF;
			energy_t wi_split = INF;
			energy_t wip_split = INF;
			energy_t wm_split = INF;
			energy_t wm2_split = INF;
			int km1 = n;
			for ( auto const [key,val, val_ml] : CL[j] ) {
				
				size_t k=key;
				if(!evaluate_restriction(i,k-1,fres)) continue;
				int sk1 = (k>1) ? S[k-1] : -1;
				bool unpairedkj = (fres[k].pair<-1 && fres[j].pair<-1);
				pairedkj = (fres[k].pair == j && fres[j].pair == k);

				bool can_pair = true;

				
				// checks to see if the unpaired bases till k can happen
				for(int m = i;m<k;++m){
					if(fres[m].pair>-1) {
						can_pair = false;
						break;
					}
				}
				
				if(pairedkj){
					
					// WM Portion
					energy_t v_kj = val_ml;
					wm_split =  WM[k-1] + v_kj;
					if(can_pair) wm_split = std::min( wm_split,static_cast<energy_t>((k-i)*params->MLbase) + v_kj );
					wm2_split = WM[k-1] + v_kj;
					km1 = k-1;
					//

					// WI portion
					v_kj = val + params->PPS_penalty;
					wi_split = WI[k] + v_kj;
					//

					// WIP portion
					v_kj = val + params->bp_penalty;
					wip_split = WIP[k]+v_kj;
					if(can_pair) wip_split = std::min(wip_split,static_cast<energy_t>((k-i)*params->cp_penalty) +v_kj);
					//


					// W portion
					v_kj = (unpairedkj || pairedkj) ? val + vrna_E_ext_stem(pair[S[k]][S[j]],sk1,sj1,params) : INF;
					w_split = W[k-1] + v_kj; 
					break;
				}else{
					//WM portion
					energy_t v_kj = E_MLStem(val,INF,INF,INF,WM,CL,S,params,k,j,n,fres);
					wm_split = std::min( wm_split, WM[k-1] + v_kj );
					if(can_pair) wm_split = std::min( wm_split,static_cast<energy_t>((k-i)*params->MLbase) + v_kj );
					wm2_split = std::min( wm2_split, WM[k-1] + v_kj );
					if(wm2_split==WM[k-1] + v_kj) km1 = k-1;
					//

					// WI portion
					v_kj = val + params->PPS_penalty;
					wi_split = std::min(wi_split,WI[k] + v_kj);
					//

					// WIP portion
					v_kj = val + params->bp_penalty;
					wip_split = std::min(wip_split,WIP[k]+v_kj);
					if(can_pair) wip_split = std::min(wip_split,static_cast<energy_t>((k-i)*params->cp_penalty) +v_kj);
					//

					// W portion
					v_kj = (unpairedkj || pairedkj) ? val + vrna_E_ext_stem(pair[S[k]][S[j]],sk1,sj1,params) : INF;
					w_split = std::min( w_split, W[k-1] + v_kj );
				}
				//
			}			 
			for (auto const [key,val] : CLWMB[j] ) {
				
				if(pairedkj) break;

				size_t k = key;
				bool can_pair = true;

				for(int m = i;m<k;++m){
					if(fres[m].pair>-1) {
						can_pair = false;
						break;
					}
				}
				
				
				// For W
				energy_t wmb_kj = val + params->PS_penalty;
				w_split = std::min( w_split, W[k-1] + wmb_kj );	
				// For WM
				wmb_kj = val + params->PSM_penalty + params->b_penalty;
				wm_split = std::min(wm_split, WM[k-1] + wmb_kj);
				if(can_pair) wm_split = std::min(wm_split,static_cast<energy_t>((k-i)*params->cp_penalty) +wmb_kj);
				wm2_split = std::min( wm2_split, WM[k-1] + wmb_kj );
				// For WI
				wmb_kj = val + params->PSP_penalty + params->PPS_penalty;
				wi_split = std::min(wi_split,WI[k] + wmb_kj);
				// For WIP
				wmb_kj = val + params->PSM_penalty + params->bp_penalty;
				wip_split = std::min(wip_split,WIP[k]+wmb_kj);
				if(can_pair) wip_split = std::min(wip_split,static_cast<energy_t>((k-i)*params->cp_penalty) +wmb_kj);
			}
			if(fres[j].pair<0){
			 	wm2_split = std::min( wm2_split, WM2[j-1] + params->MLbase );
			 	wm_split = std::min( wm_split, WM[j-1] + params->MLbase );
				w_split = std::min(w_split,W[j-1]);
				wi_split = std::min(wi_split,WI[j-1] + params->PUP_penalty);
				wip_split = std::min(wip_split,WIP[j-1] + params->cp_penalty);
			}
			

			// ------------------------------
			// WM and WM2: split cases
			
			// auto [wm_split, wm2_split] = split_cases( CL[j], WM,S, params,i,j,km1,n,fres);			
			
			// Check to see if wm and wm2 can be split
			// bool check = !(evaluate_restriction(i,km1,fres,true));
			// if(check && km1 != n) wm2_split=wm_split=INF;
			energy_t w  = w_split; // entry of W w/o contribution of V
			energy_t wm = wm_split; // entry of WM w/o contribution of V


			size_t i_mod=i%(MAXLOOP+1);

			const int ptype_closing = pair[S[i]][S[j]];
			const bool restricted = fres[i].pair == -1 || fres[j].pair == -1;
			bool unpaired = (fres[i].pair<-1 && fres[j].pair<-1);
			bool paired = (fres[i].pair == j && fres[j].pair == i);
			energy_t v = INF;
			energy_t w_v = INF;
			energy_t wm_v = INF;
			// ----------------------------------------
			// cases with base pair (i,j)
			if(ptype_closing>0 && !restricted && evaluate && j-i>3) { // if i,j form a canonical base pair

				bool canH = true;
				if((fres[i].pair>-1 && fres[i].pair != j) || (fres[j].pair>-1 && fres[j].pair != i)) canH = false;
				for(int k=i+1;k<j;k++) if(fres[k].pair>-1){canH = false;} // make more efficient later
				energy_t v_h = canH ? HairpinE(seq,S,S1,params,i,j) : INF;
				// info of best interior loop decomposition (if better than hairpin)
				size_t best_l=0;
				size_t best_k=0;
				energy_t best_e;

				energy_t v_iloop=INF;

				// constraints for interior loops
				// i<k; l<j
				// k-i+j-l-2<=MAXLOOP  ==> k <= MAXLOOP+i+1
				//            ==> l >= k+j-i-MAXLOOP-2
				// l-k>=TURN+1         ==> k <= j-TURN-2
				//            ==> l >= k+TURN+1
				// j-i>=TURN+3
				//
				size_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
				#pragma omp parallel for num_threads(6);
				for ( size_t k=i+1; k<=max_k; k++) {
					size_t k_mod=k%(MAXLOOP+1);

					size_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;

					for (size_t l=min_l; l<j; l++) {
						bool canI = true;

						
						for(int m = i+1; m<k;m++) if(fres[m].pair>-1){canI = false;}
						for(int m = l+1; m<j;m++) if(fres[m].pair>-1){canI = false;}
						if((fres[i].pair>-1 && fres[i].pair != j) || (fres[j].pair>-1 && fres[j].pair != i) || (fres[k].pair>-1 && fres[k].pair != l)) canI=false;

						assert(k-i+j-l-2<=MAXLOOP);

						const energy_t v_iloop_kl = canI ? V(k_mod,l) + ILoopE(S,S1,params,ptype_closing,i,j,k,l) : INF;
						if ( v_iloop_kl < v_iloop ) {
							v_iloop = v_iloop_kl;
							best_l=l;
							best_k=k;
							best_e=V(k_mod,l);
						}
					}
				}
				
				
				energy_t v_split = E_MbLoop(dmli1,dmli2,S,params,i,j,fres);
				// Look at case for WMB in VM
				// v_split = std::min(v_split,(dwmbi[j-1]+params->PSM_penalty+E_MLstem(ptype_closing,(i == 1) ? S[n] : S[i - 1], S[j + 1], params)));
				v = std::min(v_h,std::min(v_iloop,v_split));

				v_split = std::min(v_split,dwmbi[j-1]);

				size_t ip1_mod = (i+1)%(MAXLOOP+1);
				energy_t vk1j = V(ip1_mod,j);
				energy_t vkj1 = V(i_mod,j-1);
				energy_t vk1j1 = V(ip1_mod,j-1);

				w_v  = (unpaired || paired) ? v + vrna_E_ext_stem(ptype_closing,si1,sj1,params): INF;
				wm_v = (unpaired || paired) ? E_MLStem(v,vk1j,vkj1,vk1j1,WM,CL,S,params,i,j,n,fres): INF;
				// update w and wm by v
				if(paired){
					w = w_v;
					wm = wm_v;
				} else if(pairedkj){
					w = w_split;
					wm = wm_split;
				} else{
					w  = std::min({w_v, w_split});
					wm = std::min({wm_v, wm_split});
				}
				
				// register required trace arrows from (i,j)
				if ( v_iloop < std::min(v_h,v_split) ) {
					if ( is_candidate(CL,cand_comp,best_k,best_l) ) {
						//std::cout << "Avoid TA "<<best_k<<" "<<best_l<<std::endl;
						avoid_trace_arrow(ta);
					} else {
						//std::cout<<"Reg TA "<<i<<","<<j<<":"<<best_k<<","<<best_l<<std::endl;
						
						register_trace_arrow(ta,i,j,best_k,best_l,best_e);
					}
				}
				V(i_mod,j) = v;
			} else {
				V(i_mod,j) = INF;
			} // end if (i,j form a canonical base pair)
			

			
			int Bp_ij = getBp(fres,i,j);
			int B_ij = getB(B,fres,i,j);
			int b_ij = getb(b,fres,i,j);
			int bp_ij = getbp(fres,i,j);
			std::vector<energy_t> wiB1;
			std::vector<energy_t> wibp1;
			wiB1.resize(n+1,0);
			wibp1.resize(n+1,0);
			
			
			// Start of VP ---- Will have to change the bounds to 1 to n instead of 0 to n-1
			int weakly_closed_ij = is_weakly_closed(fres,B,b,i,j);
			if (i == j || weakly_closed_ij == 1 || fres[i].pair > -1 || fres[j].pair > -1 || ptype_closing == 0)	{
			
				VP(i_mod,j) = INF;
			}
			else{
				energy_t m1 = INF, m2 = INF, m3 = INF, m4= INF, m5 = INF, m6 = INF;
				if(fres[fres[i].last_j].pair > -1 && fres[fres[j].last_j].pair == -1 && Bp_ij >= 0 && B_ij >= 0){
					recompute_WI(wiB1,CL,CLWMB,S,params,n,B_ij+1,j-1,fres);
					// int WI_ipus1_BPminus = get_WI(i+1,Bp_i - 1) ;
					energy_t WI_ipus1_BPminus = dwibi[Bp_ij-1];
					// int WI_Bplus_jminus = get_WI(B_i + 1,j-1);
					energy_t WI_Bplus_jminus = wiB1[j-1];
					m1 =   WI_ipus1_BPminus + WI_Bplus_jminus;
				}
				// if(i==11 && j==24) std::cout << "b is " << b_ij << " and bp is " << bp_ij << std::endl; 
				if (fres[fres[i].last_j].pair > -1 && fres[fres[j].last_j].pair > -1 && b_ij>= 0 && bp_ij >= 0){
					recompute_WI(wibp1,CL,CLWMB,S,params,n,bp_ij+1,j-1,fres);
					// int WI_i_plus_b_minus = get_WI(i+1,b_i - 1);
					// if(i==11 && j==24) std::cout << i+1 << "    " << b_ij-1 << std::endl;
					energy_t WI_i_plus_b_minus = dwibi[b_ij-1];
					// int WI_bp_plus_j_minus = get_WI(bp_i + 1,j-1);
					energy_t WI_bp_plus_j_minus = wibp1[j-1];
					// if(i==10 && j==24) std::cout << WI_i_plus_b_minus << " " << WI_bp_plus_j_minus << std::endl;
					m2 = WI_i_plus_b_minus + WI_bp_plus_j_minus;
				}
				
				if(fres[fres[i].last_j].pair > -1 && fres[fres[j].last_j].pair > -1 && Bp_ij >= 0 && B_ij >= 0  && b_ij >= 0 && bp_ij>= 0){
					// int WI_i_plus_Bp_minus = get_WI(i+1,Bp_i - 1);
					energy_t WI_i_plus_Bp_minus = dwibi[Bp_ij-1];
					energy_t WI_B_plus_b_minus = wiB1[b_ij-1];
					energy_t WI_bp_plus_j_minus = wibp1[j-1];
					// int WI_B_plus_b_minus = get_WI(B_i + 1,b_i - 1);
					// int WI_bp_plus_j_minus = get_WI(bp_i +1,j - 1);
					m3 = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
				}
				if(fres[i+1].pair < -1 && fres[j-1].pair < -1){
					// m4 = get_e_stP(i,j)+ get_VP(i+1,j-1);
					int ip1_mod = (i+1)%(MAXLOOP+1);
					m4 = params->e_stP_penalty*ILoopE(S,S1,params,ptype_closing,i,j,i+1,j-1) + VP(ip1_mod,j-1);
				}
				

				int ip, jp;
				int max_borders;
				int min_borders = 1; 
				if (Bp_ij> 1 && Bp_ij < n && b_ij >1 && b_ij < n) min_borders = std::min(Bp_ij,b_ij);
				else if (b_ij > 1 && b_ij < n && (Bp_ij < 1 || Bp_ij > n)) min_borders = b_ij;
				else if (Bp_ij > 1 && Bp_ij < n && (b_ij < 1 || b_ij > n)) min_borders = Bp_ij;
				int edge_i = i+MAXLOOP+1;
				min_borders = std::min({min_borders,edge_i});
				
				for (ip = i+1; ip < min_borders; ip++){
					int empty_region_i = is_empty_region(fres,B,b,i+1,ip-1); // i+1 to ip-1
					int ip_mod = ip%(MAXLOOP+1);
					if (fres[ip].pair < -1 && (fres[fres[i].last_j].pair == fres[fres[ip].last_j].pair) && empty_region_i == 1){
						max_borders= 1;
						if (bp_ij > 1 && B_ij > 1) max_borders = std::max(bp_ij,B_ij);
						else if (B_ij > 1 && bp_ij < 1) max_borders = B_ij;
						else if (bp_ij > 1 && B_ij < 1) max_borders = bp_ij;
						int edge_j = j-31;
						max_borders = std::max({max_borders,edge_j});
						for (jp = max_borders+1; jp < j ; jp++){
							int empty_region_j = is_empty_region(fres,B,b,jp+1,j-1); // jp+1 to j-1
							// if(i==7 && j==27) std::cout << ip << " " << jp << " " << empty_region_j << std::endl;
							if (fres[jp].pair < -1 && pair[S[ip]][S[jp]]>0 && empty_region_j == 1){
								//arc to arc originally
								if (fres[j].last_j == fres[jp].last_j){
									energy_t temp = params->e_intP_penalty*ILoopE(S,S1,params,ptype_closing,i,j,ip,jp) + VP(ip_mod,jp);
									// if(i==7 && j==27 && ip ==9 && jp == 25) std::cout << VP(ip_mod,jp) << " " << params->e_intP_penalty*ILoopE(S,S1,params,ptype_closing,i,j,ip,jp) << " " << temp << std::endl;
									if (m5 > temp){
										m5 = temp;
									}
								}
							}
						}
					}
				}
				// if(i==2 && j==32) std::cout << m4 << " " << m5 << std::endl;
				// if(i==7 && j==27) std::cout << min_borders << " " << max_borders << std::endl;
				// int r;
				// int min_Bp_j = j;
				// if (Bp_ij > 0 && Bp_ij < n && Bp_ij < min_Bp_j) min_Bp_j = Bp_ij;
				// for (r = i+1; r < min_Bp_j ; r++){
				// 	if (fres[r].pair < -1){

				// 		// int tmp = get_WIP(i+1,r-1) + get_VPP(r,j-1) + ap_penalty + 2*bp_penalty;
				// 		// if (tmp < m6){
				// 		// 	m6 = tmp;
				// 		// }
				// 	}
				// }
				// std::cout << max_borders << std::endl;

				// case 6 and 7
				int left_bound = std::min(Bp_ij,b_ij);
				for ( auto const [key,val] : CLVPP[j-1] ) {
					size_t k = key;
					if(!(k>i && k<left_bound)) continue;
					// std::vector<energy_t> wipkm1;
					// wipkm1.resize(n+1,INF);
					// recompute_WIP(wipkm1,CL,CLWMB,S,params,n,i+1,k-1,fres);
					energy_t WIPVPP = dwip1[k-1] + val;
					// WIPVPP = std::min(WIPVPP,val + static_cast<energy_t>(k-i)*params->cp_penalty);
					m6 = std::min(m6, WIPVPP + params->ap_penalty + 2*params->bp_penalty);
				
				}
				
				
				VP(i_mod,j) = std::min({m1, m2, m3, m4, m5, m6});

				
				
				
				energy_t VP_j_split = INF;
				for ( auto const [key,val] : CLVP_j[j] ) {
					size_t k = key;
					VP_j_split = std::min(VP_j_split,val);
				}
				if(VP(i_mod,j)<VP_j_split){
					register_candidate(CLVP_j,i,j,VP(i_mod,j));
				}
				if(VP(i_mod,j)<VP_i_split){
					register_candidate(CLVP_i,j,i,VP(i_mod,j));
				}
				VP_i_split = std::min(VP_i_split,VP(i_mod,j));

				
			}
			// End of VP

			// Start of WMBP
			// int WMBP[n+1];
			if ((fres[i].pair >= -1 && fres[i].pair > j) || (fres[j].pair >= -1 && fres[j].pair < i) || (fres[i].pair >= -1 && fres[i].pair < i ) || (fres[j].pair >= -1 && j < fres[j].pair)) WMB[j] = INF;
			else{
				int m1 = INF, m3 = INF, m4 = INF, m5 = INF;
				if(fres[j].pair < 0 && fres[i].pair >= 0){
					int tmp = INF, l, l_min=-1;
					// for (l = i+1; l < j; l++){
					for ( auto const [key,val] : CLVP_j[j] ){
						size_t l = key;
						int bp_il = getbp(fres,i,l);
						if(bp_il >= 0 && bp_il < n && l+TURN <= j){
							energy_t BE_energy = INF;
							for ( auto const [key,val] : CLBE[bp_il] ){
								size_t k = key;
								if(i==k){
									BE_energy = val;
									break;
								}
							}
							std::vector<energy_t> wipbpl;
							wipbpl.resize(n+1,INF);
							recompute_WI(wipbpl,CL,CLWMB,S,params,n,bp_il,l-1,fres);
						// int BE_energy = get_BE(i,fres[i].pair,bp_i_l,fres[bp_i_l].pair);
						// int WI_energy = get_WI(bp_i_l +1,l-1);
						energy_t WI_energy = wipbpl[l-1];
						energy_t VP_energy = val;
						energy_t sum = BE_energy + WI_energy + VP_energy;
						if (tmp > sum){
							tmp = sum;
							l_min = l;
						}
						}
					}
					m1 = 2*params->PB_penalty + tmp;
				}

				// 3)
				if (fres[j].pair < 0){
					int temp = INF, l_min=-1;
					// for (l = i+1; l<j ; l++){
					for ( auto const [key,val] : CLVP_j[j] ){
						size_t l = key;
						int B_lj = getB(B,fres,l,j);
						int Bp_lj = getBp(fres,l,j);
						if (fres[fres[l].last_j].pair > -1 && B_lj >= 0 && B_lj < n && Bp_lj >= 0 && Bp_lj<n){
							if (b_ij >= 0 && b_ij < n && l < b_ij){
								if (i <= fres[fres[l].last_j].pair && fres[fres[l].last_j].pair < j && l+3 <=j){
									// int BE_energy = get_BE(fres[B_lj].pair,B_lj,fres[Bp_lj].pair,Bp_lj)
									energy_t BE_energy = INF;
									int b_lj = fres[B_lj].pair;
									for ( auto const [key,val] : CLBE[fres[Bp_lj].pair] ){
										size_t k = key;
										if(b_lj==k){
											BE_energy = val;
											break;
										}
									}
									energy_t WMBP_energy = WMBP[l-1];
									energy_t VP_energy = val;
									energy_t sum = BE_energy + WMBP_energy + VP_energy;
									if (temp > sum){
										temp = sum;
										l_min = l;
									}
								}
							}
						}
						m3 = 2*params->PB_penalty + temp;
					}
				}

				// 4) WMB(i,j) = VP(i,j) + P_b
				m4 = VP(i_mod,j) + params->PB_penalty;
				
				if(fres[j].pair < j){
					int l,l_min =-1;
					for(l = i+1; l<j; l++){
						if (fres[l].pair < 0 && fres[fres[l].last_j].pair > -1 && fres[fres[j].last_j].pair > -1 && fres[fres[j].last_j].pair == fres[fres[l].last_j].pair){
							energy_t WMBP_energy = WMBP[l];
							
							std::vector<energy_t> wipBpl;
							wipBpl.resize(n+1,INF);
							recompute_WI(wipBpl,CL,CLWMB,S,params,n,l+1,j,fres);
							energy_t WI_energy = wipBpl[j]; 
							// int temp = get_WMBP(i,l) + get_WI(l+1,j);
							energy_t temp = WMBP_energy + WI_energy;
							if (temp < m5){
								m5 = temp;
								l_min = l;
							}
						}
					}
				}

			// get the min for WMB
			WMBP[j] = std::min({m1,m3,m4,m5});

			// End of WMBP

			// Start of WMB

			
			int m2 = INF, mWMBP = INF;
			// 2)
			if (fres[j].pair >= 0 && j > fres[j].pair){
				int l, l_min=-1;
				int bp_j = fres[j].pair;
				int temp = INF;
				for (l = (bp_j +1); (l < j); l++){
					int Bp_lj = getBp(fres,l,j);
					if (Bp_lj >= 0 && Bp_lj<n){
						// energy_t BE_energy = get_BE(bp_j,j,fres[Bp_lj].pair,Bp_lj);
						energy_t BE_energy = INF;
						for ( auto const [key,val] : CLBE[fres[Bp_lj].pair] ){
							size_t k = key;
							if(bp_j==k){
								BE_energy = val;
								break;
							}
						}


						energy_t WMBP_energy = WMBP[l];
						std::vector<energy_t> wiplBp;
						wiplBp.resize(n+1,INF);
						recompute_WI(wiplBp,CL,CLWMB,S,params,n,l+1,Bp_lj-1,fres);
						energy_t WI_energy = wiplBp[Bp_lj-1];
						// energy_t WI_energy = get_WI(l+1,Bp_lj-1)
						energy_t sum = BE_energy + WMBP_energy + WI_energy;
						if (temp > sum){
							temp = sum;
							l_min = l;
						}

					}
				}
				m2 = params->PB_penalty + temp;
			}
			// check the WMBP_ij value
			mWMBP =  WMBP[j];

			// get the min for WMB
			WMB[j] = std::min(m2,mWMBP);
			

			// End of WMB
			}

			
			// Start of WI -- the conditions on calculating WI is the same as WIP, so we combine them
			int wi_v = INF;
			int wip_v = INF;
			int wi_wmb = INF;
			int wip_wmb = INF;
	
			if (weakly_closed_ij == 0 || fres[fres[i].last_j].pair != fres[fres[j].last_j].pair){
					WI[j] = INF;
					WIP[j] = INF;
			}
			else if(i==j){
				WI[j] = params->PUP_penalty;
			}
			else{
				if(ptype_closing>0 && !restricted && evaluate) {
					wi_v = V(i_mod,j) + params->PPS_penalty;
					wip_v = V(i_mod,j)	+ params->bp_penalty;
				}
				wi_wmb = WMB[j] + params->PSP_penalty + params->PPS_penalty;
				wip_wmb = WMB[j] + params->PSM_penalty + params->bp_penalty;
				
				energy_t WI_unp = (j-i+1)*params->PUP_penalty;
				WI[j] = std::min({wi_split,wi_v,wi_wmb,WI_unp});
				WIP[j] = std::min({wip_split,wip_v,wip_wmb});

				// if(i==12 && j==14) std::cout<<wi_split << " " << wi_v << " " << wi_wmb << " " << WI_unp << std::endl;
				
			}

			// start of VPP
			
			// This is for finding the previous VPP's for j i.e. the VPP split
			if(!is_weakly_closed(fres,B,b,i,j)){
				energy_t VPP_split = INF;
				for (auto const [key,val] : CLVPP[j] ) {
					size_t k = key;
					if(val < VPP_split) VPP_split = val;
				}
			

				// k is the j for VP but it's used for the i for WIP
				// This is for determing if this ij should be a VPP candidate - it gives the value for VPP[ij]
				energy_t VPP_ij = INF;
				int right_bound = std::max(bp_ij,B_ij);
				for (auto const [key,val] : CLVP_i[i] ) {
					size_t k = key;
					if(!(k<j && k>right_bound)) continue;
					std::vector<energy_t> wivp1;
					wivp1.resize(n+1,INF);
					recompute_WIP(wibp1,CL,CLWMB,S,params,n,k+1,j,fres);
					energy_t val_kj = val + wivp1[j]; 
					val_kj = std::min(val_kj, val + static_cast<energy_t>(params->cp_penalty*(j-k)));
					if(val_kj < VPP_ij) VPP_ij = val_kj;	 
				}

				if(VPP_ij<VPP_split){
					register_candidate(CLVPP, i, j, VPP_ij );
				}
			}
			// End of VPP


			// Start of BE
			int ip = fres[i].pair; // i's pair ip should be right side so ip = )
			int jp = fres[j].pair; // j's pair jp should be left side so jp = (
			// if (!( i >= 0 && i <= ip && ip < jp && jp <= j && j < n && fres[i].pair >= -1 && fres[j].pair >= -1 && fres[ip].pair >= -1 && fres[jp].pair >= -1 && fres[i].pair == j && fres[j].pair == i && fres[ip].pair == jp && fres[jp].pair == ip)){ //impossible cases
			
			// base case: i.j and ip.jp must be in G
			if (ip > i && j > jp && jp > i && ip > j){ // Don't need to check if they are pairs separately because it is checked by virtue of these checks

				// We are checking for the closest pair that we have already calculated to i/ip from j/jp 
				// If there is nothing, then i is the closest encompassing pair to jp
				// If it is not, then we get the energy for everything from jp to lp so that we calculate less
				energy_t BE_energy = INF;
				int lp = jp;
				for ( auto const [key,val] : CLBE[fres[jp].pair] ){
					size_t k = key;
					BE_energy = val;
					lp = k;
					
				}
				int l = fres[lp].pair; // right closing base for lp

				int m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;
				// if (fres[i+1].pair == ip-1){
				if(i+1 == lp and i-1 == ip){
					m1 = params->e_stP_penalty*ILoopE(S,S1,params,ptype_closing,i,ip,lp,l) + BE_energy;
				}

					// for (int k = i+1; k<= lp ; k++){
					// if (fres[l].pair >= -1 && j <= fres[k].pair && fres[k].pair < ip){
					
				int empty_region_ilp = is_empty_region(fres,B,b,i+1,lp-1);
				int empty_region_lip = is_empty_region(fres,B,b,l+1,ip-1);
				int weakly_closed_ilp = is_weakly_closed(fres,B,b,i+1,lp-1);
				int weakly_closed_lip = is_weakly_closed(fres,B,b,l+1,ip-1);

				std::vector<energy_t> wiilp;
				std::vector<energy_t> wilip;

				if(weakly_closed_ilp){
					wiilp.resize(n+1,INF);
					recompute_WIP(wiilp,CL,CLWMB,S,params,n,i+1,lp-1,fres);
				}
				if(weakly_closed_lip){
					wilip.resize(n+1,INF);
					recompute_WIP(wiilp,CL,CLWMB,S,params,n,l+1,ip-1,fres);
				}
					
				if (empty_region_ilp && empty_region_lip){
					m2 = params->e_intP_penalty*ILoopE(S,S1,params,ptype_closing,i,ip,lp,l) + BE_energy;
				}

					// 3)
				if (weakly_closed_ilp && weakly_closed_lip){
					// get_WIP(i+1,l-1)
					// get_WIP(lp+1,ip-1)
					m3 = wiilp[l-1] + BE_energy + wilip[ip-1]+ params->ap_penalty + 2*params->bp_penalty;
				}

					// 4)
				if (weakly_closed_ilp && empty_region_lip){
					m4 = wiilp[l-1] + BE_energy + params->cp_penalty * (ip-lp+1) + params->ap_penalty + 2*params->bp_penalty;
				}

					// 5)
				if (empty_region_ilp && weakly_closed_lip){
					m5 = params->ap_penalty + 2*params->bp_penalty + (params->cp_penalty * (l-i+1)) + BE_energy + wilip[ip-1];
				}
					// }
				register_candidate(CLBE,i,jp,std::min({m1,m2,m3,m4,m5}));
			}

				// finding the min and putting it in BE[iip]
				// BE[ip] = std::min({m1,m2,m3,m4,m5});
				
			// End of BE

			//Things that needed to happen later like W's wmb
			const energy_t w_wmb = (unpaired) ? WMB[j] + params->PS_penalty : INF;
			const energy_t wm_wmb = (unpaired) ? WMB[j] + params->PSM_penalty + params->b_penalty : INF;
			if(!pairedkj && !paired) {
				w =std::min(w,w_wmb);
				wm = std::min(wm,wm_wmb);
			}
			// check whether (i,j) is a candidate; then register
			if ( w_v < w_split || wm_v < wm_split || wi_v < wi_split || wip_v < wip_split || paired) {
		
				register_candidatetd1(CL, i, j, v, w_v);

				// always keep arrows starting from candidates
				inc_source_ref_count(ta,i,j);
			}
			if ( w_wmb < w_split || wm_wmb < wm_split || wi_wmb < wi_split || wip_wmb < wip_split) {
		
				register_candidate(CLWMB, i, j, WMB[j]);

				// always keep arrows starting from candidates
				inc_source_ref_count(ta,i,j);
			}
			W[j]       = w;
			WM[j]      = wm;
			WM2[j]     = wm2_split;
		
			
		} // end loop j
		rotate_arrays(WM2,dmli1,dmli2,WMB,dwmbi,WI,dwibi,WIP,dwip1,n);
		// Clean up trace arrows in i+MAXLOOP+1
		if (garbage_collect && i+MAXLOOP+1 <= n) {
			gc_row(ta,i + MAXLOOP + 1 );
		}

		// Reallocate candidate lists in i
		for ( auto &x: CL ) {
			if (x.capacity() > 1.5*x.size()) {
				cand_list_td1 vec(x.size());
				copy(x.begin(),x.end(),vec.begin());
				vec.swap(x);
			}
		}
		

		compactify(ta);
	}
	return W[n];
}

/**
 * @brief Fills the restriction arrays
 * p_table will contain the index of each base pair
 * X or x tells the program the base cannot pair and . sets it as unpaired but can pair
 * Pseudoknots (denoted by [ ], < >, or { } ) are filled the same way as ( )
 * That is, a structure like this (<)> is not possible.
 * @param structure Input structure
 * @param fres restriction struct
 */
void detect_restricted_pairs(auto const &structure, sparse_features *fres){
	int i, j, count = 0, length = structure.length(),last_j=length;
	std::vector<int>  pairs;
	pairs.push_back(length);

	for (i=length; i >=1; --i){
		if ((structure[i-1] == 'x') || (structure[i-1] == 'X'))
			fres[i].pair = -1;
		else if (structure[i-1] == '.')
			fres[i].pair = -2;
		if (structure[i-1] == ')' || structure[i-1] == ']' || structure[i-1] == '}' || structure[i-1] == '>'){
			pairs.push_back(i);
			count++;
		}
		fres[i].last_j = pairs[pairs.size()-1];
		fres[i].in_pair = count;
		if (structure[i-1] == '(' || structure[i-1] == '[' || structure[i-1] == '{' || structure[i-1] == '<'){
			j = pairs[pairs.size()-1];
			pairs.erase(pairs.end()-1);
			fres[i].pair = j;
			fres[j].pair = i;
			count--;
		}
	}
	pairs.pop_back();
	if (pairs.size() != 0)
	{
		fprintf (stderr, "The given structure is not valid: more left parentheses than right parentheses: \n");
		exit (1);
	}
}

/**
 * @brief Sums the number of Candidates at each index over all indices
 * 
 * @param CL_ Candidate list
 * @return total number of candidates
 */
size_t num_of_candidates(auto const& CL_)  {
	size_t c=0;
	for ( auto const &x: CL_ ) {
		c += x.size();
	}
	return c;
}
/**
 * @brief Finds the size of allocated storage capacity across all indices
 * 
 * @param CL_ Candidate List
 * @return the amount of allocated storage 
 */
size_t capacity_of_candidates(auto const& CL_) {
	size_t c=0;
	for ( auto const &x: CL_ ) {
		c += x.capacity();
	}
	return c;
}

/**
* @brief Simple driver for @see SparseMFEFold.
*
* Reads sequence from command line or stdin and calls folding and
* trace-back methods of SparseMFEFold.
*/
int
main(int argc,char **argv) {

	args_info args_info;

	// get options (call gengetopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	std::string seq;
	if (args_info.inputs_num>0) {
	seq=args_info.inputs[0];
	} else {
	std::getline(std::cin,seq);
	}
	int n = seq.length();

	std::string restricted;
    args_info.input_structure_given ? restricted = input_structure : restricted = std::string (n,'.');

	if(restricted != "" && restricted.length() != n ){
		std::cout << "input sequence and structure are not the same size" << std::endl;
		exit(0);
	}
	std::string file= "";
	args_info.paramFile_given ? file = parameter_file : file = "";
	if(file!=""){
		// FILE *fp;
    	// fp = fopen(parameter_file.c_str(),"r");
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
	}

	bool verbose;
	verbose = args_info.verbose_given;

	bool mark_candidates;
	mark_candidates = args_info.mark_candidates_given;

	SparseMFEFold sparsemfefold(seq,!args_info.noGC_given,restricted);

	if(args_info.dangles_given) sparsemfefold.params_->model_details.dangles = dangles;

	std::cout << seq << std::endl;
	
	// Psuedoknot-free setup
	detect_restricted_pairs(restricted,sparsemfefold.fres);
	// Pseudoknot setup
	setB(restricted,sparsemfefold.B);
	setb(restricted,sparsemfefold.b);
	energy_t mfe = fold(sparsemfefold.seq_,sparsemfefold.V_,sparsemfefold.cand_comp,sparsemfefold.CL_,sparsemfefold.CLWMB_,sparsemfefold.CLVP_j_,sparsemfefold.CLVP_i_,sparsemfefold.CLVPP_,sparsemfefold.CLBE_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.params_,sparsemfefold.ta_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_, sparsemfefold.dmli1_, sparsemfefold.dmli2_,sparsemfefold.VP_,sparsemfefold.WMB_,sparsemfefold.dwmbi_,sparsemfefold.WMBP_,sparsemfefold.WI_,sparsemfefold.dwib1_,sparsemfefold.WIP_,sparsemfefold.dwip1_,sparsemfefold.n_,sparsemfefold.garbage_collect_, sparsemfefold.fres,sparsemfefold.B,sparsemfefold.b);	
	
	
	std::string structure = trace_back(sparsemfefold.seq_,sparsemfefold.CL_,sparsemfefold.cand_comp,sparsemfefold.structure_,sparsemfefold.params_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.ta_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_,sparsemfefold.n_,sparsemfefold.fres, mark_candidates);
	
	
	std::ostringstream smfe;
	smfe << std::setiosflags(std::ios::fixed) << std::setprecision(2) << mfe/100.0 ;

	std::cout << structure << " ("<<smfe.str()<<")"<<std::endl;

	// float factor=1024;
	
	
	
	// const std::string unit=" kB";
	
	
	if (verbose) {
		

	std::cout <<std::endl;

	std::cout << "TA cnt:\t"<<sizeT(sparsemfefold.ta_)<<std::endl;
	std::cout << "TA max:\t"<<maxT(sparsemfefold.ta_)<<std::endl;
	std::cout << "TA av:\t"<<avoidedT(sparsemfefold.ta_)<<std::endl;
	std::cout << "TA rm:\t"<<erasedT(sparsemfefold.ta_)<<std::endl;

	std::cout <<std::endl;
	std::cout << "TAs num:\t"<<sizeT(sparsemfefold.ta_)<<std::endl;
	std::cout << "TAs cap:\t"<<capacityT(sparsemfefold.ta_)<<std::endl;
	std::cout<< std::endl;

	std::cout << "V Can num:\t"<<num_of_candidates(sparsemfefold.CL_)<<std::endl;
	std::cout << "VCan cap:\t"<<capacity_of_candidates(sparsemfefold.CL_)<<std::endl;
	std::cout << "VP_j Can num:\t"<<num_of_candidates(sparsemfefold.CLVP_j_)<<std::endl;
	std::cout << "VP_j Can cap:\t"<<capacity_of_candidates(sparsemfefold.CLVP_j_)<<std::endl;
	std::cout << "VP_i Can num:\t"<<num_of_candidates(sparsemfefold.CLVP_i_)<<std::endl;
	std::cout << "VP_i Can cap:\t"<<capacity_of_candidates(sparsemfefold.CLVP_i_)<<std::endl;
	std::cout << "VPP Can num:\t"<<num_of_candidates(sparsemfefold.CLVPP_)<<std::endl;
	std::cout << "VPP Can cap:\t"<<capacity_of_candidates(sparsemfefold.CLVPP_)<<std::endl;
	std::cout << "WMB Can num:\t"<<num_of_candidates(sparsemfefold.CLWMB_)<<std::endl;
	std::cout << "WMB Can cap:\t"<<capacity_of_candidates(sparsemfefold.CLWMB_)<<std::endl;
	std::cout << "BE Can num:\t"<<num_of_candidates(sparsemfefold.CLBE_)<<std::endl;
	std::cout << "BE Can cap:\t"<<capacity_of_candidates(sparsemfefold.CLBE_)<<std::endl;
	}

	// for(int j = 1;j<=n;++j){
	// 	for (auto const [key,val] : sparsemfefold.CLVP_j_[j] ) {
	// 	size_t k = key;
	// 	std::cout << k << " " << j << std::endl;
	// 	}
	// }
	

	return 0;
}
