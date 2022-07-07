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

#include "base.hh"
#include "trace_arrow.hh"

extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loops/all.h"
}

#include "cmdline.hh"
#include <omp.h>

#include <fstream>


typedef unsigned short int cand_pos_t;
typedef std::pair<cand_pos_t,energy_t> cand_entry_t;
typedef std::vector< cand_entry_t > cand_list_t;

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


	bool mark_candidates_;


	TraceArrows ta_;
	
	std::vector< cand_list_t > CL_;
	std::vector< cand_list_t > CLWMB_;

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
	bool operator ()(const cand_entry_t &x, size_t y) const {
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
	WI_.resize(n_+1,INF);
	dwib1_.resize(n_+1,INF);
	WIP_.resize(n_+1,INF);

	// init candidate lists
	CL_.resize(n_+1);
	CLWMB_.resize(n_+1);

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
void rotate_arrays(auto &WM, auto &WM2, auto &dmli1, auto &dmli2, auto &WMB, auto &dwmbi, auto &WI, auto dwib1, auto n){
	

	for (int j = 1; j <= n; j++){
		dmli2[j] = dmli1[j];
		dmli1[j] = WM2[j];
		dwmbi[j] = WMB[j];
		dwib1[j] = WI[j];
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
			tt  = pair[S[j]][S[i]];
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
		// case 3:
		// 	if ((p_table[i] <-1 && p_table[j] <-1) || (p_table[i] == j and p_table[j] == i)) {
		// 		e = dmli1[j - 1];

		// 		if (e != INF) {
		// 			e += E_MLstem(tt, -1, -1, params) + params->MLclosing;
		// 		}
		// 	}
		// 	break; 
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
		if ((fres[i+1].pair < -1 && fres[j].pair < -1) || (fres[i+1].pair == j && fres[j].pair == i+1 && fres[i].pair <-1)) {
      		en = vk1j;
      		if (en != INF) {
        		en += params->MLbase;

            	type = pair[S[i+1]][S[j]];
            	en += E_MLstem(type, mm5, -1, params);

        		e = MIN2(e, en);
      		}
    	}

		if ((fres[i].pair < -1 && fres[j-1].pair < -1) || (fres[i].pair == j-1 && fres[j-1].pair == i && fres[j].pair <-1)) {
      		en = vkj1; 
      		if (en != INF) {
       			en += params->MLbase;

            	type = pair[S[i]][S[j-1]];
            	en += E_MLstem(type, -1, mm3, params);
 
        		e = MIN2(e, en);
      		}
    	}

    	if ((fres[i+1].pair < -1 && fres[j-1].pair < -1) || (fres[i+1].pair == j-1 && fres[j-1].pair == i+1 && fres[j].pair < -1 && fres[i].pair < -1)) {
      		en = vk1j1; // i+1 j-1
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


auto const recompute_WIP(auto &WI, auto const &CL, auto const &CLWMB, auto const& S, auto const &params, auto const& n, size_t i, size_t max_j, sparse_features *fres) {
	

	assert(i>=1);
	assert(max_j<=n);

	
	
	for ( size_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wip = INF;
		bool paired;
		#pragma omp parallel for num_threads(6);
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			size_t k = it->first;
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
			wip = std::min( wip, WI[k-1]  + v_kj );
			if(paired) break;
		}

		for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i ; ++it ) {
			if(paired) continue;
			size_t k = it->first;
			paired = (fres[k].pair == j && fres[j].pair == k);
			const energy_t wmb_kj = it->second + params->bp_penalty + params->PPS_penalty;
			wip = std::min( wip, static_cast<energy_t>(params->MLbase*(k-i)) + wmb_kj );
			wip = std::min( wip, WI[k-1]  + wmb_kj );	
		}
		if(fres[j].pair<0) wip = std::min(wip, WI[j-1] + params->MLbase);
		WI[j] = wip;
	}
	return WI;
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
			const energy_t v_kj = E_MLStem(it->second,INF,INF,INF,WM,CL,S,params,k,j,n,fres);
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
			energy_t v_kl = E_MLStem(it->second,INF,INF,INF,WM,CL,S,params,k,j,n,fres);
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
	const cand_list_t &list = CL[j];

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
		const energy_t v_kj = E_MLStem(it->second,INF,INF,INF,WM,CL,S,params,k,j,n,fres);
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
		const energy_t v_kj = E_MLStem(it->second,INF,INF,INF,WM,CL,S,params,k,j,n,fres);
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
*/
void register_candidate(auto &CL, size_t i, size_t j, energy_t e) {
	assert(i<=j+TURN+1);
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
    if((fres[i].in_pair>fres[l].in_pair && fres[fres[l].last_j].pair < fres[fres[i].last_j].pair) || (fres[i].in_pair==fres[l].in_pair && fres[l].last_j >= fres[i].last_j && (fres[i].pair < 0 || i> fres[i].pair))) return -1;
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
    if((fres[l].in_pair<fres[j].in_pair && fres[l].last_j > fres[j].last_j) || (fres[l].in_pair==fres[j].in_pair  && fres[l].last_j >= fres[j].last_j && (fres[j].pair < 0 || j< fres[j].pair))) return -1;
    if(j<=l) return 0;
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

std::pair< energy_t, energy_t > split_cases( auto const& CL, auto const& WM, auto const& S, auto const& params, int i, int j, auto &km1, int n, sparse_features *fres) {
	energy_t wm_split = INF;
	energy_t wm2_split = INF;
	int mm3 = S[j-1];

	for ( auto const [key,val] : CL) {
		size_t k = key;
		int mm5 = S[k+1];
		bool paired = (fres[k].pair == j && fres[j].pair == k);
		energy_t v_kj = E_MLStem(val,INF,INF,INF,WM,CL,S,params,k,j,n,fres);
		wm_split = std::min( wm_split, WM[k-1] + v_kj );
		bool can_pair = true;
		// checks to see if the unpaired bases till k can happen
		for(int m = i;m<k;++m){
			if(fres[m].pair>-1) {
				can_pair = false;
				break;
			}
		}
		if(can_pair) wm_split = std::min( wm_split,static_cast<energy_t>((k-i)*params->MLbase) + v_kj );
		wm2_split = std::min( wm2_split, WM[k-1] + v_kj );
		if(wm2_split==WM[k-1] + v_kj) km1 = k-1;
		
		if(paired) return std::make_pair( wm_split, wm2_split );
	}	
	return std::make_pair( wm_split, wm2_split );

}
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
bool evaluate_restriction(int i, int j, sparse_features *fres, bool multiloop){
	bool evaluate = 1;
	if(fres[i].in_pair>fres[j].in_pair) evaluate = 0;

	if(fres[i].in_pair<fres[j].in_pair) evaluate = 0;

	if(fres[i].in_pair==fres[j].in_pair){
		if(j>fres[i].last_j) evaluate = 0;
	}
	// Resolves the cases where k-1 is the end of a restricted pair but i is less than the beginning of the k-1 pair
	// And where i is the beginning of the restricted pair but k-1 is past the end of the pair 
	if(multiloop){
		if((fres[j].pair >0 && i<fres[j].pair && j>fres[j].pair) || (fres[i].pair>0 && j > fres[i].pair && i<fres[i].pair)) evaluate = 1;
	}
	return evaluate;
}

energy_t fold(auto const& seq, auto &V, auto const& cand_comp, auto &CL, auto &CLWMB, auto const& S, auto const& S1, auto const& params, auto &ta, auto &W, auto &WM, auto &WM2, auto &dmli1, auto &dmli2, auto &VP, auto &WMB, auto &dwmbi,auto &WMBP,auto &WI,auto &dwibi,auto &WIP, auto const& n, auto const& garbage_collect, sparse_features *fres, int *B, int *b) {
	for (size_t i=n; i>0; --i) {
		int si1 = (i>1) ? S[i-1] : -1;
		for ( size_t j=i+TURN+1; j<=n; j++ ) {

			int sj1 = (j<n) ? S[j+1] : -1;
			int mm5 = S[i+1];
			int mm3 = S[j-1];
			bool evaluate = evaluate_restriction(i,j,fres,false);
			// ------------------------------
			// W: split case
			bool pairedkj = 0;
			energy_t w_split = INF;
			for ( auto const [key,val] : CL[j] ) {
				size_t k=key;
				int sk1 = (k>1) ? S[k-1] : -1;
				bool unpairedkj = (fres[k].pair<-1 && fres[j].pair<-1);
				pairedkj = (fres[k].pair == j && fres[j].pair == k);
				energy_t v_kj = (unpairedkj || pairedkj) ? val + vrna_E_ext_stem(pair[S[k]][S[j]],sk1,sj1,params) : INF;
				if(pairedkj){
					w_split = W[k-1] + v_kj; 
					break;
				}else{
					w_split = std::min( w_split, W[k-1] + v_kj );
				}
			}
			if(fres[j].pair<0) w_split = std::min(w_split,W[j-1]);

			// ------------------------------
			// WM and WM2: split cases
			int km1 = n;
			auto [wm_split, wm2_split] = split_cases( CL[j], WM,S, params,i,j,km1,n,fres);
			

			if(fres[j].pair<0) wm2_split = std::min( wm2_split, WM2[j-1] + params->MLbase );
			if(fres[j].pair<0) wm_split = std::min( wm_split, WM[j-1] + params->MLbase );
			
			
			// Check to see if wm and wm2 can be split
			bool check = !(evaluate_restriction(i,km1,fres,true));
			if(check && km1 != n) wm2_split=wm_split=INF;
			energy_t w  = w_split; // entry of W w/o contribution of V
			energy_t wm = wm_split; // entry of WM w/o contribution of V


			size_t i_mod=i%(MAXLOOP+1);

			const int ptype_closing = pair[S[i]][S[j]];
			const bool restricted = fres[i].pair == -1 || fres[j].pair == -1;

			// ----------------------------------------
			// cases with base pair (i,j)
			if(ptype_closing>0 && !restricted && evaluate) { // if i,j form a canonical base pair

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
				bool unpaired = (fres[i].pair<-1 && fres[j].pair<-1);
				bool paired = (fres[i].pair == j && fres[j].pair == i);
				
				energy_t v_split = E_MbLoop(dmli1,dmli2,S,params,i,j,fres);
				// Look at case for WMB in VM
				// v_split = std::min(v_split,(dwmbi[j-1]+params->PSM_penalty+E_MLstem(ptype_closing,(i == 1) ? S[n] : S[i - 1], S[j + 1], params)));
				const energy_t v = std::min(v_h,std::min(v_iloop,v_split));

				const energy_t w_v  = (unpaired || paired) ? v + vrna_E_ext_stem(ptype_closing,si1,sj1,params): INF;
				const energy_t wm_v = (unpaired || paired) ? E_MLStem(v,INF,INF,INF,WM,CL,S,params,i,j,n,fres): INF;
				
				// update w and wm by v
				if(paired){
					w = w_v;
					wm = wm_v;
				} else if(pairedkj){
					w = w_split;
					wm = wm_split;
				} else{
					w  = std::min(w_v, w_split);
					wm = std::min(wm_v, wm_split);
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
				// check whether (i,j) is a candidate; then register
				if ( w_v < w_split || wm_v < wm_split || paired) {
			
					register_candidate(CL, i, j, v );

					// always keep arrows starting from candidates
					inc_source_ref_count(ta,i,j);
				}
				V(i_mod,j) = v;
			} else {
				V(i_mod,j) = INF;
			} // end if (i,j form a canonical base pair)
			W[j]       = w;
			WM[j]      = wm;
			WM2[j]     = wm2_split;

			
			int Bp_ij = getBp(fres,i,j);
			int B_ij = getB(B,fres,i,j);
			int b_ij = getb(b,fres,i,j);
			int bp_ij = getbp(fres,i,j);
			std::vector<energy_t> wiB1;
			std::vector<energy_t> wibp1;
			wiB1.resize(n+1,INF);
			wibp1.resize(n+1,INF);
			
			
			const int ptype_closingp1 = pair[S[i+1]][S[j-1]];
			// Start of VP ---- Will have to change the bounds to 1 to n instead of 0 to n-1
			int weakly_closed_ij = is_weakly_closed(fres,B,b,i,j);
			if (i == j || weakly_closed_ij == 1 || fres[i].pair > -1 || fres[j].pair > -1 || ptype_closing == 0)	{
			
				VP(i_mod,j) = INF;
			}
			else{
				int m1 = INF, m2 = INF, m3 = INF, m4= INF, m5 = INF, m6 = INF;
				if(fres[fres[i].last_j].pair > -1 && fres[fres[j].last_j].pair == -1 && Bp_ij >= 0 && Bp_ij< n && B_ij >= 0 && B_ij < n){
					recompute_WIP(wiB1,CL,CLWMB,S,params,n,B_ij+1,j,fres);
					// int WI_ipus1_BPminus = get_WI(i+1,Bp_i - 1) ;
					int WI_ipus1_BPminus = dwibi[Bp_ij-1];
					// int WI_Bplus_jminus = get_WI(B_i + 1,j-1);
					int WI_Bplus_jminus = wiB1[j-1];
					m1 =   WI_ipus1_BPminus + WI_Bplus_jminus;
				}
				if (fres[fres[i].last_j].pair == -1 && fres[fres[j].last_j].pair > -1 && b_ij>= 0 && b_ij < n && bp_ij >= 0 && bp_ij < n){
					recompute_WIP(wibp1,CL,CLWMB,S,params,n,bp_ij+1,j,fres);
					// int WI_i_plus_b_minus = get_WI(i+1,b_i - 1);
					int WI_i_plus_b_minus = dwibi[b_ij-1];
					// int WI_bp_plus_j_minus = get_WI(bp_i + 1,j-1);
					int WI_bp_plus_j_minus = wibp1[j-1];
					m2 = WI_i_plus_b_minus + WI_bp_plus_j_minus;
				}
				if(fres[fres[i].last_j].pair > -1 && fres[fres[j].last_j].pair > -1 && Bp_ij >= 0 && Bp_ij < n && B_ij >= 0 && B_ij < n && b_ij >= 0 && b_ij < n && bp_ij>= 0 && bp_ij < n){
					// int WI_i_plus_Bp_minus = get_WI(i+1,Bp_i - 1);
					int WI_i_plus_Bp_minus = dwibi[Bp_ij-1];
					int WI_B_plus_b_minus = wiB1[b_ij-1];
					int WI_bp_plus_j_minus = wibp1[j-1];
					// int WI_B_plus_b_minus = get_WI(B_i + 1,b_i - 1);
					// int WI_bp_plus_j_minus = get_WI(bp_i +1,j - 1);
					m3 = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
				}
				if(fres[i+1].pair < -1 && fres[j-1].pair < -1 && ptype_closingp1>0){
					// m4 = get_e_stP(i,j)+ get_VP(i+1,j-1);
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
					if (fres[ip].pair < -1 && (fres[fres[i].last_j].pair == fres[fres[ip].last_j].pair) && empty_region_i == 1){
						max_borders= 1;
						if (bp_ij > 1 && bp_ij < n && B_ij > 1 && B_ij < n) max_borders = std::max(bp_ij,B_ij);
						else if (B_ij > 1 && B_ij < n && (bp_ij < 1 || bp_ij > n)) max_borders = B_ij;
						else if (bp_ij > 1 && bp_ij < n && (B_ij < 1 || B_ij > n)) max_borders = bp_ij;
						int edge_j = j-30;
						max_borders = std::max({max_borders,edge_j});
						for (jp = max_borders+1; jp < j ; jp++){
							int empty_region_j = is_empty_region(fres,B,b,jp+1,j-1); // jp+1 to j-1
							if (fres[jp].pair < -1 && pair[S[ip]][S[jp]]>0 && empty_region_j == 1){
								//arc to arc originally
								if (fres[j].last_j == fres[jp].last_j){
									int ip_mod = ip%(MAXLOOP+1);
									int temp = params->e_intP_penalty*ILoopE(S,S1,params,ptype_closing,i,j,ip,jp) + VP(ip_mod,jp);
									if (m5 > temp){
										m5 = temp;
									}
								}
							}
						}
					}
				}

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
				int max_i_bp = i;
				if (bp_ij > 0 && bp_ij < n && bp_ij > max_i_bp) max_i_bp = bp_ij;
				for(int l = max_i_bp; l<j;++l){
					for ( auto const [key,val] : CL[l] ) {
						size_t k=key;
						if(!is_weakly_closed(fres,B,b,k,l)){
							int upik = INF;
							if(is_empty_region(fres,B,b,i,k)) params->cp_penalty*(k-i);
							int uplj = INF;
							if(is_empty_region(fres,B,b,l,j)) params->cp_penalty*(j-l); // could perhaps optimize this by having it check just the next base after the first time instead of recalculating empty region
							std::vector<energy_t> WIlj;
							WIlj.resize(n+1,INF);
							recompute_WIP(WIlj,CL,CLWMB,S,params,n,l,j,fres);
							if(upik<dwibi[k] && uplj<WIlj[j-1]) break;
							int vp_kl = std::min(dwibi[k],upik) + val + std::min(WIlj[j-1],uplj) + params->ap_penalty + 2*params->bp_penalty;
							m6 = std::min(m6,vp_kl);	
						}
						
					}
				}

				// int max_i_bp = i;
				// if (bp_ij > 0 && bp_ij < n && bp_ij > max_i_bp) max_i_bp = bp_ij;
				// for (r = max_i_bp+1; r < j ; r++){
				// 	if (fres[r].pair < -1){
				// 		// int tmp = get_VPP(i+1,r) + get_WIP(r+1,j-1)+ ap_penalty + 2* bp_penalty;
				// 		// if (tmp < m7){
				// 		// 	m7 = tmp;
				// 		// }
				// 	}
				// }
				// I would think that if we wanted to limit VP case 6 and 7 (and VPP) to just 30 on either side as well
				// then we could combine them all into one for loop. As well, if we do this, we should be able to combine
				// case 6 and 8 into one case and do the calculation for the latter WIP through the use of candidates

				VP(i_mod,j) = std::min({m1, m2, m3, m4, m5, m6});
			}
			// End of VP

			// Start of WMBP
			int WMBP[n+1];
			if ((fres[i].pair >= -1 && fres[i].pair > j) || (fres[j].pair >= -1 && fres[j].pair < i) || (fres[i].pair >= -1 && fres[i].pair < i ) || (fres[j].pair >= -1 && j < fres[j].pair)) WMB[j] = INF;
			else{
				int m1 = INF, m3 = INF, m4 = INF, m5 = INF;
				if(fres[j].pair < 0 && fres[i].pair >= 0){
					int tmp = INF, l, l_min=-1;
					for (l = i+1; l < j; l++){
						int bp_il = getbp(fres,i,l);
						if(bp_il >= 0 && bp_il < n && l+TURN <= j){
						// int BE_energy = get_BE(i,fres[i].pair,bp_i_l,fres[bp_i_l].pair);
						// int WI_energy = get_WI(bp_i_l +1,l-1);
						// int VP_energy = get_VP(l,j);
						// int sum = BE_energy + WI_energy + VP_energy;
						// if (tmp > sum){
						// 	tmp = sum;
						// 	l_min = l;
						// }
						}
					}
					m1 = 2*params->PB_penalty + tmp;
				}

				// 3)
				if (fres[j].pair < 0){
					int l, temp = INF, l_min=-1;
					for (l = i+1; l<j ; l++){
						int B_lj = getB(B,fres,l,j);
						int Bp_lj = getBp(fres,l,j);
						if (fres[fres[l].last_j].pair > -1 && B_lj >= 0 && B_lj < n && Bp_lj >= 0 && Bp_lj<n){
							if (b_ij >= 0 && b_ij < n && l < b_ij){
								if (i <= fres[fres[l].last_j].pair && fres[fres[l].last_j].pair < j && l+3 <=j){
									// int sum = get_BE(fres[B_lj].pair,B_lj,fres[Bp_lj].pair,Bp_lj)+ get_WMBP(i,l-1)+ get_VP(l,j);
									// if (temp > sum){
									// 	temp = sum;
									// 	l_min = l;
									// }
								}
							}
						}
						m3 = 2*params->PB_penalty + temp;
					}
				}

				// 4) WMB(i,j) = VP(i,j) + P_b
				int temp = VP(i_mod,j) + params->PB_penalty;
				if (temp < m4){
					m4 = temp;
				}
				if(fres[j].pair < j){
					int l,l_min =-1;
					for(l = i+1; l<j; l++){
						if (fres[l].pair < 0 && fres[fres[l].last_j].pair > -1 && fres[fres[j].last_j].pair > -1 && fres[fres[j].last_j].pair == fres[fres[l].last_j].pair){
							// int temp = get_WMBP(i,l) + get_WI(l+1,j);
							// if (temp < m5){
							// 	m5 = temp;
							// 	l_min = l;
							// }
						}
					}
				}

			// get the min for WMB
			WMBP[j] = std::min({m1,m3,m4,m5});
			}

			// End of WMBP

			// Start of WMB

			if ((fres[i].pair >= -1 && fres[i].pair > j) || (fres[j].pair >= -1 && fres[j].pair < i) || (fres[i].pair >= -1 && fres[i].pair < i ) || (fres[j].pair >= -1 && j < fres[j].pair)) WMB[j] = INF;
			else{
			int m2 = INF, mWMBP = INF;
			// 2)
			if (fres[j].pair >= 0 && j > fres[j].pair){
				int l, l_min=-1;
				int bp_j = fres[j].pair;
				int temp = INF;
				for (l = (bp_j +1); (l < j); l++){
					int Bp_lj = getBp(fres,l,j);
					if (Bp_lj >= 0 && Bp_lj<n){
						// int sum = get_BE(bp_j,j,fres[Bp_lj].pair,Bp_lj) + get_WMBP(i,l) + get_WI(l+1,Bp_lj-1);
						// if (temp > sum){
						// 	temp = sum;
						// 	l_min = l;
						// }

					}
				}
				m2 = params->PB_penalty + temp;
			}
			// check the WMBP_ij value
			mWMBP =  WMBP[j];

			// get the min for WMB
			WMB[j] = std::min(m2,mWMBP);
			}

			// End of WMB

			// Start of WI -- the conditions on calculating WI is the same as WIP, so we combine them
			int wi_split = INF;
			int wip_split = INF;
		
			if (weakly_closed_ij == 0 || fres[fres[i].last_j].pair != fres[fres[j].last_j].pair){
					WI[j] = INF;
					WIP[j] = INF;
			}
			else{
				int wi_v = INF;
				int wip_v = INF;
				int wi_wmb = INF;
				int wip_wmb = INF;
				
				for ( auto const [key,val] : CL[j] ) {
					size_t k=key;
					// Start with WI
					energy_t v_kj = val + params->PPS_penalty;
					wi_split = std::min(wi_split,WI[k] + v_kj);
					// Then do WIP
					v_kj = val + params->bp_penalty;
					wip_split = std::min(wip_split,WIP[k]+v_kj);
					wip_split = std::min(wip_split,static_cast<energy_t>((k-i)*params->cp_penalty) +v_kj);
				}

				for ( auto const [key,val] : CLWMB[j] ) {
					size_t k=key;
					// Start with WI
					energy_t wmb_kj = val + params->PSP_penalty + params->PPS_penalty;
					wi_split = std::min(wi_split,WI[k] + wmb_kj);
					// Then do WIP
					wmb_kj = val + params->PSM_penalty + params->bp_penalty;
					wip_split = std::min(wip_split,WIP[k]+wmb_kj);
					wip_split = std::min(wip_split,static_cast<energy_t>((k-i)*params->cp_penalty) +wmb_kj);
				}
				wip_split = std::min(wip_split,WIP[j-1]) + params->cp_penalty;
				// for (int t = i; t< j; t++){
					// int wi_1 = get_WI(i,t);
					// int wi_2 = get_WI(t+1,j);
					// int energy = wi_1 + wi_2;
					// m1 = (m1 > energy)? energy : m1;
				// }

				// branch 2:

				// if ((fres[i].pair == j && fres[j].pair == i) ||(fres[i].pair < -1 && fres[j].pair < -1)){
				// 	int v_ener = (i>j)? INF: V(i_mod,j);
				// 	m2 = v_ener + params->PPS_penalty;
				// }
				if(ptype_closing>0 && !restricted && evaluate) {
					wi_v = V(i_mod,j) + params->PPS_penalty;
					wip_v = V(i_mod,j)	+ params->bp_penalty;
				}
				wi_wmb = WMB[j] + params->PSP_penalty + params->PPS_penalty;
				wip_wmb = WMB[j] + params->PSM_penalty + params->bp_penalty;

				WI[j] = std::min({wi_split,wi_v,wi_wmb});
				WIP[j] = std::min({wip_split,wip_v,wip_wmb});
			}
			// End of WI
			// int WIP[n+1];
			// // Start of WIP
			// if (fres[fres[i].last_j].pair != fres[fres[j].last_j].pair || weakly_closed_ij == 0){
			// 	WIP[j] = INF;
			// }else{
			// 	int m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;
			// 	// branch 1:
			// 	if (fres[i].pair < -1){
			// 		// m1 = get_WIP(i+1,j) + params->cp_penalty;
			// 	}
			// 	// branch 2:
			// 	if (fres[j].pair < -1){
			// 		// m2 = WIP[j-1] + params->cp_penalty;
			// 	}
			// 	//branch 3:
			// 	int t;
			// 	for (t = i; t <j; t++){
			// 		// int tmp = get_WIP(i,t) + get_WIP(t+1,j);
			// 		// if (tmp < m3){
			// 		// 	m3 = tmp;
			// 		// }
			// 	}

			// 	// branch 4:
			// 	if (fres[i].pair == j || (fres[i].pair < -1 && fres[j].pair < -1 && ptype_closing>0)){
			// 		// m4 = V(i_mod,j)	+ params->bp_penalty;
			// 	}

			// 	// branch 5:
			// 	// m5 = WMB[j] + params->PSM_penalty + params->bp_penalty;

			// 	WIP[j] = std::min({m1,m2,m3,m4,m5});
			// }
			// End of WIP



			// // start of VPP
			// int VPP[n+1];
			// if(is_weakly_closed(fres,B,b,i,j)) VPP[j] = INF;
			// else{
			// 	int m1 = INF, m2 = INF, m3 = INF, m4 = INF;
			// 	int r = -1;

			// 	int max_i_bp = i;
			// 	if (bp_ij > 0 && bp_ij < n && bp_ij > max_i_bp) max_i_bp = bp_ij;
			// 	for (r = max_i_bp+1; r < j; r++ ){
			// 		if (fres[r].pair < -1){
			// 			// int tmp = get_VP(i,r) + get_WIP(r+1,j);
			// 			// if (tmp < m1){
			// 				// m1 = tmp;
			// 			// }
			// 		}
			// 	}

			// 	int min_Bp_j = j;
			// 	if (Bp_ij > 0 && Bp_ij < n && bp_ij < min_Bp_j) min_Bp_j = Bp_ij;
			// 	for (r = i+1; r < min_Bp_j; r++){
			// 		if (fres[r].pair < -1){
			// 			// int tmp = get_WIP(i,r-1) + get_VP(r,j);
			// 			// if (tmp < m2){
			// 			// 	m2 = tmp;
			// 			// }
			// 		}
			// 	}

			// 	for (r = max_i_bp+1; r < j; r++ ){
			// 		int empty_region_rj = is_empty_region(fres,B,b,r+1,j); // r+1 to j
			// 		if (fres[r].pair < -1 && empty_region_rj){
			// 			// int tmp = get_VP(i,r) + (cp_penalty *(j-r)); // check the (j-r) part
			// 			// if (tmp < m3){
			// 			// 	m3 = tmp;
			// 			// }
			// 		}
			// 	}

			// 	for (r = i+1; r < min_Bp_j; r++){
			// 		int empty_region_ir = is_empty_region(fres,B,b,i,r-1); // i to r-1
			// 		if (fres[r].pair < -1 && empty_region_ir){
			// 			// int tmp = (params->cp_penalty * (r-i)) + get_VP(r,j);
			// 			// if (tmp < m4){
			// 			// 	m4 = tmp;
			// 			// }
			// 		}
			// 	}
			// 	VPP[j] = std::min({m1,m2,m3,m4});
			// }
			// End of VPP


			// Start of BE
			int BE[n+1];
			int ip = fres[i].pair; // might be the case that j and jp should be i and ip and vice versa.
			int jp = fres[j].pair; // currently, i is paired with ip and j with jp
			// if (!( i >= 0 && i <= ip && ip < jp && jp <= j && j < n && fres[i].pair >= -1 && fres[j].pair >= -1 && fres[ip].pair >= -1 && fres[jp].pair >= -1 && fres[i].pair == j && fres[j].pair == i && fres[ip].pair == jp && fres[jp].pair == ip)){ //impossible cases
			
			// base case: i.j and ip.jp must be in G
			if (fres[i].pair != j || fres[ip].pair != jp) BE[ip] = INF;
			else{

				int m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;
				if (fres[i+1].pair == ip-1){
					// m1 = params->e_stP_penalty*ILoopE(S,S1,params,ptype_closing,i,j,i+1,j-1) + get_BE(i+1,j-1,ip,jp);
				}

				for (int l = i+1; l<= ip ; l++){
					if (fres[l].pair >= -1 && j <= fres[l].pair && fres[l].pair < ip){
						int lp = fres[l].pair;
						int empty_region_il = is_empty_region(fres,B,b,i+1,l-1);
						int empty_region_lj = is_empty_region(fres,B,b,lp+1,ip-1);
						int weakly_region_il = is_weakly_closed(fres,B,b,i+1,l-1);
						int weakly_closed_lj = is_weakly_closed(fres,B,b,lp+1,ip-1);
						
						if (empty_region_il == 1 && empty_region_lj == 1 ){
							// int temp = params->e_intP_penalty*ILoopE(S,S1,params,ptype_closing,i,ip,l,lp)+ get_BE(l,lp,jp,j);
							// if (m2 > temp){
							// 	m2 = temp;
							// }
						}

						// 3)
						if (weakly_region_il == 1 && weakly_closed_lj == 1){
							// int temp = get_WIP(i+1,l-1) + get_BE(l,lp,jp,j) + get_WIP(lp+1,ip-1)+ params->ap_penalty + 2*params->bp_penalty;
							// if (m3 > temp){
							// 	m3 = temp;
							// }
						}

						// 4)
						if (weakly_region_il == 1 && empty_region_lj == 1){
							// int temp = get_WIP(i+1,l-1) + get_BE(l,lp,jp,j) + params->cp_penalty * (ip-lp+1) + params->ap_penalty + 2*params->bp_penalty;
							// if (m4 > temp){
							// 	m4 = temp;
							// }
						}

						// 5)
						if (empty_region_il == 1 && weakly_closed_lj == 1){
							// int temp = params->ap_penalty + 2*params->bp_penalty + (params->cp_penalty * (l-i+1)) + get_BE(l,lp,jp,j) + get_WIP(lp+1,ip-1);
							// if (m5 > temp){
							// 	m5 = temp;
							// }
						}
					}
				}

				// finding the min and putting it in BE[iip]
				BE[ip] = std::min({m1,m2,m3,m4,m5});
			}
			// End of BE

			
		
		
		} // end loop j
		rotate_arrays(WM,WM2,dmli1,dmli2,WMB,dwmbi,WI,dwibi,n);
		// Clean up trace arrows in i+MAXLOOP+1
		if (garbage_collect && i+MAXLOOP+1 <= n) {
			gc_row(ta,i + MAXLOOP + 1 );
		}

		// Reallocate candidate lists in i
		for ( auto &x: CL ) {
			if (x.capacity() > 1.5*x.size()) {
				cand_list_t vec(x.size());
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
	energy_t mfe = fold(sparsemfefold.seq_,sparsemfefold.V_,sparsemfefold.cand_comp,sparsemfefold.CL_,sparsemfefold.CLWMB_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.params_,sparsemfefold.ta_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_, sparsemfefold.dmli1_, sparsemfefold.dmli2_,sparsemfefold.VP_,sparsemfefold.WMB_,sparsemfefold.dwmbi_,sparsemfefold.WMBP_,sparsemfefold.WI_,sparsemfefold.dwib1_,sparsemfefold.WIP_,sparsemfefold.n_,sparsemfefold.garbage_collect_, sparsemfefold.fres,sparsemfefold.B,sparsemfefold.b);	
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
	std::cout << "Can num:\t"<<num_of_candidates(sparsemfefold.CL_)<<std::endl;
	std::cout << "Can cap:\t"<<capacity_of_candidates(sparsemfefold.CL_)<<std::endl;
	std::cout << "TAs num:\t"<<sizeT(sparsemfefold.ta_)<<std::endl;
	std::cout << "TAs cap:\t"<<capacityT(sparsemfefold.ta_)<<std::endl;
	}
	

	return 0;
}
