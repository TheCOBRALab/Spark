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
#include <vector>
#include <iterator>
#include <cstring>
#include <string>
#include <cassert>

#include <matrix.hh>
#include "base_types.hh"
#include "trace_arrow.hh"
#include "sparse_tree.cc"

extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/params/io.h"
}

#include "cmdline.hh"

#include <fstream>

struct quatret
{
    cand_pos_t first; 
    energy_t second;
    energy_t third;
	energy_t fourth;
	quatret(){
		first = 1;
		second = 2;
		third = 3;
		fourth = 4;
	}
	quatret(cand_pos_t x, energy_t y , energy_t z,energy_t w){
		first = x;
		second = y;
		third = z;
		fourth = w;
	}
};


typedef std::pair<cand_pos_t,energy_t> cand_entry_t;
typedef std::vector< cand_entry_t > cand_list_t;

typedef quatret cand_entry_td1;
typedef std::vector< cand_entry_td1 > cand_list_td1;

class SparseMFEFold;

energy_t ILoopE(auto const& S_,auto const& S1_, auto const& params_,const int& ptype_closing, const size_t& i, const size_t& j, const size_t& k,  const size_t& l);
energy_t MbLoopE(auto const& S_, auto const& params_, int ptype_closing,size_t i, size_t j);
energy_t Mlstem(auto const& S_, auto const& params_, int ptype_closing,size_t i, size_t j);
void trace_V(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e,std::vector<Node> tree,std::vector<int> up);
void trace_W(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j,std::vector<Node> tree,std::vector<int> up);
void trace_WM(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e,std::vector<Node> tree,std::vector<int> up) ;
void trace_WM2(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,size_t i, size_t j,std::vector<Node> tree,std::vector<int> up);

/**
* Space efficient sparsification of Zuker-type RNA folding with
* trace-back. Provides methods for the evaluation of dynamic
* programming recursions and the trace-back.
*/
class SparseMFEFold {

public:
	std::string seq_;
	cand_pos_t n_;

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
	std::vector<energy_t> VPP_;
	std::vector<energy_t> WMB_;
	std::vector<energy_t> dwmbi_; // WMB from 1 iteration ago
	std::vector<energy_t> WMBP_;
	std::vector<energy_t> WI_;
	std::vector<energy_t> dwi1_; // WI from 1 iteration ago
	std::vector<energy_t> WIP_;
	std::vector<energy_t> dwip1_; // WIP from 1 iteration ago
	std::vector<energy_t> WV_;
	std::vector<energy_t> dwvp_; // WV from 1 iteration ago;

	bool mark_candidates_;


	TraceArrows ta_;

	// TraceArrows ta_dangle_;
	
	std::vector< cand_list_td1 > CL_;
	std::vector< cand_list_t > CLVP_;
	std::vector< cand_list_t > CLWMB_;
	std::vector< cand_list_t > CLBE_;

	
	

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
	: seq_(seq), n_(seq.length()), params_(scale_parameters()), ta_(n_), garbage_collect_(garbage_collect)
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
	VPP_.resize(n_+1,INF);
	WMB_.resize(n_+1,INF);
	dwmbi_.resize(n_+1,INF);
	WMBP_.resize(n_+1,INF);
	WI_.resize(n_+1,0);
	dwi1_.resize(n_+1,0);
	WIP_.resize(n_+1,INF);
	dwip1_.resize(n_+1,INF);
	WV_.resize(n_+1,INF);
	dwvp_.resize(n_+1,INF);

	// init candidate lists
	CL_.resize(n_+1);
	CLWMB_.resize(n_+1);
	CLVP_.resize(n_+1);
	CLBE_.resize(n_+1);

	resize(ta_,n_+1);

	// resize(ta_dangle_,n_+1);

	restricted_ = restricted;
	
	}

	

	~SparseMFEFold() {
	free(params_);
	free(S_);
	free(S1_);
	}
};


/**
 * @brief This code returns the hairpin energy for a given base pair.
 * @param i The left index in the base pair
 * @param j The right index in the base pair
*/
energy_t HairpinE(auto const& seq, auto const& S, auto const& S1, auto const& params, size_t i, size_t j) {
	
	const pair_type ptype_closing = pair[S[i]][S[j]];

	if (ptype_closing==0) return INF;

	return E_Hairpin(j-i-1,ptype_closing,S1[i+1],S1[j-1],&seq.c_str()[i-1], const_cast<paramT *>(params));
}

/**
 * @brief Encodes the dangle type into the energy by shifting the energy two bits to the left and putting the dangle.
 * Dangles has 4 possible values: 00, 01, 10, 11, for no dangle, 5' dangle, 3' dangle, and 53' dangle, respectively.
 * 
 * @param e The energy being shifted. In practice, W(i,j) or WM(i,j)
 * @param d The type of dangle.
*/
energy_t encode(energy_t e, Dangle d){
    return (e << 2) | d;
}

/**
 * @brief The complement to encode. This returns the energy and dangle, respectively, in a tuple format
 * 
 * @param enc The encoded energy 
*/
std::pair<energy_t,Dangle> decode(energy_t enc){
    return std::make_pair((enc >> 2), (enc & 3));
}


energy_t E_ext_Stem(auto const& vij,auto const& vi1j,auto const& vij1,auto const& vi1j1,auto const& S, auto const& params, const cand_pos_t i,const size_t j, Dangle &d, cand_pos_t n, std::vector<Node> tree){

	energy_t e = INF, en = INF;
  	pair_type tt  = pair[S[i]][S[j]];

    if ((tree[i].pair <-1 && tree[j].pair <-1) || (tree[i].pair == j && tree[j].pair == i)) {
				en = vij; // i j

                const base_type si1 = i>1 ? S[i-1] : -1;
                const base_type sj1 = j<n ? S[j+1] : -1;
				if (en != INF) {
					if (params->model_details.dangles == 2)
                        en += vrna_E_ext_stem(tt, si1, sj1, params);
                    else
                        en += vrna_E_ext_stem(tt, -1, -1, params);

                    e = MIN2(e, en);
					
				}

	}

	if(params->model_details.dangles  == 1){
        tt  = pair[S[i+1]][S[j]];
        if (((tree[i+1].pair <-1 && tree[j].pair <-1) || (tree[i+1].pair == j)) && tree[i].pair<0) {
            en = (j-i-1>TURN) ? vi1j : INF; //i+1 j

            if (en != INF) {

                const base_type si1 = S[i];
                en += vrna_E_ext_stem(tt, si1, -1, params);
            }

            e = MIN2(e,en);
            if(e == en){
                d=1;
            }

        }
        tt  = pair[S[i]][S[j-1]];
        if (((tree[i].pair <-1 && tree[j-1].pair <-1) || (tree[i].pair == j-1)) && tree[j].pair<0) {
            en = (j-1-i>TURN) ? vij1 : INF; // i j-1
            if (en != INF) {

                const base_type sj1 = S[j];

                en += vrna_E_ext_stem(tt, -1, sj1, params);
            }
            e = MIN2(e,en);
            if(e == en){
                d=2;
            }

        }
        tt  = pair[S[i+1]][S[j-1]];
        if (((tree[i+1].pair <-1 && tree[j-1].pair <-1) || (tree[i+1].pair == j-1)) && tree[i].pair < 0 && tree[j].pair<0) {
            en = (j-1-i-1>TURN) ? vi1j1 : INF; // i+1 j-1

            if (en != INF) {

                const base_type si1 = S[i];
                const base_type sj1 = S[j];

                en += vrna_E_ext_stem(tt, si1, sj1, params);
            }
            e = MIN2(e,en);
            if(e == en){
                d=3;
            }
        }
	}
	return e;
}

/**
* @brief Rotate WM2 and WI arrays to store the previous and previous previous iterations
* @param WM2 WM2 array
* @param dmli1 WM2 from one iteration ago
* @param dmli2 WM2 from two iterations ago
*/
void rotate_arrays(auto &WM2, auto &dmli1, auto &dmli2, auto &WI, auto &WIP, auto &dwi1, auto &dwip1,auto &WMB, auto &dwmbi, auto &WV, auto &dwvp){
	dmli2.swap(dmli1);
    dmli1.swap(WM2);
	dwi1.swap(WI);
	dwip1.swap(WIP);
	dwmbi.swap(WMB);
	dwvp.swap(WV);
}


/**
* @brief Computes the multiloop V contribution. This gives back essentially VM(i,j).
* 
* @param dmli1 Row of WM2 from one iteration ago
* @param dmli2 Row of WM2 from two iterations ago 
*/
energy_t E_MbLoop(auto const& dmli1, auto const& dmli2, auto const& S, auto const& params, cand_pos_t i, cand_pos_t j, std::vector<Node> tree){

	energy_t e = INF,en = INF;
  	pair_type tt  = pair[S[j]][S[i]];
	bool pairable = (tree[i].pair <-1 && tree[j].pair <-1) || (tree[i].pair == j);
	
	/* double dangles */
	switch(params->model_details.dangles){
		case 2:
			if (pairable) {
			e = dmli1[j - 1];

			if (e != INF) {

				const base_type si1 = S[i + 1];
				const base_type sj1 = S[j - 1];

				e += E_MLstem(tt, sj1, si1, params) + params->MLclosing;
			}

			}
			break;

		case 1:
			/**
			* ML pair D0
			*  new closing pair (i,j) with mb part [i+1,j-1]  
			*/
			
			if (pairable) {
        		e = dmli1[j - 1];

        		if (e != INF) {

          			e += E_MLstem(tt, -1, -1, params) + params->MLclosing;

        		}
      		}
      		/** 
			* ML pair 5
			* new closing pair (i,j) with mb part [i+2,j-1] 
			*/
      		if (pairable && tree[i+1].pair < 0) {
        		en = dmli2[j - 1];

        		if (en != INF) {

          			const base_type si1 =  S[i + 1];

          			en += E_MLstem(tt, -1, si1, params) + params->MLclosing + params->MLbase;
      
        		}
      		}
			e   = MIN2(e, en);
			

			/** 
			* ML pair 3
			* new closing pair (i,j) with mb part [i+1, j-2] 
			*/
			if (pairable && tree[j-1].pair < 0) {
				en = dmli1[j - 2];

				if (en != INF) {
					const base_type sj1 = S[j - 1];

					en += E_MLstem(tt, sj1, -1, params) + params->MLclosing + params->MLbase; 
				}
			}
			e   = MIN2(e, en);
			/** 
			* ML pair 53
			* new closing pair (i,j) with mb part [i+2.j-2]
			*/
			if (pairable && tree[i+1].pair < 0 && tree[j-1].pair <0) {
				en = dmli2[j - 2];

				if (en != INF) {

					const base_type si1 = S[i + 1];
					const base_type sj1 = S[j - 1];

					en += E_MLstem(tt, sj1, si1, params) + params->MLclosing + 2 * params->MLbase;
				}
			}
			e   = MIN2(e, en);
      		break;
		case 0:
			if (pairable) {
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
 * @brief Gives the WM(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
 * @param vij The V(i,j) energy
 * @param vi1j The V(i+1,j) energy
 * @param vij1 The V(i,j-1) energy
 * @param vi1j1 The V(i+1,j-1) energy
*/
energy_t E_MLStem(auto const& vij,auto const& vi1j,auto const& vij1,auto const& vi1j1, auto const& WM, auto const& CL,auto const& S, auto const& params,const cand_pos_t i,const cand_pos_t j,Dangle &d, const cand_pos_t& n, std::vector<Node> tree){

	energy_t e = INF,en=INF;

	pair_type type = pair[S[i]][S[j]];


	if ((tree[i].pair < -1 && tree[j].pair < -1) || (tree[i].pair == j)) {
		en = vij; // i j
		if (en != INF) {
            const base_type si1 = i>1 ? S[i-1] : -1;
            const base_type sj1 = j<n ? S[j+1] : -1;
			if (params->model_details.dangles == 2)
				en += E_MLstem(type, si1, sj1, params);
			else
				en += E_MLstem(type, -1, -1, params);

			e = MIN2(e, en);
		}
	}
	if(params->model_details.dangles == 1 && tree[i].pair != j and tree[j].pair != i){
		const base_type mm5 = S[i], mm3 = S[j];

		if (((tree[i+1].pair < -1 && tree[j].pair < -1) || (tree[i+1].pair == j)) && tree[i].pair < 0) {
      		en = (j-i-1 >TURN) ? vi1j : INF; // i+1 j
      		if (en != INF) {
        		en += params->MLbase;

            	type = pair[S[i+1]][S[j]];
            	en += E_MLstem(type, mm5, -1, params);

        		e = MIN2(e, en);
				if(e == en){
					d=1;
				}
      		}
    	}

		if (((tree[i].pair < -1 && tree[j-1].pair < -1) || (tree[i].pair == j-1)) && tree[j].pair < 0) {
      		en = (j-1-i>TURN) ? vij1 : INF; // i j-1
      		if (en != INF) {
       			en += params->MLbase;

            	type = pair[S[i]][S[j-1]];
            	en += E_MLstem(type, -1, mm3, params);
 
        		e = MIN2(e, en);
				if(e == en){
					d=2;
				}
      		}
    	}
    	if (((tree[i+1].pair < -1 && tree[j-1].pair < -1) || (tree[i+1].pair == j-1)) && tree[i].pair < 0 && tree[j].pair<0) {
      		en = (j-1-i-1>TURN) ? vi1j1 : INF; // i+1 j-1
      		if (en != INF) {
        		en += 2 * params->MLbase;

        		type = pair[S[i+1]][S[j-1]];
        		en += E_MLstem(type, mm5, mm3, params);
        
				e = MIN2(e, en);
				if(e == en){
					d=3;
				}
      		}
    	} 
		
	}


    return e;
}


auto const recompute_WI(auto &WI, auto const &CL, auto const &CLWMB, auto const& S, auto const &params, auto const& n, cand_pos_t i, cand_pos_t max_j, std::vector<Node> tree) {
	

	assert(i>=1);
	assert(max_j<=n);

	
	for ( cand_pos_t j=i; j<=max_j; j++ ) {
		energy_t wi = (j-i+1)*params->PUP_penalty;
		bool paired;
		#pragma omp parallel for num_threads(6);
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			cand_pos_t k = it->first;
			if(k<i) continue;
			paired = (tree[k].pair == j && tree[j].pair == k);
			const energy_t v_kj = it->second + params->bp_penalty;
			wi = std::min( wi, WI[k-1]  + v_kj );
			if(paired) break;
		}

		for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i ; ++it ) {
			if(paired) continue;
			cand_pos_t k = it->first;
			if(k<i) continue;
			paired = (tree[k].pair == j && tree[j].pair == k);
			const energy_t wmb_kj = it->second + params->bp_penalty + params->PPS_penalty;
			wi = std::min( wi, WI[k-1]  + wmb_kj );	
		}
		// if(wi == 0) WI[j] = (j-i)*params->PUP_penalty;
		WI[j] = wi;
		
	}
	return WI;
}

auto const recompute_WIP(auto &WIP, auto const &CL, auto const &CLWMB, auto const& S, auto const &params, auto const& n, cand_pos_t i, cand_pos_t max_j, std::vector<Node> tree, std::vector<int> up) {
	

	assert(i>=1);
	assert(max_j<=n);

	
	
	for ( cand_pos_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wip = INF;
		bool paired;
		#pragma omp parallel for num_threads(6);
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			cand_pos_t k = it->first;
			if(k<i) continue;
			paired = (tree[k].pair == j && tree[j].pair == k);
			const energy_t v_kj = it->second + params->bp_penalty;
			bool can_pair = up[k-1] > k-i-1;
			if(can_pair) wip = std::min( wip, static_cast<energy_t>(params->MLbase*(k-i)) + v_kj );
			wip = std::min( wip, WIP[k-1]  + v_kj );
			if(paired) break;
		}

		for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i ; ++it ) {
			if(paired) continue;
			cand_pos_t k = it->first;
			if(k<i) continue;
			paired = (tree[k].pair == j && tree[j].pair == k);
			const energy_t wmb_kj = it->second + params->bp_penalty + params->PPS_penalty;
			wip = std::min( wip, static_cast<energy_t>(params->MLbase*(k-i)) + wmb_kj );
			wip = std::min( wip, WIP[k-1]  + wmb_kj );	
		}
		if(tree[j].pair<0) wip = std::min(wip, WIP[j-1] + params->MLbase);
		WIP[j] = wip;
	}
	return WIP;
}


/**
* @brief Recompute row of WM. This is used in the traceback when we haved decided the current i.j pair closes a multiloop,
* and the WM energies need to be recomputed fom the candidates.
* 
* @param WM WM array
* @param CL Candidate List
* @param i Current i
* @param max_j Current j
*/
auto const recompute_WM(auto const& WM, auto const &CL, auto const& S, auto const &params, auto const& n, cand_pos_t i, cand_pos_t max_j, std::vector<Node> tree, std::vector<int> up) {
	

	assert(i>=1);
	assert(max_j<=n);

	std::vector<energy_t> temp = WM;

	for ( cand_pos_t j=i-1; j<=std::min(i+TURN,max_j); j++ ) { temp[j]=INF; }
	
	for ( cand_pos_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wm = INF;
		bool paired = 0;
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			cand_pos_t k = it->first;
			paired += (tree[k].pair == j);
			auto const [v_kj,d] = decode(it->third);
			bool can_pair = up[k-1] >= (k-i) ? true : false;
			if(can_pair) wm = std::min( wm, static_cast<energy_t>(params->MLbase*(k-i)) + v_kj );
			wm = std::min( wm, temp[k-1]  + v_kj );
		}
		if(tree[j].pair<0 && !paired) wm = std::min(wm, temp[j-1] + params->MLbase);
		temp[j] = wm;
	}
	return temp;
}

/**
* @brief Recompute row of WM2. This is used in the traceback when we haved decided the current i.j pair closes a multiloop,
* and the WM2 energies need to be recomputed fom the candidates to get the corresponding energy for it.
* 
* @param WM WM array
* @param WM2 WM2 array
* @param CL Candidate List
* @param i Current i
* @param max_j Current j
*/
auto const recompute_WM2(auto const& WM, auto const& WM2, auto const CL, auto const& S, auto const &params, auto const& n, cand_pos_t i, cand_pos_t max_j, std::vector<Node> tree) {
	

	assert(i>=1);
	assert(max_j<= n);

	std::vector<energy_t> temp = WM2;

	for ( cand_pos_t j=i-1; j<=std::min(i+2*TURN+2,max_j); j++ ) { temp[j]=INF; }

	for ( cand_pos_t j=i+2*TURN+3; j<=max_j; j++ ) {
		energy_t wm2 = INF;
		bool paired = 0;
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>i+TURN+1 ; ++it ) {
			
			cand_pos_t k = it->first;
			paired += (tree[k].pair == j && tree[j].pair == k);
			auto const [v_kj,d] = decode(it->third);
			wm2 = std::min( wm2, WM[k-1]  + v_kj );
		}
		if(tree[j].pair<0 && !paired) wm2 = std::min(wm2, temp[j-1] + params->MLbase);
		temp[j] = wm2;
	}
	return temp;
}

/**
 * @brief Test existence of candidate. Used primarily for determining whether (i,j) is candidate for W/WM splits
 * 
 * @param CL Candidate List
 * @param cand_comp Candidate Comparator
 * @param i start
 * @param j end
 * @return  
 */
bool is_candidate(auto const& CL,auto const& cand_comp,cand_pos_t i, cand_pos_t j) {
	const cand_list_td1 &list = CL[j];

	auto it = std::lower_bound(list.begin(),list.end(),i,cand_comp);

	return it!=list.end() && it->first==i;
}
/**
 * @brief Determines the type of dangle being used for a closing multiloop while in traceback.
 * 
 * @param WM2ij The WM2 energy for the region [i,j]
 * @param WM2i1j The WM2 energy for the region [i+1,j]
 * @param WM2ij1 The WM2 energy for the region [i,j-1]
 * @param WM2i1j1 The WM2 energy for the region [i+1,j-1]
*/
void find_mb_dangle(const energy_t &WM2ij,const energy_t &WM2i1j,const energy_t &WM2ij1,const energy_t &WM2i1j1,auto const &params, auto const& S, const cand_pos_t &i, const cand_pos_t &j, cand_pos_t &k, cand_pos_t &l,std::vector<Node> tree){
	if(params->model_details.dangles == 2) return;

	pair_type tt = pair[S[j]][S[i]];
	energy_t e1 = WM2ij +  E_MLstem(tt, -1, -1, params);
	energy_t e2 = WM2i1j +  E_MLstem(tt, -1, S[i+1], params);
	energy_t e3 = WM2ij1 +  E_MLstem(tt, S[j-1], -1, params);
	energy_t e4 = WM2i1j1 +  E_MLstem(tt, S[j-1], S[i+1], params);
	energy_t e = e1;
	if(e2<e && tree[i+1].pair< 0){
		e = e2;
		k = i+2;
		l = j-1;
	}
	if(e3<e && tree[j-1].pair< 0){
		e = e3;
		k = i+1; 
		l = j-2;
	}
	if(e4<e && tree[i+1].pair< 0 && tree[j-1].pair< 0){
		e = e4;
		k = i+2;
		l = j-2;
	}
 }

/**
 * @brief Traceback from W entry.
 * pre: W contains values of row i in interval i..j
 * 
 * @param seq Sequence
 * @param structure Final structure
 * @param W W array
 * @param i row index
 * @param j column index
 */
void trace_W(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, cand_pos_t i, cand_pos_t j,std::vector<Node> tree, std::vector<int> up) {
	// printf("W at %lu and %lu with %d\n",i,j,W[j]);
	if (i+TURN+1>=j) return;
	// case j unpaired
	if (W[j] == W[j-1]) {
		trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,i,j-1,tree,up);
		return;
	}
	cand_pos_t m=j+1;
	energy_t v=INF;
	energy_t w;
    Dangle dangle =3;
	energy_t vk = INF;
	for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i;++it ) {
		m = it->first;
		auto const[v_kj,d] = decode(it->fourth);
		w = W[m-1] + v_kj;
		if (W[j] == w) {
		v =it->second;
        dangle = d;
		vk = v_kj;
		break;
		}
	}
	cand_pos_t k=m;
	cand_pos_t l=j;
	pair_type ptype = 0;
    switch(dangle){
		case 0:
			ptype = pair[S[k]][S[l]];
			v = vk - E_ExtLoop(ptype,-1,-1,params);
			break;
        case 1:
            k=m+1;
			ptype = pair[S[k]][S[l]];
			v = vk - E_ExtLoop(ptype,S[m],-1,params);
            break;
        case 2:
            l=j-1;
			ptype = pair[S[k]][S[l]];
			v = vk - E_ExtLoop(ptype,-1,S[j],params);
            break;
        case 3:
            if(params->model_details.dangles == 1){
                k=m+1;
                l=j-1;
				ptype = pair[S[k]][S[l]];
				v = vk - E_ExtLoop(ptype,S[m],S[j],params);
            }
            break;

        
    }
	assert(i<=m && m<j);
	assert(v<INF);
	// don't recompute W, since i is not changed
	trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,i,m-1,tree,up);
	trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,v,tree,up);
}

/**
* @brief Traceback from V entry
* 
* @param structure Final Structure
* @param mark_candidates Whether Candidates should be [ ]
* @param i row index
* @param j column index
*/
void trace_V(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e,std::vector<Node> tree,std::vector<int> up) {
	// printf("V at %lu and %lu with %d\n",i,j,e);
	

	assert( i+TURN+1<=j );
	assert( j<=n );
	
	
	if (mark_candidates && is_candidate(CL,cand_comp,i,j)) {
		structure[i]='{';
		structure[j]='}';
	} else {
		structure[i]='(';
		structure[j]=')';
	}
	const pair_type ptype_closing = pair[S[i]][S[j]];
	if (exists_trace_arrow_from(ta,i,j)) {
		
		const TraceArrow &arrow = trace_arrow_from(ta,i,j);
		const size_t k=arrow.k(i,j);
		const size_t l=arrow.l(i,j);
		assert(i<k);
		assert(l<j);
		
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l, arrow.target_energy(),tree,up);
		return;

	}
	else {

		// try to trace back to a candidate: (still) interior loop case
		int l_min = std::max((int)i,(int) j-31);
		for ( cand_pos_t l=j-1; l>l_min; l--) {
			// Break if it's an assured dangle case
			for ( auto it=CL[l].begin(); CL[l].end()!=it && it->first>i; ++it ) {
				const cand_pos_t k=it->first;
				if(k-i > 31) continue;
				// if (  e == it->second + ILoopE(S,S1,params,ptype_closing,i,j,k,l) ) {
				if (  e == it->second + E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S[k]][S[l]]],S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params)) ) {
					trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,it->second,tree,up);
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
	auto const temp = recompute_WM(WM,CL,S,params,n,i+1,j-1,tree,up);
	WM = temp;
	auto const temp2 = recompute_WM2(WM,WM2,CL,S,params,n,i+1,j-1,tree);
	WM2 = temp2;
	
	// Dangle for Multiloop
	cand_pos_t k = i+1;
	cand_pos_t l = j-1;
	if(params->model_details.dangles == 1){
		auto const temp3 = recompute_WM(WM,CL,S,params,n,i+2,j-1,tree,up);
		auto const temp4 = recompute_WM2(temp3,WM2,CL,S,params,n,i+2,j-1,tree);
		find_mb_dangle(temp2[j-1],temp4[j-1],temp2[j-2],temp4[j-2],params,S,i,j,k,l,tree);
		if (k>i+1){
			WM = temp3;
			WM2 = temp4;
		}
	}
	
	
	trace_WM2(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,tree,up);
}

/**
* @brief Traceback from WM
*
* @param WM WM array at [i,j]
* @param WM2 WM2 array at [i,j]
* @param i row index
* @param j column index 
* @param e energy in WM(i,j) 
*/
void trace_WM(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, std::vector<Node> tree, std::vector<int> up) {
	// printf("WM at %lu and %lu with %d\n",i,j,e);

	if (i+TURN+1>j) {return;}

	if ( e == WM[j-1] + params->MLbase ) {
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,j-1,WM[j-1],tree,up);
		return;
	}
	energy_t v = INF;
    energy_t vk = INF;
    Dangle dangle = 3;
    cand_pos_t m = j+1;
	for ( auto it=CL[j].begin();CL[j].end() != it && it->first>=i;++it ) {
		m = it->first;
		auto const [v_kj,d] = decode(it->third);
		if ( e == WM[m-1] + v_kj ) {
            dangle = d;
			vk = v_kj;
			v = it->second;
			// no recomp, same i
		    break;
		} else if ( e == static_cast<energy_t>((m-i)*params->MLbase) + v_kj ) {
            dangle = d;
			vk = v_kj;
			v = it->second;
		    break;
		}
	}
	cand_pos_t k = m;
	cand_pos_t l = j;
	pair_type ptype = 0;
    switch(dangle){
		case 0:
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,-1,-1,params);
			break;
        case 1:
            k=m+1;
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,S[m],-1,params);
            break;
        case 2:
            l=j-1;
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,-1,S[j],params);
            break;
        case 3:
			if(params->model_details.dangles == 1){
				k=m+1;
				l=j-1;
				ptype= pair[S[k]][S[l]];
				v = vk - E_MLstem(ptype,S[m],S[j],params);
			}
            break;
    }
    

    if ( e == WM[m-1] + vk ) {
		// no recomp, same i
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,m-1,WM[m-1],tree,up);
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,v,tree,up);
		return;
	} else if ( e == static_cast<energy_t>((k-i)*params->MLbase) + vk ) {
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,v,tree,up);
		return;
	}
	assert(false);
}

/**
* @brief Traceback from WM2
* 
* @param WM WM array at [i,j]
* @param WM2 Wm2 array at [i,j]
* @param i row index
* @param j column index
 */
void trace_WM2(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,cand_pos_t i, cand_pos_t j,std::vector<Node> tree, std::vector<int> up) {
	// printf("WM2 at %lu and %lu with %d\n",i,j,WM2[j]);

	if (i+2*TURN+3>j) {return;}

	const energy_t e = WM2[j];

	// case j unpaired
	if ( e == WM2[j-1] + params->MLbase ) {
		
		// same i, no recomputation
		trace_WM2(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,j-1,tree,up);
		return;
	}

    cand_pos_t m = j+1;
    energy_t v = INF;
    energy_t vk = INF;
    Dangle dangle = 4;
	for ( auto it=CL[j].begin();CL[j].end() != it  && it->first>=i+TURN+1;++it ) {
		m = it->first;
		
		auto const [v_kj,d] = decode(it->third);
		if (e == WM[m-1] + v_kj) {
			vk = v_kj;
            dangle = d;
			v = it->second;
			break;
		}
	}
	cand_pos_t k = m;
	cand_pos_t l = j;
	pair_type ptype = 0;
    switch(dangle){
		case 0:
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,-1,-1,params);
			break;
        case 1:
            k=m+1;
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,S[m],-1,params);
            break;
        case 2:
            l=j-1;
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,-1,S[j],params);
            break;
        case 3:
			if(params->model_details.dangles == 1){
				k=m+1;
				l=j-1;
				ptype= pair[S[k]][S[l]];
				v = vk - E_MLstem(ptype,S[m],S[j],params);
			}
            break;
    }
    

    if ( e == WM[m-1] + vk ) {
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,m-1,WM[m-1],tree,up);
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,v,tree,up);
		return;
	}
	assert(false);
}
/**
* @brief Trace back
* pre: row 1 of matrix W is computed
* @return mfe structure (reference)
*/
const std::string & trace_back(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n,std::vector<Node> tree, std::vector<int> up,auto const& mark_candidates=false) {

	structure.resize(n+1,'.');

	/* Traceback */
	trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,1,n,tree,up);
	structure = structure.substr(1,n);

	return structure;
}

/**
 * @brief Returns the internal loop energy for a given i.j and k.l
 * 
*/
energy_t ILoopE(auto const& S, auto const& S1, auto const& params, const int32_t& ptype_closing,const size_t &i,const size_t &j,const size_t &k,const size_t &l)  {
	if(ptype_closing<=0) std::cout << i << " " << j << " " << k << " " << l << std::endl;
	assert(ptype_closing>0);
	assert(1<=i);
	assert(i<k);
	assert(k<l);
	assert(l<j);
	//assert(l<=len); // don't know len here

	// note: enclosed bp type 'turned around' for lib call
	const pair_type ptype_enclosed = rtype[pair[S[k]][S[l]]];

	if (ptype_enclosed==0) return INF;

	return E_IntLoop(k-i-1,j-l-1,ptype_closing,ptype_enclosed,S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params));
}


/**
* @brief Register a candidate
* @param i start
* @param j end
* @param e energy of candidate "V(i,j)"
* @param wmij energy at WM(i,j)
* @param wij energy at W(i,j)
*/
void register_candidate(auto &CL, cand_pos_t const& i, cand_pos_t const& j, energy_t const& e,energy_t const& wmij,energy_t const& wij) {
	assert(i<=j+TURN+1);
	
	CL[j].push_back( cand_entry_td1(i,e,wmij,wij) );
}

void register_candidate(auto &CL, cand_pos_t const& i, cand_pos_t const& j, energy_t const& e) {
	
	CL[j].push_back( cand_entry_t(i,e) );
}
/**
 * @brief Determines the MFE energy for a given sequence
*/
energy_t fold(auto const& seq, sparse_tree sparse_tree, auto &V, auto const& cand_comp, auto &CL, auto &CLWMB,auto &CLVP, auto &CLBE, auto const& S, auto const& S1, auto const& params, auto &ta, auto &W, auto &WM, auto &WM2, auto &dmli1, auto &dmli2, auto &VP, auto &VPP, auto &WV, auto &dwvp, auto &WMB, auto &dwmbi,auto &WMBP,auto &WI,auto &dwi1,auto &WIP, auto &dwip1, auto const& n, auto const& garbage_collect) {
	Dangle d = 3;
    if(params->model_details.dangles == 0) d = 0;

	
	std::vector<energy_t> WI_Bbp;
	WI_Bbp.resize(n+1,INF);
	std::vector<energy_t> WIP_Bbp;
	WIP_Bbp.resize(n+1,INF);

	std::vector<energy_t> WI_Bp;
	WI_Bp.resize(n+1,INF);

    
    for (cand_pos_t i=n; i>0; --i) {
		energy_t VP_i_split = INF;
		for(cand_pos_t j=i;sparse_tree.tree[j].pair<0 && j<=n ;++j){
			WI[j] = (j-i+1)*params->PUP_penalty;
		}
		bool bp_bool = false;
		bool B_bool = false;
		bool Bp_bool = false;
		if(sparse_tree.tree[i-1].pair>0 && i+4 < n){
			if(sparse_tree.tree[i+4].parent->index == i-1) bp_bool = true;
			if(sparse_tree.tree[i+4].parent->index == sparse_tree.tree[i-1].parent->index) B_bool = true;
		}

		for ( cand_pos_t j=i+TURN+1; j<=n; j++ ) {

			if(sparse_tree.tree[j].pair>0){
				bp_bool = false;
				B_bool = false;
			}
			
			if(sparse_tree.tree[j].index == sparse_tree.tree[j].parent->pair-1 && sparse_tree.tree[i].parent->index == sparse_tree.tree[j].parent->index) Bp_bool = true;
			else Bp_bool = false;
			

			bool evaluate = sparse_tree.weakly_closed(i,j);
			// ------------------------------
			// W: split case
			bool pairedkj = 0;
			energy_t w_split = INF;
			energy_t wm_split = INF;
			energy_t wm2_split = INF;
			energy_t wi_split = INF;
			energy_t wip_split = INF;
			for ( auto const [key,val,val_ml,val_w] : CL[j] ) {
				size_t k=key;
				
				auto const [v_kj,d] = decode(val_ml);
				auto const [v_kjw,dw] = decode(val_w);
				bool can_pair = sparse_tree.up[k-1] >= (k-i) ? true: false;
				// WM Portion
				wm_split = std::min( wm_split, WM[k-1] + v_kj );
				if(can_pair) wm_split = std::min( wm_split,static_cast<energy_t>((k-i)*params->MLbase) + v_kj );
				// WM2 Portion
				wm2_split = std::min( wm2_split, WM[k-1] + v_kj );
				// W Portion
				w_split = std::min( w_split, W[k-1] + v_kjw );

				//WI portion
				energy_t v_kjj = val + params->PPS_penalty;
				wi_split = WI[k] + v_kj;
				

				//WIP portion
				v_kjj = val + params->bp_penalty;
				wip_split = WIP[k]+v_kj;
				if(can_pair) wip_split = std::min(wip_split,static_cast<energy_t>((k-i)*params->cp_penalty) +v_kj);
				
		
				
			}

			for ( auto const [key,val] : CLWMB[j] ) {
				
				if(pairedkj) break;   // Not needed I believe as there shouldn't be any candidates there if paired anyways
				// Maybe this would just avoid this loop however

				cand_pos_t k = key;
				bool can_pair = sparse_tree.up[k-1] >= (k-i) ? true: false;
				
				// For W
				energy_t wmb_kj = val + params->PS_penalty;
				w_split = std::min( w_split, W[k-1] + wmb_kj );	
				// For WM -> I believe this would add a PSM penalty for every pseudoknot which would be bad
				wmb_kj = val + params->PSM_penalty + params->b_penalty;
				wm_split = std::min(wm_split, WM[k-1] + wmb_kj);
				if(can_pair) wm_split = std::min(wm_split,static_cast<energy_t>((k-i)*params->MLbase) + wmb_kj);
				wm2_split = std::min( wm2_split, WM[k-1] + wmb_kj );
				if(can_pair) wm2_split = std::min( wm2_split, static_cast<energy_t>((k-i)*params->MLbase) + wmb_kj );
				// For WI
				wmb_kj = val + params->PSP_penalty + params->PPS_penalty;
				wi_split = std::min(wi_split,WI[k] + wmb_kj);
				// For WIP
				wmb_kj = val + params->PSM_penalty + params->bp_penalty;
				wip_split = std::min(wip_split,WIP[k]+wmb_kj);
				if(can_pair) wip_split = std::min(wip_split,static_cast<energy_t>((k-i)*params->cp_penalty) +wmb_kj);
			}
			if(sparse_tree.tree[j].pair<0) w_split = std::min(w_split,W[j-1]);
			if(sparse_tree.tree[j].pair<0) wm2_split = std::min( wm2_split, WM2[j-1] + params->MLbase );
			if(sparse_tree.tree[j].pair<0) wm_split = std::min( wm_split, WM[j-1] + params->MLbase );
			if(sparse_tree.tree[j].pair<0) wi_split = std::min(wi_split,WI[j-1] + params->PUP_penalty);
			if(sparse_tree.tree[j].pair<0) wip_split = std::min(wip_split,WIP[j-1] + params->cp_penalty);
			

			
			energy_t w  = w_split; // entry of W w/o contribution of V
			energy_t wm = wm_split; // entry of WM w/o contribution of V


			size_t i_mod=i%(MAXLOOP+1);

			const pair_type ptype_closing = pair[S[i]][S[j]];
			const bool restricted = sparse_tree.tree[i].pair == -1 || sparse_tree.tree[j].pair == -1;

			const bool unpaired = (sparse_tree.tree[i].pair<-1 && sparse_tree.tree[j].pair<-1);
			const bool paired = (sparse_tree.tree[i].pair == j && sparse_tree.tree[j].pair == i);
			energy_t v = INF;
			// ----------------------------------------
			// cases with base pair (i,j)
			if(ptype_closing>0 && !restricted && evaluate) { // if i,j form a canonical base pair
				bool canH = (paired || unpaired);
				if(sparse_tree.up[j-1]<(j-i-1)) canH=false;
				
				energy_t v_h = canH ? HairpinE(seq,S,S1,params,i,j) : INF;
				// info of best interior loop decomposition (if better than hairpin)
				cand_pos_t best_l=0;
				cand_pos_t best_k=0;
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
				cand_pos_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
				// #pragma omp parallel for
				if((sparse_tree.tree[i].pair<-1 && sparse_tree.tree[j].pair < -1) || sparse_tree.tree[i].pair == j) { 
					for ( cand_pos_t k=i+1; k<=max_k; k++) {
						cand_pos_t k_mod=k%(MAXLOOP+1);
						
						bool cank = sparse_tree.up[k-1]>=(k-i-1) ? true : false;
						cand_pos_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;
						for (cand_pos_t l=j-1; l>=min_l; --l) {
							assert(k-i+j-l-2<=MAXLOOP);
							if(V(k_mod,l) == INF || ((sparse_tree.tree[k].pair > -1 || sparse_tree.tree[l].pair > -1) && sparse_tree.tree[k].pair != l )) continue;
							
							const energy_t v_iloop_kl = cank && sparse_tree.up[j-1]>=(j-l-1) ? V(k_mod,l) + E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S[k]][S[l]]],S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params)) : INF;
							if ( v_iloop_kl < v_iloop) {
								v_iloop = v_iloop_kl;
								best_l=l;
								best_k=k;
								best_e=V(k_mod,l);
							}
							if(sparse_tree.tree[l].pair > -1) break;
						}
						if(sparse_tree.tree[k].pair > -1) break;
					}
				}
				
				const energy_t v_split = E_MbLoop(dmli1,dmli2,S,params,i,j,sparse_tree.tree);

				v = std::min(v_h,std::min(v_iloop,v_split));
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
            
 			cand_pos_t ip1_mod = (i+1)%(MAXLOOP+1);
			energy_t vi1j = V(ip1_mod,j);
			energy_t vij1 = V(i_mod,j-1);
			energy_t vi1j1 = V(ip1_mod,j-1);	

			// Checking the dangle positions for W
			if(params->model_details.dangles == 1) d =0;
			energy_t w_v  = E_ext_Stem(v,vi1j,vij1,vi1j1,S,params,i,j,d,n,sparse_tree.tree);
			// Checking the dangle positions for W
			const energy_t wm_v = E_MLStem(v,vi1j,vij1,vi1j1,WM,CL,S,params,i,j,d,n,sparse_tree.tree);
			cand_pos_t k = i;
            cand_pos_t l = j;
			if(params->model_details.dangles == 1){
                if(d>0){
                    switch(d){
                        case 1:
                            k = i+1;
							break;
                        case 2:
                            l = j-1;
							break;
                        case 3: 
                            k = i+1;
                            l = j-1;
							break;
                    }
                    if(exists_trace_arrow_from(ta,k,l) && (wm_v < wm_split || w_v < w_split)) inc_source_ref_count(ta,k,l);	
                }
			}
			
			cand_pos_t Bp_ij = sparse_tree.Bp(i,j);
			cand_pos_t B_ij = sparse_tree.B(i,j);
			cand_pos_t b_ij = sparse_tree.b(i,j);
			cand_pos_t bp_ij = sparse_tree.bp(i,j);
			std::vector<energy_t> wiB1;
			std::vector<energy_t> wibp1;
			wiB1.resize(n+1,0);
			wibp1.resize(n+1,0);

			// Start of VP ---- Will have to change the bounds to 1 to n instead of 0 to n-1
			bool weakly_closed_ij = sparse_tree.weakly_closed(i,j);
			if ( weakly_closed_ij || sparse_tree.tree[i].pair > -1 || sparse_tree.tree[j].pair > -1 || ptype_closing == 0)	{
			
				VP(i_mod,j) = INF;
				
			}
			else{
				energy_t m1 = INF, m2 = INF, m3 = INF, m4= INF, m5 = INF, m6 = INF;
				if(sparse_tree.tree[i].parent->index > 0 && sparse_tree.tree[j].parent->index < sparse_tree.tree[i].parent->index && Bp_ij >= 0 && B_ij >= 0 && bp_ij == -1){
					// recompute_WI(wiB1,CL,CLWMB,S,params,n,B_ij+1,j-1,sparse_tree.tree);
					// int WI_ipus1_BPminus = get_WI(i+1,Bp_i - 1) ;
					energy_t WI_ipus1_BPminus = dwi1[Bp_ij-1];
					// int WI_Bplus_jminus = get_WI(B_i + 1,j-1);
					// energy_t WI_Bplus_jminus = wiB1[j-1];

					energy_t WI_Bplus_jminus = (j-1-(B_ij+1))> 4 ? WI_Bbp[j-1] : params->PUP_penalty*(j-1-(B_ij+1)+1);
					// if(i==23 && j==41) std::cout << B_ij+1 << " and " << j-1 << " " << sparse_tree.up[j-1] << std::endl;
					m1 = WI_ipus1_BPminus + WI_Bplus_jminus;
				}
				if (sparse_tree.tree[i].parent->index < sparse_tree.tree[j].parent->index && sparse_tree.tree[j].parent->index > 0 && b_ij>= 0 && bp_ij >= 0 && Bp_ij == -1){
					// recompute_WI(wibp1,CL,CLWMB,S,params,n,bp_ij+1,j-1,sparse_tree.tree);
					// int WI_i_plus_b_minus = get_WI(i+1,b_i - 1);
					energy_t WI_i_plus_b_minus = dwi1[b_ij-1];
					// int WI_bp_plus_j_minus = get_WI(bp_i + 1,j-1);
					// energy_t WI_bp_plus_j_minus = wibp1[j-1];
					energy_t WI_bp_plus_j_minus = (j-1-(bp_ij+1))> 4 ? WI_Bbp[j-1] : params->PUP_penalty*(j-1-(bp_ij+1)+1);
					m2 = WI_i_plus_b_minus + WI_bp_plus_j_minus;
				}
				
				if(sparse_tree.tree[i].parent->index > 0 && sparse_tree.tree[j].parent->index > 0 && Bp_ij >= 0 && B_ij >= 0  && b_ij >= 0 && bp_ij>= 0){
					// int WI_i_plus_Bp_minus = get_WI(i+1,Bp_i - 1);
					
					energy_t WI_i_plus_Bp_minus = dwi1[Bp_ij-1];
					
					// recompute_WI(wiB1,CL,CLWMB,S,params,n,B_ij+1,b_ij-1,sparse_tree.tree);
					// energy_t WI_B_plus_b_minus = wiB1[b_ij-1];
					energy_t WI_B_plus_b_minus = (b_ij-1-(B_ij+1))> 4 ? WI_Bbp[b_ij-1] : params->PUP_penalty*(b_ij-1-(B_ij+1)+1);
					// recompute_WI(wibp1,CL,CLWMB,S,params,n,bp_ij+1,j-1,sparse_tree.tree);
					// energy_t WI_bp_plus_j_minus = wibp1[j-1];
					energy_t WI_bp_plus_j_minus = (j-1-(bp_ij+1))> 4 ? WI_Bbp[j-1] : params->PUP_penalty*(j-1-(bp_ij+1)+1);

					// if(i==23 && j==41) std::cout << B_ij+1 << " and " << b_ij-1 << " " << sparse_tree.up[b_ij-1] << std::endl;
					// if(i==23 && j==41) std::cout << bp_ij+1 << " and " << j-1 << " " << sparse_tree.up[j-1] << std::endl;
					// if(i==23 && j==41) printf("i is %lu and j is %lu and wibp1[j-1] is %d and wiB1 is %d\n",i,j,wibp1[j-1],wiB1[b_ij-1]);
					// int WI_B_plus_b_minus = get_WI(B_i + 1,b_i - 1);
					// int WI_bp_plus_j_minus = get_WI(bp_i +1,j - 1);
					m3 = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
				}
				if(sparse_tree.tree[i+1].pair < -1 && sparse_tree.tree[j-1].pair < -1){
					// m4 = get_e_stP(i,j)+ get_VP(i+1,j-1);
					cand_pos_t ip1_mod = (i+1)%(MAXLOOP+1);
					m4 = params->e_stP_penalty*ILoopE(S,S1,params,ptype_closing,i,j,i+1,j-1) + VP(ip1_mod,j-1);
				}
				

				cand_pos_t ip, jp;
				cand_pos_t max_borders;
				cand_pos_t min_borders = 1; 
				if (Bp_ij> 1 && Bp_ij < n && b_ij >1 && b_ij < n) min_borders = std::min(Bp_ij,b_ij);
				else if (b_ij > 1 && b_ij < n && (Bp_ij < 1 || Bp_ij > n)) min_borders = b_ij;
				else if (Bp_ij > 1 && Bp_ij < n && (b_ij < 1 || b_ij > n)) min_borders = Bp_ij;
				cand_pos_t edge_i = i+MAXLOOP+1;
				min_borders = std::min({min_borders,edge_i});
				for (ip = i+1; ip < min_borders; ip++){
					bool empty_region_i = sparse_tree.up[ip-1]>=(ip-i-1);//sparse_tree.up[ip-2] > ip-1-(i+1)-1; // i+1 to ip-1
					cand_pos_t ip_mod = ip%(MAXLOOP+1);
					if (sparse_tree.tree[ip].pair < -1 && (sparse_tree.tree[i].parent->index == sparse_tree.tree[ip].parent->index) && empty_region_i){
						max_borders= 1;
						if (bp_ij > 1 && B_ij > 1) max_borders = std::max(bp_ij,B_ij);
						else if (B_ij > 1 && bp_ij < 1) max_borders = B_ij;
						else if (bp_ij > 1 && B_ij < 1) max_borders = bp_ij;
						cand_pos_t edge_j = j-31;
						max_borders = std::max({max_borders,edge_j});
						for (jp = j; jp > max_borders+1 ; jp--){
							bool empty_region_j =  sparse_tree.up[j-1]>=(j-jp-1);//sparse_tree.up[j-2] > j-1-(jp+1)-1; // jp+1 to j-1
							if (sparse_tree.tree[jp].pair < -1 && pair[S[ip]][S[jp]]>0 && empty_region_j == 1){
								//arc to arc originally
								if (sparse_tree.tree[j].parent->index == sparse_tree.tree[jp].parent->index){
									
									energy_t temp = (int)round(((j-jp==1 && ip-i==1) ? params->e_stP_penalty :params->e_intP_penalty )*ILoopE(S,S1,params,ptype_closing,i,j,ip,jp));
									temp = temp + VP(ip_mod,jp);
									
									
									if (m5 > temp){
										m5 = temp;
									}
								}
							}
						}
					}
				}

				// case 6 and 7
				cand_pos_t left_bound = std::min(Bp_ij,b_ij);
				m6 = dwvp[j-1] + params->ap_penalty + 2*params->bp_penalty;			
				
				VP(i_mod,j) = std::min({m1, m2, m3, m4, m5, m6});
				
				
				energy_t VP_j_split = INF;
				for ( auto const [key,val] : CLVP[j] ) {
					cand_pos_t k = key;
					VP_j_split = std::min(VP_j_split,val);
				}
				if(VP(i_mod,j)<VP_j_split){
					register_candidate(CLVP,i,j,VP(i_mod,j));
				}
			}
			// -------------------------------------------End of VP----------------------------------------------------------------
			
			// Start of WMBP
			
			if ((sparse_tree.tree[i].pair >= -1 && sparse_tree.tree[i].pair > j) || (sparse_tree.tree[j].pair >= -1 && sparse_tree.tree[j].pair < i) || (sparse_tree.tree[i].pair >= -1 && sparse_tree.tree[i].pair < i ) || (sparse_tree.tree[j].pair >= -1 && j < sparse_tree.tree[j].pair)) WMB[j] = INF;
			else{
				int m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF, m6 = INF;

				// 1) WMBP(i,j) = BE(bpg(Bp(l,j)),Bp(l,j),bpg(B(l,j)),B(l,j)) + WMBP(i,l) + VP(l+1,j)
				if (sparse_tree.tree[j].pair < 0){
					energy_t temp = INF;
					cand_pos_t l_min=-1;
					// for (l = i+1; l<j ; l++){
					for ( auto const [key,val] : CLVP[j] ){
						cand_pos_t l = key;
						cand_pos_t B_lj = sparse_tree.B(l,j);
						cand_pos_t Bp_lj = sparse_tree.Bp(l,j);
						if (sparse_tree.tree[l].parent->index > -1 && B_lj >= 0 && B_lj < n && Bp_lj >= 0 && Bp_lj<n){
							if (b_ij >= 0 && b_ij < n && l < b_ij){
								if (i <= sparse_tree.tree[l].parent->index && sparse_tree.tree[l].parent->index < j && l+3 <=j){
									
									energy_t BE_energy = INF;
									cand_pos_t b_lj = sparse_tree.tree[B_lj].pair;
									for ( auto const [key,val] : CLBE[sparse_tree.tree[Bp_lj].pair] ){
										cand_pos_t k = key;
										if(k == b_lj){
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
						m1 = 2*params->PB_penalty + temp;
					}
				}

				// 2) WMBP(i,j) = WMBP(i,l) + WI(l+1,j)
				// if(sparse_tree.tree[j].pair < j){
				// 	int l,l_min =-1;
				// 	for(l = i+1; l<j; l++){
				// 		if (sparse_tree.tree[l].pair < 0 && sparse_tree.tree[l].parent->index > -1 && sparse_tree.tree[j].parent->index > -1 && sparse_tree.tree[j].parent->index == sparse_tree.tree[l].parent->index){
				// 			energy_t WMBP_energy = WMBP[l];
							
				// 			std::vector<energy_t> wipBpl;
				// 			wipBpl.resize(n+1,INF);
				// 			recompute_WI(wipBpl,CL,CLWMB,S,params,n,l+1,j,sparse_tree.tree);
				// 			// energy_t WI_energy = (j-(l+1))> 4 ? WI_Bbp[j] : params->PUP_penalty*(j-(l+1)+1); // Check later
				// 			energy_t WI_energy = wipBpl[j]; 
				// 			// int temp = get_WMBP(i,l) + get_WI(l+1,j);
				// 			energy_t temp = WMBP_energy + WI_energy;
				// 			if (temp < m2){
				// 				m2 = temp;
				// 				l_min = l;
				// 			}
				// 		}
				// 	}
				// }

				if(sparse_tree.tree[j].pair < j){
					
					for ( auto const [key,val,val_wm,val_w] : CL[j] ){
						cand_pos_t k = key;
						m2 = std::min(m2,WMBP[k-1] + val);
					
					}

					for ( auto const [key,val] : CLWMB[j] ){
						cand_pos_t k = key;
						m3 = std::min(m3,WMBP[k-1] + val);
					
					}
				
				}

				if(sparse_tree.tree[j].pair < 0) m4 = WMBP[j-1] + params->PUP_penalty;

				// 3) WMB(i,j) = VP(i,j) + P_b
				m5 = VP(i_mod,j) + params->PB_penalty;

				// check later if <0 or <-1
				// WMBP(i,j) = BE(bp(i,l),bpg(bp(i,l)),i,bpg(i)) +WI(bp(i,l)+1,l-1) + VP(l,j)
				if(sparse_tree.tree[j].pair < 0 && sparse_tree.tree[i].pair >= 0){
					energy_t tmp = INF;
					cand_pos_t l, l_min=-1;
					// for (l = i+1; l < j; l++){
					for ( auto const [key,val] : CLVP[j] ){
						cand_pos_t l = key;
						cand_pos_t bp_il = sparse_tree.bp(i,l);
						if(bp_il >= 0 && bp_il < n && l+TURN <= j){
							energy_t BE_energy = INF;
							cand_pos_t Bp_il = sparse_tree.tree[bp_il].pair;
							if (!CLBE[Bp_il].empty()){
								auto const [k,vbe] = CLBE[Bp_il].back();
								if(i==k) BE_energy = vbe;
								// if(i==2 && j==49) std::cout << i << " " << k << " " << vbe << std::endl;
							}
							
							// for ( auto const [key,val] : CLBE[bp_il] ){
							// 	size_t k = key;
							// 	if(i==k){
							// 		BE_energy = val;
							// 		break;
							// 	}
							// }
							// std::vector<energy_t> wipbpl;
							// wipbpl.resize(n+1,INF);
							// recompute_WI(wipbpl,CL,CLWMB,S,params,n,bp_il+1,l-1,sparse_tree.tree);
						// int BE_energy = get_BE(i,fres[i].pair,bp_i_l,fres[bp_i_l].pair);
						// int WI_energy = get_WI(bp_i_l +1,l-1);
						
						energy_t WI_energy = (l-1-(bp_il+1))> 4 ? WI_Bbp[l-1] : params->PUP_penalty*(l-1-(bp_il+1)+1);
						// energy_t WI_energy = wipbpl[l-1];
						energy_t VP_energy = val;
						energy_t sum = BE_energy + WI_energy + VP_energy;
						// if(i==2 && j==49) std::cout << WI_energy << " " << VP_energy << " " << BE_energy << std::endl;
						// if(i==2 && j==59) std::cout << WI_energy << " " << VP_energy << " " << BE_energy << std::endl;
						if (tmp > sum){
							tmp = sum;
							l_min = l;
						}
						}
					}
					m6 = 2*params->PB_penalty + tmp;
				}
				

			// get the min for WMB
			WMBP[j] = std::min({m1,m2,m3,m4,m5,m6});
			// if(i==2 && j==49) printf("i is %lu and j is %lu and m1 is %d and m3 is %d and m4 is %d and m5 is %d\n",i,j,m1,m3,m4,m5);
			// if(i==2 && j==59) printf("i is %lu and j is %lu and m1 is %d and m3 is %d and m4 is %d and m5 is %d\n",i,j,m1,m3,m4,m5);

			// ----------------------------------------------End of WMBP------------------------------------------------------

			// Start of WMB

			
			m2 = INF;
			energy_t mWMBP = INF;
			// 2)
			if (sparse_tree.tree[j].pair >= 0 && j > sparse_tree.tree[j].pair){
				cand_pos_t bp_j = sparse_tree.tree[j].pair;
				
				for (int l = (bp_j +1); (l < j); l++){
					cand_pos_t Bp_lj = sparse_tree.Bp(l,j);
					if (Bp_lj >= 0 && Bp_lj<n){
						
						energy_t BE_energy = INF;
						for ( auto const [key,val] : CLBE[sparse_tree.tree[Bp_lj].pair] ){
							size_t k = key;
							if(k == bp_j){
								BE_energy = val;
								break;
							}
						}


						energy_t WMBP_energy = WMBP[l];
						energy_t WI_energy = (Bp_lj-1-(l+1)+1)> 4 ? WI_Bp[l] : params->PUP_penalty*(Bp_lj-1-(l+1)+1);
						energy_t sum = BE_energy + WMBP_energy + WI_energy;
						// if(i==2 && j==59) std::cout << l << " " << WI_energy << " " << WMBP_energy << " " << BE_energy << std::endl;
						if (m2 > sum){
							m2 = sum;
						}

					}
				}
				m2 = params->PB_penalty + m2;
			}
			// check the WMBP_ij value
			mWMBP =  WMBP[j];

			// get the min for WMB
			WMB[j] = std::min(m2,mWMBP);
			// if(i==2 && j==49) printf("i is %lu and j is %lu and m2 is %d and WMBP is %d\n",i,j,m2,mWMBP);

			// -------------------------------------------------------End of WMB------------------------------------------------------
			}

			// Start of WI -- the conditions on calculating WI is the same as WIP, so we combine them
			energy_t wi_v = INF;
			energy_t wip_v = INF;
			energy_t wi_wmb = INF;
			energy_t wip_wmb = INF;
	
			if (!weakly_closed_ij || sparse_tree.tree[i].parent->index != sparse_tree.tree[j].parent->index){
					WI[j] = INF;
					WIP[j] = INF;
			}
			else{
				if(ptype_closing>0 && !restricted && evaluate) {
					wi_v = V(i_mod,j) + params->PPS_penalty;
					wip_v = V(i_mod,j)	+ params->bp_penalty;
				}
				wi_wmb = WMB[j] + params->PSP_penalty + params->PPS_penalty;
				wip_wmb = WMB[j] + params->PSM_penalty + params->bp_penalty;
				WI[j] = std::min({wi_split,wi_v,wi_wmb});
				WIP[j] = std::min({wip_split,wip_v,wip_wmb});
				if(B_bool){
					// std::cout << j << std::endl;
					WI_Bbp[j] = WI[j];
					WIP_Bbp[j] = WIP[j];
				}
				if(bp_bool){
					WI_Bbp[j] = WI[j];
					WIP_Bbp[j] = WIP[j];
				}
				if(Bp_bool) WI_Bp[i] = WI[j];
			}

			// ------------------------------------------------End of Wi/Wip--------------------------------------------------
			
			// start of WV
			if (!weakly_closed_ij){
				energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF,m5 = INF, m6 = INF;
				for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>std::max(bp_ij,B_ij)+1 ; ++it ) {
					cand_pos_t k = it->first;
					energy_t val = it->second;
					m1 = std::min(m1, VPP[k-1] + val);
					m5 = std::min(m5, WV[k-1] + val);					
				}

				for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>std::max(bp_ij,B_ij)+1 ; ++it ) {
					cand_pos_t k = it->first;
					energy_t val = it->second;
					m2 = std::min(m1, VPP[k-1] + val);
					m6 = std::min(m6, WV[k-1] + val);
				}

				for ( auto it = CLVP[j].begin();CLVP[j].end()!=it; ++it ) {
					cand_pos_t k = it->first;
					energy_t val = it->second;
					if(k>std::min(Bp_ij,b_ij)) continue;
					m3 = std::min(m1, dwip1[k-1] + val);
					
				}

				m4 = WV[j-1] + params->cp_penalty;

				WV[j] = std::min({m1,m2,m3,m4,m5,m6});
			}


			// ------------------------------------------------End of WV------------------------------------------------------


			// start of VPP
			
			// This is for finding the previous VPP's for j i.e. the VPP split
			if(!weakly_closed_ij){
				energy_t m1 = INF, m2 = INF;
				for ( auto it = CLVP[j].begin();CLVP[j].end()!=it; ++it ) {
					cand_pos_t k = it->first;
					energy_t val = it->second;
					if(k>std::min(Bp_ij,b_ij)) continue;
					m1 = std::min(m1, dwip1[k-1] + val);
					m2 = std::min(m2, static_cast<energy_t>((k-i)*params->cp_penalty) + val);
					
				}
				VPP[j] = std::min(m1,m2);

			}
			// -------------------------------------------------End of VPP------------------------------------------------------------

			/*
			The order for these should be        x i        lp               jp         j           l      ip   y
												 ( (        (    (   (       (          )      )  ) )      )    )
			i and ip are the outer base pair
			lp and l are the closest encompassing base pair to i/ip
			jp and j are some inner base pair;j has the be the closing due to the j=i+4 setup we have
			*/
			// Start of BE
			// Cannot use size_t due to overflow when pair is -2. We also cannot compare size_t and int for that same reason
			cand_pos_t ip = sparse_tree.tree[i].pair; // i's pair ip should be right side so ip = )
			cand_pos_t jp = sparse_tree.tree[j].pair; // j's pair jp should be left side so jp = ()
			// if(i==1 && j==24) std::cout << "pair is " << sparse_tree.tree[i].pair << std::endl;
			// if (!( i >= 0 && i <= ip && ip < jp && jp <= j && j < n && fres[i].pair >= -1 && fres[j].pair >= -1 && fres[ip].pair >= -1 && fres[jp].pair >= -1 && fres[i].pair == j && fres[j].pair == i && fres[ip].pair == jp && fres[jp].pair == ip)){ //impossible cases
			
			// if(i==2 && j==11) printf("i is %lu and j is %lu and ip is %lu and jp is %lu\n",i,j,ip,jp);
			// if(i==10 && j==24) printf("i is %lu and j is %lu and ip is %lu and jp is %lu\n",i,j,ip,jp);
			// base case: i.j and ip.jp must be in G
			if (ip > i && j > jp && jp > i && ip > j){ // Don't need to check if they are pairs separately because it is checked by virtue of these checks
			// if(i==2 && j==11) std::cout << "here" << std::endl;
				
				// We are checking for the closest pair that we have already calculated to i/ip from j/jp 
				// If there is nothing, then i is the closest encompassing pair to jp
				// If it is not, then we get the energy for everything from jp to lp so that we calculate less
				energy_t BE_energy = 0;
				cand_pos_t lp = jp;
				if (!CLBE[j].empty()){
					auto const [k,vbe] = CLBE[j].back();
					BE_energy = vbe;
					lp = k;
				}
				cand_pos_t l = sparse_tree.tree[lp].pair; // right closing base for lp
				 
				// printf("i is %lu and j is %lu and ip is %d and jp is %d and lp is %d and l is %d\n",i,j,ip,jp,lp,l);

				// for ( auto const [key,val] : CLBE[jp] ){ //It was CLBE[fres[jp].pair but I feel it should just be jp since that is how we save it
				// 	size_t k = key;
				// 	BE_energy = val;
				// 	lp = k;
					
				// }
				const pair_type ptype_closing_iip = pair[S[i]][S[ip]];
				energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF;
				// if (fres[i+1].pair == ip-1){
				if(i+1 == lp and ip-1 == l){
					m1 = (energy_t)round(params->e_stP_penalty*(double)ILoopE(S,S1,params,ptype_closing_iip,i,ip,lp,l)) + BE_energy;
				}
				

				// 	// for (int k = i+1; k<= lp ; k++){
				// 	// if (fres[l].pair >= -1 && j <= fres[k].pair && fres[k].pair < ip){
				bool empty_region_ilp = sparse_tree.up[lp-1] >= lp-i-1; //empty between i+1 and lp-1
				bool empty_region_lip = sparse_tree.up[ip-1] >= ip-l-1; // empty between l+1 and ip-1
				bool weakly_closed_ilp = sparse_tree.weakly_closed(i+1,lp-1); // weakly closed between i+1 and lp-1
				bool weakly_closed_lip = sparse_tree.weakly_closed(l+1,ip-1); // weakly closed between l+1 and ip-1
				// if(i==7) printf("empty ilp is %d and and i is %lu and lp is %d and up is %d\n",empty_region_ilp,i,lp,sparse_tree.up[lp-1]);
				// if(i==7) printf("empty ilp is %d and empty lip is %d\n",empty_region_ilp,empty_region_lip);

				// std::vector<energy_t> wiilp;
				// std::vector<energy_t> wilip;

				// if(weakly_closed_ilp){
				// 	wiilp.resize(n+1,INF);
				// 	recompute_WIP(wiilp,CL,CLWMB,S,params,n,i+1,lp-1,sparse_tree.tree,sparse_tree.up);
				// }
				// if(weakly_closed_lip){
				// 	wilip.resize(n+1,INF);
				// 	recompute_WIP(wiilp,CL,CLWMB,S,params,n,l+1,ip-1,sparse_tree.tree,sparse_tree.up);
				// }
					
				if (empty_region_ilp && empty_region_lip){
					m2 = (energy_t)round(params->e_intP_penalty*(double)ILoopE(S,S1,params,ptype_closing_iip,i,ip,lp,l)) + BE_energy;
				}

					// 3)
				if (weakly_closed_ilp && weakly_closed_lip){
					// get_WIP(i+1,l-1)
					// get_WIP(lp+1,ip-1)
					// m3 = wiilp[l-1] + BE_energy + wilip[ip-1]+ params->ap_penalty + 2*params->bp_penalty;
					m3 = WIP_Bbp[l-1] + BE_energy + WIP_Bbp[ip-1]+ params->ap_penalty + 2*params->bp_penalty;
				}

					// 4)
				if (weakly_closed_ilp && empty_region_lip){
					// m4 = wiilp[l-1] + BE_energy + params->cp_penalty * (ip-lp+1) + params->ap_penalty + 2*params->bp_penalty;
					m4 = WIP_Bbp[l-1] + BE_energy + params->cp_penalty * (ip-lp+1) + params->ap_penalty + 2*params->bp_penalty;
				}

					// 5)
				if (empty_region_ilp && weakly_closed_lip){
					// m5 = params->ap_penalty + 2*params->bp_penalty + (params->cp_penalty * (l-i+1)) + BE_energy + wilip[ip-1];
					m5 = params->ap_penalty + 2*params->bp_penalty + (params->cp_penalty * (l-i+1)) + BE_energy + WIP_Bbp[ip-1];
				}
					// }
				// printf("i is %lu and ip is %d and j is %lu and jp is %d and lp is %d and l is %d and m1 is %d and m2 is %d\n",i,ip,j,jp,lp,l,m1,m2);
				// printf("i is %lu and ip is %d and j is %lu and jp is %d and lp is %d and l is %d and m1 is %d and m2 is %d and m3 is %d and m4 is %d and m5 is %d\n",i,ip,j,jp,lp,l,m1,m2,m3,m4,m5);
				register_candidate(CLBE,i,j,std::min({m1,m2,m3,m4,m5}));
			} 

				// finding the min and putting it in BE[iip]
				// BE[ip] = std::min({m1,m2,m3,m4,m5});
				
			// // ------------------------------------------------End of BE---------------------------------------------------------
			//Things that needed to happen later like W's wmb
			const energy_t w_wmb = WMB[j] + params->PS_penalty;
			const energy_t wm_wmb = WMB[j] + params->PSM_penalty + params->b_penalty;
			// if(i==2 && j ==49) printf("i is %lu and j is %lu and w_wmb is %d\n",i,j,w_wmb);
			// if(i==2 && j ==59) printf("i is %lu and j is %lu and w_wmb is %d\n",i,j,w_wmb);
			w  = std::min(w_v, w_split);
			wm = std::min(wm_v, wm_split);
			if(!pairedkj) { // IS even pairedkj needed?
				w =std::min(w,w_wmb);
				wm = std::min(wm,wm_wmb);
			}
			if ( w_v < w_split || wm_v < wm_split || wi_v < wi_split || wip_v < wip_split || paired) { //wi_v < wi_split || wip_v < wip_split ||
				cand_pos_t k_mod = k%(MAXLOOP+1);
				register_candidate(CL, i, j,V(i_mod,j), encode((int) wm_v,d),encode((int) w_v,d));
				// always keep arrows starting from candidates
				inc_source_ref_count(ta,i,j);
			}	
			if ( w_wmb < w_split || wm_wmb < wm_split || wi_wmb < wi_split || wip_wmb < wip_split) {
		
				register_candidate(CLWMB, i, j, WMB[j]);

				// always keep arrows starting from candidates
				inc_source_ref_count(ta,i,j); // should i increment this?
			}	

			W[j]       = w;
			WM[j]      = wm;
			WM2[j]     = wm2_split;

		} // end loop j
		rotate_arrays(WM2,dmli1,dmli2,WI,WIP,dwi1,dwip1,WMB,dwmbi,WV,dwvp);

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
 * @brief Sums the number of Candidates at each index over all indices
 * 
 * @param CL_ Candidate list
 * @return total number of candidates
 */
size_t num_of_candidates(auto const& CL_)  {
	cand_pos_t c=0;
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
	cand_pos_t c=0;
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
	cand_pos_t n = seq.length();

	std::string restricted;
    args_info.input_structure_given ? restricted = input_structure : restricted = std::string (n,'.');

	if(restricted != "" && restricted.length() != n ){
		std::cout << "input sequence and structure are not the same size" << std::endl;
		exit(0);
	}

	std::string file= "";
	args_info.paramFile_given ? file = parameter_file : file = "";
	if(file!=""){
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
	}

	

	bool verbose;
	verbose = args_info.verbose_given;

	bool mark_candidates;
	mark_candidates = args_info.mark_candidates_given;
	sparse_tree tree(restricted,n);


	SparseMFEFold sparsemfefold(seq,!args_info.noGC_given,restricted);

	if(args_info.dangles_given) sparsemfefold.params_->model_details.dangles = dangles;
	
	cmdline_parser_free(&args_info);

	std::cout << seq << std::endl;
	
	energy_t mfe = fold(sparsemfefold.seq_, tree,sparsemfefold.V_,sparsemfefold.cand_comp,sparsemfefold.CL_,sparsemfefold.CLWMB_,sparsemfefold.CLVP_,sparsemfefold.CLBE_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.params_,sparsemfefold.ta_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_, sparsemfefold.dmli1_, sparsemfefold.dmli2_,sparsemfefold.VP_, sparsemfefold.VPP_, sparsemfefold.WV_, sparsemfefold.dwvp_, sparsemfefold.WMB_,sparsemfefold.dwmbi_,sparsemfefold.WMBP_,sparsemfefold.WI_,sparsemfefold.dwi1_,sparsemfefold.WIP_,sparsemfefold.dwip1_,sparsemfefold.n_,sparsemfefold.garbage_collect_);		
	// std::cout << mfe << std::endl;
	// std::string structure = trace_back(sparsemfefold.seq_,sparsemfefold.CL_,sparsemfefold.cand_comp,sparsemfefold.structure_,sparsemfefold.params_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.ta_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_,sparsemfefold.n_,tree.tree,tree.up, mark_candidates);
	
	std::ostringstream smfe;
	// smfe << std::setiosflags(std::ios::fixed) << std::setprecision(2) << mfe/100.0 ;

	// std::cout << structure << " ("<<smfe.str()<<")"<<std::endl;

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
