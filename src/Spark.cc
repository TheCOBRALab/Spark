/**
* @mainpage
*
* Space-efficient sparse variant of an RNA (loop-based) free energy
* minimization algorithm (RNA folding equivalent to the Zuker
* algorithm).
*
* The results are equivalent to RNAfold.
*/
// #define NDEBUG
#include "base_types.hh"
#include "PK_globals.hh"
#include "cmdline.hh"
#include "matrix.hh"
#include "trace_arrow.hh"
#include "sparse_tree.cc"
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <cstring>
#include <string>
#include <cassert>
#include <sys/stat.h>
#include <fstream>

extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/params/io.h"
}
static bool pseudoknot = false;

static bool pk_only = false;

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

energy_t ILoopE(const short* S, const short* S1, const paramT* params, const pair_type& ptype_closing,const cand_pos_t &i,const cand_pos_t &j,const cand_pos_t &k,const cand_pos_t &l);

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

	vrna_param_t *params_;

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
	std::vector<energy_t> WVe_;
	std::vector<energy_t> WMB_;
	std::vector<energy_t> dwmbi_; // WMB from 1 iteration ago
	std::vector<energy_t> WMBP_;
	std::vector<energy_t> WMBA_;
	std::vector<energy_t> WI_;
	std::vector<energy_t> dwi1_; // WI from 1 iteration ago
	std::vector<energy_t> WIP_;
	std::vector<energy_t> dwip1_; // WIP from 1 iteration ago
	std::vector<energy_t> WV_;
	std::vector<energy_t> dwvp_; // WV from 1 iteration ago;

	std::vector<energy_t> WI_Bbp; // WI from band borders on left
	std::vector<energy_t> WIP_Bbp; // WIP from band borders on left
	std::vector<energy_t> WI_Bp; // WI from band borders on right
	std::vector<energy_t> WIP_Bp; // WIP from band borders on the right

	bool mark_candidates_;


	TraceArrows ta_;
	TraceArrows taVP_;


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
	struct Cand_comp{
	bool operator ()(const cand_entry_t &x, cand_pos_t y) const {
		return x.first > y;
	}
	bool operator ()(const cand_entry_td1 &x, cand_pos_t y) const {
		return x.first > y;
	}
	} cand_comp;

	


	SparseMFEFold(const std::string &seq, bool garbage_collect, std::string restricted)
	: seq_(seq), n_(seq.length()), params_(scale_parameters()), garbage_collect_(garbage_collect), ta_(n_),taVP_(n_)
	{
	make_pair_matrix();

	S_ = encode_sequence(seq.c_str(),0);
	S1_ = encode_sequence(seq.c_str(),1);

	V_.resize(MAXLOOP+1,n_+1,INF);
	W_.resize(n_+1,0);
	WM_.resize(n_+1,INF);
	WM2_.resize(n_+1,INF);
	dmli1_.resize(n_+1,INF);
	dmli2_.resize(n_+1,INF);
	
	// Pseudoknot portion

	VP_.resize(MAXLOOP+1,n_+1,INF);
	WVe_.resize(n_+1,INF);
	WMB_.resize(n_+1,INF);
	dwmbi_.resize(n_+1,INF);
	WMBP_.resize(n_+1,INF);
	WMBA_.resize(n_+1,INF);
	WI_.resize(n_+1,0);
	dwi1_.resize(n_+1,0);
	WIP_.resize(n_+1,INF);
	dwip1_.resize(n_+1,INF);
	WV_.resize(n_+1,INF);
	dwvp_.resize(n_+1,INF);

	WI_Bbp.resize(n_+1,0);
	WIP_Bbp.resize(n_+1,INF);
	WI_Bp.resize(n_+1,0);
	WIP_Bp.resize(n_+1,INF);

	// init candidate lists
	CL_.resize(n_+1);
	CLWMB_.resize(n_+1);
	CLVP_.resize(n_+1);
	CLBE_.resize(n_+1);

	resize(ta_,n_+1);
	resize(taVP_,n_+1);

	// resize(ta_dangle_,n_+1);

	restricted_ = restricted;
	
	}

	

	~SparseMFEFold() {
	free(params_);
	free(S_);
	free(S1_);
	}
};

void trace_V(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >& WIP_Bbp,const std::vector<energy_t> &WIP_Bp,const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_W(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, const std::vector<energy_t> & W, std::vector< energy_t > &WM, std::vector< energy_t > &WM2,const std::vector< energy_t >&WI_Bbp,const std::vector<energy_t> &WIP_Bbp,const std::vector<energy_t> &WIP_Bp, const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates, cand_pos_t i, cand_pos_t j, sparse_tree &tree);
void trace_WM(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params,const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector<energy_t> &WIP_Bp,const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_WM2(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector<energy_t> &WIP_Bp,const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, sparse_tree &tree);
void trace_WMB(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector<energy_t> &WIP_Bp,const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_VP(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector<energy_t> &WIP_Bp,const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_WI(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector<energy_t> &WIP_Bp,const std::vector<energy_t> &WI_Bp,const std::vector< energy_t >&WI, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_WIP(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector<energy_t> &WIP_Bp,const std::vector<energy_t> &WI_Bp,const std::vector<energy_t> &WIP, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_WV(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector<energy_t> &WIP_Bp,const std::vector<energy_t> &WI_bBp,const std::vector<energy_t> &WVe,const std::vector<energy_t> &WIP, const std::vector<energy_t> &WV, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_WVe(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector<energy_t> &WIP_Bp,const std::vector<energy_t> &WI_bBp,const std::vector<energy_t> WVe, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_WMBP(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector<energy_t> &WIP_Bp,const std::vector<energy_t> &WI_Bp, std::vector<energy_t> &WMBP,std::vector<energy_t> &WMBA, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_BE(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector<energy_t> &WIP_Bp,const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t ip, energy_t e, sparse_tree &tree);
void trace_WMBA(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector<energy_t> &WIP_Bp,const std::vector<energy_t> &WI_Bp, std::vector<energy_t> &WMBP,std::vector<energy_t> &WMBA, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);


/**
* @brief Rotate WM2 and WI arrays to store the previous and previous previous iterations
* @param WM2 WM2 array
* @param dmli1 WM2 from one iteration ago
* @param dmli2 WM2 from two iterations ago
* @param WI WI array
* @param dwi1 WM2 from one iteration ago
* @param WIP WIP array
* @param dwip1 WIP from one iteration ago
* @param WV WV array
* @param dwvp WV from one iteration ago
*/
void rotate_arrays(std::vector<energy_t> &WM2, std::vector<energy_t> &dmli1, std::vector<energy_t> &dmli2, std::vector<energy_t> &WI, std::vector<energy_t> &WIP, std::vector<energy_t> &dwi1, std::vector<energy_t> &dwip1,std::vector<energy_t> &WMB, std::vector<energy_t> &dwmbi, std::vector<energy_t> &WV, std::vector<energy_t> &dwvp){
	dmli2.swap(dmli1);
    dmli1.swap(WM2);
	dwi1.swap(WI);
	dwip1.swap(WIP);
	dwmbi.swap(WMB);
	dwvp.swap(WV);
}

/**
 * @brief This code returns the hairpin energy for a given base pair.
 * @param i The left index in the base pair
 * @param j The right index in the base pair
*/
energy_t HairpinE(const std::string& seq, const short* S, const short* S1,  const paramT* params, cand_pos_t i, cand_pos_t j) {
	
	const pair_type ptype_closing = pair[S[i]][S[j]];

	if (ptype_closing==0) return INF;

	return E_Hairpin(j-i-1,ptype_closing,S1[i+1],S1[j-1],&seq.c_str()[i-1], const_cast<paramT *>(params));
}

/**
 * @brief Returns the internal loop energy for a given i.j and k.l
 * 
*/
energy_t ILoopE(const short* S, const short* S1, const paramT* params, const pair_type& ptype_closing,const cand_pos_t &i,const cand_pos_t &j,const cand_pos_t &k,const cand_pos_t &l)  {
	assert(ptype_closing>0);
	assert(1<=i);
	assert(i<k);
	assert(k<l);
	assert(l<j);

	// note: enclosed bp type 'turned around' for lib call
	const pair_type ptype_enclosed = rtype[pair[S[k]][S[l]]];

	// if (ptype_enclosed==0) return INF;

	return E_IntLoop(k-i-1,j-l-1,ptype_closing,ptype_enclosed,S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params));
}

/**
 * @brief Gives the W(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
 * @param vij The V(i,j) energy
 * @param vi1j The V(i+1,j) energy
 * @param vij1 The V(i,j-1) energy
 * @param vi1j1 The V(i+1,j-1) energy
*/
energy_t E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, paramT* params, const cand_pos_t i,const cand_pos_t j, Dangle &d, cand_pos_t n,const std::vector<Node> &tree){

	energy_t e = INF,en = INF;
  	pair_type tt  = pair[S[i]][S[j]];
	
    if ((tree[i].pair <-1 && tree[j].pair <-1) || (tree[i].pair == j && tree[j].pair == i)) {
				en = vij; // i j

				if (en != INF) {
					if (params->model_details.dangles == 2){
						base_type si1 = i>1 ? S[i-1] : -1;
                		base_type sj1 = j<n ? S[j+1] : -1;
                        en += vrna_E_ext_stem(tt, si1, sj1, params);
					}
                    else{
                        en += vrna_E_ext_stem(tt, -1, -1, params);
						d = 0;
					}

                    e = MIN2(e, en);
					
				}

	}

	if(params->model_details.dangles  == 1){
        tt  = pair[S[i+1]][S[j]];
        if (((tree[i+1].pair <-1 && tree[j].pair <-1) || (tree[i+1].pair == j)) && tree[i].pair<0) {
            en = (j-i-1>TURN) ? vi1j : INF; //i+1 j

            if (en != INF) {

                base_type si1 = S[i];
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

                base_type sj1 = S[j];

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

                base_type si1 = S[i];
                base_type sj1 = S[j];

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
* @brief Computes the multiloop V contribution. This gives back essentially VM(i,j).
* 
* @param dmli1 Row of WM2 from one iteration ago
* @param dmli2 Row of WM2 from two iterations ago 
*/
energy_t E_MbLoop( const std::vector<energy_t>& dmli1, const std::vector<energy_t>& dmli2, const short* S, paramT* params, cand_pos_t i, cand_pos_t j, const std::vector<Node> &tree){

	energy_t e = INF,en = INF;
  	pair_type tt  = pair[S[j]][S[i]];
	bool pairable = (tree[i].pair <-1 && tree[j].pair <-1) || (tree[i].pair == j);
	
	/* double dangles */
	switch(params->model_details.dangles){
		case 2:
			if (pairable) {
			e = dmli1[j - 1];

			if (e != INF) {

				base_type si1 = S[i + 1];
				base_type sj1 = S[j - 1];

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

          			base_type si1 =  S[i + 1];

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
					base_type sj1 = S[j - 1];

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

					base_type si1 = S[i + 1];
					base_type sj1 = S[j - 1];

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
energy_t E_MLStem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, paramT* params,cand_pos_t i, cand_pos_t j,Dangle &d, const  cand_pos_t& n, const std::vector<Node> &tree){

	energy_t e = INF,en=INF;

	pair_type type = pair[S[i]][S[j]];

	
	if ((tree[i].pair < -1 && tree[j].pair < -1) || (tree[i].pair == j)) {
		en = vij; // i j
		if (en != INF) {
			if (params->model_details.dangles == 2){
				base_type mm5 = i>1 ? S[i-1] : -1;
            	base_type mm3 = j<n ? S[j+1] : -1;
				en += E_MLstem(type, mm5, mm3, params);
			}
			else{
				en += E_MLstem(type, -1, -1, params);
				d = 0;
			}
			e = MIN2(e, en);
		}
	}
	if(params->model_details.dangles == 1){
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

/**
 * In cases where the band border is not found, if specific cases are met, the value is Inf(n) not -1.
 * Mateo Jan 2025: Added to Fix WMBP problem
*/
int compute_exterior_cases(cand_pos_t l, cand_pos_t j, sparse_tree &tree){

	// Case 1 -> l is not covered
	bool case1 = tree.tree[l].parent->index <= 0;
	// Case 2 -> l is paired
	bool case2 = tree.tree[l].pair > 0;
	// Case 3 -> l is part of a closed subregion
	bool case3 = 0;
	// Case 4 -> l.bp(l) i.e. l.j does not cross anything -- could I compare parents instead?
	bool case4 = j<tree.Bp(l,j);
	// By bitshifting each one, we have a more granular idea of what cases fail and is faster than branching
	return (case1 <<3) | (case2 << 2) | (case3 << 1) | case4;
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
const std::vector<energy_t> recompute_WM(const std::vector<energy_t>& WM, const std::vector< cand_list_td1 > &CL,const std::vector< cand_list_t > &CLWMB, paramT* params, cand_pos_t i, cand_pos_t max_j, const std::vector<Node> &tree, const std::vector<cand_pos_t> &up) {
	assert(i>=1);
	assert(max_j<=n);
	std::vector<energy_t> temp = WM;

	for ( cand_pos_t j=i-1; j<=std::min(i+TURN,max_j); j++ ) { temp[j]=INF; }
	
	for ( cand_pos_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wm = INF;
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			cand_pos_t k = it->first;
			energy_t v_kj = it->third >> 2;

			bool can_pair = up[k-1] >= (k-i);
			if(can_pair) wm = std::min( wm, static_cast<energy_t>(params->MLbase*(k-i)) + v_kj );
			wm = std::min( wm, temp[k-1]  + v_kj );
		}
		for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>i+TURN+1 ; ++it ) {
			
			cand_pos_t k = it->first;
			energy_t wmb_kj = it->second +PSM_penalty + b_penalty;
			bool can_pair = up[k-1] >= (k-i);
			if(can_pair) wm = std::min( wm, static_cast<energy_t>(params->MLbase*(k-i)) + wmb_kj );

			wm = std::min( wm, temp[k-1]  + wmb_kj );
		}
		if(tree[j].pair<0) wm = std::min(wm, temp[j-1] + params->MLbase);
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
const std::vector<energy_t> recompute_WM2(const std::vector<energy_t>& WM,const std::vector<energy_t>& WM2, const std::vector< cand_list_td1 > &CL,const std::vector< cand_list_t > &CLWMB, paramT* params, cand_pos_t i, cand_pos_t max_j, const std::vector<Node> &tree, const std::vector<cand_pos_t> &up) {
	assert(i>=1);
	assert(max_j<=n);
	std::vector<energy_t> temp = WM2;

	for ( cand_pos_t j=i-1; j<=std::min(i+2*TURN+2,max_j); j++ ) { temp[j]=INF; }

	for ( cand_pos_t j=i+2*TURN+3; j<=max_j; j++ ) {
		energy_t wm2 = INF;
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>i+TURN+1 ; ++it ) {
			
			cand_pos_t k = it->first;
			energy_t v_kj = it->third >> 2;
			
			wm2 = std::min( wm2, WM[k-1]  + v_kj );
		}
		for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i ; ++it ) {
			
			cand_pos_t k = it->first;
			energy_t wmb_kj = it->second +PSM_penalty + b_penalty;
			bool can_pair = up[k-1] >= (k-i);

			wm2 = std::min( wm2, WM[k-1]  + wmb_kj );
			if(can_pair) wm2 = std::min( wm2, static_cast<energy_t>(params->MLbase*(k-i)) + wmb_kj );

		}
		if(tree[j].pair<0) wm2 = std::min(wm2, temp[j-1] + params->MLbase);
		temp[j] = wm2;
	}
	return temp;
}

/**
* @brief Recompute row of WMBP. This is used in the traceback when we haved decided the current ij pair closes a psuedoknot,
* and the WMBP energies need to be recomputed such that the pseudoknot can be broken down. This happens infrequently enough such that it not too costly
* 
* @param CL V candidate list
* @param CLVP VP candidate list
* @param CLWMB WMB Candidate list
* @param CLBE BE candidate list
* @param WI_Bbp array which holds the values where the left index is either an outer closing base or an inner opening base
* @param i Current i
* @param max_j Current j
*/
const std::vector<energy_t> recompute_WMBP(const std::vector< cand_list_td1 > &CL,const std::vector< cand_list_t > &CLVP,const std::vector< cand_list_t > &CLWMB,const std::vector< cand_list_t > &CLBE, const std::vector<energy_t> &WI_Bbp,std::vector<energy_t> &WMBA, const cand_pos_t& n, cand_pos_t i, cand_pos_t max_j, sparse_tree &tree) {
	assert(i>=1);
	assert(max_j<=n);
	std::vector<energy_t> temp;
	temp.resize(n+1,INF);

	for ( cand_pos_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wmbp = INF;
		energy_t wmba = INF;

		energy_t vp_ij = INF;
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			
			cand_pos_t k = it->first;
			energy_t v_kj = it->second;
			
			wmba = std::min(wmba, temp[k-1]  + v_kj + PPS_penalty);
		}
		for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i ; ++it ) {
			
			cand_pos_t k = it->first;
			energy_t wmb_kj = it->second;
			
			wmba = std::min(wmba, temp[k-1]  + wmb_kj + PPS_penalty + PSM_penalty);
		}

		for ( auto it = CLVP[j].begin();CLVP[j].end()!=it && it->first>=i ; ++it ) {
			cand_pos_t k = it->first;
			if(k==i) vp_ij = it->second; //second?
			cand_pos_t bp_ik = tree.bp(i,k);
			cand_pos_t Bp_kj = tree.Bp(k,j);
			cand_pos_t b_ij = tree.b(i,j);
			int ext_case = compute_exterior_cases(k,j,tree);
			if ((b_ij > 0 && k < b_ij) || ext_case == 0){
				if (bp_ik >= 0 && k>bp_ik && Bp_kj > 0 && k<Bp_kj){
					energy_t BE_energy = INF;
					cand_pos_t B_kj = tree.B(k,j);
					cand_pos_t b_kj = (B_kj>0) ? tree.tree[B_kj].pair : -2;
					for ( auto it2 = CLBE[Bp_kj].begin();CLBE[Bp_kj].end()!=it2; ++it2 ) {
						cand_pos_t l = it2->first;
						if(l == b_kj){
							BE_energy = it2->second;
							break;
						}
					}
					wmbp = std::min(wmbp, BE_energy + WMBA[k-1]  + it->second + 2*PB_penalty);
				}
			}
		}

		for ( auto it = CLVP[j].begin();CLVP[j].end()!=it && it->first>=i ; ++it ) {
			cand_pos_t k = it->first;

			cand_pos_t bp_ik = tree.bp(i,k);
			if(bp_ik >= 0 && k+TURN <= j){
				cand_pos_t Bp_ik = tree.tree[bp_ik].pair;
				energy_t BE_energy = INF;
				for ( auto it2 = CLBE[Bp_ik].begin();CLBE[Bp_ik].end()!=it2; ++it2 ) {
					cand_pos_t l = it2->first;
					if(l == i){
						BE_energy = it2->second;
						break;
					}
				}
				wmbp = std::min(wmbp, BE_energy + WI_Bbp[k-1]  + it->second +2*PB_penalty);

			}
		}
		wmbp = std::min(wmbp, vp_ij + PB_penalty);
		if(tree.tree[j].pair<0) wmba = std::min(wmba, WMBA[j-1] + PUP_penalty);
		wmba = std::min(wmba,wmbp);
		WMBA[j] = wmba;
		temp[j] = wmbp;
	}

	return temp;

}

/**
* @brief Recompute row of WI. This is used in the traceback when we haved decided the current ij pair is a crossing base pair, we need to calculate the energy,
* and the WI energies need to be recomputed such that the VP energy can be broken down. This happens only when the calculations are confirmed to be needed such that it not too costly
* 
* @param CL V candidate list
* @param CLWMB WMB candidate list
* @param i Current i
* @param max_j Current j
*/
const std::vector<energy_t> recompute_WI(const std::vector< cand_list_td1 > &CL,const std::vector< cand_list_t> &CLWMB, const cand_pos_t& n, cand_pos_t i, cand_pos_t max_j, const std::vector<Node> &tree, const std::vector<cand_pos_t> &up) {
	assert(i>=1);
	assert(max_j<=n);
	std::vector<energy_t> temp;
	temp.resize(n+1,0);

	// Causes vector resize error if the ifs are not there because if i is close to n, it would go past the bounds
	if(i<n){ temp[i] = PUP_penalty;} if(i+1<n){ temp[i+1] = 2*PUP_penalty;} if(i+2<n){ temp[i+2] = 3*PUP_penalty;} if(i+3<n){ temp[i+3] = 4*PUP_penalty;}

	for ( cand_pos_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wi = INF;
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			
			cand_pos_t k = it->first;
			energy_t v_kj = it->second + PPS_penalty;
			
			wi = std::min(wi, temp[k-1]  + v_kj );
			bool can_pair = up[k-1] >= (k-i);
			if(can_pair) wi = std::min( wi, static_cast<energy_t>(PUP_penalty*(k-i)) + v_kj );
		}
		for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i ; ++it ) {
			
			cand_pos_t k = it->first;
			energy_t wmb_kj = it->second + PSM_penalty + PPS_penalty;
			
			wi = std::min(wi, temp[k-1]  + wmb_kj );
			bool can_pair = up[k-1] >= (k-i);
			if(can_pair) wi = std::min( wi, static_cast<energy_t>(PUP_penalty*(k-i)) + wmb_kj );
		}
		if(tree[j].pair<0) wi = std::min(wi, temp[j-1] + PUP_penalty);
		temp[j] = wi;
	
	}
	return temp;

}

/**
* @brief Recompute row of WIP. This is used in the traceback when we haved decided the current ij pair is a multiloop that spans a band or part of a BE, we need to calculate the energy,
* and the WIP energies need to be recomputed such that the VP or BE energy can be broken down. This happens only when the calculations are confirmed to be needed such that it not too costly
* 
* @param CL V candidate list
* @param CLWMB WMB candidate list
* @param i Current i
* @param max_j Current j
*/
const std::vector<energy_t> recompute_WIP(const std::vector< cand_list_td1 > &CL,const std::vector< cand_list_t> &CLWMB, const cand_pos_t& n, cand_pos_t i, cand_pos_t max_j, const std::vector<Node> &tree, const std::vector<cand_pos_t> &up) {
	assert(i>=1);
	assert(max_j<=n);
	std::vector<energy_t> temp;
	temp.resize(n+1,INF);
	for ( cand_pos_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wip = INF;
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i; ++it ) {
			
			cand_pos_t k = it->first;
			energy_t v_kj = it->second + bp_penalty;
			
			wip = std::min(wip, temp[k-1]  + v_kj );
			bool can_pair = up[k-1] >= (k-i);
			if(can_pair) wip = std::min( wip, static_cast<energy_t>(cp_penalty*(k-i)) + v_kj );
		}
		for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i; ++it ) {
			
			cand_pos_t k = it->first;
			energy_t wmb_kj = it->second + PSM_penalty + bp_penalty;
			
			wip = std::min(wip, temp[k-1]  + wmb_kj );
			bool can_pair = up[k-1] >= (k-i);
			if(can_pair) wip = std::min( wip, static_cast<energy_t>(cp_penalty*(k-i)) + wmb_kj );
		}
		if(tree[j].pair<0) wip = std::min(wip, temp[j-1] + cp_penalty);
		temp[j] = wip;
	
	}
	return temp;

}
/**
* @brief Recompute row of WVe. This is used in the traceback when we haved decided the current ij pair is a multiloop that spans a band, we need to calculate the energy,
* and the WVe energies need to be recomputed such that the VP energy can be broken down. This happens infrequently such that it not too costly
* 
* @param CLVP VP candidate list
* @param i Current i
* @param max_j Current j
*/
const std::vector<energy_t> recompute_WVe(const std::vector< cand_list_t > &CLVP, const cand_pos_t& n, cand_pos_t i, cand_pos_t max_j, sparse_tree &tree) {
	assert(i>=1);
	assert(max_j<=n);
	std::vector<energy_t> temp;
	temp.resize(n+1,INF);
	for ( cand_pos_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wve = INF;
		for ( auto it = CLVP[j].begin();CLVP[j].end()!=it && it->first>=i ; ++it ) {
			cand_pos_t k = it->first;
			energy_t vp_kj = it->second;

			bool can_pair = tree.up[k-1] >= (k-i);
			if(can_pair) wve = std::min( wve, static_cast<energy_t>(cp_penalty*(k-i)) + vp_kj );
		}
		if(tree.tree[j].pair<0) wve = std::min(wve, temp[j-1] + cp_penalty);
		temp[j] = wve;
	}
	return temp;
}
/**
* @brief Recompute row of WV. This is used in the traceback when we haved decided the current ij pair is a multiloop that spans a band, we need to calculate the energy,
* and the WV energies need to be recomputed such that the VP energy can be broken down. This happens infrequently such that it not too costly
* 
* @param CL V candidate list
* @param CLWMB WMB candidate list
* @param i Current i
* @param CLVP VP candidate list
* @param WIP WIP energies
* @param i Current i
* @param max_j Current j
*/
const std::vector<energy_t> recompute_WV(const std::vector<energy_t>& WVe, const std::vector< cand_list_td1 > &CL, const std::vector< cand_list_t > &CLWMB, const std::vector< cand_list_t > &CLVP,const std::vector<energy_t> &WIP, const cand_pos_t& n, cand_pos_t i, cand_pos_t max_j, sparse_tree &tree) {
	assert(i>=1);
	assert(max_j<=n);
	std::vector<energy_t> temp;
	temp.resize(n+1,INF);
	
	for ( cand_pos_t j=i+TURN+1; j<=max_j; j++ ) {
		cand_pos_t bound_right = std::max(tree.bp(i,j),tree.B(i,j));
		cand_pos_t bound_left = std::min((cand_pos_tu) tree.Bp(i,j),(cand_pos_tu) tree.b(i,j));

		energy_t wv = INF;
		if(!tree.weakly_closed(i,j)){
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>bound_right ; ++it ) {
			
			cand_pos_t k = it->first;
			// energy_t v_kj = it->third >> 2;
			energy_t v_kj = it->second;
			// Dangle d = it->third & 3;
			// cand_pos_t num = 0;
			// if(d == 1 || d==2) num = 1;
			// else if(d==3 && params->model_details.dangles == 1) num =2;
			// energy_t fix = num*cp_penalty - num*params->MLbase;
			
			wv = std::min( wv, WVe[k-1] + v_kj +bp_penalty);
			wv = std::min( wv, temp[k-1] + v_kj + bp_penalty);
		}
		// How do I do the VP one?

		for ( auto it = CLVP[j].begin();CLVP[j].end()!=it; ++it ) {
			
			cand_pos_t k = it->first;
			energy_t vp_kj = it->second;
			if(k<bound_left) wv = std::min(wv, WIP[k-1]  + vp_kj );
		}
		for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>bound_right; ++it ) {
			
			cand_pos_t k = it->first;
			energy_t wmb_kj = it->second;
			
			wv = std::min( wv, WVe[k-1]  + wmb_kj + PSM_penalty + bp_penalty);
			wv = std::min( wv, temp[k-1]  + wmb_kj + PSM_penalty + bp_penalty);
		}
		if(tree.tree[j].pair<0) wv = std::min(wv, temp[j-1] + cp_penalty);
		temp[j] = wv;
		}
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
bool is_candidate(const std::vector< cand_list_td1 >& CL,const SparseMFEFold::Cand_comp& cand_comp,cand_pos_t i, cand_pos_t j) {
	const cand_list_td1 &list = CL[j];

	auto it = std::lower_bound(list.begin(),list.end(),i,cand_comp);

	return it!=list.end() && it->first==i;
}
bool is_candidate(const std::vector< cand_list_t >& CL,const SparseMFEFold::Cand_comp& cand_comp,cand_pos_t i, cand_pos_t j) {
	const cand_list_t &list = CL[j];

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
void find_mb_dangle(const energy_t WM2ij,const energy_t WM2i1j,const energy_t WM2ij1,const energy_t WM2i1j1, paramT* params, const short* S, const cand_pos_t i, const cand_pos_t j, cand_pos_t &k, cand_pos_t &l,const std::vector<Node> &tree){

	const pair_type tt = pair[S[j]][S[i]];
	const energy_t e1 = WM2ij +  E_MLstem(tt, -1, -1, params);
	const energy_t e2 = WM2i1j +  E_MLstem(tt, -1, S[i+1], params);
	const energy_t e3 = WM2ij1 +  E_MLstem(tt, S[j-1], -1, params);
	const energy_t e4 = WM2i1j1 +  E_MLstem(tt, S[j-1], S[i+1], params);
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
void trace_W(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, const std::vector<energy_t> & W, std::vector< energy_t > &WM, std::vector< energy_t > &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector< energy_t >&WIP_Bp, const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates, cand_pos_t i, cand_pos_t j, sparse_tree &tree) {
	// printf("W at %d and %d with %d\n",i,j,W[j]);

	if (i+TURN+1>j) return;
	// case j unpaired
	if (W[j] == W[j-1]) {
		trace_W(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,W,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,i,j-1,tree);
		return;
	}

	cand_pos_t m=j+1;
	energy_t w = INF;
	for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i;++it ) {
		m = it->first;
		const energy_t wmb_kj = it->second;
		w = W[m-1] + wmb_kj;
		if(W[j] == w + PS_penalty){
			trace_W(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,W,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,i,m-1,tree);
			trace_WMB(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,m,j,wmb_kj,tree);
			return;
		}

	}
	energy_t v=INF;
	w = INF;
    Dangle dangle =3;
	energy_t vk = INF;
	for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i;++it ) {
		m = it->first;
		const energy_t v_kj = it->fourth >> 2;
		const Dangle d = it->fourth & 3;
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
	trace_W(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,W,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,i,m-1,tree);
	trace_V(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,l,v,tree);
}

/**
* @brief Traceback from V entry
* 
* @param structure Final Structure
* @param mark_candidates Whether Candidates should be [ ]
* @param i row index
* @param j column index
*/
void trace_V(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector< energy_t >&WIP_Bp,const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree){
	// printf("V at %d and %d with %d\n",i,j,e);
	

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
		const cand_pos_t k=arrow.k(i);
		const cand_pos_t l=arrow.l(j);
		assert(i<k);
		assert(l<j);
		trace_V(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,l,arrow.target_energy(),tree);
		return;

	}
	else {

		// try to trace back to a candidate: (still) interior loop case
		cand_pos_t l_min = std::max(i,j-31);
		for ( cand_pos_t l=j-1; l>l_min; l--) {
			// Break if it's an assured dangle case
			for ( auto it=CL[l].begin(); CL[l].end()!=it && it->first>i; ++it ) {
				const cand_pos_t k=it->first;
				if(k-i > 31) continue;
				if (  e == it->second + E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S[k]][S[l]]],S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params)) ) {
					trace_V(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,l,it->second,tree);
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
	auto const temp = recompute_WM(WM,CL,CLWMB,params,i+1,j-1,tree.tree,tree.up);
	WM = temp;
	auto const temp2 = recompute_WM2(WM,WM2,CL,CLWMB,params,i+1,j-1,tree.tree,tree.up);
	WM2 = temp2;
	
	// Dangle for Multiloop
	cand_pos_t k = i+1;
	cand_pos_t l = j-1;
	if(params->model_details.dangles == 1){
		auto const temp3 = recompute_WM(WM,CL,CLWMB,params,i+2,j-1,tree.tree,tree.up);
		auto const temp4 = recompute_WM2(temp3,WM2,CL,CLWMB,params,i+2,j-1,tree.tree,tree.up);
		find_mb_dangle(temp2[j-1],temp4[j-1],temp2[j-2],temp4[j-2],params,S,i,j,k,l,tree.tree);
		if (k>i+1){
			WM = temp3;
			WM2 = temp4;
		}
	}
	
	
	trace_WM2(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,l,tree);
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
void trace_WM(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params,const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector< energy_t >&WIP_Bp,const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree){
	// printf("WM at %d and %d with %d\n",i,j,e);

	if (i+TURN+1>j) {return;}

	if ( e == WM[j-1] + params->MLbase ) {
		trace_WM(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,i,j-1,WM[j-1],tree);
		return;
	}
	cand_pos_t m = j+1;
	for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i;++it ) {
		m = it->first;
		const energy_t wmb_kj = it->second +PSM_penalty;
		energy_t wmb_up = static_cast<energy_t>((m-i)*params->MLbase) + wmb_kj;
		energy_t wmb_wm = WM[m-1] + wmb_kj;
		if(e == wmb_up + PSM_penalty){
			trace_WMB(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,m,j,wmb_kj,tree);
			return;
		}
		else if(e == wmb_wm + PSM_penalty){
			trace_WM(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,i,m-1,WM[m-1],tree);
			trace_WMB(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,m,j,wmb_kj,tree);
			return;
		}

	}

	energy_t v = INF;
    energy_t vk = INF;
    Dangle dangle = 3;
	for ( auto it=CL[j].begin();CL[j].end() != it && it->first>=i;++it ) {
		m = it->first;
		const energy_t v_kj = it->third >> 2;
		const Dangle d = it->third & 3;
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
			v = vk - E_MLstem(ptype,S[m],-1,params) - params->MLbase;
            break;
        case 2:
            l=j-1;
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,-1,S[j],params) - params->MLbase;
            break;
        case 3:
			if(params->model_details.dangles == 1){
				k=m+1;
				l=j-1;
				ptype= pair[S[k]][S[l]];
				v = vk - E_MLstem(ptype,S[m],S[j],params) - 2*params->MLbase;
			}
            break;
    }
    

    if ( e == WM[m-1] + vk ) {
		// no recomp, same i
		trace_WM(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,i,m-1,WM[m-1],tree);
		trace_V(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,l,v,tree);
		return;
	} else if ( e == static_cast<energy_t>((m-i)*params->MLbase) + vk ) {
		trace_V(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,l,v,tree);
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
void trace_WM2(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector< energy_t >&WIP_Bp,const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, sparse_tree &tree){
	// printf("WM2 at %d and %d with %d\n",i,j,WM2[j]);

	if (i+2*TURN+3>j) {return;}
	const energy_t e = WM2[j];

	// case j unpaired
	if ( e == WM2[j-1] + params->MLbase ) {
		// same i, no recomputation
		trace_WM2(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,i,j-1,tree);
		return;
	}

   	cand_pos_t m = j+1;
	for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i;++it ) {
		m = it->first;
		const energy_t wmb_kj = it->second + PSM_penalty;
		energy_t wmb_up = static_cast<energy_t>((m-i)*params->MLbase) + wmb_kj;
		energy_t wmb_wm = WM[m-1] + wmb_kj;
		if(e == wmb_up){
			trace_WMB(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,m,j,wmb_kj,tree);
			return;
		}
		else if(e == wmb_wm){
			trace_WM(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,i,m-1,WM[m-1],tree);
			trace_WMB(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,m,j,wmb_kj,tree);
			return;
		}

	}

    energy_t v = INF;
    energy_t vk = INF;
    Dangle dangle = 4;
	for ( auto it=CL[j].begin();CL[j].end() != it  && it->first>=i+TURN+1;++it ) {
		m = it->first;
		
		const energy_t v_kj = it->third >> 2;
		const Dangle d = it->third & 3;
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
			v = vk - E_MLstem(ptype,S[m],-1,params)- params->MLbase;
            break;
        case 2:
            l=j-1;
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,-1,S[j],params)- params->MLbase;
            break;
        case 3:
			if(params->model_details.dangles == 1){
				k=m+1;
				l=j-1;
				ptype= pair[S[k]][S[l]];
				v = vk - E_MLstem(ptype,S[m],S[j],params)- 2*params->MLbase;
			}
            break;
    }
    
	
    if ( e == WM[m-1] + vk ) {
		trace_WM(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,i,m-1,WM[m-1],tree);
		trace_V(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,l,v,tree);
		return;
	}
	assert(false);
}

/**
* @brief Traceback from WMB
* 
* @param WI_Bp WI array where all j values are a inner closing base minus 1
* @param CLBE BE candidate list
* @param i row index
* @param j column index
 */
void trace_WMB(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector< energy_t >&WIP_Bp,const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree){
	// printf("WMB at i is %d and j is %d and e is %d\n",i,j,e);
	assert( i+TURN+1<=j );
	assert( j<=n);

	std::vector<energy_t> WMBA;
	WMBA.resize(n+1,INF);
	std::vector<energy_t> WMBP = recompute_WMBP(CL,CLVP,CLWMB,CLBE,WI_Bbp,WMBA,n,i,j,tree);

	cand_pos_t bp_j = tree.tree[j].pair;

	if (tree.tree[j].pair >= 0 && j > tree.tree[j].pair && tree.tree[j].pair > i){
	for(int l = bp_j+1; l<j; ++l){
		// if (tree.tree[j].pair >= 0 && j > tree.tree[j].pair){
			cand_pos_t Bp_lj = tree.Bp(l,j);
			energy_t BE_energy = INF;
			if (Bp_lj >= 0){
				for ( auto it = CLBE[Bp_lj].begin();CLBE[Bp_lj].end()!=it; ++it ) {
					cand_pos_t k = it->first;
					if(k == bp_j){
						BE_energy = it->second;
						break;
					}
				}
				if(e == WMBP[l] + WI_Bp[l+1] + BE_energy + PB_penalty){
					std::vector<energy_t> WI_lp1;
					WI_lp1 = recompute_WI(CL,CLWMB,n,l+1,Bp_lj-1,tree.tree,tree.up);
					trace_BE(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,bp_j,tree.tree[Bp_lj].pair,BE_energy,tree);
					trace_WMBP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WMBP,WMBA,n,mark_candidates,i,l,WMBP[l],tree);
					trace_WI(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI_lp1,n,mark_candidates,l+1,Bp_lj-1,WI_Bp[l+1],tree);
					return;
				}
			}
		// }
	}
	}
	trace_WMBP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WMBP,WMBA,n,mark_candidates,i,j,WMBP[j],tree);
	return;
	

}

/**
* @brief Traceback from VP entry
* 
* @param structure Final Structure
* @param taVP VP trace arrows
* @param CLVP VP candidate list
* @param i row index
* @param j column index
*/
void trace_VP(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector< energy_t >&WIP_Bp,const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree){
	// printf("VP at %d and %d with %d\n",i,j,e);
	structure[i] = '[';
	structure[j] = ']';
	if(e==0) return;

	const pair_type ptype_closing = pair[S[i]][S[j]];

	
	cand_pos_t B_ij = tree.B(i,j);
	cand_pos_t Bp_ij = tree.Bp(i,j);
	cand_pos_t b_ij = tree.b(i,j);
	cand_pos_t bp_ij = tree.bp(i,j);
	if(tree.tree[i].parent->index > 0 && tree.tree[j].parent->index < tree.tree[i].parent->index && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0){
		std::vector<energy_t> WI_ip1;
		std::vector<energy_t> WI_jm1;
		WI_ip1 = recompute_WI(CL,CLWMB,n,i+1,Bp_ij-1,tree.tree,tree.up);
		WI_jm1 = recompute_WI(CL,CLWMB,n,B_ij+1,j-1,tree.tree,tree.up);
		if(e== WI_ip1[Bp_ij-1] + WI_jm1[j-1]){
			trace_WI(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI_ip1,n,mark_candidates,i+1,Bp_ij-1,WI_ip1[Bp_ij-1],tree);
			trace_WI(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI_jm1,n,mark_candidates,B_ij+1,j-1,WI_jm1[j-1],tree);
			return;
		}
	}
	if (tree.tree[i].parent->index < tree.tree[j].parent->index && tree.tree[j].parent->index > 0 && b_ij>= 0 && bp_ij >= 0 && Bp_ij < 0){
		
		std::vector<energy_t> WI_ip1;
		std::vector<energy_t> WI_jm1;
		WI_ip1 = recompute_WI(CL,CLWMB,n,i+1,b_ij-1,tree.tree,tree.up);
		WI_jm1 = recompute_WI(CL,CLWMB,n,bp_ij+1,j-1,tree.tree,tree.up);
		if(e== (WI_ip1[b_ij-1] + WI_jm1[j-1])){
			trace_WI(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI_ip1,n,mark_candidates,i+1,b_ij-1,WI_ip1[b_ij-1],tree);
			trace_WI(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI_jm1,n,mark_candidates,bp_ij+1,j-1,WI_jm1[j-1],tree);
			return;
		}
	}
	if(tree.tree[i].parent->index > 0 && tree.tree[j].parent->index > 0 && Bp_ij >= 0 && B_ij >= 0  && b_ij >= 0 && bp_ij>= 0){
		std::vector<energy_t> WI_ip1;
		std::vector<energy_t> WI_Bb;
		std::vector<energy_t> WI_jm1;
		WI_ip1 = recompute_WI(CL,CLWMB,n,i+1,Bp_ij-1,tree.tree,tree.up);
		WI_Bb = recompute_WI(CL,CLWMB,n,B_ij+1,b_ij-1,tree.tree,tree.up);
		WI_jm1 = recompute_WI(CL,CLWMB,n,bp_ij+1,j-1,tree.tree,tree.up);
		
		if(e== WI_ip1[Bp_ij-1] + WI_Bb[b_ij-1] + WI_jm1[j-1]){
			trace_WI(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI_ip1,n,mark_candidates,i+1,Bp_ij-1,WI_ip1[Bp_ij+1],tree);
			trace_WI(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI_Bb,n,mark_candidates,B_ij+1,b_ij-1,WI_Bb[b_ij-1],tree);
			trace_WI(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI_jm1,n,mark_candidates,bp_ij+1,j-1,WI_jm1[j-1],tree);
			return;
		}
	}
	if (exists_trace_arrow_from(taVP,i,j)) {
		
		const TraceArrow &arrow = trace_arrow_from(taVP,i,j);
		const size_t k=arrow.k(i);
		const size_t l=arrow.l(j);
		assert(i<k);
		assert(l<j);
		trace_VP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,l,arrow.target_energy(),tree);
		return;

	}
	else {

		// try to trace back to a candidate: (still) interior loop case
		cand_pos_t l_min = std::max(i, j-31);
		for ( cand_pos_t l=j-1; l>l_min; l--) {
			// Break if it's an assured dangle case
			for ( auto it=CLVP[l].begin(); CLVP[l].end()!=it && it->first>i; ++it ) {
				const cand_pos_t k=it->first;

				if(k-i > 31) continue;
				// energy_t temp = (int)round(((j-l==1 && k-i==1) ? e_stP_penalty : e_intP_penalty )*ILoopE(S,S1,params,ptype_closing,i,j,k,l));
				energy_t temp = lrint(((j-l==1 && k-i==1) ? e_stP_penalty : e_intP_penalty )*E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S[k]][S[l]]],S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params)));
				if (  e == it->second + temp ) {
					trace_VP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,l,it->second,tree);
					return;
				}
			}
		}
	}

// 	// If not other cases, must be WV multiloop
	const std::vector<energy_t> temp = 	recompute_WVe(CLVP,n,i+1,j-1,tree);
	const std::vector<energy_t> WIP = recompute_WIP(CL,CLWMB,n,i+1,j-1,tree.tree,tree.up);
	// std::cout << WIP[186] << std::endl;
	// std::cout << WIP[30] << std::endl;
	const std::vector<energy_t> temp2 = recompute_WV(temp,CL,CLWMB,CLVP,WIP,n,i+1,j-1,tree);

	trace_WV(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,temp,WIP,temp2,n,mark_candidates,i+1,j-1,temp2[j-1],tree);


}

/**
* @brief Traceback from WVe entry
* 
* @param CLVP VP candidate list
* @param i row index
* @param j column index
*/
void trace_WVe(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector< energy_t >&WIP_Bp,const std::vector<energy_t> &WI_Bp,const std::vector<energy_t> WVe, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree){
	// printf("WVe at %d and %d with %d\n",i,j,e);

	if (i+TURN+1>=j) return;
	if(WVe[j] == WVe[j-1] + cp_penalty)  {
		trace_WVe(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WVe,n,mark_candidates,i,j-1,WVe[j-1],tree);
		return;
	}
	cand_pos_t bound_left = j;
	if(tree.b(i,j)>0) bound_left = tree.b(i,j);
	if(tree.Bp(i,j)>0) bound_left = std::min(bound_left,tree.Bp(i,j));
	// std::cout << bound_left << "  " << CLVP[j].size() << std::endl;
	for ( auto it = CLVP[j].begin();CLVP[j].end()!=it && it->first>=i; ++it ) {
		cand_pos_t k = it->first;
		// std::cout << k << std::endl;
		if(k>bound_left) continue;
		if(e == static_cast<energy_t>(cp_penalty*(k-i)) + it->second){
			trace_VP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
			return;
		}
	
	}
}

/**
* @brief Traceback from WV entry
* 
* @param CL V candidate list
* @param CLVP VP candidate list
* @param CLWMB WMB candidate list
* @param WIP WIP values in region [i,j]
* @param i row index
* @param j column index
*/
void trace_WV(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector< energy_t >&WIP_Bp,const std::vector<energy_t> &WI_Bp,const std::vector<energy_t> &WVe,const std::vector<energy_t> &WIP, const std::vector<energy_t> &WV, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree){
	// printf("WV at %d and %d with %d\n",i,j,e);

	if(WV[j] == WV[j-1] + cp_penalty)  {
		trace_WV(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WVe,WIP,WV,n,mark_candidates,i,j-1,WV[j-1],tree);
		return;
	}
	cand_pos_t bound_left = std::min(tree.b(i,j),tree.Bp(i,j));
	cand_pos_t bound_right = std::min(tree.b(i,j),tree.Bp(i,j));


	for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>bound_right; ++it ) {
		cand_pos_t k = it->first;
		// energy_t wm_v = it->third >> 2;
		energy_t v = it->second;
		// Dangle d = it->third & 3;
		// cand_pos_t num = 0;
		// if(d == 1 || d==2) num = 1;
		// else if(d==3 && params->model_details.dangles == 1) num =2;
		// energy_t fix = num*cp_penalty - num*params->MLbase - b_penalty;

		if(e == WV[k-1] + v + bp_penalty){
			cand_pos_t m = k;
			cand_pos_t l = j;
			// pair_type ptype = 0;
			energy_t v = it->second;
			// switch(d){
			// 	case 0:
			// 		ptype= pair[S[m]][S[l]];
			// 		v = wm_v - E_MLstem(ptype,-1,-1,params);
			// 		break;
			// 	case 1:
			// 		m=k+1;
			// 		ptype= pair[S[m]][S[l]];
			// 		v = wm_v - E_MLstem(ptype,S[k],-1,params) - params->MLbase;
			// 		break;
			// 	case 2:
			// 		l=j-1;
			// 		ptype= pair[S[m]][S[l]];
			// 		v = wm_v - E_MLstem(ptype,-1,S[j],params) - params->MLbase;
			// 		break;
			// 	case 3:
			// 		if(params->model_details.dangles == 1){
			// 			m=k+1;
			// 			l=j-1;
			// 			ptype= pair[S[m]][S[l]];
			// 			v = wm_v - E_MLstem(ptype,S[k],S[j],params) - 2*params->MLbase;
			// 		}
			// 		break;
			// }
			
			trace_WV(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WVe,WIP,WV,n,mark_candidates,i,k-1,WV[k-1],tree);
			trace_V(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,m,l,v,tree);			
			return;
		}
		if(e == WVe[k-1] + v +bp_penalty){
			cand_pos_t m = k;
			cand_pos_t l = j;
			// pair_type ptype = 0;
			energy_t v = it->second;
			
			// switch(d){
			// 	case 0:
			// 		ptype= pair[S[m]][S[l]];
			// 		v = wm_v - E_MLstem(ptype,-1,-1,params);
			// 		break;
			// 	case 1:
			// 		m=k+1;
			// 		ptype= pair[S[m]][S[l]];
			// 		v = wm_v - E_MLstem(ptype,S[k],-1,params) - params->MLbase;
			// 		break;
			// 	case 2:
			// 		l=j-1;
			// 		ptype= pair[S[m]][S[l]];
			// 		v = wm_v - E_MLstem(ptype,-1,S[j],params) - params->MLbase;
			// 		break;
			// 	case 3:
			// 		if(params->model_details.dangles == 1){
			// 			m=k+1;
			// 			l=j-1;
			// 			ptype= pair[S[m]][S[l]];
			// 			v = wm_v - E_MLstem(ptype,S[k],S[j],params) - 2*params->MLbase;
			// 		}
			// 		break;
			// }
			trace_WVe(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WVe,n,mark_candidates,i,k-1,WVe[k-1],tree);
			trace_V(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,m,l,v,tree);			
			return;
		}
	}

	for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>bound_right; ++it ) {
		cand_pos_t k = it->first;
		if(e == WV[k-1] + it->second + PSM_penalty + bp_penalty){
			trace_WV(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WVe,WIP,WV,n,mark_candidates,i,k-1,WV[k-1],tree);
			trace_WMB(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
			return;
		}
		if(e == WVe[k-1] + it->second + PSM_penalty + bp_penalty){
			trace_WVe(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WVe,n,mark_candidates,i,k-1,WV[k-1],tree);
			trace_WMB(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
			return;
		}
	}
	for ( auto it = CLVP[j].begin();CLVP[j].end()!=it && it->first>=i; ++it ) {
		cand_pos_t k = it->first;
		if(k>bound_left) continue;
		if(e == WIP[k-1] + it->second){
			trace_WIP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WIP,n,mark_candidates,i,k-1,WIP[k-1],tree);
			trace_VP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
			return;
		}
	
	}

}

/**
* @brief Traceback from WI entry
* 
* @param CL V candidate list
* @param CLWMB WMB candidate list
* @param WI WI values in region [i,j]
* @param i row index
* @param j column index
*/
void trace_WI(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector< energy_t >&WIP_Bp,const std::vector<energy_t> &WI_Bp,const std::vector< energy_t >&WI, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree){
	// printf("WI at %d and %d with %d\n",i,j,e);

	if (i+TURN+1>=j) return;

	// How to do one base backwards?

	if(e == WI[j-1] + PUP_penalty){
		trace_WI(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI,n,mark_candidates,i,j-1,WI[j-1],tree);
		return;
	}

	for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
		cand_pos_t k = it->first;
		if(e == WI[k-1] + it->second + PPS_penalty){
			trace_WI(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI,n,mark_candidates,i,k-1,WI[k-1],tree);
			trace_V(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
			return;
		}
	}

	for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i ; ++it ) {
		cand_pos_t k = it->first;
		if(e == WI[k-1] + it->second + PPS_penalty + PSM_penalty){
			trace_WI(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI,n,mark_candidates,i,k-1,WI[k-1],tree);
			trace_WMB(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
			return;
		}
	}
	assert(false);
}

/**
* @brief Traceback from WIP entry
* 
* @param CL V candidate list
* @param CLVP VP candidate list
* @param WIP WIP values in region [i,j]
* @param i row index
* @param j column index
*/
void trace_WIP(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector< energy_t >&WIP_Bp,const std::vector<energy_t> &WI_Bp,const std::vector< energy_t >& WIP, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree){
	// printf("WIP at %d and %d with %d\n",i,j,e);

	if (i+TURN+1>=j) return;

// 	// How to do one base backwards?

	if(e == WIP[j-1] + cp_penalty){
		trace_WIP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WIP,n,mark_candidates,i,j-1,WIP[j-1],tree);
		return;
	}
	for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
		cand_pos_t k = it->first;
		if(e == WIP[k-1] + it->second + bp_penalty){
			trace_WIP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WIP,n,mark_candidates,i,k-1,WIP[k-1],tree);
			trace_V(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
			return;
		}
		if(e == static_cast<energy_t>((k-i)*cp_penalty) + it->second + bp_penalty){
			trace_V(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
			return;
		}
	}

	for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i ; ++it ) {
		cand_pos_t k = it->first;
		if(e == WIP[k-1] + it->second + bp_penalty +PSM_penalty){
			trace_WIP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WIP,n,mark_candidates,i,k-1,WIP_Bbp[k-1],tree);
			trace_WMB(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
			return;
		}
		if(e == static_cast<energy_t>((k-i)*cp_penalty) + it->second + bp_penalty +PSM_penalty){
			trace_WMB(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
			return;
		}
	}
}

/**
* @brief Traceback from WMBP entry
* 
* @param CL V candidate list
* @param CLVP VP candidate list
* @param CLWMB WMB candidate list
* @param CLBE BE candidate list
* @param WI_Bbp WI values in regions where the left index is a outer closing base or an inner opening base
* @param i row index
* @param j column index
*/
void trace_WMBP(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector< energy_t >&WIP_Bp,const std::vector<energy_t> &WI_Bp, std::vector<energy_t> &WMBP,std::vector<energy_t> &WMBA, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree){
	// printf("WMBP at %d and %d with %d\n",i,j,e);

	energy_t VP_ij = INF;
	
	cand_pos_t b_ij = tree.b(i,j);
	for ( auto it = CLVP[j].begin();CLVP[j].end()!=it && it->first>=i ; ++it ) {
		cand_pos_t k = it->first;
		if(k==i) VP_ij = it->second; //second?

		int ext_case = compute_exterior_cases(k,j,tree);
		if ((b_ij > 0 && k < b_ij) || ext_case == 0){
			cand_pos_t bp_ik = tree.bp(i,k);
			cand_pos_t Bp_kj = tree.Bp(k,j);
			if (bp_ik >= 0 && k>bp_ik && Bp_kj > 0 && k<Bp_kj){
				energy_t BE_energy = INF;
				cand_pos_t B_kj = tree.B(k,j);
				cand_pos_t bp_kj = tree.tree[Bp_kj].pair;
				cand_pos_t b_kj = (B_kj>0) ? tree.tree[B_kj].pair : -2;
				for ( auto it2 = CLBE[Bp_kj].begin();CLBE[Bp_kj].end()!=it2; ++it2 ) {
					cand_pos_t l = it2->first;
					if(l == b_kj){
						BE_energy = it2->second;
						break;
					}
				}
				if(e == WMBA[k-1] + it->second + 2*PB_penalty + BE_energy){
					trace_BE(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,b_kj,bp_kj,BE_energy,tree);
					trace_WMBA(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WMBP,WMBA,n,mark_candidates,i,k-1,WMBP[k-1],tree);
					trace_VP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
					return;
				}
			}
		}
	}

	for ( auto it = CLVP[j].begin();CLVP[j].end()!=it && it->first>=i ; ++it ) {
		cand_pos_t k = it->first;

		cand_pos_t bp_ik = tree.bp(i,k);
		if(bp_ik >= 0 && k+TURN <= j){
			cand_pos_t Bp_ik = tree.tree[bp_ik].pair;
			energy_t BE_energy = INF;
			for ( auto it2 = CLBE[Bp_ik].begin();CLBE[Bp_ik].end()!=it2; ++it2 ) {
				cand_pos_t l = it2->first;
				if(l == i){
					BE_energy = it2->second;
					break;
				}
			}

			if(e == WI_Bbp[k-1] + it->second + 2*PB_penalty +BE_energy){
				std::vector<energy_t> WI_ip1km1;
				WI_ip1km1 = recompute_WI(CL,CLWMB,n,bp_ik+1,k-1,tree.tree,tree.up);
				trace_BE(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,i,bp_ik,BE_energy,tree);
				trace_WI(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI_ip1km1,n,mark_candidates,bp_ik+1,k-1,WI_Bbp[k-1],tree);
				trace_VP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
				return;
			}
		}
	}

	if(e == VP_ij + PB_penalty){
		trace_VP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,i,j,VP_ij,tree);
		return;
	}
}

/**
* @brief Traceback from WMBA entry
* 
* @param CL V candidate list
* @param CLWMB WMB candidate list
* @param i row index
* @param j column index
*/
void trace_WMBA(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector< energy_t >&WIP_Bp,const std::vector<energy_t> &WI_Bp, std::vector<energy_t> &WMBP,std::vector<energy_t> &WMBA, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree){
	// printf("WMBA at %d and %d with %d\n",i,j,e);
	
	if(WMBA[j] == WMBA[j-1] + PUP_penalty){
		trace_WMBA(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WMBP,WMBA,n,mark_candidates,i,j-1,WMBP[j-1],tree);
		return;
	}

	if(WMBA[j] == WMBP[j]){
		trace_WMBP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WMBP,WMBA,n,mark_candidates,i,j,WMBP[j],tree);
		return;
	}

	for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
		cand_pos_t k = it->first;
		if(e == WMBA[k-1] + it->second + PPS_penalty){
			trace_WMBA(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WMBP,WMBA,n,mark_candidates,i,k-1,WMBP[k-1],tree);
			trace_V(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
			return;
		}
	}

	for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>=i; ++it ) {
		cand_pos_t k = it->first;
		if(e == WMBA[k-1] + it->second + PPS_penalty +PSM_penalty){
			trace_WMBA(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WMBP,WMBA,n,mark_candidates,i,k-1,WMBP[k-1],tree);
			trace_WMB(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,k,j,it->second,tree);
			return;
		}
	}
}

/**
* @brief Traceback from BE entry
* 
* @param CLBE BE candidate list
* @param WIP WIP values in region [i,j]
* @param i row index
* @param j column index
*/
void trace_BE(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &WM, std::vector<energy_t> &WM2,const std::vector<energy_t> &WI_Bbp,const std::vector< energy_t >&WIP_Bbp,const std::vector<energy_t> &WIP_Bp,const std::vector<energy_t> &WI_Bp, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t ip, energy_t e, sparse_tree &tree){
	// printf("BE at %d and %d with %d\n",i,ip,e);
	structure[i] = '(';
	cand_pos_t j = tree.tree[i].pair;
	structure[j] = ')';
	const pair_type ptype_closing_ij = pair[S[i]][S[j]];
	if(i == ip) return;
	
	energy_t BE_energy = INF;
	int jp = tree.tree[ip].pair;
	int lp = jp;

	for ( auto it = CLBE[jp].begin();CLBE[jp].end()!=it && it->first>i; ++it ) {
		lp = it->first;
		BE_energy = it->second;
	}
	int l = tree.tree[lp].pair;
	// i lp ip       jp  l   j

	
	
	if(e == lrint(e_stP_penalty*E_IntLoop(lp-i-1,j-l-1,ptype_closing_ij,rtype[pair[S[lp]][S[l]]],S1[i+1],S1[j-1],S1[lp-1],S1[l+1],const_cast<paramT *>(params))) + BE_energy){
		trace_BE(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,lp,ip,BE_energy,tree);
		return;
	}
	if(e == lrint(e_intP_penalty*E_IntLoop(lp-i-1,j-l-1,ptype_closing_ij,rtype[pair[S[lp]][S[l]]],S1[i+1],S1[j-1],S1[lp-1],S1[l+1],const_cast<paramT *>(params))) + BE_energy){
		trace_BE(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,lp,ip,BE_energy,tree);
		return;
	}

	if(e == WIP_Bbp[lp-1] + BE_energy + WIP_Bp[l+1] + ap_penalty + 2*bp_penalty){
		std::vector<energy_t> WI_ip1lpm1;
		std::vector<energy_t> WI_lp1jm1;
		WI_ip1lpm1 = recompute_WIP(CL,CLWMB,n,i+1,lp-1,tree.tree,tree.up);
		WI_lp1jm1 = recompute_WIP(CL,CLWMB,n,l+1,j-1,tree.tree,tree.up);
		trace_WIP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI_ip1lpm1,n,mark_candidates,i+1,lp-1,WIP_Bbp[lp-1],tree);
		trace_BE(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,lp,ip,BE_energy,tree);
		trace_WIP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI_lp1jm1,n,mark_candidates,l+1,j-1,WIP_Bp[j-1],tree);
		return;
	}
	if(e == cp_penalty*((lp-i+1)) + BE_energy + WIP_Bp[l+1] + ap_penalty + 2*bp_penalty){
		std::vector<energy_t> WI_lp1jm1;
		WI_lp1jm1 = recompute_WIP(CL,CLWMB,n,l+1,j-1,tree.tree,tree.up);
		trace_BE(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,lp,ip,BE_energy,tree);
		trace_WIP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI_lp1jm1,n,mark_candidates,l+1,j-1,WIP_Bbp[j-1],tree);
		return;
	}
	if(e == WIP_Bbp[lp-1] + BE_energy + cp_penalty*((j-l+1)) + ap_penalty + 2*bp_penalty){
		std::vector<energy_t> WI_ip1lpm1;
		WI_ip1lpm1 = recompute_WIP(CL,CLWMB,n,i+1,lp-1,tree.tree,tree.up);
		trace_WIP(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,WI_ip1lpm1,n,mark_candidates,i+1,lp-1,WIP_Bbp[lp-1],tree);
		trace_BE(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,lp,ip,BE_energy,tree);
		return;
	}
	
}
/**
* @brief Trace back
* pre: row 1 of matrix W is computed
* @return mfe structure (reference)
*/
const std::string & trace_back(const std::string& seq, const std::vector< cand_list_td1 >& CL, const std::vector< cand_list_t >& CLWMB,const std::vector< cand_list_t >& CLBE,const std::vector< cand_list_t >& CLVP, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1,TraceArrows &ta,TraceArrows &taVP, const std::vector<energy_t>& W, std::vector<energy_t> &WM, std::vector<energy_t> &WM2, std::vector<energy_t> &WI_Bbp, std::vector< energy_t >&WIP_Bbp,std::vector<energy_t> &WIP_Bp, std::vector<energy_t> &WI_Bp, const cand_pos_t& n, sparse_tree &tree, const bool& mark_candidates=false) {

	structure.resize(n+1,'.');

	/* Traceback */
	trace_W(seq,CL,CLWMB,CLBE,CLVP,cand_comp,structure,params,S,S1,ta,taVP,W,WM,WM2,WI_Bbp,WIP_Bbp,WIP_Bp,WI_Bp,n,mark_candidates,1,n,tree);
	structure = structure.substr(1,n);

	return structure;
}

/**
* @brief Register a candidate
* @param i start
* @param j end
* @param e energy of candidate "V(i,j)"
* @param wmij energy at WM(i,j)
* @param wij energy at W(i,j)
*/
 void register_candidate(std::vector<cand_list_td1> &CL, cand_pos_t const& i, cand_pos_t const& j, energy_t const& e,energy_t const& wmij,energy_t const& wij) {
	assert(i<=j+TURN+1);
	
	CL[j].emplace_back( cand_entry_td1(i,e,wmij,wij) );
}
/**
* @brief Register a candidate
* @param i start
* @param j end
* @param e energy of candidate "V(i,j)"
*/
void register_candidate(std::vector<cand_list_t> &CL, cand_pos_t const& i, cand_pos_t const& j, energy_t const& e) {
	assert(i<=j+TURN+1);
	
	CL[j].emplace_back(cand_entry_t(i,e));
}

energy_t compute_WVe(cand_pos_t i, cand_pos_t j, std::vector<energy_t> &WVe, std::vector< cand_list_t > &CLVP, sparse_tree &tree){
	
	energy_t wve = INF;
	cand_pos_t bound = std::min((cand_pos_tu) tree.Bp(i,j),(cand_pos_tu) tree.b(i,j));
	for ( auto it = CLVP[j].begin();CLVP[j].end()!=it; ++it ) {
		cand_pos_t k = it->first;
		bool can_pair = tree.up[k-1] >= (k-i);
		if(can_pair && k<bound) wve = std::min(wve, static_cast<energy_t>(cp_penalty*(k-i)) + it->second);	
	}
	if(tree.tree[j].pair<0) wve = std::min(wve,WVe[j-1] + cp_penalty);

	return wve;
}


/**
* @brief Computes the values for the WV and WVe matrix
* @param i start
* @param j end
* @param WV WV array for region [i,j]
* @param WVe WVe array for region [i,j]
* @param CL V candidate list
* @param CLVP VP candidate list
* @param CLWMB WMB candidate list
*/
energy_t compute_WV(cand_pos_t j, cand_pos_t bound_left, cand_pos_t bound_right, std::vector<energy_t> &WV,std::vector<energy_t> &WVe, std::vector<energy_t> &dwip1, std::vector< cand_list_td1 > &CL, std::vector< cand_list_t> &CLWMB, std::vector< cand_list_t > &CLVP, sparse_tree &tree){
	energy_t m1 = INF, m2 = INF, m3 = INF,m5 = INF, m6 = INF,wv = INF;
	
	for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>bound_right ; ++it ) {
		cand_pos_t k = it->first;
		// energy_t val = it->third >> 2;
		energy_t val = it->second;
		// Dangle d = it->third & 3;
		// cand_pos_t num = 0;
		// if(d == 1 || d==2) num = 1;
		// else if(d==3 && params->model_details.dangles == 1) num =2;
		// energy_t fix = num*cp_penalty - num*params->MLbase-b_penalty;
		// m1 = std::min(m1, WVe[k-1] + val + fix + bp_penalty);
		// m5 = std::min(m5, WV[k-1] + val+ fix + bp_penalty);
		m1 = std::min(m1, WVe[k-1] + val + bp_penalty);
		m5 = std::min(m5, WV[k-1] + val + bp_penalty);
	
	}
	for ( auto it = CLWMB[j].begin();CLWMB[j].end()!=it && it->first>bound_right ; ++it ) {
		cand_pos_t k = it->first;
		energy_t val = it->second;
		m2 = std::min(m1, WVe[k-1] + val + PSM_penalty + bp_penalty);
		m6 = std::min(m6, WV[k-1] + val + PSM_penalty + bp_penalty);
	}

	
	for ( auto it = CLVP[j].begin();CLVP[j].end()!=it; ++it ) {
		cand_pos_t k = it->first;
		energy_t val = it->second;
		if(k<bound_left) m3 = std::min(m1, dwip1[k-1] + val);		
	}
	wv = std::min({m1,m2,m3,m5,m6});


	if(tree.tree[j].pair<0) wv = std::min(wv,WV[j-1] + cp_penalty);

	return wv;
}

/**
* @brief computes entries for BE
* 
* @param CLBE BE candidate list
* @param WIP WIP values in region [i,j]
* @param i row index
* @param j column index
*/
energy_t compute_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, std::vector< cand_list_t> &CLBE, sparse_tree &sparse_tree,short *S,short *S1,paramT* params, const std::vector<energy_t> dwip1, const std::vector< energy_t >& WIP_Bp){
// int compute_BE(int i, int j, int ip, int jp, int lp, int l, auto &S, auto &S1, auto &WIP_Bbp){
    // We are checking for the closest pair that we have already calculated to i/ip from j/jp 
    // If there is nothing, then i is the closest encompassing pair to jp
    // If it is not, then we get the energy for everything from jp to lp so that we calculate less

	// (.....(..(....)..).....)
	// i     lp      j  l     ip
	energy_t BE_energy = INF;
	cand_pos_t lp = jp;
	if (!CLBE[j].empty()){
		auto const [k,vbe] = CLBE[j].back();
		BE_energy = vbe;
		lp = k;
	}
	cand_pos_t l = sparse_tree.tree[lp].pair; // right closing base for lp
	
	const pair_type ptype_closing_iip = pair[S[i]][S[ip]];
    
    energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF,val=INF;
    if(i+1 == lp && ip-1 == l){
		m1 = lrint(e_stP_penalty*E_IntLoop(lp-i-1,ip-l-1,ptype_closing_iip,rtype[pair[S[lp]][S[l]]],S1[i+1],S1[ip-1],S1[lp-1],S1[l+1],const_cast<paramT *>(params))) + BE_energy;

        val = std::min(val,m1);
    }

    bool empty_region_ilp = (sparse_tree.up[lp-1] >= lp-i-1); //empty between i+1 and lp-1
    bool empty_region_lip = (sparse_tree.up[ip-1] >= ip-l-1); // empty between l+1 and ip-1
    bool weakly_closed_ilp = sparse_tree.weakly_closed(i+1,lp-1); // weakly closed between i+1 and lp-1
	bool weakly_closed_lip = sparse_tree.weakly_closed(l+1,ip-1); // weakly closed between l+1 and ip-1
        
    if (empty_region_ilp && empty_region_lip){
		m2 = lrint(e_intP_penalty*E_IntLoop(lp-i-1,ip-l-1,ptype_closing_iip,rtype[pair[S[lp]][S[l]]],S1[i+1],S1[ip-1],S1[lp-1],S1[l+1],const_cast<paramT *>(params))) + BE_energy;
        val = std::min(val,m2);
    }

        // 3
    if (weakly_closed_ilp && weakly_closed_lip){
        m3 = dwip1[lp-1] + BE_energy + WIP_Bp[l+1] + ap_penalty + 2*bp_penalty;
        val = std::min(val,m3);
    }

        // 4
    if (weakly_closed_ilp && empty_region_lip){
        m4 = dwip1[lp-1] + BE_energy + cp_penalty * (ip-l+1) + ap_penalty + 2*bp_penalty;
        val = std::min(val,m4);
    }

        // 5
    if (empty_region_ilp && weakly_closed_lip){
		
        m5 = ap_penalty + 2*bp_penalty + (cp_penalty * (lp-i+1)) + BE_energy + WIP_Bp[l+1];
        val = std::min(val,m5);
    }


    return(val);
        
}
/**
* @brief Computes entries for WMBP
* 
* @param CLVP VP candidate list
* @param CLBE BE candidate list
* @param WI WI values in region [i,j]
* @param i row index
* @param j column index
*/
energy_t compute_WMBP(cand_pos_t i, cand_pos_t j, sparse_tree &sparse_tree, std::vector<cand_list_t> &CLVP, std::vector<cand_list_t> &CLBE, LocARNA::Matrix<energy_t> &VP, std::vector<energy_t> &WI_Bbp, std::vector<energy_t> &WMBA){
	energy_t m1 = INF, m2 = INF, m3 = INF, wmbp = INF;
	// 1) WMBP(i,j) = BE(bpg(Bp(l,j)),Bp(l,j),bpg(B(l,j)),B(l,j)) + WMBP(i,l) + VP(l+1,j)
	if (sparse_tree.tree[j].pair < 0){
		energy_t tmp = INF;
		cand_pos_t b_ij = sparse_tree.b(i,j);
		for ( auto it = CLVP[j].begin();CLVP[j].end()!=it; ++it ) {
			cand_pos_t l = it->first;
			// Mateo Jan 2025 Added exterior cases to consider when looking at band borders. Solved case of [.(.].[.).]
			int ext_case = compute_exterior_cases(l,j,sparse_tree);
			if ((b_ij > 0 && l < b_ij) || ext_case == 0){
				cand_pos_t bp_il = sparse_tree.bp(i,l);
				cand_pos_t Bp_lj = sparse_tree.Bp(l,j);
				if (bp_il >= 0 && l>bp_il && Bp_lj > 0 && l<Bp_lj){//if(sparse_tree.b(i,j)>=0 && l <sparse_tree.b(i,j)){//
					cand_pos_t B_lj = sparse_tree.B(l,j);
					if (i <= sparse_tree.tree[l].parent->index && sparse_tree.tree[l].parent->index < j && l+3 <=j){
						energy_t BE_energy = INF;
						cand_pos_t b_lj = sparse_tree.tree[B_lj].pair;
						for ( auto it2 = CLBE[Bp_lj].begin();CLBE[Bp_lj].end()!=it2; ++it2 ) {
							cand_pos_t k = it2->first;
							if(k == b_lj){
								BE_energy = it2->second;
								break;
							}
						}
						energy_t WMBA_energy = WMBA[l-1];
						energy_t VP_energy = it->second;
						energy_t sum = BE_energy + WMBA_energy + VP_energy;

						tmp = std::min(tmp,sum);
					}
				}
			}
			
			m1 = 2*PB_penalty + tmp;
		}
	}
	

	// 2) WMBP(i,j) = VP(i,j) + P_b
	int i_mod = i % (MAXLOOP+1);
	m2 = VP(i_mod,j) + PB_penalty;

	// check later if <0 or <-1

	// WMBP(i,j) = BE(i,,,) _ WI(bp(i,k),k-1) + VP(k,j)
	if(sparse_tree.tree[j].pair < 0 && sparse_tree.tree[i].pair >= 0){
		energy_t tmp = INF;
		for ( auto it = CLVP[j].begin();CLVP[j].end()!=it; ++it ) {
			cand_pos_t l = it->first;
			cand_pos_t bp_il = sparse_tree.bp(i,l);
			if(bp_il >= 0 && l+TURN <= j){
				energy_t BE_energy = INF;
				cand_pos_t Bp_il = sparse_tree.tree[bp_il].pair;
				if (!CLBE[Bp_il].empty()){
					auto const [k,vbe] = CLBE[Bp_il].back();
					if(i==k) BE_energy = vbe;
				}
			
				energy_t WI_energy = (l-1-(bp_il+1))> 4 ? WI_Bbp[l-1] : PUP_penalty*(l-1-(bp_il+1)+1);
				energy_t VP_energy = it->second;
				energy_t sum = BE_energy + WI_energy + VP_energy;

				tmp = std::min(tmp,sum);
			}
		}
		
		m3 = 2*PB_penalty + tmp;
	}
	// get the min for WMB
	wmbp = std::min({m1,m2,m3});

	return(wmbp);

}


energy_t compute_WMBA(cand_pos_t j, sparse_tree &sparse_tree, std::vector<energy_t> &WMBP, std::vector<energy_t> &WMBA,std::vector<cand_list_td1> &CL, std::vector<cand_list_t> &CLWMB){

	// WMBA criteria
	energy_t wmba = INF;
	if(sparse_tree.tree[j].parent->index>0){

		for ( auto it = CL[j].begin(); CL[j].end()!=it; ++it ) {
			cand_pos_t k = it->first;
			if (sparse_tree.tree[k].pair < 0 && sparse_tree.tree[k].parent->index > -1 && sparse_tree.tree[j].parent->index > -1 && sparse_tree.tree[j].parent->index == sparse_tree.tree[k].parent->index){
				wmba = std::min(wmba,WMBA[k-1] + it->second + PPS_penalty);
			}
		}

		for ( auto it = CLWMB[j].begin(); CLWMB[j].end()!=it; ++it ) {
			cand_pos_t k = it->first;
			if (sparse_tree.tree[k].pair < 0 && sparse_tree.tree[k].parent->index > -1 && sparse_tree.tree[j].parent->index > -1 && sparse_tree.tree[j].parent->index == sparse_tree.tree[k].parent->index){
			wmba = std::min(wmba,WMBA[k-1] + it->second + PPS_penalty + PSM_penalty);

			}
		}
		if(sparse_tree.tree[j].pair < 0) wmba = std::min(wmba,WMBA[j-1] + PUP_penalty);
	}
	else{
		wmba = INF;
	}
	wmba = std::min(wmba,WMBP[j]);
	return wmba;

}
/**
* @brief Computes entries for WMB
* 
* @param CLVP VP candidate list
* @param CLBE BE candidate list
* @param WI WI values in region [i,j]
* @param i row index
* @param j column index
*/
energy_t compute_WMB(cand_pos_t i, cand_pos_t j, sparse_tree &sparse_tree, std::vector<cand_list_t> &CLBE, std::vector<energy_t> &WMBP, std::vector<energy_t> &WMBA, std::vector<energy_t> &WI_Bp){
	energy_t m1 = INF, m2 = INF, wmb = INF;

	// 2)
	if (sparse_tree.tree[j].pair >= 0 && j > sparse_tree.tree[j].pair && sparse_tree.tree[j].pair > i){
		cand_pos_t bp_j = sparse_tree.tree[j].pair;
		
		for (int l = (bp_j +1);l<j; ++l){
			cand_pos_t Bp_lj = sparse_tree.Bp(l,j);
			if (Bp_lj > 0){
				
				energy_t BE_energy = INF;
				for ( auto it = CLBE[Bp_lj].begin();CLBE[Bp_lj].end()!=it; ++it ) {
					cand_pos_t k = it->first;
					if(k == bp_j){
						BE_energy = it->second;
						break;
					}
				}

				energy_t WMBP_energy = WMBP[l];
				energy_t WI_energy = (Bp_lj-1-(l+1)+1)> 4 ? WI_Bp[l+1] : PUP_penalty*(Bp_lj-1-(l+1)+1);
				energy_t sum = BE_energy + WMBP_energy + WI_energy;

				m1 = std::min(m1,sum);

			}
		}
		
		m1 = PB_penalty + m1;
		
	}
	// check the WMBP_ij value
	m2 =  WMBA[j];

	wmb = std::min(m1,m2);
	return wmb;
}

/**
* @brief Computes VP internal loop value for case 5 of VP
* 
* @param VP VP array
* @param i row index
* @param j column index
*/
energy_t compute_VP_internal(cand_pos_t i, cand_pos_t j, cand_pos_t b_ij, cand_pos_t bp_ij, cand_pos_t Bp_ij, cand_pos_t B_ij,cand_pos_t &best_k,cand_pos_t &best_l, energy_t &best_e, sparse_tree& sparse_tree, short* S, short* S1,LocARNA::Matrix<energy_t> &VP, paramT* params){
	
	energy_t m5 = INF;

	cand_pos_t min_borders = std::min((cand_pos_tu) Bp_ij, (cand_pos_tu) b_ij);
	cand_pos_t edge_i = std::min(i+MAXLOOP+1,j-TURN-1);
	min_borders = std::min({min_borders,edge_i});
	const pair_type ptype_closing = pair[S[i]][S[j]];
		for ( cand_pos_t k=i+1; k<=min_borders; k++) {
			cand_pos_t k_mod= k%(MAXLOOP+1);

			energy_t cank = ((sparse_tree.up[k-1]>=(k-i-1))-1);
			cand_pos_t max_borders = std::max(bp_ij,B_ij)+1;
			cand_pos_t edge_j = k+j-i-MAXLOOP-2;
			max_borders = std::max({max_borders,edge_j});
			
			for (cand_pos_t l=j-1; l>=max_borders; --l) {
				assert(k-i+j-l-2<=MAXLOOP);

				energy_t canl = (((sparse_tree.up[j-1]>=(j-l-1))-1) | cank);
				energy_t v_iloop_kl = INF & canl;
				v_iloop_kl = v_iloop_kl + VP(k_mod,l) + lrint(e_intP_penalty*E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S[k]][S[l]]],S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params)));

				if ( v_iloop_kl < m5) {
					m5 = v_iloop_kl;
					best_l=l;
					best_k=k;
					best_e=VP(k_mod,l);
				}
				
			}
			
			
		}
	return m5;
}

/**
* @brief Computes the energy values for VP
* 
* @param CLVP VP candidate list
* @param taVP VP trace arrows
* @param WI WI values in region [i,j]
* @param i row index
* @param j column index
*/
energy_t compute_VP(cand_pos_t i, cand_pos_t j, cand_pos_t b_ij, cand_pos_t bp_ij, cand_pos_t Bp_ij, cand_pos_t B_ij, sparse_tree& sparse_tree, std::vector<energy_t> &dwi1, std::vector<energy_t> &dwvp, std::vector<energy_t> &WI_Bbp, short* S, short* S1,LocARNA::Matrix<energy_t> &VP, std::vector<cand_list_t> &CLVP, TraceArrows &taVP, const SparseMFEFold::Cand_comp& cand_comp , paramT* params){
	energy_t m1 = INF, m2 = INF, m3 = INF, m4= INF, m5 = INF, m6 = INF, vp = INF;
	const pair_type ptype_closing = pair[S[i]][S[j]];
	if(sparse_tree.tree[i].parent->index > 0 && sparse_tree.tree[j].parent->index < sparse_tree.tree[i].parent->index && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0){
		energy_t WI_ipus1_BPminus = dwi1[Bp_ij-1];
		energy_t WI_Bplus_jminus = (j-1-(B_ij+1))> 4 ? WI_Bbp[j-1] : PUP_penalty*(j-1-(B_ij+1)+1);
		
		m1 = WI_ipus1_BPminus + WI_Bplus_jminus;
	}
	if (sparse_tree.tree[i].parent->index < sparse_tree.tree[j].parent->index && sparse_tree.tree[j].parent->index > 0 && b_ij>= 0 && bp_ij >= 0 && Bp_ij < 0){
		energy_t WI_i_plus_b_minus = dwi1[b_ij-1];
		energy_t WI_bp_plus_j_minus = (j-1-(bp_ij+1))> 4 ? WI_Bbp[j-1] : PUP_penalty*(j-1-(bp_ij+1)+1);
		
		m2 = WI_i_plus_b_minus + WI_bp_plus_j_minus;
	}
	
	if(sparse_tree.tree[i].parent->index > 0 && sparse_tree.tree[j].parent->index > 0 && Bp_ij >= 0 && B_ij >= 0  && b_ij >= 0 && bp_ij>= 0){						
		energy_t WI_i_plus_Bp_minus = dwi1[Bp_ij-1];
		energy_t WI_B_plus_b_minus = (b_ij-1-(B_ij+1))> 4 ? WI_Bbp[b_ij-1] : PUP_penalty*(b_ij-1-(B_ij+1)+1);
		energy_t WI_bp_plus_j_minus = (j-1-(bp_ij+1))> 4 ? WI_Bbp[j-1] : PUP_penalty*(j-1-(bp_ij+1)+1);
		
		m3 = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
	}
	if(sparse_tree.tree[i+1].pair < -1 && sparse_tree.tree[j-1].pair < -1){
		cand_pos_t ip1_mod = (i+1)%(MAXLOOP+1);
		
		m4 = lrint(e_stP_penalty*ILoopE(S,S1,params,ptype_closing,i,j,i+1,j-1)) + VP(ip1_mod,j-1);
	}

	
	cand_pos_t best_k = 0;
	cand_pos_t best_l = 0;
	energy_t best_e = 0;

	m5 = compute_VP_internal(i,j,b_ij,bp_ij,Bp_ij,B_ij,best_k,best_l,best_e,sparse_tree,S,S1,VP,params);
	
	// case 6 and 7
	m6 = dwvp[j-1] + ap_penalty + 2*bp_penalty;

	energy_t vp_h = std::min({m1,m2,m3});
	energy_t vp_iloop = std::min(m4,m5);
	if(m4<m5){
		best_k = i+1;
		best_l = j-1;
		cand_pos_t ip1_mod = (i+1)%(MAXLOOP+1);
		best_e = VP(ip1_mod,j-1);
	}
	energy_t vp_split = m6;
	vp = std::min({vp_h,vp_iloop,vp_split});

	if ( vp_iloop < std::min(vp_h,vp_split) ) {
		if ( is_candidate(CLVP,cand_comp,best_k,best_l) ) {
			//std::cout << "Avoid TA "<<best_k<<" "<<best_l<<std::endl;
			avoid_trace_arrow(taVP);
		} else {
			//std::cout<<"Reg TA "<<i<<","<<j<<":"<<best_k<<","<<best_l<<std::endl;
			
			register_trace_arrow(taVP,i,j,best_k,best_l,best_e);
		}
	}
	
	return vp;

}

/**
* @brief Computes the internal loop value for V
* 
* @param V V array
* @param i row index
* @param j column index
*/
energy_t compute_internal(cand_pos_t i, cand_pos_t j,cand_pos_t &best_k, cand_pos_t &best_l, energy_t &best_e, sparse_tree &sparse_tree, short* S, short* S1,LocARNA::Matrix<energy_t> &V, paramT *params){
	energy_t v_iloop = INF;
	cand_pos_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
	const pair_type ptype_closing = pair[S[i]][S[j]];
	for ( cand_pos_t k=i+1; k<=max_k; k++) {
		cand_pos_t k_mod=k%(MAXLOOP+1);
		
		energy_t cank = ((sparse_tree.up[k-1]>=(k-i-1))-1);
		cand_pos_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;
		for (cand_pos_t l=j-1; l>=min_l; --l) {
			assert(k-i+j-l-2<=MAXLOOP);
			energy_t canl = (((sparse_tree.up[j-1]>=(j-l-1))-1) | cank);
			energy_t v_iloop_kl = INF & canl;
			
			v_iloop_kl = v_iloop_kl + V(k_mod,l) + E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S[k]][S[l]]],S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params));
			if ( v_iloop_kl < v_iloop) {
				v_iloop = v_iloop_kl;
				best_l=l;
				best_k=k;
				best_e=V(k_mod,l);
			}
		}
	}
	return v_iloop;
}

/**
 * @brief Determines the MFE energy for a given sequence
*/
energy_t fold(const std::string& seq, sparse_tree &sparse_tree, LocARNA::Matrix<energy_t> &V, const SparseMFEFold::Cand_comp& cand_comp, std::vector<cand_list_td1> &CL, std::vector<cand_list_t> &CLWMB,std::vector<cand_list_t> &CLVP, std::vector<cand_list_t> &CLBE, short* S, short* S1, paramT* params, TraceArrows &ta,TraceArrows &taVP, std::vector<energy_t> &W, std::vector<energy_t> &WM, std::vector<energy_t> &WM2, std::vector<energy_t> &dmli1, std::vector<energy_t> &dmli2, LocARNA::Matrix<energy_t> &VP, std::vector<energy_t> &WVe, std::vector<energy_t> &WV, std::vector<energy_t> &dwvp, std::vector<energy_t> &WMB, std::vector<energy_t> &dwmbi,std::vector<energy_t> &WMBP,std::vector<energy_t> &WMBA,std::vector<energy_t> &WI,std::vector<energy_t> &dwi1,std::vector<energy_t> &WIP, std::vector<energy_t> &dwip1, std::vector<energy_t> &WI_Bbp, std::vector<energy_t> &WIP_Bbp, std::vector< energy_t >&WIP_Bp, std::vector<energy_t> &WI_Bp , const cand_pos_t n, const bool garbage_collect) {
	Dangle d = 3;
	// int BE_avoided = 0;
    if(params->model_details.dangles == 0 || params->model_details.dangles == 1) d = 0;
    
    for (cand_pos_t i=n; i>0; --i) {
		if(pseudoknot){
			for(cand_pos_t j=i;j<=n && sparse_tree.tree[j].pair<0;++j){
				WI[j] = (j-i+1)*PUP_penalty;
			}
		}

		for ( cand_pos_t j=i+TURN+1; j<=n; j++ ) {

			bool evaluate = sparse_tree.weakly_closed(i,j);
			// ------------------------------
			// W: split case
			bool pairedkj = 0;
			energy_t w_split = INF;
			energy_t wm_split = INF;
			energy_t wm2_split = INF;
			energy_t wi_split = INF;
			energy_t wip_split = INF;
			for ( auto it=CL[j].begin();CL[j].end() != it;++it ) {
				cand_pos_t k=it->first;

				const energy_t v_kj = it->third >> 2;
				const energy_t v_kjw = it ->fourth >> 2;
				bool can_pair = sparse_tree.up[k-1] >= (k-i);
				// WM Portion
				wm_split = std::min( wm_split, WM[k-1] + v_kj );
				if(can_pair) wm_split = std::min( wm_split,static_cast<energy_t>((k-i)*params->MLbase) + v_kj );
				// WM2 Portion
				wm2_split = std::min( wm2_split, WM[k-1] + v_kj );
				// W Portion
				w_split = std::min( w_split, W[k-1] + v_kjw );

				//WI portion
				energy_t v_kjj = it->second + PPS_penalty;
				wi_split = std::min(wi_split,WI[k-1] + v_kjj);
				//WIP portion
				v_kjj = it->second + bp_penalty;
				wip_split = std::min(wip_split,WIP[k-1]+v_kjj);
				if(can_pair) wip_split = std::min(wip_split,static_cast<energy_t>((k-i)*cp_penalty) +v_kjj);
				
			}

			if(sparse_tree.weakly_closed(i,j)){
				for ( auto it=CLWMB[j].begin();CLWMB[j].end() != it;++it ) {
					
					if(pairedkj) break;   // Not needed I believe as there shouldn't be any candidates there if paired anyways
					// Maybe this would just avoid this loop however

					cand_pos_t k = it->first;
					bool can_pair = sparse_tree.up[k-1] >= (k-i);

					// For W
					energy_t wmb_kj = it->second + PS_penalty;
					w_split = std::min( w_split, W[k-1] + wmb_kj );
					// For WM -> I believe this would add a PSM penalty for every pseudoknot which would be bad
					wmb_kj = it->second + PSM_penalty + b_penalty;
					wm_split = std::min(wm_split, WM[k-1] + wmb_kj);
					if(can_pair) wm_split = std::min(wm_split,static_cast<energy_t>((k-i)*params->MLbase) + wmb_kj);
					wm2_split = std::min( wm2_split, WM[k-1] + wmb_kj );
					if(can_pair) wm2_split = std::min( wm2_split, static_cast<energy_t>((k-i)*params->MLbase) + wmb_kj );
					// For WI
					wmb_kj = it->second + PSM_penalty + PPS_penalty;
					wi_split = std::min(wi_split,WI[k-1] + wmb_kj);
					// if(i==50 && j==279) printf("i is %d and k is %d and j is %d and wi is %d and wmbkj is %d\n",i,k,j,WI[k-1],wmb_kj);

					// For WIP
					wmb_kj = it->second + PSM_penalty + bp_penalty;
					wip_split = std::min(wip_split,WIP[k-1]+wmb_kj);
					if(can_pair) wip_split = std::min(wip_split,static_cast<energy_t>((k-i)*cp_penalty) +wmb_kj);
				}
			}

			if(sparse_tree.tree[j].pair<0) w_split = std::min(w_split,W[j-1]);
			if(sparse_tree.tree[j].pair<0) wm2_split = std::min( wm2_split, WM2[j-1] + params->MLbase );
			if(sparse_tree.tree[j].pair<0) wm_split = std::min( wm_split, WM[j-1] + params->MLbase );
			if(sparse_tree.tree[j].pair<0) wi_split = std::min(wi_split,WI[j-1] + PUP_penalty);
			if(sparse_tree.tree[j].pair<0) wip_split = std::min(wip_split,WIP[j-1] + cp_penalty);

			
			energy_t w  = w_split; // entry of W w/o contribution of V
			energy_t wm = wm_split; // entry of WM w/o contribution of V

			size_t i_mod=i%(MAXLOOP+1);

			const pair_type ptype_closing = pair[S[i]][S[j]];
			const bool restricted = sparse_tree.tree[i].pair == -1 || sparse_tree.tree[j].pair == -1;

			const bool unpaired = (sparse_tree.tree[i].pair<-1 && sparse_tree.tree[j].pair<-1);
			const bool paired = (sparse_tree.tree[i].pair == j && sparse_tree.tree[j].pair == i);
			const bool pkonly = (!pk_only || paired);
			energy_t v = INF;
			// ----------------------------------------
			// cases with base pair (i,j)
			// if(ptype_closing>0 && !restricted && evaluate) { // if i,j form a canonical base pair
			if(ptype_closing>0 && !restricted && evaluate && pkonly) { 
				bool canH = (paired || unpaired);
				if(sparse_tree.up[j-1]<(j-i-1)) canH=false;
				
				energy_t v_h = canH ? HairpinE(seq,S,S1,params,i,j) : INF;
				// info of best interior loop decomposition (if better than hairpin)
				cand_pos_t best_l=0;
				cand_pos_t best_k=0;
				energy_t best_e = INF;

				energy_t v_iloop=INF;

				// constraints for interior loops
				// i<k; l<j
				// k-i+j-l-2<=MAXLOOP  ==> k <= MAXLOOP+i+1
				//            ==> l >= k+j-i-MAXLOOP-2
				// l-k>=TURN+1         ==> k <= j-TURN-2
				//            ==> l >= k+TURN+1
				// j-i>=TURN+3
				//
				if((sparse_tree.tree[i].pair<-1 && sparse_tree.tree[j].pair < -1) || sparse_tree.tree[i].pair == j) {
					v_iloop = compute_internal(i,j,best_k,best_l,best_e,sparse_tree,S,S1,V,params); 
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
			energy_t w_v  = E_ext_Stem(v,vi1j,vij1,vi1j1,S,params,i,j,d,n,sparse_tree.tree);
			// Checking the dangle positions for W
			const energy_t wm_v = E_MLStem(v,vi1j,vij1,vi1j1,S,params,i,j,d,n,sparse_tree.tree);

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
			energy_t wi_v = INF;
			energy_t wip_v = INF;
			energy_t wi_wmb = INF;
			energy_t wip_wmb = INF;
			energy_t w_wmb = INF,wm_wmb = INF;
			if(pseudoknot){
				cand_pos_t Bp_ij = sparse_tree.Bp(i,j);
				cand_pos_t B_ij = sparse_tree.B(i,j);
				cand_pos_t b_ij = sparse_tree.b(i,j);
				cand_pos_t bp_ij = sparse_tree.bp(i,j);

				// Start of VP ---- Will have to change the bounds to 1 to n instead of 0 to n-1
				bool weakly_closed_ij = sparse_tree.weakly_closed(i,j);
				if ( weakly_closed_ij || sparse_tree.tree[i].pair >= -1 || sparse_tree.tree[j].pair >= -1 || ptype_closing == 0)	{
				
					VP(i_mod,j) = INF;
					
				}
				else{
					const energy_t vp = compute_VP(i,j,b_ij,bp_ij,Bp_ij,B_ij,sparse_tree,dwi1,dwvp,WI_Bbp,S,S1,VP,CLVP,taVP,cand_comp,params);

					VP(i_mod,j) = vp;
				}

				// -------------------------------------------End of VP----------------------------------------------------------------

				// Start of WMBP
				if ((sparse_tree.tree[i].pair >= -1 && sparse_tree.tree[i].pair > j) || (sparse_tree.tree[j].pair >= -1 && sparse_tree.tree[j].pair < i) || (sparse_tree.tree[i].pair >= -1 && sparse_tree.tree[i].pair < i ) || (sparse_tree.tree[j].pair >= -1 && j < sparse_tree.tree[j].pair)){
					WMB[j] = INF;
					WMBP[j] = INF;
					WMBA[j] = INF;

				}
				else{
					
					const energy_t wmbp = compute_WMBP(i,j,sparse_tree,CLVP,CLBE,VP,WI_Bbp,WMBA);
					WMBP[j] = wmbp;

					const energy_t wmba = compute_WMBA(j,sparse_tree,WMBP,WMBA,CL,CLWMB);
					WMBA[j] = wmba;

					const energy_t wmb = compute_WMB(i,j,sparse_tree,CLBE,WMBP,WMBA,WI_Bp);
					WMB[j] = wmb;

					// const energy_t wmb_vp = VP(i_mod,j) + PB_penalty;
				}

				// -------------------------------------------------------End of WMB------------------------------------------------------
				

				// Start of WI -- the conditions on calculating WI is the same as WIP, so we combine them
		
				if (!weakly_closed_ij){
					WI[j] = INF;
					WIP[j] = INF;
				}
				else{
					
					wi_v = V(i_mod,j) + PPS_penalty;
					wip_v = V(i_mod,j)	+ bp_penalty;
					
					wi_wmb = WMB[j] + PSM_penalty + PPS_penalty;
					wip_wmb = WMB[j] + PSM_penalty + bp_penalty;


					WI[j] = std::min({wi_v,wi_wmb,wi_split});
					WIP[j] = std::min({wip_v,wip_wmb,wip_split});
					
					if((sparse_tree.tree[i-1].pair > (i-1) && sparse_tree.tree[i-1].pair > j) || sparse_tree.tree[i-1].pair < (i-1)){
							WI_Bbp[j] = WI[j];
							WIP_Bbp[j] = WIP[j];
					}
					if(j+1<n && sparse_tree.tree[j+1].pair < i){
						WIP_Bp[i] = WIP[j];
						WI_Bp[i] = WI[j];
					}

				}
				
				// ------------------------------------------------End of Wi/Wip--------------------------------------------------
				
				// start of WV and WVe
				if (!weakly_closed_ij){
					cand_pos_t bound_right = std::max(bp_ij,B_ij);
					cand_pos_t bound_left = std::min((cand_pos_tu) sparse_tree.Bp(i,j),(cand_pos_tu) sparse_tree.b(i,j));
					
					const energy_t wve = compute_WVe(i,j,WVe,CLVP,sparse_tree);
					const energy_t wv = compute_WV(j,bound_left,bound_right,WV,WVe,dwip1,CL,CLWMB,CLVP,sparse_tree);
					
					WVe[j] = wve;

					WV[j] = wv;
				}
				else{
					WV[j] = INF;
					WVe[j] = INF;
				}


				// ------------------------------------------------End of WV------------------------------------------------------


				/*
				The order for these should be        x i        lp               jp         j           l      ip   y
													( (        (    (   (       (          )      )  ) )      )    )
				i and ip are the outer base pair
				lp and l are the closest encompassing base pair to i/ip
														//    jp>i   j>jp  ip>i ip>j
				jp and j are some inner base pair;j has the be the closing due to the j=i+4 setup we have
				*/
				// Start of BE
				// Cannot use size_t due to overflow when pair is -2. We also cannot compare size_t and int for that same reason
				cand_pos_t ip = sparse_tree.tree[i].pair; // i's pair ip should be right side so ip = )
				cand_pos_t jp = sparse_tree.tree[j].pair; // j's pair jp should be left side so jp = (
				
				// base case: i.j and ip.jp must be in G
				
				if (jp > i && j > jp && ip > j && ip > i){ // Don't need to check if they are pairs separately because it is checked by virtue of these checks
					if(sparse_tree.tree[jp+1].pair == j-1){
						// BE_avoided++;
					}
					else{
						energy_t BE = compute_BE(i,j,ip,jp,CLBE,sparse_tree,S,S1,params,dwip1,WIP_Bp);
						register_candidate(CLBE,i,j,BE);
					}
					
				}
				else if(i == jp && ip == j){
					if(sparse_tree.tree[jp+1].pair == j-1){
						// BE_avoided++;
					}
					else{
						energy_t BE = 0;
						register_candidate(CLBE,i,j,BE);
					}
				}
				
					
				// // ------------------------------------------------End of BE---------------------------------------------------------
			}

			energy_t vp_min1 = INF;
            energy_t vp_min2 = INF;
            for ( auto it=CLVP[j].begin();CLVP[j].end() != it;++it ) {
                const cand_pos_t k = it->first;
                energy_t WIk = WI[k-1];
                energy_t WIPk = WIP[k-1];

                    vp_min1 = std::min(vp_min1,WIk + it->second);
                    vp_min2 = std::min(vp_min2,WIPk + it->second);
            }
			if((VP(i_mod,j) < INF/2) && (VP(i_mod,j) < vp_min1 || VP(i_mod,j) < vp_min2)){
				register_candidate(CLVP,i,j,VP(i_mod,j));
				inc_source_ref_count(taVP,i,j);
			}
        

			
			//Things that needed to happen later like W's wmb
			w_wmb = sparse_tree.weakly_closed(i,j) ? WMB[j] + PS_penalty : INF;
			wm_wmb = sparse_tree.weakly_closed(i,j) ? WMB[j] + PSM_penalty + b_penalty : INF;
			w  = std::min({w_v, w_split,w_wmb});
			wm = std::min({wm_v, wm_split,wm_wmb});

			// Some case of WI candidate splits will show INf as a vkj value as we are unable to do dangle versions yet but there are dangle candidates
			if (w_v < w_split || wm_v < wm_split || wi_v < wi_split || wip_v < wip_split || paired) {
				// cand_pos_t k_mod = k%(MAXLOOP+1);
				// Encode the dangles into the energies
				energy_t w_enc = (w_v << 2) | d;
				energy_t wm_enc = (wm_v << 2) | d;
				register_candidate(CL, i, j,V(i_mod,j), wm_enc,w_enc);
				// always keep arrows starting from candidates
				inc_source_ref_count(ta,i,j);
			}	
			if ((WMB[j] < INF/2) && (w_wmb < w_split || wm_wmb < wm_split || wi_wmb < wi_split || wip_wmb < wip_split)) {
		
				register_candidate(CLWMB, i, j, WMB[j]);

			}

			W[j]       = w;
			WM[j]      = wm;
			WM2[j]     = std::min(wm2_split,WMB[j] + PSM_penalty + b_penalty);

		} // end loop j
		rotate_arrays(WM2,dmli1,dmli2,WI,WIP,dwi1,dwip1,WMB,dwmbi,WV,dwvp);

		// Clean up trace arrows in i+MAXLOOP+1
		if (garbage_collect && i+MAXLOOP+1 <= n) {
            gc_row(ta,i + MAXLOOP + 1 );
			gc_row(taVP,i + MAXLOOP + 1 );
		}
		// Reallocate candidate lists in i
		for ( auto &x: CL ) {
			if (x.capacity() > 1.5*x.size()) {
				cand_list_td1 vec(x.size());
				copy(x.begin(),x.end(),vec.begin());
				vec.swap(x);
			}
		}
		for ( auto &x: CLVP ) {
			if (x.capacity() > 1.5*x.size()) {
				cand_list_t vec(x.size());
				copy(x.begin(),x.end(),vec.begin());
				vec.swap(x);
			}
		}
		for ( auto &x: CLWMB ) {
			if (x.capacity() > 1.5*x.size()) {
				cand_list_t vec(x.size());
				copy(x.begin(),x.end(),vec.begin());
				vec.swap(x);
			}
		}
		for ( auto &x: CLBE ) {
			if (x.capacity() > 1.5*x.size()) {
				cand_list_t vec(x.size());
				copy(x.begin(),x.end(),vec.begin());
				vec.swap(x);
			}
		}

		compactify(ta);
		compactify(taVP);
	}
	// std::cout << "BE avoided: " << BE_avoided << std::endl;
	return W[n];
}

/**
 * @brief Sums the number of Candidates at each index over all indices
 * 
 * @param CL_ Candidate list
 * @return total number of candidates
 */
cand_pos_t num_of_candidates(const std::vector<cand_list_td1>& CL_)  {
	cand_pos_t c=0;
	for ( const cand_list_td1 &x: CL_ ) {
		c += x.size();
	}
	return c;
}
cand_pos_t num_of_candidates(const std::vector<cand_list_t>& CL_)  {
	cand_pos_t c=0;
	for ( const cand_list_t &x: CL_ ) {
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
cand_pos_t capacity_of_candidates(const std::vector<cand_list_td1>& CL_) {
	cand_pos_t c=0;
	for ( const cand_list_td1 &x: CL_ ) {
		c += x.capacity();
	}
	return c;
}
cand_pos_t capacity_of_candidates(const std::vector<cand_list_t>& CL_) {
	cand_pos_t c=0;
	for ( const cand_list_t &x: CL_ ) {
		c += x.capacity();
	}
	return c;
}

void seqtoRNA(std::string &sequence){
    for (char &c : sequence) {
      	if (c == 'T') c = 'U';
    }
}

void validate_structure(std::string structure){
	int count = 0;
	int n = structure.length();
	for(int i = 0; i<n;++i){
		if(structure[i] == '(') count++;
		if(structure[i] == ')') count--;
	}
	if(count !=0){
		std::cout << "incorrect input as count is " << count << std::endl;
		exit(0);
	}
}

bool exists (const std::string path) {
  struct stat buffer;   
  return (stat (path.c_str(), &buffer) == 0); 
}

void get_input(std::string file, std::string &sequence, std::string &structure){
	if(!exists(file)){
		std::cout << "Input file does not exist" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::ifstream in(file.c_str());
	std::string str;
	int i = 0;
	while(getline(in,str)){
		if(str[0] == '>') continue;
		if(i==0) sequence = str;
		if(i==1) structure = str;
		++i;
	}
	in.close();
}

/**
* @brief Simple driver for @see Spark.
*
* Reads sequence from command line or stdin and calls folding and
* trace-back methods of SparseMFEFold.
*/
int main(int argc,char **argv) {

	args_info args_info;

	// get options (call gengetopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	std::string seq;
	if (args_info.inputs_num>0) {
	seq=args_info.inputs[0];
	} else {
		if(!args_info.input_file_given) std::getline(std::cin,seq);
	}

	std::string restricted = "";
	if(args_info.input_structure_given) restricted = input_structure;

	std::string fileI;
    args_info.input_file_given ? fileI = input_file : fileI = "";

	if(fileI != ""){
		
		if(exists(fileI)){
			get_input(fileI,seq,restricted);
		}
		if(seq == ""){
			std::cout << "sequence is missing from file" << std::endl; 
		}
		
	}
	cand_pos_t n = seq.length();
	std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
	if(!args_info.noConv_given) seqtoRNA(seq);
	if(restricted == "") restricted = std::string(n,'.');

	if(restricted.length() != (cand_pos_tu) n ){
		std::cout << "input sequence and structure are not the same size" << std::endl;
		std::cout << seq << std::endl;
		std::cout << restricted << std::endl;
		exit(0);
	}
	validate_structure(restricted);

	std::string file= args_info.paramFile_given ? parameter_file : "params/rna_DirksPierce09.par";
	if(exists(file)){
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
	}
	else if (seq.find('T') != std::string::npos){
		vrna_params_load_DNA_Mathews2004();
	}

	noGU = args_info.noGU_given;
	// seqtoRNA(seq);
	

	bool verbose;
	verbose = args_info.verbose_given;

	bool mark_candidates;
	mark_candidates = args_info.mark_candidates_given;
	sparse_tree tree(restricted,n);

	SparseMFEFold sparsemfefold(seq,!args_info.noGC_given,restricted);

	if(args_info.dangles_given) sparsemfefold.params_->model_details.dangles = dangle_model;
	pseudoknot = ~args_info.pseudoknot_free_given;
	pk_only = args_info.pk_only_given;  
	
	cmdline_parser_free(&args_info);

	int count = 0;
	for(int i = 1;i<=n;++i){
		if(tree.tree[i].pair > i || (tree.tree[i].pair < i && tree.tree[i].pair > 0)) count = 4;
		if(tree.tree[i].pair < 0 && count > 0){
			sparsemfefold.WI_Bbp[i] = (5 - count)*PUP_penalty;
			count--;
		}
	}
	count = 0;
	for(int i = n;i>0;--i){
		if(tree.tree[i].pair > i || (tree.tree[i].pair < i && tree.tree[i].pair > 0)) count = 4;
		if(tree.tree[i].pair < 0 && count > 0){
			sparsemfefold.WI_Bp[i] = (5 - count)*PUP_penalty;
			count--;
		}
	}
	
	energy_t mfe = fold(sparsemfefold.seq_, tree,sparsemfefold.V_,sparsemfefold.cand_comp,sparsemfefold.CL_,sparsemfefold.CLWMB_,sparsemfefold.CLVP_,sparsemfefold.CLBE_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.params_,sparsemfefold.ta_,sparsemfefold.taVP_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_, sparsemfefold.dmli1_, sparsemfefold.dmli2_,sparsemfefold.VP_, sparsemfefold.WVe_, sparsemfefold.WV_, sparsemfefold.dwvp_, sparsemfefold.WMB_,sparsemfefold.dwmbi_,sparsemfefold.WMBP_,sparsemfefold.WMBA_,sparsemfefold.WI_,sparsemfefold.dwi1_,sparsemfefold.WIP_,sparsemfefold.dwip1_,sparsemfefold.WI_Bbp,sparsemfefold.WIP_Bbp,sparsemfefold.WIP_Bp,sparsemfefold.WI_Bp,sparsemfefold.n_,sparsemfefold.garbage_collect_);		
	// std::cout << mfe << std::endl;
	std::string structure = trace_back(sparsemfefold.seq_,sparsemfefold.CL_,sparsemfefold.CLWMB_,sparsemfefold.CLBE_,sparsemfefold.CLVP_,sparsemfefold.cand_comp,sparsemfefold.structure_,sparsemfefold.params_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.ta_,sparsemfefold.taVP_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_,sparsemfefold.WI_Bbp,sparsemfefold.WIP_Bbp,sparsemfefold.WIP_Bp,sparsemfefold.WI_Bp,sparsemfefold.n_,tree, mark_candidates);
	
	std::ostringstream smfe;
	smfe << std::setiosflags(std::ios::fixed) << std::setprecision(2) << mfe/100.0 ;
	std::cout << seq << std::endl;

	std::cout << structure << " ("<<smfe.str()<<")"<<std::endl;
	// std::cout < " ("<<smfe.str()<<")"<<std::endl;
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

	std::cout << "\nPsuedoknotted\n" << std::endl;
	std::cout << "TAs num:\t"<<sizeT(sparsemfefold.taVP_)<<std::endl;
	std::cout << "TAs cap:\t"<<capacityT(sparsemfefold.taVP_)<<std::endl;
	std::cout << "TA av:\t"<<avoidedT(sparsemfefold.taVP_)<<std::endl;
	std::cout << "TA rm:\t"<<erasedT(sparsemfefold.taVP_)<<std::endl;
	std::cout <<std::endl;
	std::cout << "WMB Can num:\t"<<num_of_candidates(sparsemfefold.CLWMB_)<<std::endl;
	std::cout << "WMB Can cap:\t"<<capacity_of_candidates(sparsemfefold.CLWMB_)<<std::endl;
	std::cout << "VP Can num:\t"<<num_of_candidates(sparsemfefold.CLVP_)<<std::endl;
	std::cout << "VP Can cap:\t"<<capacity_of_candidates(sparsemfefold.CLVP_)<<std::endl;
	std::cout << "BE Can num:\t"<<num_of_candidates(sparsemfefold.CLBE_)<<std::endl;
	std::cout << "BE Can cap:\t"<<capacity_of_candidates(sparsemfefold.CLBE_)<<std::endl;
	}

	return 0;
}
