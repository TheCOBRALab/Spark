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
#define debug 0
#include "PK_globals.hh"
#include "base_types.hh"
#include "cmdline.hh"
#include "matrix.hh"
#include "sparse_tree.cc"
#include "trace_arrow.hh"
#include <cassert>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

extern "C" {
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/params/io.h"
}
static bool pseudoknot = false;

static bool pk_only = false;

#define INFover2 5000000 /* (INT_MAX/20) */

struct quatret {
    cand_pos_t first;
    energy_t second;
    energy_t third;
    energy_t fourth;
    quatret() {
        first = 1;
        second = 2;
        third = 3;
        fourth = 4;
    }
    quatret(cand_pos_t x, energy_t y, energy_t z, energy_t w) {
        first = x;
        second = y;
        third = z;
        fourth = w;
    }
};

typedef std::pair<cand_pos_t, energy_t> cand_entry_t;
typedef std::vector<cand_entry_t> cand_list_t;

typedef quatret cand_entry_td1;
typedef std::vector<cand_entry_td1> cand_list_td1;

class Spark;

energy_t ILoopE(const short *S, const short *S1, const paramT *params, const pair_type &ptype_closing, const cand_pos_t &i, const cand_pos_t &j,
                const cand_pos_t &k, const cand_pos_t &l);

/**
 * Space efficient sparsification of Zuker-type RNA folding with
 * trace-back. Provides methods for the evaluation of dynamic
 * programming recursions and the trace-back.
 */
class Spark {

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
    std::vector<energy_t> WMBP_;
    std::vector<energy_t> WMBA_;
    std::vector<energy_t> WI_;
    std::vector<energy_t> dwi1_; // WI from 1 iteration ago
    std::vector<energy_t> WIP_;
    std::vector<energy_t> dwip1_; // WIP from 1 iteration ago
    std::vector<energy_t> WV_;
    std::vector<energy_t> dwvp_; // WV from 1 iteration ago;

    std::vector<energy_t> WI_Bbp;  // WI from band borders on left
    std::vector<energy_t> WIP_Bbp; // WIP from band borders on left
    std::vector<energy_t> WIP_Bp;  // WIP from band borders on the right

    bool mark_candidates_;

    TraceArrows ta_;
    TraceArrows taVP_;

    // TraceArrows ta_dangle_;

    std::vector<cand_list_td1> CL_;
    std::vector<cand_list_t> CLVP_;
    std::vector<cand_list_t> CLWMB_;
    std::vector<cand_list_t> CLBE_;
    std::vector<cand_list_t> CLBEO_;

    /**
    candidate list for decomposition in W or WM

    @note Avoid separate candidate lists CLW and CLWM for split cases in W and
    WM to save even more space; here, this works after
    reformulating the recursions such that both split-cases recurse to
    V-entries. (compare OCTs)
    */

    // compare candidate list entries by keys (left index i) in descending order
    struct Cand_comp {
        bool operator()(const cand_entry_t &x, cand_pos_t y) const { return x.first > y; }
        bool operator()(const cand_entry_td1 &x, cand_pos_t y) const { return x.first > y; }
    } cand_comp;

    Spark(const std::string &seq, bool garbage_collect, std::string restricted)
        : seq_(seq), n_(seq.length()), params_(scale_parameters()), garbage_collect_(garbage_collect), ta_(n_), taVP_(n_) {
        make_pair_matrix();

        S_ = encode_sequence(seq.c_str(), 0);
        S1_ = encode_sequence(seq.c_str(), 1);

        V_.resize(MAXLOOP + 1, n_ + 1, INF);
        W_.resize(n_ + 1, 0);
        WM_.resize(n_ + 1, INF);
        WM2_.resize(n_ + 1, INF);
        dmli1_.resize(n_ + 1, INF);
        dmli2_.resize(n_ + 1, INF);

        // Pseudoknot portion

        VP_.resize(MAXLOOP + 1, n_ + 1, INF);
        WVe_.resize(n_ + 1, INF);
        WMB_.resize(n_ + 1, INF);
        WMBP_.resize(n_ + 1, INF);
        WMBA_.resize(n_ + 1, INF);
        WI_.resize(n_ + 1, 0);
        dwi1_.resize(n_ + 1, 0);
        WIP_.resize(n_ + 1, INF);
        dwip1_.resize(n_ + 1, INF);
        WV_.resize(n_ + 1, INF);
        dwvp_.resize(n_ + 1, INF);

        WI_Bbp.resize(n_ + 1, 0);
        WIP_Bbp.resize(n_ + 1, INF);
        WIP_Bp.resize(n_ + 1, INF);

        // init candidate lists
        CL_.resize(n_ + 1);
        CLWMB_.resize(n_ + 1);
        CLVP_.resize(n_ + 1);
        CLBE_.resize(n_ + 1);
        CLBEO_.resize(n_ + 1);

        resize(ta_, n_ + 1);
        resize(taVP_, n_ + 1);

        // resize(ta_dangle_,n_+1);

        restricted_ = restricted;
    }

    ~Spark() {
        free(params_);
        free(S_);
        free(S1_);
    }
};

void trace_V(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_W(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, sparse_tree &tree);
void trace_WM(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_WM2(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, sparse_tree &tree);
void trace_WMB(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_VP(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_WI(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_WIP(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_WV(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_WVe(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_WMBP(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);
void trace_BE(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t ip, energy_t e, sparse_tree &tree);
void trace_WMBA(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree);

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
void rotate_arrays(Spark &spark) {
    spark.dmli2_.swap(spark.dmli1_);
    spark.dmli1_.swap(spark.WM2_);
    spark.dwi1_.swap(spark.WI_);
    spark.dwip1_.swap(spark.WIP_);
    spark.dwvp_.swap(spark.WV_);
}

/**
 * @brief This code returns the hairpin energy for a given base pair.
 * @param i The left index in the base pair
 * @param j The right index in the base pair
 */
energy_t HairpinE(const std::string &seq, const short *S, const short *S1, const paramT *params, cand_pos_t i, cand_pos_t j) {

    const pair_type ptype_closing = pair[S[i]][S[j]];

    if (ptype_closing == 0) return INF;

    return E_Hairpin(j - i - 1, ptype_closing, S1[i + 1], S1[j - 1], &seq.c_str()[i - 1], const_cast<paramT *>(params));
}

/**
 * @brief Returns the internal loop energy for a given i.j and k.l
 *
 */
energy_t ILoopE(const short *S, const short *S1, const paramT *params, const pair_type &ptype_closing, const cand_pos_t &i, const cand_pos_t &j,
                const cand_pos_t &k, const cand_pos_t &l) {
    assert(ptype_closing > 0);
    assert(1 <= i);
    assert(i < k);
    assert(k < l);
    assert(l < j);

    // note: enclosed bp type 'turned around' for lib call
    const pair_type ptype_enclosed = rtype[pair[S[k]][S[l]]];

    if (ptype_enclosed == 0) return INF;

    return E_IntLoop(k - i - 1, j - l - 1, ptype_closing, ptype_enclosed, S1[i + 1], S1[j - 1], S1[k - 1], S1[l + 1], const_cast<paramT *>(params));
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
energy_t E_ext_Stem(const energy_t &vij, const energy_t &vi1j, const energy_t &vij1, const energy_t &vi1j1, const short *S, paramT *params,
                    const cand_pos_t i, const cand_pos_t j, Dangle &d, cand_pos_t n, const std::vector<Node> &tree) {

    energy_t e = INF, en = INF;
    pair_type tt = pair[S[i]][S[j]];

    if ((tree[i].pair < -1 && tree[j].pair < -1) || (tree[i].pair == j && tree[j].pair == i)) {
        en = vij; // i j

        if (en != INF) {
            if (params->model_details.dangles == 2) {
                base_type si1 = i > 1 ? S[i - 1] : -1;
                base_type sj1 = j < n ? S[j + 1] : -1;
                en += vrna_E_ext_stem(tt, si1, sj1, params);
            } else {
                en += vrna_E_ext_stem(tt, -1, -1, params);
                d = 0;
            }

            e = std::min(e, en);
        }
    }

    if (params->model_details.dangles == 1) {
        tt = pair[S[i + 1]][S[j]];
        if (((tree[i + 1].pair < -1 && tree[j].pair < -1) || (tree[i + 1].pair == j)) && tree[i].pair < 0) {
            en = (j - i - 1 > TURN) ? vi1j : INF; // i+1 j

            if (en != INF) {

                base_type si1 = S[i];
                en += vrna_E_ext_stem(tt, si1, -1, params);
            }

            e = std::min(e, en);
            if (e == en) {
                d = 1;
            }
        }
        tt = pair[S[i]][S[j - 1]];
        if (((tree[i].pair < -1 && tree[j - 1].pair < -1) || (tree[i].pair == j - 1)) && tree[j].pair < 0) {
            en = (j - 1 - i > TURN) ? vij1 : INF; // i j-1
            if (en != INF) {

                base_type sj1 = S[j];

                en += vrna_E_ext_stem(tt, -1, sj1, params);
            }
            e = std::min(e, en);
            if (e == en) {
                d = 2;
            }
        }
        tt = pair[S[i + 1]][S[j - 1]];
        if (((tree[i + 1].pair < -1 && tree[j - 1].pair < -1) || (tree[i + 1].pair == j - 1)) && tree[i].pair < 0 && tree[j].pair < 0) {
            en = (j - 1 - i - 1 > TURN) ? vi1j1 : INF; // i+1 j-1

            if (en != INF) {

                base_type si1 = S[i];
                base_type sj1 = S[j];

                en += vrna_E_ext_stem(tt, si1, sj1, params);
            }
            e = std::min(e, en);
            if (e == en) {
                d = 3;
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
energy_t E_MbLoop(const std::vector<energy_t> &dmli1, const std::vector<energy_t> &dmli2, const short *S, paramT *params, cand_pos_t i, cand_pos_t j,
                  const std::vector<Node> &tree) {

    energy_t e = INF, en = INF;
    pair_type tt = pair[S[j]][S[i]];
    bool pairable = (tree[i].pair < -1 && tree[j].pair < -1) || (tree[i].pair == j);

    /* double dangles */
    switch (params->model_details.dangles) {
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

        if (pairable && tree[i + 1].pair < 0) {
            en = dmli2[j - 1];

            if (en != INF) {

                base_type si1 = S[i + 1];

                en += E_MLstem(tt, -1, si1, params) + params->MLclosing + params->MLbase;
            }
        }
        e = std::min(e, en);

        /**
         * ML pair 3
         * new closing pair (i,j) with mb part [i+1, j-2]
         */
        if (pairable && tree[j - 1].pair < 0) {
            en = dmli1[j - 2];

            if (en != INF) {
                base_type sj1 = S[j - 1];

                en += E_MLstem(tt, sj1, -1, params) + params->MLclosing + params->MLbase;
            }
        }
        e = std::min(e, en);
        /**
         * ML pair 53
         * new closing pair (i,j) with mb part [i+2.j-2]
         */
        if (pairable && tree[i + 1].pair < 0 && tree[j - 1].pair < 0) {
            en = dmli2[j - 2];

            if (en != INF) {

                base_type si1 = S[i + 1];
                base_type sj1 = S[j - 1];

                en += E_MLstem(tt, sj1, si1, params) + params->MLclosing + 2 * params->MLbase;
            }
        }
        e = std::min(e, en);
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
energy_t E_MLStem(const energy_t &vij, const energy_t &vi1j, const energy_t &vij1, const energy_t &vi1j1, const short *S, paramT *params,
                  cand_pos_t i, cand_pos_t j, Dangle &d, const cand_pos_t &n, const std::vector<Node> &tree) {

    energy_t e = INF, en = INF;

    pair_type type = pair[S[i]][S[j]];

    if ((tree[i].pair < -1 && tree[j].pair < -1) || (tree[i].pair == j)) {
        en = vij; // i j
        if (en != INF) {
            if (params->model_details.dangles == 2) {
                base_type mm5 = i > 1 ? S[i - 1] : -1;
                base_type mm3 = j < n ? S[j + 1] : -1;
                en += E_MLstem(type, mm5, mm3, params);
            } else {
                en += E_MLstem(type, -1, -1, params);
                d = 0;
            }
            e = std::min(e, en);
        }
    }
    if (params->model_details.dangles == 1) {
        const base_type mm5 = S[i], mm3 = S[j];

        if (((tree[i + 1].pair < -1 && tree[j].pair < -1) || (tree[i + 1].pair == j)) && tree[i].pair < 0) {
            en = (j - i - 1 > TURN) ? vi1j : INF; // i+1 j
            if (en != INF) {
                en += params->MLbase;

                type = pair[S[i + 1]][S[j]];
                en += E_MLstem(type, mm5, -1, params);

                e = std::min(e, en);
                if (e == en) {
                    d = 1;
                }
            }
        }

        if (((tree[i].pair < -1 && tree[j - 1].pair < -1) || (tree[i].pair == j - 1)) && tree[j].pair < 0) {
            en = (j - 1 - i > TURN) ? vij1 : INF; // i j-1
            if (en != INF) {
                en += params->MLbase;

                type = pair[S[i]][S[j - 1]];
                en += E_MLstem(type, -1, mm3, params);

                e = std::min(e, en);
                if (e == en) {
                    d = 2;
                }
            }
        }
        if (((tree[i + 1].pair < -1 && tree[j - 1].pair < -1) || (tree[i + 1].pair == j - 1)) && tree[i].pair < 0 && tree[j].pair < 0) {
            en = (j - 1 - i - 1 > TURN) ? vi1j1 : INF; // i+1 j-1
            if (en != INF) {
                en += 2 * params->MLbase;

                type = pair[S[i + 1]][S[j - 1]];
                en += E_MLstem(type, mm5, mm3, params);

                e = std::min(e, en);
                if (e == en) {
                    d = 3;
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
int compute_exterior_cases(cand_pos_t l, cand_pos_t j, sparse_tree &tree) {

    // Case 1 -> l is not covered
    bool case1 = tree.tree[l].parent->index <= 0;
    // Case 2 -> l is paired
    bool case2 = tree.tree[l].pair > 0;
    // Case 3 -> l is part of a closed subregion
    bool case3 = 0;
    // Case 4 -> l.bp(l) i.e. l.j does not cross anything -- could I compare parents instead?
    bool case4 = j < tree.Bp(l, j);
    // By bitshifting each one, we have a more granular idea of what cases fail and is faster than branching
    return (case1 << 3) | (case2 << 2) | (case3 << 1) | case4;
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
const std::vector<energy_t> recompute_WM(Spark &spark, cand_pos_t i, cand_pos_t max_j, const std::vector<Node> &tree,
                                         const std::vector<cand_pos_t> &up) {
    assert(i >= 1);
    assert(max_j <= spark.n_);
    std::vector<energy_t> temp = spark.WM_;

    for (cand_pos_t j = i - 1; j <= std::min(i + TURN, max_j); j++) {
        temp[j] = INF;
    }

    for (cand_pos_t j = i + TURN + 1; j <= max_j; j++) {
        energy_t wm = INF;
        for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first >= i; ++it) {
            cand_pos_t k = it->first;
            energy_t v_kj = it->third >> 2;

            bool can_pair = up[k - 1] >= (k - i);
            if (can_pair) wm = std::min(wm, static_cast<energy_t>(spark.params_->MLbase * (k - i)) + v_kj);
            wm = std::min(wm, temp[k - 1] + v_kj);
        }
        for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first > i + TURN + 1; ++it) {

            cand_pos_t k = it->first;
            energy_t wmb_kj = it->second + PSM_penalty + b_penalty;
            bool can_pair = up[k - 1] >= (k - i);
            if (can_pair) wm = std::min(wm, static_cast<energy_t>(spark.params_->MLbase * (k - i)) + wmb_kj);

            wm = std::min(wm, temp[k - 1] + wmb_kj);
        }
        if (tree[j].pair < 0) wm = std::min(wm, temp[j - 1] + spark.params_->MLbase);
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
const std::vector<energy_t> recompute_WM2(Spark &spark, cand_pos_t i, cand_pos_t max_j,
                                          const std::vector<Node> &tree, const std::vector<cand_pos_t> &up) {
    assert(i >= 1);
    assert(max_j <= spark.n_);
    std::vector<energy_t> temp = spark.WM2_;

    for (cand_pos_t j = i - 1; j <= std::min(i + 2 * TURN + 2, max_j); j++) {
        temp[j] = INF;
    }

    for (cand_pos_t j = i + 2 * TURN + 3; j <= max_j; j++) {
        energy_t wm2 = INF;
        for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first > i + TURN + 1; ++it) {

            cand_pos_t k = it->first;
            energy_t v_kj = it->third >> 2;

            wm2 = std::min(wm2, spark.WM_[k - 1] + v_kj);
        }
        for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t wmb_kj = it->second + PSM_penalty + b_penalty;
            bool can_pair = up[k - 1] >= (k - i);

            wm2 = std::min(wm2, spark.WM_[k - 1] + wmb_kj);
            if (can_pair) wm2 = std::min(wm2, static_cast<energy_t>(spark.params_->MLbase * (k - i)) + wmb_kj);
        }
        if (tree[j].pair < 0) wm2 = std::min(wm2, temp[j - 1] + spark.params_->MLbase);
        temp[j] = wm2;
    }
    return temp;
}

/**
 * @brief Recompute row of WMBP. This is used in the traceback when we haved decided the current ij pair closes a psuedoknot,
 * and the WMBP energies need to be recomputed such that the pseudoknot can be broken down. This happens infrequently enough such that it not too
 * costly
 *
 * @param CL V candidate list
 * @param CLVP VP candidate list
 * @param CLWMB WMB Candidate list
 * @param CLBE BE candidate list
 * @param WI_Bbp array which holds the values where the left index is either an outer closing base or an inner opening base
 * @param i Current i
 * @param max_j Current j
 */
void recompute_WMBP(Spark &spark, cand_pos_t i, cand_pos_t max_j, sparse_tree &tree) {
    assert(i >= 1);
    assert(max_j <= spark.n_);

    for (cand_pos_t j = i + TURN + 1; j <= max_j; j++) {
        energy_t wmbp = INF;
        energy_t wmba = INF;

        energy_t vp_ij = INF;
        for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t v_kj = it->second;

            wmba = std::min(wmba, spark.WMBP_[k - 1] + v_kj + PPS_penalty);
        }
        for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t wmb_kj = it->second;

            wmba = std::min(wmba, spark.WMBP_[k - 1] + wmb_kj + PPS_penalty + PSM_penalty);
        }
        // m2
        if (tree.tree[j].pair < 0) {
            for (auto it = spark.CLVP_[j].begin(); spark.CLVP_[j].end() != it && it->first >= i; ++it) {
                cand_pos_t k = it->first;
                if (k == i) vp_ij = it->second; // second?
                cand_pos_t bp_ik = tree.bp(i, k);
                cand_pos_t Bp_kj = tree.Bp(k, j);
                cand_pos_t b_ij = tree.b(i, j);
                int ext_case = compute_exterior_cases(k, j, tree);
                if ((b_ij > 0 && k < b_ij) || (b_ij < 0 && ext_case == 0)) {
                    if (bp_ik >= 0 && k > bp_ik && Bp_kj > 0 && k < Bp_kj) {
                        energy_t BE_energy = INF;
                        cand_pos_t B_kj = tree.B(k, j);
                        cand_pos_t b_kj = (B_kj > 0) ? tree.tree[B_kj].pair : -2;
                        for (auto it2 = spark.CLBE_[Bp_kj].begin(); spark.CLBE_[Bp_kj].end() != it2; ++it2) {
                            cand_pos_t l = it2->first;
                            if (l == b_kj) {
                                BE_energy = it2->second;
                                break;
                            }
                        }
                        wmbp = std::min(wmbp, BE_energy + spark.WMBA_[k - 1] + it->second + 2 * PB_penalty);
                    }
                }
            }
        }
        // m3
        if (tree.tree[j].pair < 0 && tree.tree[i].pair >= 0 && tree.tree[i].pair < j) {
            for (auto it = spark.CLVP_[j].begin(); spark.CLVP_[j].end() != it && it->first >= i; ++it) {
                cand_pos_t k = it->first;

                cand_pos_t bp_ik = tree.bp(i, k);
                if (bp_ik >= 0 && k + TURN <= j) {
                    cand_pos_t Bp_ik = tree.tree[bp_ik].pair;
                    energy_t BE_energy = INF;
                    for (auto it2 = spark.CLBE_[Bp_ik].begin(); spark.CLBE_[Bp_ik].end() != it2; ++it2) {
                        cand_pos_t l = it2->first;
                        if (l == i) {
                            BE_energy = it2->second;
                            break;
                        }
                    }
                    wmbp = std::min(wmbp, BE_energy + spark.WI_Bbp[k - 1] + it->second + 2 * PB_penalty);
                }
            }
        }
        wmbp = std::min(wmbp, vp_ij + PB_penalty);
        if (tree.tree[j].pair < 0) wmba = std::min(wmba, spark.WMBA_[j - 1] + PUP_penalty);
        wmba = std::min(wmba, wmbp);
        spark.WMBA_[j] = wmba;
        spark.WMBP_[j] = wmbp;
    }
}

/**
 * @brief Recompute row of WI. This is used in the traceback when we haved decided the current ij pair is a crossing base pair, we need to calculate
 * the energy, and the WI energies need to be recomputed such that the VP energy can be broken down. This happens only when the calculations are
 * confirmed to be needed such that it not too costly
 *
 * @param CL V candidate list
 * @param CLWMB WMB candidate list
 * @param i Current i
 * @param max_j Current j
 */
void recompute_WI(Spark &spark, cand_pos_t i, cand_pos_t max_j, const std::vector<Node> &tree, const std::vector<cand_pos_t> &up) {
    assert(i >= 1);
    assert(max_j <= spark.n_);

    // Causes vector resize error if the ifs are not there because if i is close to n, it would go past the bounds
    for (cand_pos_t j = 0; j < 4; ++j) {
        if (i + j < spark.n_) {
            spark.WI_[i + j] = (j + 1) * PUP_penalty;
        }
    }

    for (cand_pos_t j = i + TURN + 1; j <= max_j; j++) {
        energy_t wi = INF;
        for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t v_kj = it->second + PPS_penalty;

            wi = std::min(wi, spark.WI_[k - 1] + v_kj);
            bool can_pair = up[k - 1] >= (k - i);
            if (can_pair) wi = std::min(wi, static_cast<energy_t>(PUP_penalty * (k - i)) + v_kj);
        }
        for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t wmb_kj = it->second + PSM_penalty + PPS_penalty;

            wi = std::min(wi, spark.WI_[k - 1] + wmb_kj);
            bool can_pair = up[k - 1] >= (k - i);
            if (can_pair) wi = std::min(wi, static_cast<energy_t>(PUP_penalty * (k - i)) + wmb_kj);
        }
        if (tree[j].pair < 0) wi = std::min(wi, spark.WI_[j - 1] + PUP_penalty);
        spark.WI_[j] = wi;
    }
}

/**
 * @brief Recompute row of WIP. This is used in the traceback when we haved decided the current ij pair is a multiloop that spans a band or part of a
 * BE, we need to calculate the energy, and the WIP energies need to be recomputed such that the VP or BE energy can be broken down. This happens only
 * when the calculations are confirmed to be needed such that it not too costly
 *
 * @param CL V candidate list
 * @param CLWMB WMB candidate list
 * @param i Current i
 * @param max_j Current j
 */
void recompute_WIP(Spark &spark, cand_pos_t i, cand_pos_t max_j, const std::vector<Node> &tree, const std::vector<cand_pos_t> &up) {
    assert(i >= 1);
    assert(max_j <= spark.n_);
    for (cand_pos_t j = i + TURN + 1; j <= max_j; j++) {
        energy_t wip = INF;
        for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t v_kj = it->second + bp_penalty;

            wip = std::min(wip, spark.WIP_[k - 1] + v_kj);
            bool can_pair = up[k - 1] >= (k - i);
            if (can_pair) wip = std::min(wip, static_cast<energy_t>(cp_penalty * (k - i)) + v_kj);
        }
        for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first >= i; ++it) {

            cand_pos_t k = it->first;
            energy_t wmb_kj = it->second + PSM_penalty + bp_penalty;

            wip = std::min(wip, spark.WIP_[k - 1] + wmb_kj);
            bool can_pair = up[k - 1] >= (k - i);
            if (can_pair) wip = std::min(wip, static_cast<energy_t>(cp_penalty * (k - i)) + wmb_kj);
        }
        if (tree[j].pair < 0) wip = std::min(wip, spark.WIP_[j - 1] + cp_penalty);
        spark.WIP_[j] = wip;
    }
}
/**
 * @brief Recompute row of WVe. This is used in the traceback when we haved decided the current ij pair is a multiloop that spans a band, we need to
 * calculate the energy, and the WVe energies need to be recomputed such that the VP energy can be broken down. This happens infrequently such that it
 * not too costly
 *
 * @param CLVP VP candidate list
 * @param i Current i
 * @param max_j Current j
 */
void recompute_WVe(Spark &spark, cand_pos_t i, cand_pos_t max_j, sparse_tree &tree) {
    assert(i >= 1);
    assert(max_j <= spark.n_);
    for (cand_pos_t j = i + TURN + 1; j <= max_j; j++) {
        energy_t wve = INF;
        for (auto it = spark.CLVP_[j].begin(); spark.CLVP_[j].end() != it && it->first >= i; ++it) {
            cand_pos_t k = it->first;
            energy_t vp_kj = it->second;

            bool can_pair = tree.up[k - 1] >= (k - i);
            if (can_pair) wve = std::min(wve, static_cast<energy_t>(cp_penalty * (k - i)) + vp_kj);
        }
        if (tree.tree[j].pair < 0) wve = std::min(wve, spark.WVe_[j - 1] + cp_penalty);
        spark.WVe_[j] = wve;
    }
}
/**
 * @brief Recompute row of WV. This is used in the traceback when we haved decided the current ij pair is a multiloop that spans a band, we need to
 * calculate the energy, and the WV energies need to be recomputed such that the VP energy can be broken down. This happens infrequently such that it
 * not too costly
 *
 * @param CL V candidate list
 * @param CLWMB WMB candidate list
 * @param i Current i
 * @param CLVP VP candidate list
 * @param WIP WIP energies
 * @param i Current i
 * @param max_j Current j
 */
void recompute_WV(Spark &spark, cand_pos_t i, cand_pos_t max_j, sparse_tree &tree) {
    assert(i >= 1);
    assert(max_j <= spark.n_);

    for (cand_pos_t j = i + TURN + 1; j <= max_j; j++) {
        cand_pos_t bound_right = std::max(tree.bp(i, j), tree.B(i, j));
        cand_pos_t bound_left = std::min((cand_pos_tu)tree.Bp(i, j), (cand_pos_tu)tree.b(i, j));

        energy_t wv = INF;
        if (!tree.weakly_closed(i, j)) {
            for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first > bound_right; ++it) {

                cand_pos_t k = it->first;
                // energy_t v_kj = it->third >> 2;
                energy_t v_kj = it->second;
                // Dangle d = it->third & 3;
                // cand_pos_t num = 0;
                // if(d == 1 || d==2) num = 1;
                // else if(d==3 && params->model_details.dangles == 1) num =2;
                // energy_t fix = num*cp_penalty - num*params->MLbase;

                wv = std::min(wv, spark.WVe_[k - 1] + v_kj + bp_penalty);
                wv = std::min(wv, spark.WV_[k - 1] + v_kj + bp_penalty);
            }
            // How do I do the VP one?

            for (auto it = spark.CLVP_[j].begin(); spark.CLVP_[j].end() != it; ++it) {

                cand_pos_t k = it->first;
                energy_t vp_kj = it->second;
                if (k < bound_left) wv = std::min(wv, spark.WIP_[k - 1] + vp_kj);
            }
            for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first > bound_right; ++it) {

                cand_pos_t k = it->first;
                energy_t wmb_kj = it->second;

                wv = std::min(wv, spark.WVe_[k - 1] + wmb_kj + PSM_penalty + bp_penalty);
                wv = std::min(wv, spark.WV_[k - 1] + wmb_kj + PSM_penalty + bp_penalty);
            }
            if (tree.tree[j].pair < 0) wv = std::min(wv, spark.WV_[j - 1] + cp_penalty);
            spark.WV_[j] = wv;
        }
    }
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
bool is_candidate(const std::vector<cand_list_td1> &CL, const Spark::Cand_comp &cand_comp, cand_pos_t i, cand_pos_t j) {
    const cand_list_td1 &list = CL[j];

    auto it = std::lower_bound(list.begin(), list.end(), i, cand_comp);

    return it != list.end() && it->first == i;
}
bool is_candidate(const std::vector<cand_list_t> &CL, const Spark::Cand_comp &cand_comp, cand_pos_t i, cand_pos_t j) {
    const cand_list_t &list = CL[j];

    auto it = std::lower_bound(list.begin(), list.end(), i, cand_comp);

    return it != list.end() && it->first == i;
}

/**
 * @brief Sparse function to move through WMB case 1 without linear check
 *
 * @param j right closing base pair
 */
void compute_WMB_case1(cand_pos_t j, energy_t &m1, energy_t &BE_en, cand_pos_t &best_border, const sparse_tree &tree,
                       const std::vector<cand_list_t> &CLBEO, const std::vector<energy_t> &WMBA) {
    // We are moving through the elements of BE/ the stucture in G
    // We are taking advantage of the fact that WMBA holds both the actual pseudoknotted base pair from WMB' and anything to the right of it
    // If there exists candidates in BEO, then they encompass the areas to be sparsely decomposed
    // We can calculate the BE from the candidate energy and the WMBA encompasses all values within and can be sparsely decomposed
    cand_pos_t bp_j = tree.tree[j].pair; // j' or opening base pair
    for (auto it = CLBEO[bp_j].begin(); CLBEO[bp_j].end() != it; ++it) {
        cand_pos_t candidate_index = it->first; // j or an inner closing base pair
        energy_t BE_energy = it->second;
        energy_t WMBA_energy = WMBA[candidate_index - 1];
        if (BE_energy + WMBA_energy < m1) {
            m1 = BE_energy + WMBA_energy;
            BE_en = BE_energy;
            best_border = candidate_index;
        }
    }
}
/**
 * @brief Determines the type of dangle being used for a closing multiloop while in traceback.
 *
 * @param WM2ij The WM2 energy for the region [i,j]
 * @param WM2i1j The WM2 energy for the region [i+1,j]
 * @param WM2ij1 The WM2 energy for the region [i,j-1]
 * @param WM2i1j1 The WM2 energy for the region [i+1,j-1]
 */
void find_mb_dangle(const energy_t WM2ij, const energy_t WM2i1j, const energy_t WM2ij1, const energy_t WM2i1j1, paramT *params, const short *S,
                    const cand_pos_t i, const cand_pos_t j, cand_pos_t &k, cand_pos_t &l, const std::vector<Node> &tree) {

    const pair_type tt = pair[S[j]][S[i]];
    const energy_t e1 = WM2ij + E_MLstem(tt, -1, -1, params);
    const energy_t e2 = WM2i1j + E_MLstem(tt, -1, S[i + 1], params);
    const energy_t e3 = WM2ij1 + E_MLstem(tt, S[j - 1], -1, params);
    const energy_t e4 = WM2i1j1 + E_MLstem(tt, S[j - 1], S[i + 1], params);
    energy_t e = e1;

    if (e2 < e && tree[i + 1].pair < 0) {
        e = e2;
        k = i + 2;
        l = j - 1;
    }
    if (e3 < e && tree[j - 1].pair < 0) {
        e = e3;
        k = i + 1;
        l = j - 2;
    }
    if (e4 < e && tree[i + 1].pair < 0 && tree[j - 1].pair < 0) {
        e = e4;
        k = i + 2;
        l = j - 2;
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
void trace_W(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, sparse_tree &tree) {
    if (debug) printf("W at %d and %d with %d\n", i, j, spark.W_[j]);

    if (i + TURN + 1 > j) return;
    // case j unpaired
    if (spark.W_[j] == spark.W_[j - 1]) {
        trace_W(spark, mark_candidates, i, j - 1, tree);
        return;
    }

    cand_pos_t m = j + 1;
    energy_t w = INF;
    for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first >= i; ++it) {
        m = it->first;
        const energy_t wmb_kj = it->second;
        w = spark.W_[m - 1] + wmb_kj;
        if (spark.W_[j] == w + PS_penalty) {
            trace_W(spark, mark_candidates, i, m - 1, tree);
            trace_WMB(spark, mark_candidates, m, j, wmb_kj, tree);
            return;
        }
    }
    energy_t v = INF;
    w = INF;
    Dangle dangle = 3;
    energy_t vk = INF;
    for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first >= i; ++it) {
        m = it->first;
        const energy_t v_kj = it->fourth >> 2;
        const Dangle d = it->fourth & 3;
        w = spark.W_[m - 1] + v_kj;
        if (spark.W_[j] == w) {
            v = it->second;
            dangle = d;
            vk = v_kj;
            break;
        }
    }
    cand_pos_t k = m;
    cand_pos_t l = j;
    pair_type ptype = 0;
    switch (dangle) {
    case 0:
        ptype = pair[spark.S_[k]][spark.S_[l]];
        v = vk - E_ExtLoop(ptype, -1, -1, spark.params_);
        break;
    case 1:
        k = m + 1;
        ptype = pair[spark.S_[k]][spark.S_[l]];
        v = vk - E_ExtLoop(ptype, spark.S_[m], -1, spark.params_);
        break;
    case 2:
        l = j - 1;
        ptype = pair[spark.S_[k]][spark.S_[l]];
        v = vk - E_ExtLoop(ptype, -1, spark.S_[j], spark.params_);
        break;
    case 3:
        if (spark.params_->model_details.dangles == 1) {
            k = m + 1;
            l = j - 1;
            ptype = pair[spark.S_[k]][spark.S_[l]];
            v = vk - E_ExtLoop(ptype, spark.S_[m], spark.S_[j], spark.params_);
        }
        break;
    }
    assert(i <= m && m < j);
    assert(v < INF);
    // don't recompute W, since i is not changed
    trace_W(spark, mark_candidates, i, m - 1, tree);
    trace_V(spark, mark_candidates, k, l, v, tree);
}

/**
 * @brief Traceback from V entry
 *
 * @param structure Final Structure
 * @param mark_candidates Whether Candidates should be [ ]
 * @param i row index
 * @param j column index
 */
void trace_V(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree) {
    if (debug) printf("V at %d and %d with %d\n", i, j, e);

    assert(i + TURN + 1 <= j);
    assert(j <= spark.n_);

    if (mark_candidates && is_candidate(spark.CL_, spark.cand_comp, i, j)) {
        spark.structure_[i] = '{';
        spark.structure_[j] = '}';
    } else {
        spark.structure_[i] = '(';
        spark.structure_[j] = ')';
    }
    const pair_type ptype_closing = pair[spark.S_[i]][spark.S_[j]];
    if (exists_trace_arrow_from(spark.ta_, i, j)) {

        const TraceArrow &arrow = trace_arrow_from(spark.ta_, i, j);
        const cand_pos_t k = arrow.k(i);
        const cand_pos_t l = arrow.l(j);
        assert(i < k);
        assert(l < j);
        trace_V(spark, mark_candidates, k, l, arrow.target_energy(), tree);
        return;

    } else {

        // try to trace back to a candidate: (still) interior loop case
        cand_pos_t l_min = std::max(i, j - 31);
        for (cand_pos_t l = j - 1; l > l_min; l--) {
            // Break if it's an assured dangle case
            for (auto it = spark.CL_[l].begin(); spark.CL_[l].end() != it && it->first > i; ++it) {
                const cand_pos_t k = it->first;
                if (k - i > 31) continue;
                if (e
                    == it->second
                           + E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[spark.S_[k]][spark.S_[l]]], spark.S1_[i + 1], spark.S1_[j - 1],
                                       spark.S1_[k - 1], spark.S1_[l + 1], const_cast<paramT *>(spark.params_))) {
                    trace_V(spark, mark_candidates, k, l, it->second, tree);
                    return;
                }
            }
        }
    }
    // is this a hairpin?
    if (e == HairpinE(spark.seq_, spark.S_, spark.S1_, spark.params_, i, j)) {
        return;
    }

    // if we are still here, trace to wm2 (split case);
    // in this case, we know the 'trace arrow'; the next row has to be recomputed
    std::vector<energy_t> temp;
    if (spark.params_->model_details.dangles == 1) {
        temp = recompute_WM(spark, i + 2, j - 1, tree.tree, tree.up);
        spark.WM_ = temp;
        spark.dmli1_ = recompute_WM2(spark, i + 2, j - 1, tree.tree, tree.up);
    }
    spark.WM_  = recompute_WM(spark, i + 1, j - 1, tree.tree, tree.up);
    spark.WM2_  = recompute_WM2(spark, i + 1, j - 1, tree.tree, tree.up);

    // Dangle for Multiloop
    cand_pos_t k = i + 1;
    cand_pos_t l = j - 1;
    if (spark.params_->model_details.dangles == 1) {
        find_mb_dangle(spark.WM2_[j - 1], spark.dmli1_[j - 1], spark.WM2_[j - 2], spark.dmli1_[j - 2], spark.params_, spark.S_, i, j, k, l, tree.tree);
        if (k > i + 1) {
            spark.WM_.swap(temp);
            spark.WM2_.swap(spark.dmli1_);
        }
    }

    trace_WM2(spark, mark_candidates, k, l, tree);
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
void trace_WM(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree) {
    if (debug) printf("WM at %d and %d with %d\n", i, j, e);

    if (i + TURN + 1 > j) {
        return;
    }

    if (e == spark.WM_[j - 1] + spark.params_->MLbase) {
        trace_WM(spark, mark_candidates, i, j - 1, spark.WM_[j - 1], tree);
        return;
    }
    cand_pos_t m = j + 1;
    for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first >= i; ++it) {
        m = it->first;
        const energy_t wmb_kj = it->second + PSM_penalty;
        energy_t wmb_up = static_cast<energy_t>((m - i) * spark.params_->MLbase) + wmb_kj;
        energy_t wmb_wm = spark.WM_[m - 1] + wmb_kj;
        if (e == wmb_up + PSM_penalty) {
            trace_WMB(spark, mark_candidates, m, j, wmb_kj, tree);
            return;
        } else if (e == wmb_wm + PSM_penalty) {
            trace_WM(spark, mark_candidates, i, m - 1, spark.WM_[m - 1], tree);
            trace_WMB(spark, mark_candidates, m, j, wmb_kj, tree);
            return;
        }
    }

    energy_t v = INF;
    energy_t vk = INF;
    Dangle dangle = 3;
    for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first >= i; ++it) {
        m = it->first;
        const energy_t v_kj = it->third >> 2;
        const Dangle d = it->third & 3;
        if (e == spark.WM_[m - 1] + v_kj) {
            dangle = d;
            vk = v_kj;
            v = it->second;
            // no recomp, same i
            break;
        } else if (e == static_cast<energy_t>((m - i) * spark.params_->MLbase) + v_kj) {
            dangle = d;
            vk = v_kj;
            v = it->second;
            break;
        }
    }
    cand_pos_t k = m;
    cand_pos_t l = j;
    pair_type ptype = 0;
    switch (dangle) {
    case 0:
        ptype = pair[spark.S_[k]][spark.S_[l]];
        v = vk - E_MLstem(ptype, -1, -1, spark.params_);
        break;
    case 1:
        k = m + 1;
        ptype = pair[spark.S_[k]][spark.S_[l]];
        v = vk - E_MLstem(ptype, spark.S_[m], -1, spark.params_) - spark.params_->MLbase;
        break;
    case 2:
        l = j - 1;
        ptype = pair[spark.S_[k]][spark.S_[l]];
        v = vk - E_MLstem(ptype, -1, spark.S_[j], spark.params_) - spark.params_->MLbase;
        break;
    case 3:
        if (spark.params_->model_details.dangles == 1) {
            k = m + 1;
            l = j - 1;
            ptype = pair[spark.S_[k]][spark.S_[l]];
            v = vk - E_MLstem(ptype, spark.S_[m], spark.S_[j], spark.params_) - 2 * spark.params_->MLbase;
        }
        break;
    }

    if (e == spark.WM_[m - 1] + vk) {
        // no recomp, same i
        trace_WM(spark, mark_candidates, i, m - 1, spark.WM_[m - 1], tree);
        trace_V(spark, mark_candidates, k, l, v, tree);
        return;
    } else if (e == static_cast<energy_t>((m - i) * spark.params_->MLbase) + vk) {
        trace_V(spark, mark_candidates, k, l, v, tree);
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
void trace_WM2(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, sparse_tree &tree) {
    if (debug) printf("WM2 at %d and %d with %d\n", i, j, spark.WM2_[j]);

    if (i + 2 * TURN + 3 > j) {
        return;
    }
    const energy_t e = spark.WM2_[j];

    // case j unpaired
    if (e == spark.WM2_[j - 1] + spark.params_->MLbase) {
        // same i, no recomputation
        trace_WM2(spark, mark_candidates, i, j - 1, tree);
        return;
    }

    cand_pos_t m = j + 1;
    for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first >= i; ++it) {
        m = it->first;
        const energy_t wmb_kj = it->second + PSM_penalty;
        energy_t wmb_up = static_cast<energy_t>((m - i) * spark.params_->MLbase) + wmb_kj;
        energy_t wmb_wm = spark.WM_[m - 1] + wmb_kj;
        if (e == wmb_up) {
            trace_WMB(spark, mark_candidates, m, j, wmb_kj, tree);
            return;
        } else if (e == wmb_wm) {
            trace_WM(spark, mark_candidates, i, m - 1, spark.WM_[m - 1], tree);
            trace_WMB(spark, mark_candidates, m, j, wmb_kj, tree);
            return;
        }
    }

    energy_t v = INF;
    energy_t vk = INF;
    Dangle dangle = 4;
    for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first >= i + TURN + 1; ++it) {
        m = it->first;

        const energy_t v_kj = it->third >> 2;
        const Dangle d = it->third & 3;
        if (e == spark.WM_[m - 1] + v_kj) {
            vk = v_kj;
            dangle = d;
            v = it->second;
            break;
        }
    }
    cand_pos_t k = m;
    cand_pos_t l = j;
    pair_type ptype = 0;
    switch (dangle) {
    case 0:
        ptype = pair[spark.S_[k]][spark.S_[l]];
        v = vk - E_MLstem(ptype, -1, -1, spark.params_);
        break;
    case 1:
        k = m + 1;
        ptype = pair[spark.S_[k]][spark.S_[l]];
        v = vk - E_MLstem(ptype, spark.S_[m], -1, spark.params_) - spark.params_->MLbase;
        break;
    case 2:
        l = j - 1;
        ptype = pair[spark.S_[k]][spark.S_[l]];
        v = vk - E_MLstem(ptype, -1, spark.S_[j], spark.params_) - spark.params_->MLbase;
        break;
    case 3:
        if (spark.params_->model_details.dangles == 1) {
            k = m + 1;
            l = j - 1;
            ptype = pair[spark.S_[k]][spark.S_[l]];
            v = vk - E_MLstem(ptype, spark.S_[m], spark.S_[j], spark.params_) - 2 * spark.params_->MLbase;
        }
        break;
    }

    if (e == spark.WM_[m - 1] + vk) {
        trace_WM(spark, mark_candidates, i, m - 1, spark.WM_[m - 1], tree);
        trace_V(spark, mark_candidates, k, l, v, tree);
        return;
    }
    assert(false);
}

/**
 * @brief Traceback from WMB
 *
 * @param CLBEO BE candidate list outer->in
 * @param i row index
 * @param j column index
 */
void trace_WMB(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree) {
    if (debug) printf("WMB at i is %d and j is %d and e is %d\n", i, j, e);
    assert(i + TURN + 1 <= j);
    assert(j <= spark.n_);

    recompute_WMBP(spark, i, j, tree);

    cand_pos_t bp_j = tree.tree[j].pair;

    if (tree.tree[j].pair >= 0 && j > tree.tree[j].pair && tree.tree[j].pair > i) {
        energy_t en = INF;
        energy_t BE_energy = INF;
        cand_pos_t best_border = j;
        compute_WMB_case1(j, en, BE_energy, best_border, tree, spark.CLBEO_, spark.WMBA_);
        if (e == en + PB_penalty) {
            trace_BE(spark, mark_candidates, bp_j, tree.tree[best_border].pair, BE_energy, tree);
            trace_WMBA(spark, mark_candidates, i, best_border - 1, en - BE_energy, tree);
        }

        return;
    }
    trace_WMBP(spark, mark_candidates, i, j, spark.WMBP_[j], tree);
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
void trace_VP(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree) {
    if (debug) printf("VP at %d and %d with %d\n", i, j, e);
    spark.structure_[i] = '[';
    spark.structure_[j] = ']';
    if (e == 0) return;

    const pair_type ptype_closing = pair[spark.S_[i]][spark.S_[j]];

    cand_pos_t B_ij = tree.B(i, j);
    cand_pos_t Bp_ij = tree.Bp(i, j);
    cand_pos_t b_ij = tree.b(i, j);
    cand_pos_t bp_ij = tree.bp(i, j);
    if (tree.tree[i].parent->index > 0 && tree.tree[j].parent->index < tree.tree[i].parent->index && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0) {
        recompute_WI(spark, i + 1, Bp_ij - 1, tree.tree, tree.up);
        recompute_WI(spark, B_ij + 1, j - 1, tree.tree, tree.up);
        if (e == spark.WI_[Bp_ij - 1] + spark.WI_[j - 1]) {
            trace_WI(spark, mark_candidates, i + 1, Bp_ij - 1, spark.WI_[Bp_ij - 1], tree);
            trace_WI(spark, mark_candidates, B_ij + 1, j - 1, spark.WI_[j - 1], tree);
            return;
        }
    }
    if (tree.tree[i].parent->index < tree.tree[j].parent->index && tree.tree[j].parent->index > 0 && b_ij >= 0 && bp_ij >= 0 && Bp_ij < 0) {
        recompute_WI(spark, i + 1, b_ij - 1, tree.tree, tree.up);
        recompute_WI(spark, bp_ij + 1, j - 1, tree.tree, tree.up);
        if (e == (spark.WI_[b_ij - 1] + spark.WI_[j - 1])) {
            trace_WI(spark, mark_candidates, i + 1, b_ij - 1, spark.WI_[b_ij - 1], tree);
            trace_WI(spark, mark_candidates, bp_ij + 1, j - 1, spark.WI_[j - 1], tree);
            return;
        }
    }
    if (tree.tree[i].parent->index > 0 && tree.tree[j].parent->index > 0 && Bp_ij >= 0 && B_ij >= 0 && b_ij >= 0 && bp_ij >= 0) {
        recompute_WI(spark, i + 1, Bp_ij - 1, tree.tree, tree.up);
        recompute_WI(spark, B_ij + 1, b_ij - 1, tree.tree, tree.up);
        recompute_WI(spark, bp_ij + 1, j - 1, tree.tree, tree.up);

        if (e == spark.WI_[Bp_ij - 1] + spark.WI_[b_ij - 1] + spark.WI_[j - 1]) {
            trace_WI(spark, mark_candidates, i + 1, Bp_ij - 1, spark.WI_[Bp_ij + 1], tree);
            trace_WI(spark, mark_candidates, B_ij + 1, b_ij - 1, spark.WI_[b_ij - 1], tree);
            trace_WI(spark, mark_candidates, bp_ij + 1, j - 1, spark.WI_[j - 1], tree);
            return;
        }
    }
    if (exists_trace_arrow_from(spark.taVP_, i, j)) {

        const TraceArrow &arrow = trace_arrow_from(spark.taVP_, i, j);
        const size_t k = arrow.k(i);
        const size_t l = arrow.l(j);
        assert(i < k);
        assert(l < j);
        trace_VP(spark, mark_candidates, k, l, arrow.target_energy(), tree);
        return;

    } else {

        // try to trace back to a candidate: (still) interior loop case
        cand_pos_t l_min = std::max(i, j - 31);
        for (cand_pos_t l = j - 1; l > l_min; l--) {
            // Break if it's an assured dangle case
            for (auto it = spark.CLVP_[l].begin(); spark.CLVP_[l].end() != it && it->first > i; ++it) {
                const cand_pos_t k = it->first;

                if (k - i > 31) continue;
                energy_t temp = lrint(((j - l == 1 && k - i == 1) ? e_stP_penalty : e_intP_penalty)
                                      * E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[spark.S_[k]][spark.S_[l]]], spark.S1_[i + 1],
                                                  spark.S1_[j - 1], spark.S1_[k - 1], spark.S1_[l + 1], const_cast<paramT *>(spark.params_)));
                if (e == it->second + temp) {
                    trace_VP(spark, mark_candidates, k, l, it->second, tree);
                    return;
                }
            }
        }
    }
    // 	// If not other cases, must be WV multiloop
    recompute_WVe(spark, i + 1, j - 1, tree);
    recompute_WIP(spark, i + 1, j - 1, tree.tree, tree.up);
    recompute_WV(spark, i + 1, j - 1, tree);

    trace_WV(spark, mark_candidates, i + 1, j - 1, spark.WV_[j - 1], tree);
}

/**
 * @brief Traceback from WVe entry
 *
 * @param CLVP VP candidate list
 * @param i row index
 * @param j column index
 */
void trace_WVe(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree) {
    if (debug) printf("WVe at %d and %d with %d\n", i, j, e);

    if (i + TURN + 1 >= j) return;
    if (spark.WVe_[j] == spark.WVe_[j - 1] + cp_penalty) {
        trace_WVe(spark, mark_candidates, i, j - 1, spark.WVe_[j - 1], tree);
        return;
    }
    cand_pos_t bound_left = j;
    if (tree.b(i, j) > 0) bound_left = tree.b(i, j);
    if (tree.Bp(i, j) > 0) bound_left = std::min((cand_pos_tu)bound_left, (cand_pos_tu)tree.Bp(i, j));
    for (auto it = spark.CLVP_[j].begin(); spark.CLVP_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (k > bound_left) continue;
        if (e == static_cast<energy_t>(cp_penalty * (k - i)) + it->second) {
            trace_VP(spark, mark_candidates, k, j, it->second, tree);
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
void trace_WV(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree) {
    if (debug) printf("WV at %d and %d with %d\n", i, j, e);

    if (spark.WV_[j] == spark.WV_[j - 1] + cp_penalty) {
        trace_WV(spark, mark_candidates, i, j - 1, spark.WV_[j - 1], tree);
        return;
    }
    cand_pos_t bound_left = std::min(tree.b(i, j), tree.Bp(i, j));
    cand_pos_t bound_right = std::min(tree.b(i, j), tree.Bp(i, j));

    for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first > bound_right; ++it) {
        cand_pos_t k = it->first;
        // energy_t wm_v = it->third >> 2;
        energy_t v = it->second;
        // Dangle d = it->third & 3;
        // cand_pos_t num = 0;
        // if(d == 1 || d==2) num = 1;
        // else if(d==3 && params->model_details.dangles == 1) num =2;
        // energy_t fix = num*cp_penalty - num*params->MLbase - b_penalty;

        if (e == spark.WV_[k - 1] + v + bp_penalty) {
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

            trace_WV(spark, mark_candidates, i, k - 1, spark.WV_[k - 1], tree);
            trace_V(spark, mark_candidates, m, l, v, tree);
            return;
        }
        if (e == spark.WVe_[k - 1] + v + bp_penalty) {
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
            trace_WVe(spark, mark_candidates, i, k - 1, spark.WVe_[k - 1], tree);
            trace_V(spark, mark_candidates, m, l, v, tree);
            return;
        }
    }

    for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first > bound_right; ++it) {
        cand_pos_t k = it->first;
        if (e == spark.WV_[k - 1] + it->second + PSM_penalty + bp_penalty) {
            trace_WV(spark, mark_candidates, i, k - 1, spark.WV_[k - 1], tree);
            trace_WMB(spark, mark_candidates, k, j, it->second, tree);
            return;
        }
        if (e == spark.WVe_[k - 1] + it->second + PSM_penalty + bp_penalty) {
            trace_WVe(spark, mark_candidates, i, k - 1, spark.WV_[k - 1], tree);
            trace_WMB(spark, mark_candidates, k, j, it->second, tree);
            return;
        }
    }
    for (auto it = spark.CLVP_[j].begin(); spark.CLVP_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (k > bound_left) continue;
        if (e == spark.WIP_[k - 1] + it->second) {
            trace_WIP(spark, mark_candidates, i, k - 1, spark.WIP_[k - 1], tree);
            trace_VP(spark, mark_candidates, k, j, it->second, tree);
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
void trace_WI(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree) {
    if (debug) printf("WI at %d and %d with %d\n", i, j, e);

    if (i + TURN + 1 >= j) return;

    // How to do one base backwards?

    if (e == spark.WI_[j - 1] + PUP_penalty) {
        trace_WI(spark, mark_candidates, i, j - 1, spark.WI_[j - 1], tree);
        return;
    }

    for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (e == spark.WI_[k - 1] + it->second + PPS_penalty) {
            trace_WI(spark, mark_candidates, i, k - 1, spark.WI_[k - 1], tree);
            trace_V(spark, mark_candidates, k, j, it->second, tree);
            return;
        }
    }

    for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (e == spark.WI_[k - 1] + it->second + PPS_penalty + PSM_penalty) {
            trace_WI(spark, mark_candidates, i, k - 1, spark.WI_[k - 1], tree);
            trace_WMB(spark, mark_candidates, k, j, it->second, tree);
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
void trace_WIP(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree) {
    if (debug) printf("WIP at %d and %d with %d\n", i, j, e);

    if (i + TURN + 1 >= j) return;

    // 	// How to do one base backwards?

    if (e == spark.WIP_[j - 1] + cp_penalty) {
        trace_WIP(spark, mark_candidates, i, j - 1, spark.WIP_[j - 1], tree);
        return;
    }
    for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (e == spark.WIP_[k - 1] + it->second + bp_penalty) {
            trace_WIP(spark, mark_candidates, i, k - 1, spark.WIP_[k - 1], tree);
            trace_V(spark, mark_candidates, k, j, it->second, tree);
            return;
        }
        if (e == static_cast<energy_t>((k - i) * cp_penalty) + it->second + bp_penalty) {
            trace_V(spark, mark_candidates, k, j, it->second, tree);
            return;
        }
    }

    for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (e == spark.WIP_[k - 1] + it->second + bp_penalty + PSM_penalty) {
            // Why do I pick two different variables for WIP
            trace_WIP(spark, mark_candidates, i, k - 1, spark.WIP_Bbp[k - 1], tree);
            trace_WMB(spark, mark_candidates, k, j, it->second, tree);
            return;
        }
        if (e == static_cast<energy_t>((k - i) * cp_penalty) + it->second + bp_penalty + PSM_penalty) {
            trace_WMB(spark, mark_candidates, k, j, it->second, tree);
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
void trace_WMBP(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree) {
    if (debug) printf("WMBP at %d and %d with %d\n", i, j, e);

    energy_t VP_ij = INF;

    if (tree.tree[j].pair < 0) {
        cand_pos_t b_ij = tree.b(i, j);
        for (auto it = spark.CLVP_[j].begin(); spark.CLVP_[j].end() != it && it->first >= i; ++it) {
            cand_pos_t k = it->first;
            if (k == i) VP_ij = it->second; // second?

            int ext_case = compute_exterior_cases(k, j, tree);
            if ((b_ij > 0 && k < b_ij) || (b_ij < 0 && ext_case == 0)) {
                cand_pos_t bp_ik = tree.bp(i, k);
                cand_pos_t Bp_kj = tree.Bp(k, j);
                if (bp_ik >= 0 && k > bp_ik && Bp_kj > 0 && k < Bp_kj) {
                    energy_t BE_energy = INF;
                    cand_pos_t B_kj = tree.B(k, j);
                    cand_pos_t bp_kj = tree.tree[Bp_kj].pair;
                    cand_pos_t b_kj = (B_kj > 0) ? tree.tree[B_kj].pair : -2;
                    for (auto it2 = spark.CLBE_[Bp_kj].begin(); spark.CLBE_[Bp_kj].end() != it2; ++it2) {
                        cand_pos_t l = it2->first;
                        if (l == b_kj) {
                            BE_energy = it2->second;
                            break;
                        }
                    }
                    if (e == spark.WMBA_[k - 1] + it->second + 2 * PB_penalty + BE_energy) {
                        trace_BE(spark, mark_candidates, b_kj, bp_kj, BE_energy, tree);
                        trace_WMBA(spark, mark_candidates, i, k - 1, spark.WMBP_[k - 1], tree);
                        trace_VP(spark, mark_candidates, k, j, it->second, tree);
                        return;
                    }
                }
            }
        }
    }

    if (tree.tree[j].pair < 0 && tree.tree[i].pair >= 0 && tree.tree[i].pair < j) {
        for (auto it = spark.CLVP_[j].begin(); spark.CLVP_[j].end() != it && it->first >= i; ++it) {
            cand_pos_t k = it->first;
            cand_pos_t bp_ik = tree.bp(i, k);
            if (bp_ik >= 0 && k + TURN <= j) {
                cand_pos_t Bp_ik = tree.tree[bp_ik].pair;
                energy_t BE_energy = INF;
                for (auto it2 = spark.CLBE_[Bp_ik].begin(); spark.CLBE_[Bp_ik].end() != it2; ++it2) {
                    cand_pos_t l = it2->first;
                    if (l == i) {
                        BE_energy = it2->second;
                        break;
                    }
                }

                if (e == spark.WI_Bbp[k - 1] + it->second + 2 * PB_penalty + BE_energy) {
                    recompute_WI(spark, bp_ik + 1, k - 1, tree.tree, tree.up);
                    trace_BE(spark, mark_candidates, i, bp_ik, BE_energy, tree);
                    trace_WI(spark, mark_candidates, bp_ik + 1, k - 1, spark.WI_[k - 1], tree);
                    trace_VP(spark, mark_candidates, k, j, it->second, tree);
                    return;
                }
            }
        }
    }

    if (e == VP_ij + PB_penalty) {
        trace_VP(spark, mark_candidates, i, j, VP_ij, tree);
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
void trace_WMBA(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e, sparse_tree &tree) {
    if (debug) printf("WMBA at %d and %d with %d\n", i, j, e);

    if (spark.WMBA_[j] == spark.WMBA_[j - 1] + PUP_penalty) {
        trace_WMBA(spark, mark_candidates, i, j - 1, spark.WMBA_[j - 1], tree);
        return;
    }

    if (spark.WMBA_[j] == spark.WMBP_[j]) {
        trace_WMBP(spark, mark_candidates, i, j, spark.WMBP_[j], tree);
        return;
    }

    for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (e == spark.WMBA_[k - 1] + it->second + PPS_penalty) {
            trace_WMBA(spark, mark_candidates, i, k - 1, spark.WMBP_[k - 1], tree);
            trace_V(spark, mark_candidates, k, j, it->second, tree);
            return;
        }
    }

    for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first >= i; ++it) {
        cand_pos_t k = it->first;
        if (e == spark.WMBA_[k - 1] + it->second + PPS_penalty + PSM_penalty) {
            trace_WMBA(spark, mark_candidates, i, k - 1, spark.WMBP_[k - 1], tree);
            trace_WMB(spark, mark_candidates, k, j, it->second, tree);
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
void trace_BE(Spark &spark, const bool &mark_candidates, cand_pos_t i, cand_pos_t ip, energy_t e, sparse_tree &tree) {
    cand_pos_t j = tree.tree[i].pair;
    cand_pos_t jp = tree.tree[ip].pair;
    cand_pos_t lp = jp;
    // currently this is i lp ip  l j and it should be i lp jp j l ip

    if (debug) printf("BE at [%d,%d] U [%d,%d] with %d\n", i, ip, jp, j, e);

    spark.structure_[i] = '(';
    spark.structure_[j] = ')';
    const pair_type ptype_closing_ij = pair[spark.S_[i]][spark.S_[j]];
    if (i == ip) return;

    energy_t BE_energy = INF;

    for (auto it = spark.CLBE_[jp].begin(); spark.CLBE_[jp].end() != it && it->first > i; ++it) {
        lp = it->first;
        BE_energy = it->second;
    }
    cand_pos_t l = tree.tree[lp].pair;
    // i lp ip       jp  l   j

    if (e == lrint(e_stP_penalty * ILoopE(spark.S_,spark.S1_,spark.params_,ptype_closing_ij,i,j,lp,l)) + BE_energy) {
        trace_BE(spark, mark_candidates, lp, ip, BE_energy, tree);
        return;
    }
    if (e == lrint(e_intP_penalty * ILoopE(spark.S_,spark.S1_,spark.params_,ptype_closing_ij,i,j,lp,l)) + BE_energy) {
        trace_BE(spark, mark_candidates, lp, ip, BE_energy, tree);
        return;
    }

    if (e == spark.WIP_Bbp[lp - 1] + BE_energy + spark.WIP_Bp[l + 1] + ap_penalty + 2 * bp_penalty) {
        recompute_WIP(spark, i + 1, lp - 1, tree.tree, tree.up);
        recompute_WIP(spark, l + 1, j - 1, tree.tree, tree.up);
        trace_WIP(spark, mark_candidates, i + 1, lp - 1, spark.WIP_Bbp[lp - 1], tree);
        trace_BE(spark, mark_candidates, lp, ip, BE_energy, tree);
        trace_WIP(spark, mark_candidates, l + 1, j - 1, spark.WIP_Bp[j - 1], tree);
        return;
    }
    if (e == cp_penalty * ((lp - i - 1)) + BE_energy + spark.WIP_Bp[l + 1] + ap_penalty + 2 * bp_penalty) {
        recompute_WIP(spark, l + 1, j - 1, tree.tree, tree.up);
        trace_BE(spark, mark_candidates, lp, ip, BE_energy, tree);
        trace_WIP(spark, mark_candidates, l + 1, j - 1, spark.WIP_Bbp[j - 1], tree);
        return;
    }
    if (e == spark.WIP_Bbp[lp - 1] + BE_energy + cp_penalty * ((j - l - 1)) + ap_penalty + 2 * bp_penalty) {
        recompute_WIP(spark, i + 1, lp - 1, tree.tree, tree.up);
        trace_WIP(spark, mark_candidates, i + 1, lp - 1, spark.WIP_Bbp[lp - 1], tree);
        trace_BE(spark, mark_candidates, lp, ip, BE_energy, tree);
        return;
    }
}
/**
 * @brief Trace back
 * pre: row 1 of matrix W is computed
 * @return mfe structure (reference)
 */
const std::string &trace_back(Spark &spark, sparse_tree &tree, const bool &mark_candidates = false) {

    spark.structure_.resize(spark.n_ + 1, '.');

    /* Traceback */
    trace_W(spark, mark_candidates, 1, spark.n_, tree);
    spark.structure_ = spark.structure_.substr(1, spark.n_);

    return spark.structure_;
}

/**
 * @brief Register a candidate
 * @param i start
 * @param j end
 * @param e energy of candidate "V(i,j)"
 * @param wmij energy at WM(i,j)
 * @param wij energy at W(i,j)
 */
void register_candidate(std::vector<cand_list_td1> &CL, cand_pos_t const &i, cand_pos_t const &j, energy_t const &e, energy_t const &wmij,
                        energy_t const &wij) {
    assert(i <= j + TURN + 1);

    CL[j].emplace_back(cand_entry_td1(i, e, wmij, wij));
}
/**
 * @brief Register a candidate
 * @param i start
 * @param j end
 * @param e energy of candidate "V(i,j)"
 */
void register_candidate(std::vector<cand_list_t> &CL, cand_pos_t const &i, cand_pos_t const &j, energy_t const &e) {
    assert(i <= j + TURN + 1);

    CL[j].emplace_back(cand_entry_t(i, e));
}

/**
 * @brief Computes the values for the WVe matrix
 * @param i start
 * @param j end
 * @param spark datastructure
 * @param tree tree for boundary determination
 */
energy_t compute_WVe(cand_pos_t i, cand_pos_t j, Spark &spark, sparse_tree &tree) {

    energy_t wve = INF;
    cand_pos_t bound = std::min((cand_pos_tu)tree.Bp(i, j), (cand_pos_tu)tree.b(i, j));
    for (auto it = spark.CLVP_[j].begin(); spark.CLVP_[j].end() != it; ++it) {
        cand_pos_t k = it->first;
        bool can_pair = tree.up[k - 1] >= (k - i);
        if (can_pair && k < bound) wve = std::min(wve, static_cast<energy_t>(cp_penalty * (k - i)) + it->second);
    }
    if (tree.tree[j].pair < 0) wve = std::min(wve, spark.WVe_[j - 1] + cp_penalty);

    return wve;
}

/**
 * @brief Computes the values for the WV matrix
 * @param j end
 * @param bound_left boundary for non-closed region
 * @param bound_right boundary for closed region
 * @param spark datastructure
 * @param tree tree for boundary determination
 */
energy_t compute_WV(cand_pos_t j, cand_pos_t bound_left, cand_pos_t bound_right, Spark &spark, sparse_tree &tree) {
    energy_t m1 = INF, m2 = INF, m3 = INF, m5 = INF, m6 = INF, wv = INF;

    for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it && it->first > bound_right; ++it) {
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
        m1 = std::min(m1, spark.WVe_[k - 1] + val + bp_penalty);
        m5 = std::min(m5, spark.WV_[k - 1] + val + bp_penalty);
    }
    for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it && it->first > bound_right; ++it) {
        cand_pos_t k = it->first;
        energy_t val = it->second;
        m2 = std::min(m1, spark.WVe_[k - 1] + val + PSM_penalty + bp_penalty);
        m6 = std::min(m6, spark.WV_[k - 1] + val + PSM_penalty + bp_penalty);
    }

    for (auto it = spark.CLVP_[j].begin(); spark.CLVP_[j].end() != it; ++it) {
        cand_pos_t k = it->first;
        energy_t val = it->second;
        if (k < bound_left) m3 = std::min(m1, spark.dwip1_[k - 1] + val);
    }
    wv = std::min({m1, m2, m3, m5, m6});

    if (tree.tree[j].pair < 0) wv = std::min(wv, spark.WV_[j - 1] + cp_penalty);

    return wv;
}

/**
 * @brief computes entries for BE
 * @param i outer opening base pair
 * @param j inner closing base pair
 * @param ip outer closing base pair
 * @param jp inner opening base pair
 * @param spark datastructure
 * @param tree tree for boundary estimation
 */
energy_t compute_BE(cand_pos_t i, cand_pos_t j, cand_pos_t ip, cand_pos_t jp, Spark &spark, sparse_tree &tree) {
    // We are checking for the closest pair that we have already calculated to i/ip from j/jp
    // If there is nothing, then i is the closest encompassing pair to jp
    // If it is not, then we get the energy for everything from jp to lp so that we calculate less

    // (.....(..(....)..).....)
    // i     lp      j  l     ip
    energy_t BE_energy = INF;
    cand_pos_t lp = jp;
    if (!spark.CLBE_[j].empty()) {
        auto const [k, vbe] = spark.CLBE_[j].back();
        BE_energy = vbe;
        lp = k;
    }
    cand_pos_t l = tree.tree[lp].pair; // right closing base for lp

    const pair_type ptype_closing_iip = pair[spark.S_[i]][spark.S_[ip]];

    energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF, val = INF;
    // 1
    if (i + 1 == lp && ip - 1 == l) {
        m1 = lrint(e_stP_penalty * ILoopE(spark.S_, spark.S1_, spark.params_, ptype_closing_iip, i, j, lp, l)) + BE_energy;
        val = std::min(val, m1);
    }

    bool empty_region_ilp = (tree.up[lp - 1] >= lp - i - 1);    // empty between i+1 and lp-1
    bool empty_region_lip = (tree.up[ip - 1] >= ip - l - 1);    // empty between l+1 and ip-1
    bool weakly_closed_ilp = tree.weakly_closed(i + 1, lp - 1); // weakly closed between i+1 and lp-1
    bool weakly_closed_lip = tree.weakly_closed(l + 1, ip - 1); // weakly closed between l+1 and ip-1

    // 2
    if (empty_region_ilp && empty_region_lip) {
        m2 = lrint(e_intP_penalty * ILoopE(spark.S_, spark.S1_, spark.params_, ptype_closing_iip, i, j, lp, l)) + BE_energy;
        val = std::min(val, m2);
    }

    // 3
    if (weakly_closed_ilp && weakly_closed_lip) {
        m3 = spark.dwip1_[lp - 1] + BE_energy + spark.WIP_Bp[l + 1] + ap_penalty + 2 * bp_penalty;
        val = std::min(val, m3);
    }

    // 4
    if (weakly_closed_ilp && empty_region_lip) {
        m4 = spark.dwip1_[lp - 1] + BE_energy + cp_penalty * (ip - l - 1) + ap_penalty + 2 * bp_penalty;
        val = std::min(val, m4);
    }

    // 5
    if (empty_region_ilp && weakly_closed_lip) {

        m5 = ap_penalty + 2 * bp_penalty + (cp_penalty * (lp - i - 1)) + BE_energy + spark.WIP_Bp[l + 1];
        val = std::min(val, m5);
    }

    return val;
}
/**
 * @brief Computes entries for WMBP
 * @param i start
 * @param j end
 * @param spark datastructure
 * @param tree tree for boundary estimation
 */
energy_t compute_WMBP(cand_pos_t i, cand_pos_t j, Spark &spark, sparse_tree &tree) {
    energy_t m1 = INF, m2 = INF, m3 = INF, wmbp = INF;
    // 1) WMBP(i,j) = BE(bpg(Bp(l,j)),Bp(l,j),bpg(B(l,j)),B(l,j)) + WMBP(i,l) + VP(l+1,j)
    if (tree.tree[j].pair < 0) {
        energy_t tmp = INF;
        cand_pos_t b_ij = tree.b(i, j);
        for (auto it = spark.CLVP_[j].begin(); spark.CLVP_[j].end() != it; ++it) {
            cand_pos_t k = it->first;
            // Mateo Jan 2025 Added exterior cases to consider when looking at band borders. Solved case of [.(.].[.).]
            int ext_case = compute_exterior_cases(k, j, tree);
            if ((b_ij > 0 && k < b_ij) || (b_ij < 0 && ext_case == 0)) {
                cand_pos_t bp_ik = tree.bp(i, k);
                cand_pos_t Bp_kj = tree.Bp(k, j);
                if (bp_ik >= 0 && k > bp_ik && Bp_kj > 0 && k < Bp_kj) { // if(sparse_tree.b(i,j)>=0 && l <sparse_tree.b(i,j)){//
                    cand_pos_t B_kj = tree.B(k, j);
                    if (i <= tree.tree[k].parent->index && tree.tree[k].parent->index < j && k + 3 <= j) {
                        energy_t BE_energy = INF;
                        cand_pos_t b_kj = tree.tree[B_kj].pair;
                        for (auto it2 = spark.CLBE_[Bp_kj].begin(); spark.CLBE_[Bp_kj].end() != it2; ++it2) {
                            cand_pos_t l = it2->first;
                            if (l == b_kj) {
                                BE_energy = it2->second;
                                break;
                            }
                        }
                        energy_t WMBA_energy = spark.WMBA_[k - 1];
                        energy_t VP_energy = it->second;
                        energy_t sum = BE_energy + WMBA_energy + VP_energy;

                        tmp = std::min(tmp, sum);
                    }
                }
            }

            m1 = 2 * PB_penalty + tmp;
        }
    }

    // 2) WMBP(i,j) = VP(i,j) + P_b
    cand_pos_t i_mod = i % (MAXLOOP + 1);
    m2 = spark.VP_(i_mod, j) + PB_penalty;

    // check later if <0 or <-1

    // WMBP(i,j) = BE(i,,,) _ WI(bp(i,k),k-1) + VP(k,j)
    if (tree.tree[j].pair < 0 && tree.tree[i].pair >= 0 && tree.tree[i].pair < j) {
        energy_t tmp = INF;
        for (auto it = spark.CLVP_[j].begin(); spark.CLVP_[j].end() != it; ++it) {
            cand_pos_t k = it->first;
            cand_pos_t bp_ik = tree.bp(i, k);
            if (bp_ik >= 0 && k + TURN <= j) {
                energy_t BE_energy = INF;
                cand_pos_t Bp_ik = tree.tree[bp_ik].pair;
                if (!spark.CLBE_[Bp_ik].empty()) {
                    auto const [l, vbe] = spark.CLBE_[Bp_ik].back();
                    if (i == l) BE_energy = vbe;
                }

                energy_t WI_energy = (k - 1 - (bp_ik + 1)) > 4 ? spark.WI_Bbp[k - 1] : PUP_penalty * (k - 1 - (bp_ik + 1) + 1);
                energy_t VP_energy = it->second;
                energy_t sum = BE_energy + WI_energy + VP_energy;

                tmp = std::min(tmp, sum);
            }
        }

        m3 = 2 * PB_penalty + tmp;
    }
    // get the min for WMB
    wmbp = std::min({m1, m2, m3});

    return (wmbp);
}
/**
 * @brief Computes entries for WMBA
 * @param j end
 * @param spark datastructure
 * @param tree tree for boundary estimation
 */
energy_t compute_WMBA(cand_pos_t j, Spark &spark, sparse_tree &tree) {

    // WMBA criteria
    energy_t wmba = INF;
    if (tree.tree[j].parent->index > 0) {

        for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it; ++it) {
            cand_pos_t k = it->first;
            if (tree.tree[k].pair < 0 && tree.tree[k].parent->index > -1 && tree.tree[j].parent->index > -1
                && tree.tree[j].parent->index == tree.tree[k].parent->index) {
                wmba = std::min(wmba, spark.WMBA_[k - 1] + it->second + PPS_penalty);
            }
        }

        for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it; ++it) {
            cand_pos_t k = it->first;
            if (tree.tree[k].pair < 0 && tree.tree[k].parent->index > -1 && tree.tree[j].parent->index > -1
                && tree.tree[j].parent->index == tree.tree[k].parent->index) {
                wmba = std::min(wmba, spark.WMBA_[k - 1] + it->second + PPS_penalty + PSM_penalty);
            }
        }
        if (tree.tree[j].pair < 0) wmba = std::min(wmba, spark.WMBA_[j - 1] + PUP_penalty);
    } else {
        wmba = INF;
    }
    wmba = std::min(wmba, spark.WMBP_[j]);
    return wmba;
}
/**
 * @brief Computes entries for WMB
 * @param i start
 * @param j end
 * @param spark datastructure
 * @param tree tree for boundary estimation
 */
energy_t compute_WMB(cand_pos_t i, cand_pos_t j, Spark &spark, sparse_tree &tree) {
    energy_t m1 = INF, m2 = INF, wmb = INF;

    // 2)

    if (tree.tree[j].pair >= 0 && j > tree.tree[j].pair && tree.tree[j].pair > i) {
        cand_pos_t best_border = j - 1;
        energy_t BE_energy = INF;
        compute_WMB_case1(j, m1, BE_energy, best_border, tree, spark.CLBEO_, spark.WMBA_);
        m1 = PB_penalty + m1;
    }
    // check the WMBP_ij value
    m2 = spark.WMBP_[j];

    wmb = std::min(m1, m2);
    return wmb;
}

/**
 * @brief Computes VP internal loop value for case 5 of VP
 *
 * @param VP VP array
 * @param i row index
 * @param j column index
 */
energy_t compute_VP_internal(cand_pos_t i, cand_pos_t j, cand_pos_t b_ij, cand_pos_t bp_ij, cand_pos_t Bp_ij, cand_pos_t B_ij, cand_pos_t &best_k,
                             cand_pos_t &best_l, energy_t &best_e, sparse_tree &sparse_tree, short *S, short *S1, LocARNA::Matrix<energy_t> &VP,
                             paramT *params) {

    energy_t m5 = INF;
    // By doing uint, we make sure it can't be negative
    cand_pos_t min_borders = std::min((cand_pos_tu)Bp_ij, (cand_pos_tu)b_ij);
    cand_pos_t edge_i = std::min(i + MAXLOOP + 1, j - TURN - 1);
    min_borders = std::min({min_borders, edge_i});
    const pair_type ptype_closing = pair[S[i]][S[j]];
    for (cand_pos_t k = i + 1; k <= min_borders; k++) {
        cand_pos_t k_mod = k % (MAXLOOP + 1);

        energy_t cank = ((sparse_tree.up[k - 1] >= (k - i - 1)) - 1);
        cand_pos_t max_borders = std::max(bp_ij, B_ij) + 1;
        cand_pos_t edge_j = k + j - i - MAXLOOP - 2;
        max_borders = std::max({max_borders, edge_j});

        for (cand_pos_t l = j - 1; l >= max_borders; --l) {
            assert(k - i + j - l - 2 <= MAXLOOP);

            energy_t canl = (((sparse_tree.up[j - 1] >= (j - l - 1)) - 1) | cank);
            energy_t v_iloop_kl = INF & canl;
            v_iloop_kl = v_iloop_kl + VP(k_mod, l)
                         + lrint(e_intP_penalty
                                 * E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S[k]][S[l]]], S1[i + 1], S1[j - 1], S1[k - 1], S1[l + 1],
                                             const_cast<paramT *>(params)));

            if (v_iloop_kl < m5) {
                m5 = v_iloop_kl;
                best_l = l;
                best_k = k;
                best_e = VP(k_mod, l);
            }
        }
    }
    return m5;
}

/**
 * @brief Computes the energy values for VP
 * @param i row index
 * @param j column index
 * @param b_ij left outer boundary
 * @param bp_ij left inner boundary
 * @param BP_ij right inner boundary
 * @param B_ij right outer boundary
 * @param spark datastructure
 * @param tree tree for boundary determination
 */
energy_t compute_VP(cand_pos_t i, cand_pos_t j, cand_pos_t b_ij, cand_pos_t bp_ij, cand_pos_t Bp_ij, cand_pos_t B_ij, Spark &spark, sparse_tree &tree) {
    energy_t m1 = INF, m2 = INF, m3 = INF, m4 = INF, m5 = INF, m6 = INF, vp = INF;
    const pair_type ptype_closing = pair[spark.S_[i]][spark.S_[j]];

    if (tree.tree[i].parent->index > 0 && tree.tree[j].parent->index < tree.tree[i].parent->index && Bp_ij >= 0 && B_ij >= 0 && bp_ij < 0) {
        energy_t WI_ipus1_BPminus = spark.dwi1_[Bp_ij - 1];
        energy_t WI_Bplus_jminus = (j - 1 - (B_ij + 1)) > 4 ? spark.WI_Bbp[j - 1] : PUP_penalty * (j - 1 - (B_ij + 1) + 1);

        m1 = WI_ipus1_BPminus + WI_Bplus_jminus;
    }
    if (tree.tree[i].parent->index < tree.tree[j].parent->index && tree.tree[j].parent->index > 0 && b_ij >= 0 && bp_ij >= 0 && Bp_ij < 0) {
        energy_t WI_i_plus_b_minus = spark.dwi1_[b_ij - 1];
        energy_t WI_bp_plus_j_minus = (j - 1 - (bp_ij + 1)) > 4 ? spark.WI_Bbp[j - 1] : PUP_penalty * (j - 1 - (bp_ij + 1) + 1);

        m2 = WI_i_plus_b_minus + WI_bp_plus_j_minus;
    }

    if (tree.tree[i].parent->index > 0 && tree.tree[j].parent->index > 0 && Bp_ij >= 0 && B_ij >= 0 && b_ij >= 0 && bp_ij >= 0) {
        energy_t WI_i_plus_Bp_minus = spark.dwi1_[Bp_ij - 1];
        energy_t WI_B_plus_b_minus = (b_ij - 1 - (B_ij + 1)) > 4 ? spark.WI_Bbp[b_ij - 1] : PUP_penalty * (b_ij - 1 - (B_ij + 1) + 1);
        energy_t WI_bp_plus_j_minus = (j - 1 - (bp_ij + 1)) > 4 ? spark.WI_Bbp[j - 1] : PUP_penalty * (j - 1 - (bp_ij + 1) + 1);

        m3 = WI_i_plus_Bp_minus + WI_B_plus_b_minus + WI_bp_plus_j_minus;
    }
    if (tree.tree[i + 1].pair < -1 && tree.tree[j - 1].pair < -1) {
        cand_pos_t ip1_mod = (i + 1) % (MAXLOOP + 1);

        m4 = lrint(e_stP_penalty * ILoopE(spark.S_, spark.S1_, spark.params_, ptype_closing, i, j, i + 1, j - 1)) + spark.VP_(ip1_mod, j - 1);
    }

    cand_pos_t best_k = 0;
    cand_pos_t best_l = 0;
    energy_t best_e = 0;

    m5 = compute_VP_internal(i, j, b_ij, bp_ij, Bp_ij, B_ij, best_k, best_l, best_e, tree, spark.S_, spark.S1_, spark.VP_, spark.params_);

    // case 6 and 7
    m6 = spark.dwvp_[j - 1] + ap_penalty + 2 * bp_penalty;

    energy_t vp_h = std::min({m1, m2, m3});
    energy_t vp_iloop = std::min(m4, m5);
    if (m4 < m5) {
        best_k = i + 1;
        best_l = j - 1;
        cand_pos_t ip1_mod = (i + 1) % (MAXLOOP + 1);
        best_e = spark.VP_(ip1_mod, j - 1);
    }
    energy_t vp_split = m6;
    vp = std::min({vp_h, vp_iloop, vp_split});

    if (vp_iloop < std::min(vp_h, vp_split)) {
        if (is_candidate(spark.CLVP_, spark.cand_comp, best_k, best_l)) {
            avoid_trace_arrow(spark.taVP_);
        } else {
            register_trace_arrow(spark.taVP_, i, j, best_k, best_l, best_e);
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
energy_t compute_internal(cand_pos_t i, cand_pos_t j, cand_pos_t &best_k, cand_pos_t &best_l, energy_t &best_e, sparse_tree &sparse_tree, short *S,
                          short *S1, LocARNA::Matrix<energy_t> &V, paramT *params) {
    energy_t v_iloop = INF;
    cand_pos_t max_k = std::min(j - TURN - 2, i + MAXLOOP + 1);
    const pair_type ptype_closing = pair[S[i]][S[j]];
    for (cand_pos_t k = i + 1; k <= max_k; k++) {
        cand_pos_t k_mod = k % (MAXLOOP + 1);

        energy_t cank = ((sparse_tree.up[k - 1] >= (k - i - 1)) - 1);
        cand_pos_t min_l = std::max(k + TURN + 1 + MAXLOOP + 2, k + j - i) - MAXLOOP - 2;
        // cand_pos_t ind = k_mod*V.ydim_;
        for (cand_pos_t l = j - 1; l >= min_l; --l) {
            assert(k - i + j - l - 2 <= MAXLOOP);
            energy_t canl = (((sparse_tree.up[j - 1] >= (j - l - 1)) - 1) | cank);
            energy_t v_iloop_kl = INF & canl;

            v_iloop_kl = v_iloop_kl + V(k_mod, l)
                         + E_IntLoop(k - i - 1, j - l - 1, ptype_closing, rtype[pair[S[k]][S[l]]], S1[i + 1], S1[j - 1], S1[k - 1], S1[l + 1],
                                     const_cast<paramT *>(params));
            if (v_iloop_kl < v_iloop) {
                v_iloop = v_iloop_kl;
                best_l = l;
                best_k = k;
                best_e = V(k_mod, l);
            }
        }
    }
    return v_iloop;
}

/**
 * @brief Determines the MFE energy for a given sequence
 * @param spark datastructure
 * @param tree tree for boundaary determination
 * @param n length of sequence
 * @param garbage_collect whether to garbage collect for trace_arrows
 */
energy_t fold(Spark &spark, sparse_tree &tree, const cand_pos_t n, const bool garbage_collect) {
    Dangle d = 3;
    if (spark.params_->model_details.dangles == 0 || spark.params_->model_details.dangles == 1) d = 0;

    for (cand_pos_t i = n; i > 0; --i) {
        if (pseudoknot) {
            for (cand_pos_t j = i; j <= n && tree.tree[j].pair < 0; ++j) {
                spark.WI_[j] = (j - i + 1) * PUP_penalty;
            }
        }

        for (cand_pos_t j = i + TURN + 1; j <= n; j++) {

            bool evaluate = tree.weakly_closed(i, j);
            // ------------------------------
            // W: split case
            bool pairedkj = 0;
            energy_t w_split = INF;
            energy_t wm_split = INF;
            energy_t wm2_split = INF;
            energy_t wi_split = INF;
            energy_t wip_split = INF;
            for (auto it = spark.CL_[j].begin(); spark.CL_[j].end() != it; ++it) {
                cand_pos_t k = it->first;

                const energy_t v_kj = it->third >> 2;
                const energy_t v_kjw = it->fourth >> 2;
                bool can_pair = tree.up[k - 1] >= (k - i);
                // WM Portion
                wm_split = std::min(wm_split, spark.WM_[k - 1] + v_kj);
                if (can_pair) wm_split = std::min(wm_split, static_cast<energy_t>((k - i) * spark.params_->MLbase) + v_kj);
                // WM2 Portion
                wm2_split = std::min(wm2_split, spark.WM_[k - 1] + v_kj);
                // W Portion
                w_split = std::min(w_split, spark.W_[k - 1] + v_kjw);

                // WI portion
                energy_t v_kjj = it->second + PPS_penalty;
                wi_split = std::min(wi_split, spark.WI_[k - 1] + v_kjj);
                // WIP portion
                v_kjj = it->second + bp_penalty;
                wip_split = std::min(wip_split, spark.WIP_[k - 1] + v_kjj);
                if (can_pair) wip_split = std::min(wip_split, static_cast<energy_t>((k - i) * cp_penalty) + v_kjj);
            }

            if (tree.weakly_closed(i, j)) {
                for (auto it = spark.CLWMB_[j].begin(); spark.CLWMB_[j].end() != it; ++it) {

                    if (pairedkj) break; // Not needed I believe as there shouldn't be any candidates there if paired anyways
                    // Maybe this would just avoid this loop however

                    cand_pos_t k = it->first;
                    bool can_pair = tree.up[k - 1] >= (k - i);

                    // For W
                    energy_t wmb_kj = it->second + PS_penalty;
                    w_split = std::min(w_split, spark.W_[k - 1] + wmb_kj);
                    // For WM -> I believe this would add a PSM penalty for every pseudoknot which would be bad
                    wmb_kj = it->second + PSM_penalty + b_penalty;
                    wm_split = std::min(wm_split, spark.WM_[k - 1] + wmb_kj);
                    if (can_pair) wm_split = std::min(wm_split, static_cast<energy_t>((k - i) * spark.params_->MLbase) + wmb_kj);
                    wm2_split = std::min(wm2_split, spark.WM_[k - 1] + wmb_kj);
                    if (can_pair) wm2_split = std::min(wm2_split, static_cast<energy_t>((k - i) * spark.params_->MLbase) + wmb_kj);
                    // For WI
                    wmb_kj = it->second + PSM_penalty + PPS_penalty;
                    wi_split = std::min(wi_split, spark.WI_[k - 1] + wmb_kj);

                    // For WIP
                    wmb_kj = it->second + PSM_penalty + bp_penalty;
                    wip_split = std::min(wip_split, spark.WIP_[k - 1] + wmb_kj);
                    if (can_pair) wip_split = std::min(wip_split, static_cast<energy_t>((k - i) * cp_penalty) + wmb_kj);
                }
            }

            if (tree.tree[j].pair < 0) w_split = std::min(w_split, spark.W_[j - 1]);
            if (tree.tree[j].pair < 0) wm2_split = std::min(wm2_split, spark.WM2_[j - 1] + spark.params_->MLbase);
            if (tree.tree[j].pair < 0) wm_split = std::min(wm_split, spark.WM_[j - 1] + spark.params_->MLbase);
            if (tree.tree[j].pair < 0) wi_split = std::min(wi_split, spark.WI_[j - 1] + PUP_penalty);
            if (tree.tree[j].pair < 0) wip_split = std::min(wip_split, spark.WIP_[j - 1] + cp_penalty);

            energy_t w = w_split;   // entry of W w/o contribution of V
            energy_t wm = wm_split; // entry of WM w/o contribution of V

            size_t i_mod = i % (MAXLOOP + 1);

            const pair_type ptype_closing = pair[spark.S_[i]][spark.S_[j]];
            const bool restricted = tree.tree[i].pair == -1 || tree.tree[j].pair == -1;

            const bool unpaired = (tree.tree[i].pair < -1 && tree.tree[j].pair < -1);
            const bool paired = (tree.tree[i].pair == j && tree.tree[j].pair == i);
            const bool pkonly = (!pk_only || paired);
            energy_t v = INF;
            // ----------------------------------------
            // cases with base pair (i,j)
            // if(ptype_closing>0 && !restricted && evaluate) { // if i,j form a canonical base pair
            if (ptype_closing > 0 && !restricted && evaluate && pkonly) {
                bool canH = (paired || unpaired);
                if (tree.up[j - 1] < (j - i - 1)) canH = false;

                energy_t v_h = canH ? HairpinE(spark.seq_, spark.S_, spark.S1_, spark.params_, i, j) : INF;
                // info of best interior loop decomposition (if better than hairpin)
                cand_pos_t best_l = 0;
                cand_pos_t best_k = 0;
                energy_t best_e = INF;

                energy_t v_iloop = INF;

                // constraints for interior loops
                // i<k; l<j
                // k-i+j-l-2<=MAXLOOP  ==> k <= MAXLOOP+i+1
                //            ==> l >= k+j-i-MAXLOOP-2
                // l-k>=TURN+1         ==> k <= j-TURN-2
                //            ==> l >= k+TURN+1
                // j-i>=TURN+3
                //
                if ((tree.tree[i].pair < -1 && tree.tree[j].pair < -1) || tree.tree[i].pair == j) {
                    v_iloop = compute_internal(i, j, best_k, best_l, best_e, tree, spark.S_, spark.S1_, spark.V_, spark.params_);
                }
                const energy_t v_split = E_MbLoop(spark.dmli1_, spark.dmli2_, spark.S_, spark.params_, i, j, tree.tree);

                v = std::min(v_h, std::min(v_iloop, v_split));
                // register required trace arrows from (i,j)
                if (v_iloop < std::min(v_h, v_split)) {
                    if (is_candidate(spark.CL_, spark.cand_comp, best_k, best_l)) {
                        avoid_trace_arrow(spark.ta_);
                    } else {
                        register_trace_arrow(spark.ta_, i, j, best_k, best_l, best_e);
                    }
                }

                spark.V_(i_mod, j) = v;
            } else {
                spark.V_(i_mod, j) = INF;
            } // end if (i,j form a canonical base pair)

            cand_pos_t ip1_mod = (i + 1) % (MAXLOOP + 1);
            energy_t vi1j = spark.V_(ip1_mod, j);
            energy_t vij1 = spark.V_(i_mod, j - 1);
            energy_t vi1j1 = spark.V_(ip1_mod, j - 1);

            // Checking the dangle positions for W
            energy_t w_v = E_ext_Stem(v, vi1j, vij1, vi1j1, spark.S_, spark.params_, i, j, d, n, tree.tree);
            // Checking the dangle positions for W
            const energy_t wm_v = E_MLStem(v, vi1j, vij1, vi1j1, spark.S_, spark.params_, i, j, d, n, tree.tree);

            cand_pos_t k = i;
            cand_pos_t l = j;
            if (spark.params_->model_details.dangles == 1) {
                if (d > 0) {
                    switch (d) {
                    case 1:
                        k = i + 1;
                        break;
                    case 2:
                        l = j - 1;
                        break;
                    case 3:
                        k = i + 1;
                        l = j - 1;
                        break;
                    }
                    if (exists_trace_arrow_from(spark.ta_, k, l) && (wm_v < wm_split || w_v < w_split)) inc_source_ref_count(spark.ta_, k, l);
                }
            }
            energy_t wi_v = INF;
            energy_t wip_v = INF;
            energy_t wi_wmb = INF;
            energy_t wip_wmb = INF;
            energy_t w_wmb = INF, wm_wmb = INF;
            if (pseudoknot) {
                cand_pos_t Bp_ij = tree.Bp(i, j);
                cand_pos_t B_ij = tree.B(i, j);
                cand_pos_t b_ij = tree.b(i, j);
                cand_pos_t bp_ij = tree.bp(i, j);

                // Start of VP ---- Will have to change the bounds to 1 to n instead of 0 to n-1
                bool weakly_closed_ij = tree.weakly_closed(i, j);
                if (weakly_closed_ij || tree.tree[i].pair >= -1 || tree.tree[j].pair >= -1 || ptype_closing == 0) {

                    spark.VP_(i_mod, j) = INF;

                } else {
                    const energy_t vp = compute_VP(i, j, b_ij, bp_ij, Bp_ij, B_ij, spark, tree);

                    spark.VP_(i_mod, j) = vp;
                }

                // -------------------------------------------End of VP----------------------------------------------------------------

                // Start of WMBP
                if ((tree.tree[i].pair >= -1 && tree.tree[i].pair > j) || (tree.tree[j].pair >= -1 && tree.tree[j].pair < i)
                    || (tree.tree[i].pair >= -1 && tree.tree[i].pair < i) || (tree.tree[j].pair >= -1 && j < tree.tree[j].pair)) {
                    spark.WMB_[j] = INF;
                    spark.WMBP_[j] = INF;
                    spark.WMBA_[j] = INF;

                } else {

                    const energy_t wmbp = compute_WMBP(i, j, spark, tree);
                    spark.WMBP_[j] = wmbp;

                    const energy_t wmba = compute_WMBA(j, spark, tree);
                    spark.WMBA_[j] = wmba;

                    const energy_t wmb = compute_WMB(i, j, spark, tree);
                    spark.WMB_[j] = wmb;
                }

                // -------------------------------------------------------End of WMB------------------------------------------------------

                // Start of WI -- the conditions on calculating WI is the same as WIP, so we combine them

                if (!weakly_closed_ij) {
                    spark.WI_[j] = INF;
                    spark.WIP_[j] = INF;
                } else {

                    wi_v = spark.V_(i_mod, j) + PPS_penalty;
                    wip_v = spark.V_(i_mod, j) + bp_penalty;

                    wi_wmb = spark.WMB_[j] + PSM_penalty + PPS_penalty;
                    wip_wmb = spark.WMB_[j] + PSM_penalty + bp_penalty;

                    spark.WI_[j] = std::min({wi_v, wi_wmb, wi_split});
                    spark.WIP_[j] = std::min({wip_v, wip_wmb, wip_split});

                    if ((tree.tree[i - 1].pair > (i - 1) && tree.tree[i - 1].pair > j) || tree.tree[i - 1].pair < (i - 1)) {
                        spark.WI_Bbp[j] = spark.WI_[j];
                        spark.WIP_Bbp[j] = spark.WIP_[j];
                    }
                    if (j + 1 < n && tree.tree[j + 1].pair < i) {
                        spark.WIP_Bp[i] = spark.WIP_[j];
                    }
                }

                // ------------------------------------------------End of Wi/Wip--------------------------------------------------

                // start of WV and WVe
                if (!weakly_closed_ij) {
                    cand_pos_t bound_right = std::max(bp_ij, B_ij);
                    // Since bound_left is a min, by doing an uint, I make it so that it's a large number if negative (i.e. not possible)
                    cand_pos_t bound_left = std::min((cand_pos_tu)tree.Bp(i, j), (cand_pos_tu)tree.b(i, j));

                    const energy_t wve = compute_WVe(i, j, spark, tree);
                    const energy_t wv = compute_WV(j, bound_left, bound_right, spark, tree);

                    spark.WVe_[j] = wve;

                    spark.WV_[j] = wv;
                } else {
                    spark.WV_[j] = INF;
                    spark.WVe_[j] = INF;
                }

                // ------------------------------------------------End of WV------------------------------------------------------

                /*
                The order for these should be          i        lp               jp         j           l      ip
                                                     ( (        (    (   (       (          )      )  ) )      ) )
                i and ip are the outer base pair
                lp and l are the closest encompassing base pair to i/ip
                                                                                                //    jp>i   j>jp  ip>i ip>j
                jp and j are some inner base pair;j has the be the closing due to the j=i+4 setup we have
                */
                // Start of BE
                cand_pos_t ip = tree.tree[i].pair; // i's pair ip should be right side so ip = )
                cand_pos_t jp = tree.tree[j].pair; // j's pair jp should be left side so jp = (

                // base case: i.j and ip.jp must be in G

                if (jp > i && j > jp && ip > j && ip > i) { // Don't need to check if they are pairs separately because it is checked this
                    if (tree.tree[jp + 1].pair == j - 1) {
                        // BE_avoided++;
                    } else {
                        energy_t BE = compute_BE(i, j, ip, jp, spark, tree);
                        register_candidate(spark.CLBE_, i, j, BE);
                        register_candidate(spark.CLBEO_, j, i, BE);
                    }

                } else if (i == jp && ip == j) {
                    if (tree.tree[jp + 1].pair == j - 1) { // huh?
                        // BE_avoided++;
                    } else {
                        register_candidate(spark.CLBE_, i, j, 0);
                        register_candidate(spark.CLBEO_, j, i, 0);
                    }
                }

                // // ------------------------------------------------End of BE---------------------------------------------------------
            }

            energy_t vp_min1 = INF;
            energy_t vp_min2 = INF;
            for (auto it = spark.CLVP_[j].begin(); spark.CLVP_[j].end() != it; ++it) {
                const cand_pos_t k = it->first;
                bool can_pair = tree.up[k - 1] >= (k - i);
                energy_t WIk = spark.WI_[k - 1];
                energy_t WIPk = spark.WIP_[k - 1];
                vp_min1 = std::min(vp_min1, WIk + it->second);
                vp_min2 = std::min(vp_min2, WIPk + it->second);
                if (can_pair) vp_min2 = std::min(vp_min2, static_cast<energy_t>((k - i) * cp_penalty) + it->second);
            }
            if ((spark.VP_(i_mod, j) < INFover2)) {
                if ((spark.VP_(i_mod, j) < vp_min1 || spark.VP_(i_mod, j) < vp_min2)) {
                    register_candidate(spark.CLVP_, i, j, spark.VP_(i_mod, j));
                    inc_source_ref_count(spark.taVP_, i, j);
                }
            }

            // Things that needed to happen later like W's wmb
            w_wmb = tree.weakly_closed(i, j) ? spark.WMB_[j] + PS_penalty : INF;
            wm_wmb = tree.weakly_closed(i, j) ? spark.WMB_[j] + PSM_penalty + b_penalty : INF;
            w = std::min({w_v, w_split, w_wmb});
            wm = std::min({wm_v, wm_split, wm_wmb});

            // Some case of WI candidate splits will show INf as a vkj value as we are unable to do dangle versions yet but there are dangle
            // candidates
            if (w_v < w_split || wm_v < wm_split || wi_v < wi_split || wip_v < wip_split || paired) {
                // cand_pos_t k_mod = k%(MAXLOOP+1);
                // Encode the dangles into the energies
                energy_t w_enc = (w_v << 2) | d;
                energy_t wm_enc = (wm_v << 2) | d;
                register_candidate(spark.CL_, i, j, spark.V_(i_mod, j), wm_enc, w_enc);
                // always keep arrows starting from candidates
                inc_source_ref_count(spark.ta_, i, j);
            }
            if ((spark.WMB_[j] < INFover2) && (w_wmb < w_split || wm_wmb < wm_split || wi_wmb < wi_split || wip_wmb < wip_split)) {

                register_candidate(spark.CLWMB_, i, j, spark.WMB_[j]);
            }

            spark.W_[j] = w;
            spark.WM_[j] = wm;
            spark.WM2_[j] = std::min(wm2_split, spark.WMB_[j] + PSM_penalty + b_penalty);

        } // end loop j
        rotate_arrays(spark);

        // Clean up trace arrows in i+MAXLOOP+1
        if (garbage_collect && i + MAXLOOP + 1 <= n) {
            gc_row(spark.ta_, i + MAXLOOP + 1);
            gc_row(spark.taVP_, i + MAXLOOP + 1);
        }
        // Reallocate candidate lists in i
        for (auto &x : spark.CL_) {
            if (x.capacity() > 1.5 * x.size()) {
                cand_list_td1 vec(x.size());
                copy(x.begin(), x.end(), vec.begin());
                vec.swap(x);
            }
        }
        for (auto &x : spark.CLVP_) {
            if (x.capacity() > 1.5 * x.size()) {
                cand_list_t vec(x.size());
                copy(x.begin(), x.end(), vec.begin());
                vec.swap(x);
            }
        }
        for (auto &x : spark.CLWMB_) {
            if (x.capacity() > 1.5 * x.size()) {
                cand_list_t vec(x.size());
                copy(x.begin(), x.end(), vec.begin());
                vec.swap(x);
            }
        }
        for (auto &x : spark.CLBE_) {
            if (x.capacity() > 1.5 * x.size()) {
                cand_list_t vec(x.size());
                copy(x.begin(), x.end(), vec.begin());
                vec.swap(x);
            }
        }

        compactify(spark.ta_);
        compactify(spark.taVP_);
    }
    return spark.W_[n];
}

/**
 * @brief Sums the number of Candidates at each index over all indices
 *
 * @param CL_ Candidate list
 * @return total number of candidates
 */
cand_pos_t num_of_candidates(const std::vector<cand_list_td1> &CL_) {
    cand_pos_t c = 0;
    for (const cand_list_td1 &x : CL_) {
        c += x.size();
    }
    return c;
}
cand_pos_t num_of_candidates(const std::vector<cand_list_t> &CL_) {
    cand_pos_t c = 0;
    for (const cand_list_t &x : CL_) {
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
cand_pos_t capacity_of_candidates(const std::vector<cand_list_td1> &CL_) {
    cand_pos_t c = 0;
    for (const cand_list_td1 &x : CL_) {
        c += x.capacity();
    }
    return c;
}
cand_pos_t capacity_of_candidates(const std::vector<cand_list_t> &CL_) {
    cand_pos_t c = 0;
    for (const cand_list_t &x : CL_) {
        c += x.capacity();
    }
    return c;
}

void seqtoRNA(std::string &sequence) {
    for (char &c : sequence) {
        if (c == 'T') c = 'U';
    }
}

void validate_structure(std::string &seq, std::string &structure) {
    cand_pos_t n = structure.length();
    std::vector<cand_pos_t> pairs;
    for (cand_pos_t j = 0; j < n; ++j) {
        if (structure[j] == '(') pairs.push_back(j);
        if (structure[j] == ')') {
            if (pairs.empty()) {
                std::cout << "Incorrect input: More left parentheses than right" << std::endl;
                exit(0);
            } else {
                cand_pos_t i = pairs.back();
                pairs.pop_back();
                if (seq[i] == 'A' && seq[j] == 'U') {
                } else if (seq[i] == 'C' && seq[j] == 'G') {
                } else if ((seq[i] == 'G' && seq[j] == 'C') || (seq[i] == 'G' && seq[j] == 'U')) {
                } else if ((seq[i] == 'U' && seq[j] == 'G') || (seq[i] == 'U' && seq[j] == 'A')) {
                } else {
                    std::cout << "Incorrect input: " << seq[i] << " does not pair with " << seq[j] << std::endl;
                    exit(0);
                }
            }
        }
    }
    if (!pairs.empty()) {
        std::cout << "Incorrect input: More left parentheses than right" << std::endl;
        exit(0);
    }
}

bool exists(const std::string path) {
    struct stat buffer;
    return (stat(path.c_str(), &buffer) == 0);
}

void get_input(std::string file, std::string &sequence, std::string &structure) {
    if (!exists(file)) {
        std::cout << "Input file does not exist" << std::endl;
        exit(EXIT_FAILURE);
    }
    std::ifstream in(file.c_str());
    std::string str;
    cand_pos_t i = 0;
    while (getline(in, str)) {
        if (str[0] == '>') continue;
        if (i == 0) sequence = str;
        if (i == 1) structure = str;
        ++i;
    }
    in.close();
}

/**
 * @brief Simple driver for @see Spark.
 *
 * Reads sequence from command line or stdin and calls folding and
 * trace-back methods of Spark.
 */
int main(int argc, char **argv) {

    args_info args_info;

    // get options (call gengetopt command line parser)
    if (cmdline_parser(argc, argv, &args_info) != 0) {
        exit(1);
    }

    std::string seq;
    if (args_info.inputs_num > 0) {
        seq = args_info.inputs[0];
    } else {
        if (!args_info.input_file_given) std::getline(std::cin, seq);
    }

    std::string restricted = "";
    if (args_info.input_structure_given) restricted = input_structure;

    std::string fileI;
    args_info.input_file_given ? fileI = input_file : fileI = "";

    if (fileI != "") {

        if (exists(fileI)) {
            get_input(fileI, seq, restricted);
        }
        if (seq == "") {
            std::cout << "sequence is missing from file" << std::endl;
        }
    }
    cand_pos_t n = seq.length();
    std::transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
    if (!args_info.noConv_given) seqtoRNA(seq);
    if (restricted == "") restricted = std::string(n, '.');

    if (restricted.length() != (cand_pos_tu)n) {
        std::cout << "input sequence and structure are not the same size" << std::endl;
        std::cout << seq << std::endl;
        std::cout << restricted << std::endl;
        exit(0);
    }

    std::string file = args_info.paramFile_given ? parameter_file : "params/rna_DirksPierce09.par";
    if (exists(file)) {
        vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
    } else if (seq.find('T') != std::string::npos) {
        vrna_params_load_DNA_Mathews2004();
    }
    bool verbose;
    verbose = args_info.verbose_given;

    bool mark_candidates;
    mark_candidates = args_info.mark_candidates_given;

    noGU = args_info.noGU_given;
    seqtoRNA(seq);
    validate_structure(seq, restricted);

    sparse_tree tree(restricted, n);

    Spark spark(seq, !args_info.noGC_given, restricted);

    if (args_info.dangles_given) spark.params_->model_details.dangles = dangle_model;
    pseudoknot = !args_info.pseudoknot_free_given;
    pk_only = args_info.pk_only_given;

    cmdline_parser_free(&args_info);

    int count = 0;
    for (cand_pos_t i = 1; i <= n; ++i) {
        if (tree.tree[i].pair > i || (tree.tree[i].pair < i && tree.tree[i].pair > 0)) count = 4;
        if (tree.tree[i].pair < 0 && count > 0) {
            spark.WI_Bbp[i] = (5 - count) * PUP_penalty;
            count--;
        }
    }

    energy_t mfe = fold(spark, tree, spark.n_, spark.garbage_collect_);
    std::string structure = trace_back(spark, tree, mark_candidates);

    std::ostringstream smfe;
    smfe << std::setiosflags(std::ios::fixed) << std::setprecision(2) << mfe / 100.0;
    std::cout << seq << std::endl;
    std::cout << structure << " (" << smfe.str() << ")" << std::endl;

    if (verbose) {

        std::cout << std::endl;

        std::cout << "TA cnt:\t" << sizeT(spark.ta_) << std::endl;
        std::cout << "TA max:\t" << maxT(spark.ta_) << std::endl;
        std::cout << "TA av:\t" << avoidedT(spark.ta_) << std::endl;
        std::cout << "TA rm:\t" << erasedT(spark.ta_) << std::endl;

        std::cout << std::endl;
        std::cout << "Can num:\t" << num_of_candidates(spark.CL_) << std::endl;
        std::cout << "Can cap:\t" << capacity_of_candidates(spark.CL_) << std::endl;
        std::cout << "TAs num:\t" << sizeT(spark.ta_) << std::endl;
        std::cout << "TAs cap:\t" << capacityT(spark.ta_) << std::endl;

        std::cout << "\nPsuedoknotted\n" << std::endl;
        std::cout << "TAs num:\t" << sizeT(spark.taVP_) << std::endl;
        std::cout << "TAs cap:\t" << capacityT(spark.taVP_) << std::endl;
        std::cout << "TA av:\t" << avoidedT(spark.taVP_) << std::endl;
        std::cout << "TA rm:\t" << erasedT(spark.taVP_) << std::endl;
        std::cout << std::endl;
        std::cout << "WMB Can num:\t" << num_of_candidates(spark.CLWMB_) << std::endl;
        std::cout << "WMB Can cap:\t" << capacity_of_candidates(spark.CLWMB_) << std::endl;
        std::cout << "VP Can num:\t" << num_of_candidates(spark.CLVP_) << std::endl;
        std::cout << "VP Can cap:\t" << capacity_of_candidates(spark.CLVP_) << std::endl;
        std::cout << "BE Can num:\t" << num_of_candidates(spark.CLBE_) << std::endl;
        std::cout << "BE Can cap:\t" << capacity_of_candidates(spark.CLBE_) << std::endl;
    }

    return 0;
}
