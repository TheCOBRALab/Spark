#include "PK_globals.hh"
#include "base_types.hh"
#include "sparse_tree.cc"



std::string remove_structure_intersection(std::string restricted, std::string structure){
	cand_pos_t length = structure.length();
	for(cand_pos_t i=0; i< length; ++i){
		if(restricted[i] == '(' || restricted[i] == ')') structure[i] = '.';
		
		if (output[i] == '['){
			output[i] = '(';
		}
		if (output[i] == ']'){
			output[i] = ')';
		}
	}
	return output;
}

std::string find_disjoint_substructure(std::string structure){
	cand_pos_t length = structure.length();
	string restricted= std::string (n,'.')
	for(cand_pos_t i = 0; i< length;++i){
		// if()
	}
}

std::string obtainRelaxedStems(std::string restricted, std::string pkfree_structure, sparse_tree tree){
	cand_pos_t length = restricted.length();

	//Gresult <- G1
	std::string relaxed = pkfree_structure;

	cand_pos_t i = 0;
	cand_pos_t j = 0;

	for(cand_pos_t k=0;k<length;k++){
		if(tree.tree[k] > -1){
			i = k;
			j = G2_pair[k];
			if(i < j){ //for each ij in G2
				if( (G1[i] != G2[i]) && (G1[j] != G2[j]) ){//if ij not in G1
					//include bulges of size 1
					if(paired_structure(i-1,j+1,G1_pair,length) || paired_structure(i+1,j-1,G1_pair,length) ){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
					//include loops of size 1x1
					}else if( paired_structure(i-2,j+1,G1_pair,length) || paired_structure(i-1,j+2,G1_pair,length) || \
							paired_structure(i+1,j-2,G1_pair,length) || paired_structure(i+2,j-1,G1_pair,length) ){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
					//include loops of size 1x2 or 2x1
					}else if( paired_structure(i-2,j+2,G1_pair,length) || paired_structure(i+2,j-2,G1_pair,length) ){
						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
					}else if( paired_structure(i-3,j+2,G1_pair,length) || paired_structure(i-2,j+3,G1_pair,length) || \
							paired_structure(i+2,j-3,G1_pair,length) || paired_structure(i+3,j-2,G1_pair,length) ){

						Gresult[i] = G2[i];
						Gresult[j] = G2[j];
					}
				}
			}
		}
	}
}