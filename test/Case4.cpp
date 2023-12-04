#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <tuple>

#include <iomanip>

#include <sstream>

#include <limits>
#include <iterator>
#include <cstring>
#include <cassert>
#include <numeric>
#include <map>
//#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace std;
// Used to change files from rnastrand into standard fasta format
// changes sequence and structure into a single line
// Changes extra pseudoknot characters <>,{},Aa into just []
// Outputs into two files for seqs and structures 
int main(int argc,char **argv) {
    
    string file = "/home/mgray7/Spark/rnastrand.dp";

    // string fileO = "/home/mgray7/output2/fasta/rnasoftO.txt";

    vector<string> names;
    vector<string> seqs;
    vector<string> structures;
    

    ifstream in(file);
    // ofstream out(fileO);
    string str;
    string seq ="";
    string structure = "";
    string name = "";
    bool struc = false;
    while(getline(in,str)){
        if(str[0] == '#'){
            name = str.substr(7,str.length()-7);
            names.push_back(name);
            getline(in,str);
            getline(in,str);
            getline(in,str);
            getline(in,str);
        }
        
        if(seq == "" && !struc){
            // if(str[0] == 'A' || str[0] == 'C' || str[0] == 'G' || str[0] == 'U' || str[0] == 'N'){
            seq = str;
            // }
               
        }
        else if (seq != "" && !struc){
            seq = seq + str;
        }

        if(structure == "" && struc){
                // if(str[0] == '(' || str[0] == ')' || str[0] == '[' || str[0] == ']' || str[0] == '{' || str[0] == '}' || str[0] == '{'){
            structure = str;
            // }
        }
        else if (structure != "" && struc){
            structure = structure + str;
        }
        if(str == "" || str.length() < 50){
            if(seq != ""){
                struc = true;
                seqs.push_back(seq);
                seq = "";
            }
            if(structure != ""){
                struc = false;
                structures.push_back(structure);
                structure = "";
            }
        }
       
    }
    // std::cout << names.size() << " " << seqs.size() << " " << structures.size() << std::endl;
    // std::cout << names[0] << "\n" << seqs[0]  << std::endl;
    ofstream out("/home/mgray7/Spark/RNAstrand.txt");
    for(int i = 0;i<names.size();++i){
        out << ">" << names[i] << endl;
        out << seqs[i] << endl;
    }
    out.close();

    for(int i = 0;i<structures.size();++i){
        for(int j = 0; j<structures[i].length();++j){
            if(structures[i][j] == '<' || structures[i][j] == '{' || structures[i][j] == 'A') structures[i][j] = '[';
            if(structures[i][j] == '>' || structures[i][j] == '}' || structures[i][j] == 'a') structures[i][j] = ']';
        }
    
    }

    ofstream out1("/home/mgray7/Spark/RNAstrandstructures.txt");
    for(int i = 0;i<names.size();++i){
        out1 << ">" << names[i] << endl;
        out1 << structures[i] << endl;
    }
    out1.close();


    return 0;
}