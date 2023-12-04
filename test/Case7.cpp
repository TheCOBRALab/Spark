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

bool canMatch(char x, char y){

if((x == 'A' && y == 'T') || (x == 'T' && y == 'A')) {return true;}
else if((x == 'C' && y == 'G') || (x == 'G' && y == 'C')) {return true;}
else if((x == 'A' && y == 'U') || (x == 'G' && y == 'U') || (x == 'U' && y == 'G') || (x == 'U' && y == 'A')) {return true;}
else{return false;}

}

using namespace std;
// Used to separate a reference structure into it's constituent pieces, so G_big and G_small
int main(int argc,char **argv) {

    string file = "/home/mgray7/Spark/RNAstrand.txt";

    string file1 = "/home/mgray7/Spark/RNAstrandstructures.txt";


    vector<string> names;
    vector<string> seqs;

    vector<string> strucs;

    int k = 0;
    ifstream in(file);
    string str;
    while(getline(in,str)){
        if(k%2 == 1){
            seqs.push_back(str);
        }
        else{
            names.push_back(str);
        }
        ++k;
    }

    in.close();

    k=0;
    ifstream in1(file1);
    while(getline(in1,str)){
        if(k%2 == 1){
            strucs.push_back(str);
        }
         ++k;
    }
    in1.close();

    


    /*
    * Remove seqs that have iupac characters
    */
    vector<string> namesAfter;
    vector<string> seqsAfter;
    vector<string> strucsafter;
    std::cout << names.size() << " " << seqs.size() << " " << strucs.size() << endl;
    for(int i =0; i<seqs.size();++i){
        string seq = seqs[i];
        transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
        if(seq.length()>1000) continue;
        if(seq.find('X') != std::string::npos || seq.find('P') != std::string::npos || seq.find('N') != std::string::npos || seq.find('Y') != std::string::npos || seq.find('K') != std::string::npos || seq.find('R') != std::string::npos || seq.find('S') != std::string::npos || seq.find('W') != std::string::npos || seq.find('~') != std::string::npos || seq.find('V') != std::string::npos || seq.find('H') != std::string::npos || seq.find('D') != std::string::npos || seq.find('B') != std::string::npos){
            continue;
        }
        else{
            namesAfter.push_back(names[i]);
            seqsAfter.push_back(seq);
            strucsafter.push_back(strucs[i]);
            if(seqsAfter.size() == 69) std::cout << seqs[i] << std::endl;
        }
    }
    std::cout << namesAfter.size() << endl;
    /*
    * Remove structures that have crossing pseudoknots
    */
    vector<string> seqsx;
    vector<string> strucsx;
    vector<string> namesx;
    for(int i = 0; i < strucsafter.size();++i){
        int length = strucsafter[i].length();
        string structure = strucsafter[i];
        bool x = true;
        for(int j = 0;j<length;++j){
            if(isupper(structure[j])){
                x = false;
                break;
            }
        }
        if(x){
            strucsx.push_back(structure);
            seqsx.push_back(seqsAfter[i]);
            namesx.push_back(namesAfter[i]);
        } 
    }
    std::cout << strucsx.size() << endl;


    /*
    * Split structures into big and small
    */
    vector<string> strucs1;
    vector<string> strucs2;
    for(int i = 0; i < strucsx.size();++i){
        string structure1 = strucsx[i];
        string structure2 = strucsx[i];
        int length = strucsx[i].length();
        for(int j = 0;j<length;++j){
            if(structure1[j] != '.' && structure1[j] != '(' && structure1[j] != ')'){
                structure1[j] = '.';
            }
            if(structure2[j] == '(' || structure2[j] == ')'){
                structure2[j] = '.';
            }
            if(structure2[j] != '.' && structure2[j] != '(' && structure2[j] != ')'){
                if(structure2[j] == '{' || structure2[j] == '<' || structure2[j] == '['){
                    structure2[j] = '(';
                }
                if(structure2[j] == '}' || structure2[j] == '>' || structure2[j] == ']'){
                    structure2[j] = ')';
                }
                
            }
        }
        strucs1.push_back(structure1);
        strucs2.push_back(structure2);
    }
    std::cout << strucs1.size() << "  " << strucs2.size() << std::endl;
    //  exit(0);
    /*
    * If non-canonical, change to .
    */
    vector<string> gbig;
    vector<string> gsmall;
   for(int i = 0; i<strucs1.size();++i){
        vector<int> paren;
        string sequence = seqsx[i];
        string structure = strucs1[i];
        int length = structure.length();

        for(int j=0;j<length;j++){

            if(structure[j] == '(') {
                paren.push_back(j);
                continue;
            }

            int x;
            bool close = false;
            if (structure[j] == ')' && !paren.empty()){
                x = paren[paren.size()-1];
                paren.erase(paren.begin()+(paren.size()-1));
                close = true;
            }

            if(close){
                if(canMatch(sequence[j],sequence[x]) && j-x > 3){
                    structure[j] = ')';
                    structure[x] = '(';
                }
                else{
                    structure[j] = '.';
                    structure[x] = '.';
                }
            }
            
        }
        gbig.push_back(structure);

    }

     for(int i = 0; i<strucs2.size();++i){
        vector<int> paren;
        string sequence = seqsx[i];
        string structure = strucs2[i];
        int length = structure.length();

        for(int j=0;j<length;j++){

            if(structure[j] == '(') {
                paren.push_back(j);
                continue;
            }

            int x;
            bool close = false;
            if (structure[j] == ')' && !paren.empty()){
                x = paren[paren.size()-1];
                paren.erase(paren.begin()+(paren.size()-1));
                close = true;
            }

            if(close){
                if(canMatch(sequence[j],sequence[x])&& j-x > 3){
                    structure[j] = ')';
                    structure[x] = '(';
                }
                else{
                    structure[j] = '.';
                    structure[x] = '.';
                }
            }
            
        }
        gsmall.push_back(structure);

    }
    std::cout << seqsx.size() << " " << gbig.size() << std::endl;

   /*
   * Output
   */
    string fileO = "RNAseq.txt";
    string fileO1 = "Gbig.txt";
    string fileO2 = "Gsmall.txt";

    ofstream out(fileO);
    for(int i = 0; i< namesx.size();++i){
        out << namesx[i] << endl;
        out << seqsx[i] << endl;
    }

    ofstream out1(fileO1);
    for(int i = 0; i< namesx.size();++i){
        out1 << namesx[i] << endl;
        out1 << gbig[i] << endl;
    }

    ofstream out2(fileO2);
    for(int i = 0; i< namesx.size();++i){
        out2 << namesx[i] << endl;
        out2 << gsmall[i] << endl;
    }

    out.close();
    out1.close();
    out2.close();



    return 0;
}