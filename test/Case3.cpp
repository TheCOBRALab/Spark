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

string make_structure(string seq);
bool check_Pseudoknot(vector<tuple<int,int> > used, int i, int j);

static int pairs[8][8] =
  /* _  A  C  G  U  X  K  I */
{ { 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0, 0, 0, 0, 5, 0, 0, 5 },
  { 0, 0, 0, 1, 0, 0, 0, 0 },
  { 0, 0, 2, 0, 3, 0, 0, 0 },
  { 0, 6, 0, 4, 0, 0, 0, 6 },
  { 0, 0, 0, 0, 0, 0, 2, 0 },
  { 0, 0, 0, 0, 0, 1, 0, 0 },
  { 0, 6, 0, 0, 5, 0, 0, 0 } };

int main(int argc,char **argv) {
    vector<string> names;
    vector<string> seqs;
    string family;

    string file = "/Users/mateo2/Documents/Code/output/fasta/";
    cout << "Put Family: ";
    cin >> family;
    file = file + family + ".txt";
    ifstream in(file);

    string str;
    int i = 0;
    while(getline(in,str)){
    // cout << str << endl;
        if(i%2 == 0){
            names.push_back(str);
        }
        else if(i%2==1){
            seqs.push_back(str);
        }
        ++i;
    }  
    in.close();

    string fileO = "/Users/mateo2/Documents/Code/output/structures/" + family + ".txt";
    ofstream out(fileO);
    for(int i= 0;i<seqs.size();++i){
        string seq = seqs[i];
        string structure = make_structure(seq);
        out << names[i] << endl;
        // out << seqs[i] << endl;
        out << structure << endl;

    }

 
    return 0;
}

string make_structure(string seq){
    std::map<char,uint> base;
    base['A']=1;
    base['C']=2;
    base['G']=3;
    base['U']=4;

    int len = seq.length();
    string structure (len,'.');
    int k = 0;
    vector<tuple<int,int> > used;
    for(int i =0;i<len-25;++i){

        for(int j =i+25;j<len;++j){

        
            if(pairs[base[seq[i]]][base[seq[j]]] > 0 && pairs[base[seq[i+1]]][base[seq[j-1]]] > 0 && pairs[base[seq[i+2]]][base[seq[j-2]]] > 0 && pairs[base[seq[i+3]]][base[seq[j-3]]] > 0){
                structure[i] = '(';
                structure[i+1] = '(';
                structure[i+2] = '(';
                structure[i+3] = '(';
                structure[j] = ')';
                structure[j-1] = ')';
                structure[j-2] = ')';
                structure[j-3] = ')';
                ++k;
                i=j+30;
                break;
                

            }
        }
        if(k>1) break;
        
    }

    return structure;

}

bool check_Pseudoknot(vector<tuple<int,int> > used, int i, int j){
  for(int m = 0; m<used.size();++m){
        int k = get<0>(used[m]);
        int l = get<1>(used[m]);
        // cout << i << "\t" << k << "\t" << l << "\t" << j << endl;
        if(i==k || j==l || i==l || j==k) return true;
        if((i < k  && j > k && j < l) || (i < l  && j > l && i > k)) return true;
    }
  return false;
}