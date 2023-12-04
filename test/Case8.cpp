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
// Used to separate a reference structure into it's constituent pieces, so G_big and G_small
int main(int argc,char **argv) {

    string file1 = "/home/mgray7/output3/proof/Spark/rnastrandBig.txt";
    string file2 = "/home/mgray7/output3/proof/Spark/rnastrandSmall.txt";

    string fileH1 = "/home/mgray7/output3/proof/HFold/rnastrandBig.txt";
    string fileH2 = "/home/mgray7/output3/proof/HFold/rnastrandSmall.txt";


    vector<string> structures1;
    vector<double> energies1;

    vector<string> structures2;
    vector<double> energies2;

    vector<string> structuresH1;
    vector<double> energiesH1;

    vector<string> structuresH2;
    vector<double> energiesH2;


    // Get sequences, structures and energies
    ifstream in(file1);
    string str;
    string seq;
    int k = 0;
    vector<string> seqs;
    while(getline(in,str)){
        if(k%2 == 0){
            seqs.push_back(str);
        }
        else{
            string structure = str.substr(0,seq.length());
            double energy = stod(str.substr(seq.length()+2,str.length()-1));
            structures1.push_back(structure);
            energies1.push_back(energy);
        }
        seq = str;
        ++k;
    }
    k=0;

    ifstream in1(fileH1);
    while(getline(in1,str)){
        if(k%2 == 0){
        }
        else{
            string structure = str.substr(5,seq.length()-5);
            double energy = stod(str.substr(seq.length()+2,str.length()));
            structuresH1.push_back(structure);
            energiesH1.push_back(energy);
        }
        seq = str;
        ++k;
    }

    k=0;
    ifstream in2(file2);
    while(getline(in2,str)){
        if(k%2 == 0){
        }
        else{
            string structure = str.substr(0,seq.length());
            double energy = stod(str.substr(seq.length()+2,str.length()-1));
            structures2.push_back(structure);
            energies2.push_back(energy);
        }
        seq = str;
        ++k;
    }

    k=0;
    ifstream in3(fileH2);
    while(getline(in3,str)){
        if(k%2 == 0){
        }
        else{
            string structure = str.substr(5,seq.length()-5);
            double energy = stod(str.substr(seq.length()+2,str.length()));

            structuresH2.push_back(structure);
            energiesH2.push_back(energy);
        }
        seq = str;
        ++k;
    }



    // check energies and structures
    int n = seqs.size();
    int count = 0;
    int countS = 0;
    vector<int> needsToBeChecked;
    for(int i = 0; i < n;++i){
        if(energies1[0] == energiesH1[0]) count++;
        // if(structures1[0] == structuresH1[0]) countS++;
    }
    for(int i = 0; i < n;++i){
        int length = seqs[i].length();
        bool check = true;
        string structure1 = structures1[i];;
        string structure2 = structuresH1[i];
        for(int j = 0; j < length; ++j){
            if(structure1[j] != structure2[j]){
                if((structure1[j] == '(' && structure2[j] == '[') || (structure1[j] == ')' && structure2[j] == ']')){

                }
                else{
                    check = false;
                    
                }
            }
        }
        if(check) countS++;
        else{
            needsToBeChecked.push_back(i);
        }
    }


    int count2 = 0;
    int countS2 = 0;
    vector<int> needsToBeChecked2;
    for(int i = 0; i < n;++i){
        if(energies2[0] == energiesH2[0]) count2++;
        // if(structures1[0] == structuresH1[0]) countS++;
    }
    for(int i = 0; i < n;++i){
        int length = seqs[i].length();
        bool check = true;
        string structure1 = structures2[i];;
        string structure2 = structuresH2[i];
        for(int j = 0; j < length; ++j){
            if(structure1[j] != structure2[j]){
                if((structure1[j] == '(' && structure2[j] == '[') || (structure1[j] == ')' && structure2[j] == ']')){

                }
                else{
                    check = false;
                    
                }
            }
        }
        if(check) countS2++;
        else{
            needsToBeChecked2.push_back(i);
        }
    }
    
    std::cout << count << "/" << n << "  " << countS << "/" << n << endl;
    std::cout << count2 << "/" << n << "  " << countS2 << "/" << n << endl;




    std::cout << needsToBeChecked2[0] << endl;



    // for(int i = 0; i < needsToBeChecked.size();++i){
    //     string command = "/home/mgray7/hotknots/bin/computeEnergy -m \"DP\" -p params/parameters_DP09.txt -s " + seqs[needsToBeChecked[i]] + " \"";

    //     string commandS = command + structures1[needsToBeChecked[i]] + "\" > out.txt";
    //     string commandH = command + structuresH1[needsToBeChecked[i]] + "\" > out.txt";

    //     std::cout << commandS << std::endl;
    //     system(commandS.c_str());
    //     string strr;
    //     ifstream in("out.txt");

    //     while(getline(in,strr)){

    //     }
    //     // std::cout << strr << std::endl;


    // }




}