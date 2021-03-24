#include <array>
#include <iostream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include "KeepN.hpp"
#include <vector>
#include <sstream>
#include <string>
#include <unordered_map>

struct USRCAT_SL_Descriptor {
     unsigned long linenum;
     float usrcatvals[60];
};

inline float score_query_candidate(USRCAT_SL_Descriptor &query, USRCAT_SL_Descriptor &candidate ){
	  float v,diff = 0;
	  for (unsigned int i = 0; i < 60; i++) {
          v=(query.usrcatvals[i] - candidate.usrcatvals[i]);
		  diff += v*v;
	  }
	  return diff;
  }

int main(int argc, char* argv[])
{
    //Arguments must be 
    //1: Binary file location without last .bin extension, so that .bin and .smi file locations can be derived
    //2: Number of best to keep
    //3-63: USRCAT descriptors of query

    std::string binfilename = std::string(argv[1]) + ".bin";
    std::string smifilename = std::string(argv[1]) + ".smi";

    int num_to_keep=atoi(argv[2]);
    KeepNAscending<unsigned long> keepn(num_to_keep);

    // Set up USRCAT query molecule
    USRCAT_SL_Descriptor query;
    for(int i=0;i<60;++i){
        query.usrcatvals[i]=strtof(argv[i+3], NULL);
    }

    // Find how many mols are in the descriptor file
    uint64_t begin_pos, end_pos;
    std::ifstream binary_stream;
    binary_stream.open(binfilename.c_str());
    binary_stream.seekg(0);
    begin_pos = binary_stream.tellg();
    binary_stream.seekg(0, std::ios::end);
    end_pos = binary_stream.tellg();
    unsigned int numobjects = (end_pos - begin_pos) / sizeof(query);
    binary_stream.seekg(0);
    std::cerr << "Usrcat\tReading " << numobjects << " candidate molecules\n";
    USRCAT_SL_Descriptor candidate;
    float score;
    for (unsigned int i = 0; i < numobjects; i++) {
        binary_stream.read((char*)&candidate, sizeof(candidate));
        score=score_query_candidate(query, candidate);
        keepn.insert(candidate.linenum, score);
        }
    std::vector<unsigned long> best_line_numbers;

    // Insert best line numbers into best_line_numbers vector, and also make a sorted version
    best_line_numbers.resize(num_to_keep);
    int bestcounter=0;
    for (auto &&i : keepn.best) {
        best_line_numbers[bestcounter]=i.first;
        ++bestcounter;
    }
    std::vector<unsigned long> sorted_best_line_numbers = best_line_numbers;
    std::sort(sorted_best_line_numbers.begin(), sorted_best_line_numbers.end());

    std::unordered_map<unsigned long, std::string> linenumber_to_string_map;
    std::string line;
    std::ifstream infile;
    infile.open(smifilename.c_str());
    unsigned long nreads, cur_line=0;
    int sorted_index=0;
    unsigned long i;

    while(sorted_index<sorted_best_line_numbers.size()){
        nreads=sorted_best_line_numbers[sorted_index]-cur_line;
        for(i=0;i<nreads;++i){
            std::getline(infile, line);
        }
        cur_line+=nreads;
        linenumber_to_string_map[sorted_best_line_numbers[sorted_index]]=line;
        ++sorted_index;
    }
    for (auto const& pair: keepn.best){
        std::cout<<linenumber_to_string_map[pair.first]<<","<<pair.second<<"\n";
    }
}

