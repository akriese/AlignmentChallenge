#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <math.h>
#include <omp.h>
#include "AlignmentChallenge.h"

int main(int argc, char* argv[]) {
  if (argc != 6 && argc != 8) {
    std::cerr << "Wrong number of inputs!" << std::endl;
    return 1;
  }

  const std::string source = argv[1];
  const int threads = std::stoi(argv[2]);
  const int match = std::stoi(argv[3]);
  const int mismatch = std::stoi(argv[4]);
  const int gap = std::stoi(argv[5]);
  

  std::vector<std::string> db{};
  std::ifstream ifs{source};
  std::string line{};

  if (ifs.good()) {
    while (getline(ifs, line)) {
      if (line[0] != '>') db.push_back(line); 
    }
  } else {
    std::cerr << "Could not read file!" << std::endl;
    return 1;    
  }
  
  std::cout << db.size() << " Sequenzen" << std::endl;  
  
  Alignment al{match, mismatch, gap};
  
  for (size_t i = 0; i < 5; ++i) {    
    std::cout << "Computing summed score with SSE method..." << std::endl;
    auto start = omp_get_wtime();
    int score = al.compute(db, threads);
    auto end = omp_get_wtime();
    std::cout << "Score: " << score << std::endl;
    std::cout << end - start << " secs needed! With " << threads << " thread(s)!" << "\n\n";
  }
  
  
/**/  
/**/
  return 0;
}
