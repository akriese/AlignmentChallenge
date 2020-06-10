#include <string>
#include <vector>
#include <smmintrin.h>


class Alignment
{
public:
  // DO NOT CHANGE THE PUBLIC INTERFACE!
  
  /// This makes the default constructor private
  /// i.e. an object of the class can only be created with sequences (see below)
  Alignment() = delete;  // no NOT implement this function. Just leave it as is.
  
  /// Constructor with two sequences
  /// Makes an internal copy of the sequences.
  Alignment(const int match, const int mismatch, const int gap);
  
  /// compute the global aligment of all vs. all sequences and return the sum of scores
  int compute(const std::vector<std::string>& seqs, const int threads);

private:
  const int16_t _match;
  const int16_t _mismatch;
  const int16_t _gap;
  const int arraysz = 1208;
  __m128i maskSSE[151];
  __m128i maskSSE64[1208];
  int16_t maskNormal[151];
  const __m128i _oben_mask;
  const __m128i _ms;
  const __m128i _mms;
  const __m128i _gaps;
  const __m128i _shuff_mask;
  int computeSSE(const std::string& s0, const int16_t(& v)[1200], const uint32_t number_comp);
  int computeSSE64Blocks(const std::string& s0, __m128i(& v)[150], const __m128i(& w)[150]); 
  int needleWunschTemp(const std::string& seq1, const std::string& seq2);

  int16_t getVal(const char a, const char b);
  // add your private functions and member variables here
  ///....
  

};
