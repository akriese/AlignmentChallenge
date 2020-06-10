#include <vector>
#include <string>
#include <iostream>
#include <omp.h>
#include "AlignmentChallenge.h"
#include <smmintrin.h>
#include <tmmintrin.h>
#include <algorithm>
#include <iterator>

Alignment::Alignment(const int match, const int mismatch, const int gap) :
                    _match(match),
                    _mismatch(mismatch),
                    _gap(gap),
                    _oben_mask(_mm_set1_epi16(150*gap)),
                    _ms(_mm_set1_epi16(match)),
                    _mms(_mm_set1_epi16(mismatch)),
                    _gaps(_mm_set1_epi16(gap)),
                    _shuff_mask(_mm_set_epi8(13,12,11,10,9,8,7,6,5,4,3,2,1,0,15,14)) {
  for (size_t j = 0; j < 151; ++j)
    maskSSE[j] = _mm_set1_epi16(j*gap);
  for (size_t j = 0; j < 1208; ++j)
    maskSSE64[j] = _mm_set1_epi16((j/8)*gap);
  for (size_t j = 0; j < 151; ++j)
    maskNormal[j] = j*gap;
}

void merge8Seqs(int16_t(& v)[1200], const std::string& s1, const std::string& s2,
              const std::string& s3, const std::string& s4, const std::string& s5,
              const std::string& s6, const std::string& s7, const std::string& s8) {
  for (size_t i = 0; i < 150; ++i) {
    v[(i*8)] = s1[i];
    v[(i*8)+1] = s2[i];
    v[(i*8)+2] = s3[i];
    v[(i*8)+3] = s4[i];
    v[(i*8)+4] = s5[i];
    v[(i*8)+5] = s6[i];
    v[(i*8)+6] = s7[i];
    v[(i*8)+7] = s8[i];
  }
}

void merge8SeqsBlocks(__m128i(& v)[150], const std::string& s1, const std::string& s2,
              const std::string& s3, const std::string& s4, const std::string& s5,
              const std::string& s6, const std::string& s7, const std::string& s8) {
  for (size_t i = 0; i < 150; ++i) {
    v[i] = _mm_set_epi16(s1[i],s2[i],s3[i],s4[i],s5[i],s6[i],s7[i],s8[i]);
  }
}

// eine Sequenz wird acht Mal verglichen
int Alignment::computeSSE(const std::string& s0, const int16_t(& v)[1200], const uint32_t number_cmp) {
  __m128i linear[151];
  std::copy(std::begin(maskSSE), std::end(maskSSE), std::begin(linear));

  const __m128i ms = _ms;                                                       // xmm mit 8x match score gefüllt
  const __m128i mms = _mms;                                                     // xmm mit 8x mismatch score gefüllt
  const __m128i gaps = _gaps;                                                   // xmm mit 8x gap score gefüllt

  __m128i oben {};
  __m128i diag {};

  for (size_t i = 1; i < 151; ++i)
  {
    oben = _mm_set1_epi16(i * _gap);                                            // temp mit gap gefüllt
    const __m128i kmask = _mm_set1_epi16(s0[i-1]);                              // buchstaben aus s0 in register gefüllt


    for (__m128i *ptrl = (linear+1), *ptrd = linear, *ptrv = (__m128i*) v,
        *end = (__m128i*) (v+1200); ptrv < end; ++ptrd, ++ptrl, ++ptrv) {
      diag = *ptrd;
      *ptrd = oben;

      __m128i *ptr_oben = &oben, *ptr_diag = &diag;                             // pointer zu oben/diag

      __m128i score_left = _mm_add_epi16(*ptrl, gaps);                          // links+gap
      __m128i score_upper = _mm_add_epi16(*ptr_oben, gaps);                     // oben+gaps
      __m128i comparison = _mm_cmpeq_epi16(kmask, *ptrv);                       // cmp(s0,seqs)
      __m128i scores = _mm_max_epi16(score_left, score_upper);                  //max(oben, unten)
      __m128i diag_p_match = comparison & ms;                                   // ms an matchenden-Stellen
      __m128i diag_p_mismatch = ~comparison & mms;                              // mms an mismatch-Stellen
      __m128i score_diag = _mm_add_epi16(*ptr_diag,                             // diag. score mit matchscore add.
          diag_p_match | diag_p_mismatch);                                      // match und mismatch verodert
      *ptr_oben = _mm_max_epi16(score_diag, scores);                            // max(diag, seiten) in temp geschoben
    }

    *(linear+150) = oben;                                                       // letzte 8pack gefüllt
  }

  int16_t score = 0;
  int16_t erg[8];
  *(__m128i*) erg = oben;                                                       //ganzes register wird kopiert
  for (size_t i = 0; i < number_cmp; ++i)                                       //aber nur #comp werden addiert
    score += erg[i];

  return score;
}

//acht Seqeunzen werden mit acht Sequenzen verglichen
int Alignment::computeSSE64Blocks(const std::string& s0, __m128i (& v)[150], const __m128i (& w)[150]) {
  __m128i linear[1208];
  std::copy(std::begin(maskSSE64), std::end(maskSSE64), std::begin(linear));

  const __m128i ms = _ms;                                                       // xmm mit 8x match score gefüllt
  const __m128i mms = _mms;                                                     // xmm mit 8x mismatch score gefüllt
  const __m128i gaps = _gaps;                                                   // xmm mit 8x gap score gefüllt
  const __m128i shuff_mask = _shuff_mask;

  __m128i oben[8] {};
  __m128i diag[8] {};

  for (size_t i = 1; i < 151; ++i)
  {
    for (__m128i *ptr = oben; ptr != (__m128i*) (oben+8); ++ptr) {
      *ptr = _mm_set1_epi16(i * _gap);                                          // temp mit gap gefüllt
    }

    __m128i kmask = *(w+(i-1));                                                 // buchstaben aus s0 in register gefüllt

    // iteratoren (pointer) werden initialisiert
    // *ptrd wird auch nach der folgenden Schleife noch benötigt
    __m128i *ptrl = (linear+8), *ptrd = linear;

    for (__m128i *ptrv = v, *end = (v+150); ptrv < end; ++ptrv)
    {
      for (size_t m = 0; m<8; ++m, ++ptrd, ++ptrl)
      {
        diag[m] = *(ptrd);                                                      //temps werden benutzt
        *(ptrd) = oben[m];

        __m128i *ptr_oben = &(oben[m]), *ptr_diag = &(diag[m]);                 // pointer zu oben/diag

        __m128i scores = _mm_max_epi16(*ptrl, *ptr_oben);
        scores = _mm_add_epi16(scores, gaps); 
        __m128i comparison = _mm_cmpeq_epi16(kmask, *ptrv);                     // cmp(s0,seqs)
        __m128i diag_p_match = comparison & ms;                                 // ms an matchenden-Stellen
        __m128i diag_p_mismatch = ~comparison & mms;                            // mms an mismatch-Stellen
        __m128i score_diag = _mm_add_epi16(*ptr_diag,                           // diag. score mit matchscore add.
            diag_p_match | diag_p_mismatch);                                    // match und mismatch verodert
        *ptr_oben = _mm_max_epi16(score_diag, scores);                          // max(diag, seiten) in temp geschoben

        //shuffle
        // erste zwei 8bit-ints werden hinten dran gehängt
        // und der rest um 16bits geshiftet
        kmask = _mm_shuffle_epi8(kmask, shuff_mask);
      }
    }
                                                                                //letzter Block jeder berechnung wird gefüllt
    for (size_t m = 0; m<8; ++m, ++ptrd)
      *ptrd = oben[m];
  }

  int score = 0;
  int16_t erg[64];
  for (size_t m = 0; m<8; ++m)
    *(__m128i*) (erg+m*8) = oben[m];

  for (size_t i = 0; i < 64; ++i)
    score += erg[i];

  return score;
}

int Alignment::compute(const std::vector<std::string>& seqs, const int threads) {
  int score = 0;

  omp_set_num_threads(threads);
  //folgendes geht nur, weil 1000 durch 8 teilbar ist. ansonsten müssten die Randfälle unten
  //abgearbeitet werden. Diese habe ich zur Übersichtlichkeit rausgenommen. Randfälle würden bspws
  //mit NWTemp (oben auskommentiert) computed werden
  size_t aufachter = (seqs.size() >> 3) << 3;                                   //(size/8)*8
  #pragma omp parallel for schedule(dynamic,3) reduction(+:score)
  for (size_t i = 0; i < aufachter; i+=8) {
    __m128i v[150]{};
    merge8SeqsBlocks(v, seqs[i], seqs[i+1], seqs[i+2], seqs[i+3],seqs[i+4],seqs[i+5],seqs[i+6],seqs[i+7]);
    for (size_t j = i+8; j < seqs.size(); j+=8) {                               //i+8 ist letzter Fall, wo alle acht seq benutzt werden
      __m128i w[150]{};
      merge8SeqsBlocks(w, seqs[j], seqs[j+1], seqs[j+2], seqs[j+3],seqs[j+4],seqs[j+5],seqs[j+6],seqs[j+7]);
      score += computeSSE64Blocks(seqs[j], v, w);                               //ausführung der sse-fkt mit 8x8 sequenzen (shufflen)
    }

    int16_t v2[1200]{};
    merge8Seqs(v2, seqs[i], seqs[i+1], seqs[i+2], seqs[i+3],seqs[i+4],seqs[i+5],seqs[i+6],seqs[i+7]);
    for (size_t j = i+1; j < i+8; ++j) {                                         // von einer bis sieben Sequenzen SSE
      const uint32_t comp = j-i;
      score += computeSSE(seqs[j], v2, comp);
    }
  }

  return score;
}
