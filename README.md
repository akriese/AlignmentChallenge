# AlignmentChallenge

This is a project from my third semester in my Bioinformatics bachelor.

## Brief
The challenge is to compute the sum of scores of all pairwise
alignments of 1000 DNA sequences of length 150 each as fast as possible.

The score shall be computed with the [Needleman-Wunsch-Algorithm](https://en.wikipedia.org/wiki/Needleman%E2%80%93Wunsch_algorithm).  
Using multithreading and SSE-Instructions (up to sse4.2) was allowed.  
More information on the challenge can be found in the original task "challenge_task.txt".

## Files
* **AlignmentChallenge.h**: header for the class definition
* **AlignmentChallenge.cpp**: implementation of the Algorithm
* **bench_AC.cpp**: main function, reading of the file and calling of the Algorithm
* **mason_illummina_hg38_chr10_1K_150bp.fasta**: input file with 1K sequences, each 150bp long
* **Makefile**: prepared makefile to build the binary bench_AC.

## Call
    make bench_AC  
    ./bench_AC mason_illumina_hg38_chr10_1K_150bp.fa 16 3 -1 -2  
    # ./bench_AC <input_file> <num_threads> <match_score> <mismatch_score> <gap_score>

## Result
This has been the first time for me using SIMD instructions.
It took a lot of time understanding the idea and even more time to implement and tweak the solution.  
Note that the solution is optimized for this particular problem. Knowing the fixed length of the sequences made things easier.

In the end, my implementation was actually faster than the program from the [SeqAn](https://github.com/seqan/seqan3) library,
that was "the end boss" of the challenge.  

(To not fill up this repo too much, I haven't included the binary of the SeqAn-program.
Feel free to ask me to send it to you, or try to build it from the library on your own.)
