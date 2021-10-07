## Summary: 
    Given a query string (length n) and a pattern (length m), output all matches including inexact matches. This extends the original Ukkonen's algorithm (O(min(m, n)*k) complexity) but allows partial matches at the 5' or 3' end and mismatches under given error rate. 

    Compared with https://github.com/rust-bio/rust-bio/blob/master/src/pattern_matching/ukkonen.rs, it uses error rate instead
    of fixed number of errors and also removes duplicate hits by retaining the best one; 
    Compared with https://github.com/marcelm/cutadapt/blob/main/src/cutadapt/_align.pyx, it outputs all hits instead of the 
    5'-most one.

## Usage:
    {} --query text --pattern pattern [--min_overlap 6] [--max_err_rate 0.1] [--debug info] [--help] [--version]

## Output:
    a vector of triplets: start, end, cost (i.e. mismatch/insertion/deletion)
                 0         10        20        30        40        50        60        70
    text         TCTAAGAAGTAGGCGCTCTCCCTCTACGAAGTACTCTAATGAACCCTCTACAGAAGTAGAGTATGTGACCCTCTAA
                 |||||||||||        |||||||X||||||          |||||||I|||||||          ||||||||
    pattern   CCCTCTAAGAAGTA--------CCCTCTAAGAAGTA----------CCCTCTA AGAAGTA----------CCCTCTAAGAAGTA

    hits         [(0, 10, 0), (19, 32, 1), (43, 57, 1), (68, 75, 0)]
    
## Example:

    $ cargo run -- --query TCTAAGAAGTAGGCGCTCTCCCTCTACGAAGTACTCTAATGAACCCTCTACAGAAGTAGAGTATGTGACCCTCTAA --pattern CCCTCTAAGAAGTA --min_overlap 6
    hits = [(0, 10, 0), (19, 32, 1), (43, 57, 1), (68, 75, 0)]
    $ cargo run -- --query TCTAAGAAGTAGGCGCTCTCCCTCTACGAAGTACTCTAATGAACCCTCTACAGAAGTAGAGTATGTGACCCTCTAA --pattern CCCTCTAAGAAGTA --min_overlap 11
    hits = [(0, 10, 0), (19, 32, 1), (43, 57, 1)]
    $ cargo run -- --query TCTAAGAAGTAGGCGCTCTCCCTCTACGAAGTACTCTAATGAACCCTCTACAGAAGTAGAGTATGTGACCCTCTAA --pattern CCCTCTAAGAAGTA --min_overlap 14
    hits = [(19, 32, 1), (43, 57, 1)]
    $ cargo run -- --query TCTAAGAAGTAGGCGCTCTCCCTCTACGAAGTACTCTAATGAACCCTCTACAGAAGTAGAGTATGTGACCCTCTAA --pattern CCCTCTAAGAAGTA --min_overlap 15
    hits = [(43, 57, 1)]
