// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
//
// Youtao Lu@Kim Lab, 2016-2020

use std::borrow::Borrow;
use std::cmp::{min, max};
use std::iter;
use std::iter::repeat;
use chrono::Local;
use env_logger::{self, Builder};
use getopts::Options;
use log::{debug, error, info, trace, warn, LevelFilter};
use std::env;
use std::io::{Write};
use std::process::exit;

static VERSION: &str = "0.1.1";

fn usage(arg0: &str, opts: Options) {
    let s = format!("\
Summary: 
    Given a query string (length n) and a pattern (length m), output all matches including inexact matches. This extends the
    original Ukkonen's algorithm (O(min(m, n)*k) complexity) but allows partial matches at the 5' or 3' end and mismatches
    under given error rate. 

    Compared with https://github.com/rust-bio/rust-bio/blob/master/src/pattern_matching/ukkonen.rs, it uses error rate instead
    of fixed number of errors and also removes duplicate hits by retaining the best one; 
    Compared with https://github.com/marcelm/cutadapt/blob/main/src/cutadapt/_align.pyx, it outputs all hits instead of the 
    5'-most one.

Usage:
    {} --query text --pattern pattern [--min_overlap 6] [--max_err_rate 0.1] [--debug info] [--help] [--version]

Output:
    a vector of triplets: start, end, cost (i.e. mismatch/insertion/deletion)

Example:
                 0         10        20        30        40        50        60        70
    text         TCTAAGAAGTAGGCGCTCTCCCTCTACGAAGTACTCTAATGAACCCTCTACAGAAGTAGAGTATGTGACCCTCTAA
                 |||||||||||        |||||||X||||||          |||||||I|||||||          ||||||||
    pattern   CCCTCTAAGAAGTA--------CCCTCTAAGAAGTA----------CCCTCTA AGAAGTA----------CCCTCTAAGAAGTA

    hits         [(0, 10, 0), (19, 32, 1), (43, 57, 1), (68, 75, 0)]
", arg0);
    eprintln!("{}", opts.usage(&s));
}

fn init_logger(debug: &str) {
    Builder::new().format(|buf, record| {
        writeln!(buf,
            "[{} {}] {}",
            Local::now().format("%Y-%m-%d %H:%M:%S%.3f %z"),
            record.level(),
            record.args(),
        )
    }).filter(None,
        match debug {
            "error" => LevelFilter::Error, 
            "warn"  => LevelFilter::Warn,
            "info"  => LevelFilter::Info,
            "debug" => LevelFilter::Debug,
            "trace" => LevelFilter::Trace,
            _ => panic!("unknown debug level!"),
        }).init();
}

#[derive(Debug)]
struct Params {
    query: String,
    pattern: String,
    min_overlap: usize,
    max_err_rate: f32,
    debug: String,
}

fn parse_args(args: &Vec<String>, mut opts: Options) -> Params {
    opts.optopt("q", "query", "query text", "STRING");
    opts.optopt("p", "pattern", "pattern", "STRING");
    opts.optopt("o", "min_overlap", "minimum length for a match (default: 6)", "INTEGER");
    opts.optopt("e", "max_err_rate", "maximum error rate (cost over overlap)", "FLOAT");
    opts.optopt("", "debug", "level of debugging info, choose from 'error', 'warn', 'info', 'debug', 'trace'", "");
    opts.optflag("h", "help", "print usage");
    opts.optflag("v", "version", "print version");

    let m = opts.parse(&args[1..]).expect("failed to parse arguments!");
    if m.opt_present("h") {
        usage(&args[0], opts);
        exit(0);
    }

    if m.opt_present("v") {
        println!("v{}", VERSION);
        exit(0);
    }

    let query = match m.opt_str("query") {
        Some(x) => x,
        None => panic!("--query missing!"),
    };
    
    let pattern = match m.opt_str("pattern") {
        Some(x) => x,
        None => panic!("--pattern missing!"),
    };

    let min_overlap: usize = m.opt_get_default("min_overlap", 6).expect("invalid --min_overlap!");

    let max_err_rate: f32 = m.opt_get_default("max_err_rate", 0.1).expect("invalid --max_err_rate!");

    let debug = m.opt_get_default("debug", String::from("info"))
                .expect("invalid --debug, choose from 'info', 'warn', 'error', 'debug', 'trace'");

    Params { query, pattern, min_overlap, max_err_rate, debug }
}
/// Default cost function (unit costs).
pub fn is_match(a: u8, b: u8, is_n_match: bool) -> bool {
    if is_n_match && (a == 'N' as u8 || b == 'N' as u8) { return true }
    else { return a == b }
}

/// a hybrid of
/// https://github.com/rust-bio/rust-bio/blob/master/src/pattern_matching/ukkonen.rs
/// and https://github.com/marcelm/cutadapt/blob/main/src/cutadapt/_align.pyx
pub struct Ukkonen<F>
where
    F: Fn(u8, u8, bool) -> bool,
{
    cost: [Vec<usize>; 2],
    origin: [Vec<isize>; 2],
    matches: [Vec<usize>; 2],
    is_match: F,
    penalty: [usize; 3],
}

impl<F> Ukkonen<F>
where
    F: Fn(u8, u8, bool) -> bool,
{
    /// Initialize algorithm with given capacity and cost function.
    pub fn new(m: usize, is_match: F, penalty: [usize;3]) -> Self {
        let mut ukkonen = Ukkonen {
            cost: [Vec::<usize>::with_capacity(m + 1), Vec::<usize>::with_capacity(m + 1)],
            origin: [Vec::<isize>::with_capacity(m + 1), Vec::<isize>::with_capacity(m + 1)],
            matches: [Vec::<usize>::with_capacity(m + 1), Vec::<usize>::with_capacity(m + 1)],
            is_match,
            penalty,
        };
        ukkonen.cost[0].clear();
        ukkonen.cost[0].extend(repeat(0).take(m + 1));
        ukkonen.origin[0].clear();
        ukkonen.origin[0].extend((-(m as isize)..=0).rev());
        ukkonen.matches[0].clear();
        ukkonen.matches[0].extend(repeat(0).take(m + 1));       
        ukkonen.cost[1].clear();
        ukkonen.cost[1].extend(repeat(0).take(m + 1));
        ukkonen.origin[1].clear();
        ukkonen.origin[1].extend((-(m as isize)..=0).rev());
        ukkonen.matches[1].clear();
        ukkonen.matches[1].extend(repeat(0).take(m + 1));
        return ukkonen;
    }

    /// Find all matches between pattern and text with up to k errors.
    /// Matches are returned as an iterator over pairs of end position and distance.
    pub fn search<'a, C, T>(
        &'a mut self,
        pattern: &'a [u8],
        text: T,
        n: usize,
        k: usize,
        min_overlap: usize,
        max_err_rate: f32,
    ) -> Matches<'_, F, C, T::IntoIter>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        let m = pattern.len();
        Matches {
            ukkonen: self,
            pattern,
            text: text.into_iter().enumerate(),
            lastk: m, 
            m,
            n,
            k,
            min_overlap,
            max_err_rate,
        }
    }
}

/// Iterator over pairs of end positions and distance of matches.
pub struct Matches<'a, F, C, T>
where
    F: Fn(u8, u8, bool) -> bool,
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    ukkonen: &'a mut Ukkonen<F>,
    pattern: &'a [u8],
    text: iter::Enumerate<T>,
    lastk: usize,
    m: usize,
    n: usize,
    k: usize,
    min_overlap: usize,
    max_err_rate: f32,   
}

impl<'a, F, C, T> Iterator for Matches<'a, F, C, T>
where
    F: 'a + Fn(u8, u8, bool) -> bool,
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    type Item = (usize, usize, usize);

    fn next(&mut self) -> Option<(usize, usize, usize)> {
        for (i, c) in &mut self.text {
            let col = i % 2; 
            let prev = 1 - col;
            self.ukkonen.cost[col][0] = 0;
            self.ukkonen.origin[col][0] = i as isize + 1;
            // in the next round of iteration, (lastk+1)-th cost could be the same or even worse
            self.lastk = min(self.lastk + 1, self.m);
            for j in 1..=self.lastk {
                if (self.ukkonen.is_match)(self.pattern[j - 1], *c.borrow(), false) { 
                    // It's a match; go through the diagnoal
                    self.ukkonen.cost[col][j] = self.ukkonen.cost[prev][j-1];
                    self.ukkonen.origin[col][j] = self.ukkonen.origin[prev][j-1];
                    self.ukkonen.matches[col][j] = self.ukkonen.matches[prev][j-1] + 1;
                } else { // not a match; we need to determine whether it is a mismatch, or insertion, or deletion
                    let cost_mismatch = self.ukkonen.cost[prev][j-1] + self.ukkonen.penalty[0];
                    let cost_insertion = self.ukkonen.cost[prev][j] + self.ukkonen.penalty[1]; 
                    let cost_deletion = self.ukkonen.cost[col][j-1] + self.ukkonen.penalty[2];      
                    let cost = min(cost_insertion, min(cost_deletion, cost_mismatch));
                    if cost == cost_mismatch { 
                        // It's a mismatch; go through the diagnoal
                        self.ukkonen.matches[col][j] = self.ukkonen.matches[prev][j-1];
                        self.ukkonen.origin[col][j] = self.ukkonen.origin[prev][j-1];
                    } else if cost == cost_insertion { 
                        // it is an insertion; go from the left horizontal
                        self.ukkonen.matches[col][j] = self.ukkonen.matches[prev][j];
                        self.ukkonen.origin[col][j] = self.ukkonen.origin[prev][j];
                    } else { 
                        // it is a deletion; go from the top vertical
                        self.ukkonen.matches[col][j] = self.ukkonen.matches[col][j-1];
                        self.ukkonen.origin[col][j] = self.ukkonen.origin[col][j-1];
                    }
                    self.ukkonen.cost[col][j] = cost;
                }
            }

            debug!("i = {}, self.k = {}, self.lastk = {}", i, self.k, self.lastk);
            debug!("i = {}, self.ukkonen.cost[prev]   = {:?}", i, self.ukkonen.cost[prev]);
            debug!("i = {}, self.ukkonen.cost[col]    = {:?}", i, self.ukkonen.cost[col]);
            debug!("i = {}, self.ukkonen.origin[prev] = {:?}", i, self.ukkonen.origin[prev]);
            debug!("i = {}, self.ukkonen.origin[col]  = {:?}", i, self.ukkonen.origin[col]);

            // Ukkonen's trick: only record up to lastk rows, where lastk has the max. allowable cost 
            while self.ukkonen.cost[col][self.lastk] > self.k {
                self.lastk -= 1;
            }

            if i == self.n - 1 { 
                // if this is the 3' end (partial pattern at the end of the query string), we need to search upwards for the best (minimal cost, length-qualifying) item
                loop {
                    let l = max(0isize, self.ukkonen.origin[col][self.lastk]) as usize;
                    let qmatches = i + 1 - l;
                    let err_rate = self.ukkonen.cost[col][self.lastk] as f32 / qmatches as f32;
                    debug!("i = {}, c = {}, lastk = {}, l = {}, qmatches = {}, prev.origin = {}, curr.origin = {}, prev.cost = {}, curr.cost = {}, err_rate = {}", i, *c.borrow() as char, self.lastk, l, qmatches, self.ukkonen.origin[prev][self.lastk], self.ukkonen.origin[col][self.lastk], self.ukkonen.cost[prev][self.lastk], self.ukkonen.cost[col][self.lastk], err_rate);
                    if qmatches >= self.min_overlap && err_rate <= self.max_err_rate {
                        return Some((l, i, self.ukkonen.cost[col][self.lastk]));
                    }
                    self.lastk -= 1;
                    if self.ukkonen.cost[col][self.lastk] > self.k || self.lastk == 0 {
                        break;
                    }
                }
            }

            let l = max(0isize, self.ukkonen.origin[col][self.lastk]) as usize;
            let qmatches = i + 1 - l;
            let err_rate = self.ukkonen.cost[col][self.lastk] as f32 / qmatches as f32;
            debug!("i = {}, c = {}, lastk = {}, l = {}, qmatches = {}, prev.cost = {}, curr.cost = {}, err_rate = {}", i, *c.borrow() as char, self.lastk, l, qmatches, self.ukkonen.cost[prev][self.lastk], self.ukkonen.cost[col][self.lastk], err_rate);
 
            // partial pattern at the beginning of the query string or the complete pattern in the middle of the query string
            if self.lastk == self.m && qmatches >= self.min_overlap && err_rate <= self.max_err_rate {
                // we cannot do the following step to reset the cost matrix, because the first match is usually not the best one, e.g. the match can start at 9 and we have a deletion; it can start at 10 and we have a perfect match; it can start at 11 and we have a deletion. If we reset the matrix, we will always miss the other matches than the first one; 
                // self.ukkonen.cost[col].clear();
                // self.ukkonen.cost[col].extend(0..self.m+1);
                // self.ukkonen.cost[prev].clear();
                // self.ukkonen.cost[prev].extend(0..self.m+1);
                // self.ukkonen.origin[col].clear();
                // self.ukkonen.origin[col].extend(repeat(i as isize).take(self.m + 1));
                // self.ukkonen.origin[prev].clear();
                // self.ukkonen.origin[prev].extend(repeat(i as isize).take(self.m + 1));
                // self.lastk = self.k;
                return Some((l, i, self.ukkonen.cost[col][self.lastk]));
            }
        }
        return None;
    }
}

fn dedup_hits(hits: Vec<(usize, usize, usize)>) -> Vec<(usize, usize, usize)> {
    if hits.len() == 0 {
        return hits.clone();
    }
    let mut s: usize = 0;
    let mut c: usize = 0;
    let mut best_i: usize = 0;
    let mut v: Vec<usize> = vec![];
    for (i, hit) in hits.iter().enumerate() {
        if i == 0 {
            s = hit.0;
            c = hit.2;
            continue;
        }
        if hit.0 == s {
            if hit.2 < c {
                c = hit.2;
                best_i = i;
            }
        } else {
            v.push(best_i);
            s = hit.0;
            c = hit.2;
            best_i = i;
        }
    }
    v.push(best_i);
    return v.into_iter().map(|i| hits[i].clone()).collect::<Vec<(usize, usize, usize)>>();
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let params = parse_args(&args, Options::new());

    let query = params.query;
    let pattern = params.pattern;
    let min_overlap = params.min_overlap;
    let max_err_rate = params.max_err_rate;
    let debug = params.debug;
    init_logger(&debug);

    let m = pattern.len();
    let n = query.len();
    let k = (m as f32 * max_err_rate).round() as usize; 
    let mut ukkonen = Ukkonen::new(m, is_match, [1, 1, 1]);

    let hits: Vec<(usize, usize, usize)> = ukkonen.search(pattern.as_bytes(), query.as_bytes(), n, k, min_overlap, max_err_rate).collect();
    let hits = dedup_hits(hits); 
    println!("hits = {:?}", hits);
}
