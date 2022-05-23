#ifndef UNIQ_HPP
#define UNIQ_HPP
#include <iostream>
#include <string>
#include <fstream>
#include <map> 
#include <vector>
#include <tuple>
#include <algorithm>
#include <sstream>
#include <set>
#include <iterator>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <bits/stdc++.h>
#include "pdqsort.h"
#include "robin_hood.h"
using namespace std;
using namespace std::chrono;
typedef robin_hood::unordered_map<unsigned long long, uint32_t> counter_kmer;
void always_call_this();
void rm_nonprinting (std::string& str);
void reverse_c(string& r);

void find_kmers(string& read, vector<string>& read_kmers, int k);

void find_mins(string& read, vector<uint32_t>& cnts, vector<int>& line_mins,  int& k, int& m, int& cov_t);

void update_kmerdb(string& kmer, unordered_map<string, vector<uint32_t>>& kmer_db, int read_number);

void make_kmer_db(vector<string>& seqs, unordered_map<string, vector<uint32_t>>& kmer_db, int k);

void get_counts(string& read, unordered_map<string, vector<uint32_t>>& kmer_db, vector <uint32_t>& cnts_q, int k);

vector<string> split(string line, char delim);

unsigned long long kmerize(const string& s, int x, int k);
void kmerize_all(const string& s, int k, std::vector<unsigned long long>& ret);
void kmerize_all_small(const string& s, int k, std::vector<uint>& ret);

uint kmer_revcomp(uint kmer, int k);

void minmerize(const string& s, int x, int k, long long& kmer);

uint64_t ReverseComp(const uint64_t mer, uint8_t kmerSize);

void turn_kmer_string_to_int(string& kmer, int& kmer_int, int k);

struct modified_hash;

void make_local_kmer_db(int read_size, counter_kmer &kmer_hood, robin_hood::unordered_map<uint32_t, std::vector<unsigned long long>>& seq_kmers, robin_hood::unordered_set<uint32_t>& shared_reads);
void make_local_kmer_db(int read_size, counter_kmer &kmer_hood, robin_hood::unordered_map<uint32_t, std::string>& seqs, robin_hood::unordered_set<uint32_t>& shared_reads, int& k);

void get_local_counts_precomputed_kmers(std::vector<unsigned long long>& seq_kmers, counter_kmer &kmer_hood, vector<uint32_t>& v);

int get_local_counts_with_error_rate(counter_kmer &kmer_hood, string& read, vector<uint32_t>& v, int k);

void get_local_counts(counter_kmer &kmer_hood, string& read, vector <uint32_t>& cnts_q, int k);

void find_kmers_local(string& read, vector<string>& read_kmers, int k, std::unordered_map<string, uint32_t>& kmer_db);

int find_start(vector <string> s);

int find_stop(vector <string> s);

void find_pos_controls(int read_number, vector<string>& que_names, vector <string>& seqs, int ref_start, int ref_stop, unordered_map<string, vector<uint32_t>>& kmer_db, hash<string>& h, vector<string>& pos_names, int& k, int& m, int& cov_t, vector<set <int> >& pos_mins, vector<int>& pos_overlaps, vector<int>& read_starts, vector<int>& read_stops);

uint32_t find_sims(vector <uint32_t>& s1, vector <uint32_t>& s2);

void find_intersection(const vector<uint32_t>& kmer, vector <int>& shared_reads, int& cnt);

//void find_shared_reads(int t, std::vector<int>& shared_reads, std::vector<std::string>& seqs, std::vector<uint32_t>& ref_counts,std::vector<uint32_t>& rev_counts , int k, std::hash<std::string>& h, int m, int cov_t, int read_number, std::vector <vector <int>>& min_db, std::string& read, std::string& read_rev);

#endif
