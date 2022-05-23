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
#include "uniq.hpp"
#include <unordered_map>
#include <unordered_set>
#include <bits/stdc++.h>
#include "pdqsort.h"
#include "robin_hood.h"

using namespace std; 
typedef robin_hood::unordered_map<unsigned long long, uint32_t> counter_kmer;
using namespace std::chrono;

int vals[256];
//char ALPHABET[4];

void always_call_this()   
{
    for (int i=0; i<256; i++) vals[i]=0;
    vals['A']=0;
    vals['C']=1;
    vals['G']=2;
    vals['T']=3;
	//ALPHABET['A'] = 'T';
	//ALPHABET['T'] = 'A';
	//ALPHABET['G'] = 'C';
	//ALPHABET['C'] = 'G';
	
}

//char ALPHABET[4] 

void reverse_c(string& r)
{
	//string r_string;
	reverse(r.begin(), r.end());
	for (std::size_t i = 0; i < r.length(); ++i){
        	switch (r[i]){
        	case 'A':
            		r[i] = 'T';
            		break;    
        	case 'C':
            		r[i] = 'G';
            		break;
        	case 'G':
            		r[i] = 'C';
            		break;
        	case 'T':
            		r[i] = 'A';
           		break;
        	}
		//r[i] = ALPHABET[r[i]];
    }
	//return s; 
}

string to_numbers(string& s)
{	
	//int the_number = 0;
	string tmp_s = s;
	for (int i = 0; i < s.size(); i++) {
		switch (s[i]){
                	case 'A':
                        	tmp_s[i] = '0';
                        	break;
                	case 'C':
                        	tmp_s[i] = '1';
                        	break;
                	case 'G':
                        	tmp_s[i] = '2';
                        	break;
                	case 'T':
                        	tmp_s[i] = '3';
                        	break;
			default :
				tmp_s[i] = '4';
				break;
		}
	}
	return tmp_s;
}

void find_kmers(string& read, vector<string>& read_kmers, int k)
{
	//vector<string> read_kmers;
	//std::hash<std::string> h;	
	for (int x=0; x < (read.length()-k+1); x++)
	{
		string kmer = read.substr(x, k);
		//string kmer; kmer = to_numbers(kmer_tmp);
		//kmer_tmp.clear();
		//vector<string> tmp_arr = 
		read_kmers.push_back(kmer);
		//Come back here and add the RC of the string. Right now that seems too complicated 
		//string r_kmer = reverse_c(kmer);
		//read_kmers.push_back(r_kmer);
		
	}

	//return read_kmers;
}

void find_kmers_local(string& read, vector<string>& read_kmers, int k, std::unordered_map<string, uint32_t>& kmer_db)
{
        //vector<string> read_kmers;

        for (int x=0; x < (read.length()-k+1); x++)
        {
                string kmer; kmer = read.substr(x, k);
                read_kmers.push_back(kmer);
                //Come back here and add the RC of the string. Right now that seems too complicated 
                string r_kmer = kmer;
	       	reverse_c(r_kmer);
                read_kmers.push_back(r_kmer);
                kmer_db[r_kmer] = 0;
                kmer_db[kmer] = 0;
        }

        //return read_kmers;
}

void find_mins(string& read, vector<uint32_t>& cnts, vector<int>& line_mins, int& k, int& m, int& cov_t)
{
	std::hash<std::string> h;
	//cout << "Looking for minimizers\n";
	//vector<string> line_mins;
	//line_mins.push_back(0);
	for (size_t x=0; x <= (read.length() - k) ; x+=15)
	{
		string forward; 
		forward = read.substr(x, k);
		//static_cast<std::uint32_t>(min);
		int min = 1000000000; //h("Vvvvvoyager");
		//static_cast<std::uint32_t>(min);
		for (int j=0; j < k-m ; j++)
		{
			int sub; sub=static_cast<int>(h(forward.substr(j, m))) ;
			if (sub < min)
			{
				//std::cout << "At least the hashed value is smaller :) \n" ;
				//std::cout << cnts[x] << "\n";
				//std::cout.flush();
				if (cnts[x] < cov_t)
				{
					if (cnts[x] > 3) {
						//std::cout << "At least the hashed value is smaller :) \n" ;
                                		//std::cout << cnts[x] << "\n";
                                		//std::cout.flush();
						min = sub;
					}
			
				}
			}
			
		}
		if (min != 1000000000 ) { //static_cast<std::uint32_t>(h("Vvvvvoyager"))){
			//std::cout << "We changed the min value! \n";
			//std::cout.flush();
			
			//if (std::find(line_mins.begin(), line_mins.end(), min)==line_mins.end()){
				
				//std::cout << "And we addedthe min value to the line mins\n" ;
				//std::cout.flush();
				line_mins.push_back(min);
			//}
			//line_mins.push_back(min); //Changed insert to pushback
		}
	}
	//std::cout << line_mins.size() << "\n" ;
	//std::cout.flush();

	//sort(line_mins.begin(), line_mins.end());
	/*if (line_mins.size() < 2) {
		std::cout << "THERE ARE NO MINIMIZERS :( \n";
		std::cout.flush();
	}*/

	//return line_mins;
}

void update_kmerdb(string& kmer, unordered_map<string, vector<uint32_t>>& kmer_db, int read_number)
{
	if (kmer_db.count(kmer))
	{
		vector<uint32_t> tmp = kmer_db[kmer];
		tmp.push_back(read_number);
		kmer_db[kmer] = tmp;
	}
	else 
	{
		vector<uint32_t> tmp = {read_number};
		kmer_db[kmer] = tmp;
	}
	
}

void make_kmer_db(vector<string>& seqs, unordered_map<string, vector<uint32_t>>& kmer_db, int k)
{
	//map<string, int> kmer_db;
	//What about finding all of the kmers, then sorting, then countint up all of  the kmers? 	
	/*for (int x=0; x < seqs.size(); x++){
		vector<int> read_kmers;
		find_kmers(seqs[x], read_kmers, k);
		//cout << "Found the kmers!\n";
	
		for (int j=0; j < read_kmers.size() ; j++){
			//cout << j << "\n";
			update_kmerdb(read_kmers[j], kmer_db, x);
		}
	}*/
	/*vector<string> read_kmers;
	for  (int x=0; x < seqs.size(); x++){
		for (int y=0; y<(seqs[x].size()-k+1), y++){
			//
			string kmer; kmer = seqs[x].substr(y, k);
                	read_kmers.push_back(kmer);
                	//Come back here and add the RC of the string. Right now that seems too complicated
                	string r_kmer = reverse_c(kmer);
               		read_kmers.push_back(r_kmer);
		}
	}
	//Sort
	sort(read_kmers.begin(), read_kmers.end())	
	//Go through and make the counts 
	string current_kmer;
	string last_kmer = read_kmers[0];
	int cntr=0;
	for (int j=0; j<read_kmers.size(); j++){
		if (read_kmers[j] == last_kmer)
	}*/
	//return kmer_db;
}

void get_counts(string& read, unordered_map<string, vector<uint32_t>>& kmer_db, vector <uint32_t>& cnts_q, int k)
{
	//vector <int> cnts_q;
	for (int x=0; x < (read.length() - k + 1  ) ; x++){
		string kmer;
		kmer = read.substr(x, k);

		int cnt = kmer_db[kmer].size();
		cnts_q.push_back(cnt);
	}
	/*rev_read = reverse_c(read);
	for (int x=0; x < (rev_read.length() - k + 1) ; x++){
                string kmer;
                kmer = rev_read.substr(x, k);

                int cnt = kmer_db[kmer].size();
                cnts_q.push_back(cnt);
        }*/
	//return cnts_q;
}

void find_intersection(const vector<uint32_t>& kmer, vector <int>& shared_reads, int& cnt){

	cnt = 0;
	for (uint32_t ele: shared_reads){
		if (std::binary_search (kmer.begin(), kmer.end(), ele)){
                        ++cnt;
                }
	}
	/*for (auto ele: shared_reads){
		std::cout << ele << "\n";
	}*/
	



}

uint64_t ReverseComp(const uint64_t mer, uint8_t kmerSize)
{
    uint64_t res = ~mer;

    res = ((res >> 2 & 0x3333333333333333) | (res & 0x3333333333333333) << 2);
    res = ((res >> 4 & 0x0F0F0F0F0F0F0F0F) | (res & 0x0F0F0F0F0F0F0F0F) << 4);
    res = ((res >> 8 & 0x00FF00FF00FF00FF) | (res & 0x00FF00FF00FF00FF) << 8);
    res = ((res >> 16 & 0x0000FFFF0000FFFF) | (res & 0x0000FFFF0000FFFF) << 16);
    res = ((res >> 32 & 0x00000000FFFFFFFF) | (res & 0x00000000FFFFFFFF) << 32);

    return (res >> (2 * (32 - kmerSize)));
}

void characterT0Bits(char character)
{
    int vals[256];
    for (int i=0; i<256; i++) vals[i]=0;
    vals['A']=0;
    vals['C']=1;
    vals['G']=2;
    vals['T']=3;
}

//long long kmerize(const string& s, int start, int k)
//{
	/*
  long long kmer = 0;
	for(size_t i = start; i<(size_t)(k) + start; i++)
  {
    kmer = (kmer << 2) | vals[s[i]];
  }
	return kmer;*/

//	int alpha=2;
//	long long kmer = 0;
//    	for(int i = 0; i<k; i++) kmer = (kmer << alpha) | vals[(size_t)s[i]];
//    	return kmer;
//}

void minmerize(const string& s, int x, int k, long long& kmer)
{
        //k = 31;
        int alpha = 2;
        //int kmer = 0;
        //for(int i = 0; i<k; i++) kmer = (kmer << 2) | vals[s[i]];

        //std::cout << kmer << "\n";
        //std::cout.flush();
        //return kmer;

        kmer = 0;
        //for(int i = 0; i<k; i++) kmer = (kmer << alpha) | vals[(size_t)s[i]];
        for(int i = 0; i<k; i++) kmer = (kmer << alpha) | vals[(size_t)s[i+x]];
        
	//return kmer;

}


//based on https://www.biostars.org/p/113640/
unsigned long long kmer_revcomp(unsigned long long kmer, int k) {
    kmer = ~kmer;
    kmer = ( (kmer >> 2  & 0x3333333333333333) | (kmer & 0x3333333333333333) << 2 );
    kmer = ( (kmer >> 4  & 0x0F0F0F0F0F0F0F0F) | (kmer & 0x0F0F0F0F0F0F0F0F) << 4 );
    kmer = ( (kmer >> 8  & 0x00FF00FF00FF00FF) | (kmer & 0x00FF00FF00FF00FF) << 8 );
    kmer = ( (kmer >> 16 & 0x0000FFFF0000FFFF) | (kmer & 0x0000FFFF0000FFFF) << 16 );
    kmer = ( (kmer >> 32 & 0x00000000FFFFFFFF) | (kmer & 0x00000000FFFFFFFF) << 32 );
    return kmer >> (2 * (32 - k));
}
uint kmer_revcomp(uint kmer, int k) {
    kmer = ~kmer;
    kmer = ( (kmer >> 2  & 0x3333333333333333) | (kmer & 0x3333333333333333) << 2 );
    kmer = ( (kmer >> 4  & 0x0F0F0F0F0F0F0F0F) | (kmer & 0x0F0F0F0F0F0F0F0F) << 4 );
    kmer = ( (kmer >> 8  & 0x00FF00FF00FF00FF) | (kmer & 0x00FF00FF00FF00FF) << 8 );
    kmer = ( (kmer >> 16 & 0x0000FFFF0000FFFF) | (kmer & 0x0000FFFF0000FFFF) << 16 );
    //kmer = ( (kmer >> 32 & 0x00000000FFFFFFFF) | (kmer & 0x00000000FFFFFFFF) << 32 );
    return kmer >> (2 * (16 - k));
}

void kmerize_all(const string& s, int k, std::vector<unsigned long long>& ret) {
    //std::vector<unsigned long long> ret;
    unsigned long long KMER_MASK = (static_cast<unsigned long long>(1) << (2*k)) - 1;
    int alpha = 2;
    unsigned long long kmer = 0;
    for(int i = 0; i < s.size(); i++) {
        kmer = ((kmer << alpha) & KMER_MASK) | vals[(size_t) s[i]];

        if (i >= k-1) {
            auto rev = kmer_revcomp(kmer, k);
            if (kmer < rev) {
                ret.push_back(kmer);
            } else {
                ret.push_back(rev);
            }
        }
    }
	//return ret;
}

void kmerize_all_small(const string& s, int k, std::vector<uint>& ret) {
    //std::vector<unsigned long long> ret;
    uint KMER_MASK = (static_cast<uint>(1) << (2*k)) - 1;
    int alpha = 2;
    uint kmer = 0;
    for(int i = 0; i < s.size(); i++) {
        kmer = ((kmer << alpha) & KMER_MASK) | vals[(size_t) s[i]];

        if (i >= k-1) {
            auto rev = kmer_revcomp(kmer, k);
            if (kmer < rev) {
                ret.push_back(kmer);
            } else {
                ret.push_back(rev);
            }
        }
    }
        //return ret;
}


unsigned long long kmerize(const string& s, int x, int k) {
    int alpha = 2;
    unsigned long long kmer = 0;
    for(int i = 0; i < k; i++) {
        kmer = (kmer << alpha) | vals[(size_t)s[i+x]];
    }
	return kmer;
}


void turn_kmer_string_to_int(string& kmer, int& kmer_int, int k)
{
    //int kmer_int = 0;
    for (int y=0; y< kmer.length(); ++y) {
        char character = kmer[y];
        int bits = 0;//kmerize(character, k);

        kmer_int = kmer_int << 2;
        kmer_int = kmer_int|bits;
    }
}

void get_local_counts_precomputed_kmers(std::vector<unsigned long long>& seq_kmers, counter_kmer &kmer_hood, vector<uint32_t>& v) 
{
	
	//if (read.length() > 10){
        //        v.reserve(read.size());
                for (auto ele: seq_kmers) {
                        v.push_back(kmer_hood[ele]);
                }


        //} else {
         //       v.push_back(0);
       // }
}

void rm_nonprinting (std::string& str)
{
    str.erase (std::remove_if (str.begin(), str.end(),
                                [](unsigned char c){
                                    return !std::isprint(c);
                                }),
                                str.end());
}

//adjacent_kmer, local_kmer_db, v, k, shared_reads
void get_local_counts(counter_kmer &kmer_hood, string& read, vector<uint32_t>& v, int k)
{		
	//std::cout << "About to get counts \n";
        //std::cout.flush();
	uint32_t cnt1 = 0;
	v.clear();
	rm_nonprinting(read);
	if (read.size() > (k-1)){
		v.reserve(read.size());	
		std::vector<unsigned long long> read_counts;
		kmerize_all(read, k, read_counts);
		for (int i = 0; i < (read.length()-k+1); i++) {
			v.push_back(kmer_hood[read_counts[i]]);
		}

		
	} else {
		v.push_back(cnt1);
	}
}

int get_local_counts_with_error_rate(counter_kmer &kmer_hood, string& read, vector<uint32_t>& v, int k)
{
        //std::cout << "About to get counts \n";
        //std::cout.flush();
        uint32_t cnt1 = 0;
        v.clear();
	int error_rate = 0;
        rm_nonprinting(read);
        if (read.size() > (k-1)){
                v.reserve(read.size());
                std::vector<unsigned long long> read_counts;
                kmerize_all(read, k, read_counts);
                for (int i = 0; i < (read.length()-k+1); i++) {
                        v.push_back(kmer_hood[read_counts[i]]);
			if (kmer_hood[read_counts[i]] < 5 ){
				error_rate = error_rate + 1;
			}
                }


        } else {
                v.push_back(cnt1);
        }
	return error_rate;
}


void make_local_kmer_db(int read_size, counter_kmer &kmer_hood, robin_hood::unordered_map<uint32_t, std::vector<unsigned long long>>& seq_kmers, robin_hood::unordered_set<uint32_t>& shared_reads)
{
	//std::cout << shared_reads.size() << "\n";
	
	kmer_hood.reserve(read_size*2);
	for (int y: shared_reads){
		for (auto x: seq_kmers[y]) {
			++kmer_hood[x];
		}
	}
	
}

void make_local_kmer_db(int read_size, counter_kmer &kmer_hood, robin_hood::unordered_map<uint32_t, std::string>& seqs, robin_hood::unordered_set<uint32_t>& shared_reads, int& k)
{
        //std::cout << shared_reads.size() << "\n";

        kmer_hood.reserve(read_size*2);
        for (int y: shared_reads){
                std::vector<unsigned long long> read_counts;
		kmerize_all(seqs[y], k, read_counts);
		for (auto x: read_counts) {
                        ++kmer_hood[x];
                }
        }

}

vector<string> split(string line, char delim)
{
	//https://stackoverflow.com/questions/20755140/split-string-by-a-character
	//cout << line << "\n";
	vector<string> arr_split;
	//replace(line.begin(), line.end(), delim, ' '); 
	//cout << line << "\n";
	stringstream ss(line);
	string item;

	while (getline (ss, item, delim)){
		arr_split.push_back(item);
		//cout << item << "\n";
	}
	return arr_split;
}


struct modified_hash {
 
    static uint64_t splitmix64(uint64_t x)
    {
        x += 0x9e3779b97f4a7c15;
        x = (x ^ (x >> 30))
            * 0xbf58476d1ce4e5b9;
        x = (x ^ (x >> 27))
            * 0x94d049bb133111eb;
        return x ^ (x >> 31);
    }
 
    int operator()(uint64_t x) const
    {
        static const uint64_t random
            = steady_clock::now()
                  .time_since_epoch()
                  .count();
        return splitmix64(x + random);
    }
};


/*int find_start(vector <string> s)
{
	int start;

	//vector <string> s; s = split(header, '|');
        string strand; strand = s[2];
        if (strand == "forward") {
		string tmp; tmp = s[3];
                start = std::stoi(tmp);
        } else {
		string tmp2; tmp2 = s[5];
                int len; len = std::stoi(tmp2);
		string tmp3; tmp3 = s[3];
                int tmp4; tmp4 = std::stoi(tmp3);
                int stop; stop = 1000000 - tmp4;
		start = stop - len;
        }

	return start;
}

int find_stop(vector <string> s)
{
        int stop;
	//vector <string> s; s = split(header, '|'); 
	string strand; strand = s[2];
	if (strand == "forward") {
		string tmp1; tmp1 = s[3];
		int start; start = std::stoi(tmp1);
		string tmp2; tmp2 = s[5];
		int len; len = std::stoi(tmp2); //TODO 
		stop = start + len;
	} else { 
		string tmp3; tmp3 = s[5];
		int len; len = std::stoi(tmp3);
		string tmp4; tmp4 = s[3] ;
		int start; start = std::stoi(tmp4);
		stop = 1000000 - start;
	}

        return stop;
}

void find_pos_controls(int read_number, vector<string>& que_names, vector <string>& seqs, int ref_start, int ref_stop, unordered_map<string, vector<uint32_t>>& kmer_db, hash<string>& h, vector<string>& pos_names, int& k, int& m, int& cov_t, vector<set <int> >& pos_mins, vector<int>& pos_overlaps, vector<int>& read_starts, vector<int>& read_stops) 
{
	//cout << que_names[0] << "\n";		int this_stop;
		//cout << "About to check the strand \n";
		if (this_strand == "forward") { 
			//
			this_start = std::stoi(rd[3]);
			this_stop = this_start + std::stoi(this_len);
		} else {
			//
			this_stop = 1000000 - std::stoi(rd[3]);
			this_start = this_stop - std::stoi(this_len);
		}
		//cout << "Found the strand \n";
		//Calculate the overlap 
		if (ref_start <= this_start) {
			if (this_start <= ref_stop) {
				//
				int overlap; overlap = ref_stop - this_start;
				pos_overlaps.push_back(overlap);
				pos_names.push_back(r);
				vector<uint32_t> cnts_q;
				get_counts(seqs[i], kmer_db, cnts_q, k);
				vector<uint32_t> this_min; 
				find_mins(seqs[i], cnts_q, this_min, h, k, m, cov_t);
				//pos_mins.push_back(this_min);
				//read_stops.push_back(this_stop);
				//read_starts.push_back(this_start);
			}
		} else if (ref_start <= this_stop) {
			if (this_stop <= ref_stop) {
				//
				int overlap; overlap = this_stop - ref_start;
				pos_overlaps.push_back(overlap);
				pos_names.push_back(r);
				vector<uint32_t> cnts_q;
				get_counts(seqs[i], kmer_db, cnts_q, k);
				vector<uint32_t> this_min;
				find_mins(seqs[i], cnts_q, this_min, h, k, m, cov_t);
				//pos_mins.push_back(this_min);
				//read_stops.push_back(this_stop);
				//read_starts.push_back(this_start);
			}
		}
		//cout << "Overlap calculated \n";
		
	}
	//cout << "Size of Pos's: " << pos_names.size() <<"\n";
	return void();
}
*/
uint32_t find_sims(vector <uint32_t>& s1, vector <uint32_t>& s2){
	
	//int total_sim; 
	//set<string> s1(v1.begin(), v1.end());
	//set<string> s2(v2.begin(), v2.end());
	//std::multiset<int> s3;
	//std::vector<int> s3;
	//set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(), std::inserter(s3, s3.begin())); //std::back_inserter(s3));
	//set_intersection(s1.begin(),s1.end(),s2.begin(),s2.end(), std::back_inserter(s3));
	uint32_t cnt = 0;
        for (uint32_t ele: s2){
		//if (std::find(s2.begin(), s2.end(), ele) != s2.end()) { //(shared_reads.count(ele)){
                //        ++cnt;
                //}
		//Binary search, this should be the fastest? 
        	if (std::binary_search (s1.begin(), s1.end(), ele)){	
			++cnt;
		}
	}
	
	
	//total_sim = s3.size();
	
	return cnt;
}














