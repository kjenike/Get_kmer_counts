#include "./../kmc_api/kmc_file.h"
#include "./../kmc_api/kmer_api.h"
#include "nc_utils.h"
#include "stdafx.h"
#include <algorithm>
#include <bits/stdc++.h>
#include <ctime>
#include <fstream>
#include <iostream>
#include <list>
#include <math.h>
#include <mutex>
#include <numeric>
#include <stdlib.h>
#include <thread>
#include <unistd.h>
#include <vector>
#include "uniq.hpp"
#include <unordered_map>
#include <chrono>
#include "pdqsort.h"
#include "robin_hood.h"

using namespace std::chrono;
typedef robin_hood::unordered_flat_map<unsigned long long, std::vector<uint32_t>> counter_t;
typedef robin_hood::unordered_flat_map<uint32_t, std::vector<unsigned long long>> read_min_map;
//typedef std::unordered_map<long long, std::vector<uint32_t>> counter_t;
typedef robin_hood::unordered_flat_map<unsigned long long,uint32_t> counter_kmer;
typedef robin_hood::unordered_flat_map<uint32_t, std::vector<uint32_t>> allvall;
typedef robin_hood::unordered_flat_map<uint32_t, std::vector<uint32_t>> all_v_all_type;

#define DEBUG_OUT

//smooth kmcdb reads.fastq errremoved.fasta erredits.fasta errpaths.fasta hetremoved.fasta hetedits.fasta hetpath1.fasta hetpath2.fasta hetpath3.fasta hetpath4.fasta hetpath5.fasta hetpath6.fasta error_threshold het_threshold unique_threshold anchor_threshold allowed_err_fraction allowed_rep_fraction max_nodes_to_search distance_multiplier strict > counts.txt
robin_hood::unordered_flat_map<uint32_t, std::string> seqs;
std::vector<std::string> headers;
//robin_hood::unordered_flat_map<uint32_t, std::vector<unsigned long long>> seq_kmers;
std::vector<std::string> seqs_r;
std::vector<all_v_all_type> all_v_all;
//int k = 0;
int shared_read_t = 10 ; //25;
int k_arr_size = 92500;

auto s_test = chrono::high_resolution_clock::now();
auto e_test = chrono::high_resolution_clock::now();
auto section1 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section2 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section3 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section4 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section5 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section6 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section7 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section8 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section9 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section10 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section11 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section12 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section13 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section14 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section15 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section16 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section17 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section18 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section19 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section20 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section21 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);
auto section22 = chrono::duration_cast<chrono::nanoseconds>(e_test-s_test);

int get_type(uint32_t& coverage, int& error_threshold, int& het_threshold, int& unique_threshold)
{
    if (coverage <= error_threshold)
    {
        return 0;
    }
    else if (coverage <= het_threshold)
    {
        return 1;
    }
    else if (coverage <= unique_threshold)
    {
        return 2;
    }
    else
    {
        return 3;
    }

}

std::vector<std::string> get_adjacent(int& k, counter_kmer &kmer_hood, int shared_reads_l, CKMCFile& file, std::string& kmer, int& error_threshold, int& het_threshold, int& unique_threshold, bool& going_right)
{
    //std::cout << "Getting adjacent\n";
    //std::cout.flush();

    std::vector<uint32_t> v;
    std::vector<std::string> adjacent_kmers;
    //for all possible nucleotide extensions from kmer
    //std::cout << kmer << "\n";
    //std::cout.flush();
    char nucs[] = "ACGT";
    for (int p=0; p < 5; p++) //(char const &c: nucs)
    {
        char c = nucs[p];
        //std::cout << c << "\n";
        //std::cout.flush();    
        
        std::string adjacent_kmer;
        if (going_right)
        {
            
            adjacent_kmer = kmer.substr(1)+c;
            //std::cout << "KMER: " << adjacent_kmer << "\n";
                        //std::cout.flush();
        }
        else
        {
            //cout << "Going left\n";
            //cout << kmer << "\n";
            //cout << kmer.length() << "\n";
            //cout << c << "\n";
            adjacent_kmer = c+kmer.substr(0, kmer.length()-1);
            //std::cout << "KMER: " << adjacent_kmer << "\n";
            //std::cout.flush();
        }
        //std::cout << kmer << "\n" ; 
        //std::cout.flush();
        v.clear();
        //v.shrink_to_fit();
        if (shared_reads_l == 1) { //(shared_reads.size() > shared_read_t) {        
            //std::cout << adjacent_kmer << "\n";
            //std::cout.flush();
            //get_local_counts(adjacent_kmer, local_kmer_db, v, k);//TODO
            //file.GetCountersForRead(adjacent_kmer, v);
            //std::cout << v[0] << "\n";
            //v.clear();
            get_local_counts(kmer_hood, adjacent_kmer, v, k);
            //std::cout << v[0] << "\n";
            //std::cout << "****************************\n";
            //std::cout.flush();
        } else {
            file.GetCountersForRead(adjacent_kmer, v);
        }
        if (v.size() == 0){
            v.push_back(0);
        }
        //std::cout << v.size() << "\n" ;
        //std::cout.flush();
        //std::cout << kmer << "\t" << v[0] << "\n" ; 
        /*for (uint32_t something: v) {
            std::cout << something << "\n" ;
        }*/
        int current_type = get_type(v[0], error_threshold, het_threshold, unique_threshold);
        
        //if the adjacent kmer is not an error
        if ((current_type > 0)) //&& (current_type < 3)) //Added the 2nd condition 
        {
          adjacent_kmers.push_back(adjacent_kmer);
        }
    }
    /*for (auto mer: adjacent_kmers){
        std::cout << mer << "\n";
        std::cout.flush();
    }*/
    return adjacent_kmers;
}

bool is_left_anchor(int& k, counter_kmer &kmer_hood, int shared_reads_l, CKMCFile& file, std::string& previous_kmer, int& previous_count, int& error_threshold, int& het_threshold, int& unique_threshold, int& anchor_threshold)
{
    if ((previous_kmer.length() != k) || (previous_count > unique_threshold)) // > unique_threshold or anchor_threshold
    //if (previous_kmer.length() != k)
    {
        return false;
    }
    bool going_right = true;
    std::vector<std::string> adjacent_kmers = get_adjacent(k, kmer_hood, shared_reads_l, file, previous_kmer, error_threshold, het_threshold, unique_threshold, going_right);
    if (adjacent_kmers.size() == 2)
    {
        std::vector<uint32_t> v1;
        //file.GetCountersForRead(adjacent_kmers[0], v1);
        if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
            get_local_counts(kmer_hood, adjacent_kmers[0], v1, k);
            //file.GetCountersForRead(adjacent_kmers[0], v1);
        } else {
                        file.GetCountersForRead(adjacent_kmers[0], v1);
                }
        int count1 = v1[0];
        std::vector<uint32_t> v2;
        //file.GetCountersForRead(adjacent_kmers[1], v2);
        if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
            //file.GetCountersForRead(adjacent_kmers[1], v2);
            get_local_counts(kmer_hood, adjacent_kmers[1], v2, k);
        } else {
                        file.GetCountersForRead(adjacent_kmers[1], v2);
                }
        int count2 = v2[0];
        //If the sum of the coverages of the two branches is within 3 of the homozygous portion
        if ((previous_count - 3 <= count1 + count2) && (count1 + count2 <= previous_count + 3))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

bool is_right_anchor(int& k, counter_kmer &kmer_hood, int shared_reads_l, CKMCFile& file, std::string& current_kmer, int& current_count, int& error_threshold, int& het_threshold, int& unique_threshold, int& anchor_threshold)
{
    if ((current_kmer.length() != k) || (current_count > unique_threshold)) // > unique_threshold or anchor_threshold?
    //if (current_kmer.length() != k)
    {
        return false;
    }
    bool going_right = false;
    std::vector<std::string> adjacent_kmers = get_adjacent(k, kmer_hood, shared_reads_l, file, current_kmer, error_threshold, het_threshold, unique_threshold, going_right);
    if (adjacent_kmers.size() == 2)
    {
        std::vector<uint32_t> v1;
        //file.GetCountersForRead(adjacent_kmers[0], v1);
        if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
            //file.GetCountersForRead(adjacent_kmers[0], v1);
            get_local_counts(kmer_hood, adjacent_kmers[0], v1, k);
        } else {
                        file.GetCountersForRead(adjacent_kmers[0], v1);
                }
        int count1 = v1[0];
        std::vector<uint32_t> v2;
        //file.GetCountersForRead(adjacent_kmers[1], v2);
        if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
            //file.GetCountersForRead(adjacent_kmers[1], v2);
            get_local_counts(kmer_hood, adjacent_kmers[1], v2, k);
        } else {
                        file.GetCountersForRead(adjacent_kmers[1], v2);
                }
        int count2 = v2[0];
        //If the sum of the coverages of the two branches is within 3 of the homozygous portion
        if ((current_count - 3 <= count1 + count2) && (count1 + count2 <= current_count + 3))
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

int get_type_het(int& k, counter_kmer &kmer_hood, int shared_reads_l, CKMCFile& file, int& previous_type, std::string& previous_kmer, std::string& current_kmer, int& previous_count, int& current_count, int& error_threshold, int& het_threshold, int& unique_threshold, int& anchor_threshold, std::string& anchor_found)
{
    //If "previously" at the beginning of the read
    if (previous_type == -1)
    {
        //if current kmer is a right anchor
        if ((current_count > (previous_count + error_threshold)) && is_right_anchor(k, kmer_hood, shared_reads_l, file, current_kmer, current_count, error_threshold, het_threshold, unique_threshold, anchor_threshold))
        //((current_count > (previous_count + error_threshold)) && is_right_anchor(shared_reads, current_kmer, current_count, k, local_kmer_db, error_threshold, het_threshold, unique_threshold, anchor_threshold))
        {
            
            anchor_found = "right";
            return 2;
        }
        //if current kmer is a left anchor
        else if (is_left_anchor(k, kmer_hood, shared_reads_l, file, current_kmer, current_count, error_threshold, het_threshold, unique_threshold, anchor_threshold))
        {
            
            anchor_found = "left";
            return 2;
        }
        //To keep it consistent with previous get_type, we need type to be -1 only before we start the read
        //So, in this case we will do what was done before, using the coverage and coverage thresholds to
        //determine whether the kmer is homozygous or not
        else if ((het_threshold < current_count) && (current_count <= unique_threshold))
        {
            anchor_found = "none";
            return 2;
        }
        else
        {
            anchor_found = "none";
            return 1;
        }
    }
    //If previously in a nonhomozygous region
    else if (previous_type == 1)
    {
        //If current kmer is a right anchor
        if ((current_count > (previous_count + error_threshold )) && is_right_anchor(k, kmer_hood, shared_reads_l, file, current_kmer, current_count, error_threshold, het_threshold, unique_threshold, anchor_threshold))
        //((current_count > (previous_count + error_threshold )) && is_right_anchor(shared_reads, current_kmer, current_count, k, local_kmer_db, error_threshold, het_threshold, unique_threshold, anchor_threshold))
        {
            //we have left the nonhom region
            
            anchor_found = "right";
            return 2;
        }
        //If current kmer is a left anchor
        else if (is_left_anchor(k, kmer_hood, shared_reads_l, file, current_kmer, current_count, error_threshold, het_threshold, unique_threshold, anchor_threshold))
        {
            //this is a weird case where we find two left anchors in a row before a right anchor
            
            anchor_found = "left";
            return 2;
        }
        //else current kmer is not an anchor
        else
        {
            //we continue the nonhom region
            anchor_found = "none";
            return 1;
        }
    }
    //If previously in a homozygous region
    //else if (previous_type == 2)
    else
    {
        //If current kmer is a right anchor
        //if ((current_count > (previous_count + error_threshold)) && is_right_anchor(current_kmer, current_count, k, file, error_threshold, het_threshold, unique_threshold, anchor_threshold))
        //{
            //this is a weird case where we find two right anchors in a row before a left anchor
            //we don't need to set anchor_found or return 2 because we really care
            //if current kmer is a left anchor, ending the homozygous region
            
        //}
        //If current kmer is a left anchor
        if (is_left_anchor(k, kmer_hood, shared_reads_l, file, current_kmer, current_count, error_threshold, het_threshold, unique_threshold, anchor_threshold))
        {
            
            anchor_found = "left";
            return 2;
        }
        //If previous kmer is a left anchor
        else if ((current_count < (previous_count - error_threshold )) && is_left_anchor(k, kmer_hood, shared_reads_l, file, previous_kmer, previous_count, error_threshold, het_threshold, unique_threshold, anchor_threshold))
        //((current_count < (previous_count - error_threshold )) && is_left_anchor(shared_reads, previous_kmer, previous_count, k, local_kmer_db, er_k, local_kmer_db_cror_threshold, het_threshold, unique_threshold, anchor_threshold))
        {
            //we have left the hom region
            
            anchor_found = "left";
            return 1;
        }
        //else previous kmer is not a left anchor (and current kmer is not a left anchor)
        else
        {
            //we continue the hom region
            anchor_found = "none";
            return 2;
        }
    }
}

std::vector<std::string> get_paths(int& k, counter_kmer &kmer_hood, int shared_reads_l, CKMCFile& file, int& error_threshold, int& het_threshold, int& unique_threshold, std::string& left_anchor_kmer, std::string& right_anchor_kmer, int& min_distance_of_path, int& max_distance_of_path, int& max_nodes_to_search, bool& queue_broken)
{
    //std::cout << "Getting paths\n";
    //std::cout.flush();

    std::string starting_anchor_kmer;
    bool going_right;
    
    if (left_anchor_kmer.empty())
    {
        //We are at the beginning of the read.
        //This function will find paths starting from right_anchor_kmer continuing left to the
        //beginning of the read where each kmer in the path is a nonerror kmer.
        starting_anchor_kmer = right_anchor_kmer;
        going_right = false;
    }
    else if (right_anchor_kmer.empty())
    {
        //We are at the end of the read.
        //This function will find paths starting from left_anchor_kmer continuing right to the
        //end of the read where each kmer in the path is a nonerror kmer.
        starting_anchor_kmer = left_anchor_kmer;
        going_right = true;
    }
    else
    {
        //We are in the middle of the read.
        //This function will find paths starting from left_anchor_kmer continuing right until
        //right_anchor_kmer where each kmer in the path is a nonerror kmer.
        starting_anchor_kmer = left_anchor_kmer;
        going_right = true;
    }
    //std::cout << "We are either going left or right\n";
    //std::cout.flush();
    
    //In every case, we follow all paths, but cut the depth of any path to max_distance_of_path.
    //We also ensure a minimum depth equal to min_distance_of_path.
    std::list<std::string> queue;
    queue.push_back(starting_anchor_kmer);
    
    //Initialize paths to store all the paths that are found.
    std::vector<std::string> paths;
    //We use i as a counter for how many nodes have been visited in the search.
    //If we haven't finished the search within max_nodes_to_search nodes, 
    //we break the search.
    //This drastically speeds up the run time for some regions.
    //Thankfully, it doesn't seem to impact effectiveness, since most searches
    //complete before this threshold.
    int i = 0;
    //This flag keeps track of whether we had to stop the search early.
    queue_broken = false;
    while(!queue.empty())
    {
        //std::cout << "Queue not empty\n";
        //std::cout.flush();    
        i++;
        std::string current_path = queue.front();
        std::string current_kmer;
        if (going_right)
        {
            
            current_kmer = current_path.substr(current_path.length()-k);
            
        }
        else
        {
            
            current_kmer = current_path.substr(0, k);
        }
        queue.pop_front();
        int current_depth = current_path.length()-k;
        //If we have to terminate search early
        if (i > max_nodes_to_search)
        {
            
            queue_broken = true;
            break;
        }
        //cout << "Checking the node depth\n";
        //If the depth of this node hasn't exceeded the max distance of the path
        if (current_depth <= max_distance_of_path)
        {
            //std::cout << "About to get adjacent\n";
            //std::cout.flush();
            //Extend the path by one nucleotide, keep the ones that are not error kmers
            std::vector<std::string> adjacent_kmers = get_adjacent(k, kmer_hood, shared_reads_l, file, current_kmer, error_threshold, het_threshold, unique_threshold, going_right);
            //std::cout << "Extending the path by one nuc\n";
            //std::cout.flush();
            for (auto adjacent_kmer : adjacent_kmers)
            {
                std::string path;
                if (going_right)
                {
                    path = current_path + adjacent_kmer.back();
                }
                else
                {
                    path = adjacent_kmer.front() + current_path;
                }
                
                bool end_condition;
                //If we are in the middle of the read, we end when we have found a path
                //of nonerror kmers which bridges the anchor kmers
                //and doesn't terminate too early (i.e. before min_distance_of_path)
                if (!left_anchor_kmer.empty() && !right_anchor_kmer.empty())
                {
                    end_condition = ((adjacent_kmer == right_anchor_kmer) && (current_depth + 1 >= min_distance_of_path));
                }
                //If we are at either end of the read, we end when we have found a path
                //of nonerror kmers which continues until the end of the read
                //and doesn't terminate too early (i.e. before min_distance_of_path)
                else
                {
                    end_condition = ((current_depth + 1 == max_distance_of_path) && (current_depth + 1 >= min_distance_of_path));
                }
                if (end_condition)
                {
                    //added this case for if the right and left anchors overlap
                    if (!left_anchor_kmer.empty() && !right_anchor_kmer.empty() && (path.size() < 2*k))
                    {
                        int total_overlaps = 2*k - path.size();
                        path.clear();
                        for (int number_overlaps = 0; number_overlaps < total_overlaps; number_overlaps++)
                        {
                            path += "-";
                        }
                        paths.push_back(path);
                    }
                    else
                    {
                        if (!left_anchor_kmer.empty())
                        {
                            //remove left_anchor_kmer from path
                            path.erase(path.begin(), path.begin()+k);
                        }
                        if (!right_anchor_kmer.empty())
                        {
                            //remove right_anchor_kmer from path
                            path.erase(path.end()-k, path.end());
                        }
                        paths.push_back(path);
                    }
                }
                // "No path yet\n";
                //Else we haven't found a path yet
                else
                {
                    // "No path yet\n";
                    queue.push_back(path);
                }
            }
        }
    }
    return paths;
}

int minDis(std::string s1, std::string s2, int n, int m, std::vector<std::vector<int>> &dp)
{
    // If any string is empty,
    // return the remaining characters of other string
    if (n==0)
    {
        return m;
    }
    if (m==0)
    {
        return n;
    }
    // To check if the recursive tree
    // for given n & m has already been executed    
    if (dp[n][m]!=-1)
    {
        return dp[n][m];
    }
    // If characters are equal, execute
    // recursive function for n-1, m-1
    
    if (s1[n-1] == s2[m-1])
    {
        if (dp[n-1][m-1] == -1)
        {
            return dp[n][m] = minDis(s1, s2, n-1, m-1, dp);
        }
        else
        {
            return dp[n][m] = dp[n-1][m-1];
        }
    }
    // If characters are nt equal, we need to
    // find the minimum cost out of all 3 operations.
    else
    {
        int m1, m2, m3; // temp variables
        
        if (dp[n-1][m]!=-1)
        {
            m1 = dp[n-1][m];
        }
        else
        {
            m1 = minDis(s1, s2, n-1, m, dp);
        }
        
        if (dp[n][m-1]!=-1)
        {
            m2 = dp[n][m-1];
        }
        else
        {
            m2 = minDis(s1, s2, n, m-1, dp);
        }
        
        if (dp[n-1][m-1]!=-1)
        {
            m3 = dp[n-1][m-1];
        }
        else
        {
            m3 = minDis(s1, s2, n-1, m-1, dp);
        }
        return dp[n][m] = 1 + std::min(m1, std::min(m2, m3));
    }
}

void write_error_paths(int& k, counter_kmer &kmer_hood, int shared_reads_l, CKMCFile& file, bool& queue_broken, std::vector<std::string>& edited_error_portions, std::ofstream& errpaths_queuecomplete0_numpaths0_output_file, std::ofstream& errpaths_queuecomplete0_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& errpaths_queuecomplete1_numpaths0_output_file, std::ofstream& errpaths_queuecomplete1_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& erredits_output_file, std::string& edited_read, int& read_number, int& first_error_idx, int& last_error_idx, std::string& before_first_error_kmer, std::string& original_error_portion, std::string& after_last_error_kmer)
    
{
    //std::cout << "Writing error paths\n";
    //std::cout.flush();
    //#ifdef DEBUG_OUT
    std::string left_part;
    if (!before_first_error_kmer.empty())
    {
        left_part = before_first_error_kmer.substr(1);
    }
    else
    {
        left_part = before_first_error_kmer;
    }
    std::string right_part;
    if (!after_last_error_kmer.empty())
    {
        right_part = after_last_error_kmer.substr(0, after_last_error_kmer.length()-1);
    }
    else
    {
        right_part = after_last_error_kmer;
    }
    std::size_t found = original_error_portion.find("-");
    std::string original;
    if (found!=std::string::npos)
    {
        original = left_part + right_part.substr(original_error_portion.size());
    }
    else
    {
        original = left_part + original_error_portion + right_part;
    }
    auto is_less_than = [&](std::string edited_error_portion1, std::string edited_error_portion2)
    {
        int n = original.length();
        std::size_t found1 = edited_error_portion1.find("-");
        std::string portion1;
        if (found1!=std::string::npos)
        {
            portion1 = left_part + right_part.substr(edited_error_portion1.size());
        }
        else
        {
            portion1 = left_part + edited_error_portion1 + right_part;
        }
        int m = portion1.length();
        std::vector<std::vector<int>> dp(n+1, std::vector<int>(m+1, -1));
        int dist1 = minDis(original, portion1, n, m, dp);
        std::size_t found2 = edited_error_portion2.find("-");
        std::string portion2;
        if (found2!=std::string::npos)
        {
            portion2 = left_part + right_part.substr(edited_error_portion2.size());
        }
        else
        {
            portion2 = left_part + edited_error_portion2 + right_part;
        }
        m = portion2.length();
        std::vector<std::vector<int>> dp2(n+1, std::vector<int>(m+1, -1));
        int dist2 = minDis(original, portion2, n, m, dp2);
        return dist1 < dist2;
    };
    
    pdqsort(edited_error_portions.begin(), edited_error_portions.end(), is_less_than);
    
    std::ofstream* errwrite_output_file;
    bool was_edited = false;
    //we finished the search and presumably we have found one homozygous path or two heterozygous paths
    //if ((!queue_broken) && ((edited_error_portions.size() == 1) or (edited_error_portions.size() == 2)))
    if (edited_error_portions.size() >= 1)
    {
        //errwrite_output_file = &erredits_output_file;
        //std::cout << "Finished the search\n";
            //std::cout.flush();
        found = edited_error_portions[0].find("-");
        if (found!=std::string::npos)
        {
            edited_read += left_part.substr(0, left_part.size() - edited_error_portions[0].size());
        }
        else
        {
            edited_read += left_part + edited_error_portions[0];
        }
        if (edited_error_portions[0] != original_error_portion)
        {
            was_edited = true;
        }
        //if (!before_first_error_kmer.empty())
        //{
        //    edited_read += before_first_error_kmer.substr(1) + edited_error_portions[0];
        //}
        //else
        //{
        //    edited_read += edited_error_portions[0];
        //}
    }
    //there are no paths, or there are more than two paths, or the search wasn't finished.
    //We are currently not editing.
    else
    {
        //errwrite_output_file = &errpaths_output_file;
        found = original_error_portion.find("-");
        if (found!=std::string::npos)
        {
            edited_read += left_part.substr(0, left_part.size() - original_error_portion.size());
        }
        else
        {
            edited_read += left_part + original_error_portion;
        }
        //if (!before_first_error_kmer.empty())
        //{
        //    edited_read += before_first_error_kmer.substr(1) + original_error_portion;
        //}
        //else
        //{
        //    edited_read += original_error_portion;
        //}
    }
    if ((!queue_broken) && ((edited_error_portions.size() == 1) or (edited_error_portions.size() == 2)))
    {
        errwrite_output_file = &errpaths_queuecomplete1_numpaths1to2_output_file;
    }
    if ((!queue_broken) && (edited_error_portions.size() < 1))
    {
        errwrite_output_file = &errpaths_queuecomplete1_numpaths0_output_file;
    }
    if ((!queue_broken) && (edited_error_portions.size() > 2))
    {
        errwrite_output_file = &errpaths_queuecomplete1_numpaths3plus_output_file;
    }
    if ((queue_broken) && ((edited_error_portions.size() == 1) or (edited_error_portions.size() == 2)))
    {
        errwrite_output_file = &errpaths_queuecomplete0_numpaths1to2_output_file;
    }
    if ((queue_broken) && (edited_error_portions.size() < 1))
    {
        errwrite_output_file = &errpaths_queuecomplete0_numpaths0_output_file;
    }
    if ((queue_broken) && (edited_error_portions.size() > 2))
    {
        errwrite_output_file = &errpaths_queuecomplete0_numpaths3plus_output_file;
    }
    //outputfileMutex.lock();
    std::string original_error_block;
    found = original_error_portion.find("-");
    if (found!=std::string::npos)
    {
        original_error_block = before_first_error_kmer + after_last_error_kmer.substr(original_error_portion.size());
    }
    else
    {
        original_error_block = before_first_error_kmer + original_error_portion + after_last_error_kmer;
    }
    *errwrite_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
    *errwrite_output_file << before_first_error_kmer << " " << original_error_portion << " " << after_last_error_kmer << '\n';
    if (was_edited)
    {
        erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_original" << '\n';
        erredits_output_file << before_first_error_kmer << " " << original_error_portion << " " << after_last_error_kmer << '\n';
    }
    std::vector<uint32_t> w;
    //file.GetCountersForRead(original_error_block, w);
    if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
        //file.GetCountersForRead(original_error_block, w);
        get_local_counts(kmer_hood, original_error_block, w, k);
    } else {
        file.GetCountersForRead(original_error_block, w);
        }

    for (int j=0; j < w.size(); j++)
    {
        *errwrite_output_file << w.at(j) << " ";
        if (was_edited)
        {
            erredits_output_file << w.at(j) << " ";
        }
    }
    *errwrite_output_file << '\n';
    if (was_edited)
    {
        erredits_output_file << '\n';
    }
    for (int l = 0; l < edited_error_portions.size(); l++)
    {
        std::string edited_error_portion = edited_error_portions[l];
        std::string edited_error_block;
        found = edited_error_portion.find("-");
        if (found!=std::string::npos)
        {
            edited_error_block = before_first_error_kmer + after_last_error_kmer.substr(edited_error_portion.size());
        }
        else
        {
            edited_error_block = before_first_error_kmer + edited_error_portion + after_last_error_kmer;
        }
        *errwrite_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_path" << l << '\n';
        *errwrite_output_file << before_first_error_kmer << " " << edited_error_portion << " " << after_last_error_kmer << '\n';
        if ((l==0) && (was_edited))
        {
            erredits_output_file << ">read" << read_number << "_firsterrorkmer" << first_error_idx << "_lasterrorkmer" << last_error_idx << "_edited" << '\n';
            erredits_output_file << before_first_error_kmer << " " << edited_error_portion << " " << after_last_error_kmer << '\n';
        }
        //file.GetCountersForRead(edited_error_block, w);
        w.clear();
        //w.shrink_to_fit();//clear();

        if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
            //file.GetCountersForRead(edited_error_block, w);
                        //std::cout << w[0] << "\n";
                        //w.clear();
            get_local_counts(kmer_hood, edited_error_block, w, k);
            //std::cout << w[0] << "\n";
            //std::cout << "**********************************\n";
            //std::cout.flush();
        } else {
                    file.GetCountersForRead(edited_error_block, w);
            }


        for (int j=0; j < w.size(); j++)
        {
            *errwrite_output_file << w.at(j) << " ";
            if ((l==0) && (was_edited))
            {
                erredits_output_file << w.at(j) << " ";
            }
        }
        *errwrite_output_file << '\n';
        if ((l==0) && (was_edited))
        {
            erredits_output_file << '\n';
        }
    }
    //outputfileMutex.unlock();
    //#endif
}

void write_nonhom_paths(int& k, counter_kmer &kmer_hood, int shared_reads_l, CKMCFile& file, bool& queue_broken, std::vector<std::string>& smoothed_nonhom_portions, std::ofstream& hetpaths_queuecomplete0_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& hetedits_output_file, std::string& smoothed_read, int& read_number, int& first_nonhom_idx, int& last_nonhom_idx, std::string& before_first_nonhom_kmer, std::string& original_nonhom_portion, std::string& after_last_nonhom_kmer, int& strict)
{    
    //#ifdef DEBUG_OUT 
    std::string left_part;
    if (!before_first_nonhom_kmer.empty())
    {
        left_part = before_first_nonhom_kmer.substr(1);
    }
    else
    {
        left_part = before_first_nonhom_kmer;
    }
    std::string right_part;
    if (!after_last_nonhom_kmer.empty())
    {
        right_part = after_last_nonhom_kmer.substr(0, after_last_nonhom_kmer.length()-1);
    }
    else
    {
        right_part = after_last_nonhom_kmer;
    }
    auto is_greater_than = [&](std::string smoothed_nonhom_portion1, std::string smoothed_nonhom_portion2)
    {
        std::vector<uint32_t> w;
        std::size_t found1 = smoothed_nonhom_portion1.find("-");
        //added this case to account for when anchors overlap
        if (found1!=std::string::npos)
        {
            string left_right = left_part + right_part.substr(smoothed_nonhom_portion1.size());
            if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
                //file.GetCountersForRead(left_right, w);
                get_local_counts(kmer_hood, left_right, w, k);
            } else {
                            file.GetCountersForRead(left_right, w);
                    }
            //file.GetCountersForRead(left_part + right_part.substr(smoothed_nonhom_portion1.size()), w);
        }
        else
        {
            string left_right = left_part + smoothed_nonhom_portion1 + right_part;
            if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
                //file.GetCountersForRead(left_right, w);
                get_local_counts(kmer_hood,  left_right, w, k);
            } else {
                                file.GetCountersForRead(left_right, w);
                        }
            //file.GetCountersForRead(left_part + smoothed_nonhom_portion1 + right_part, w);
        }
        //std::sort(w.begin(), w.end());
        float average1 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
        std::size_t found2 = smoothed_nonhom_portion2.find("-");
        //added this case to account for when anchors overlap
        if (found2!=std::string::npos)
        {
            string left_right = left_part + right_part.substr(smoothed_nonhom_portion2.size());
            if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
                //file.GetCountersForRead(left_right, w);
                get_local_counts(kmer_hood, left_right, w, k);
            } else {
                                file.GetCountersForRead(left_right, w);
                        }
            //file.GetCountersForRead(left_part + right_part.substr(smoothed_nonhom_portion2.size()), w);
        }
        else
        {
            string left_right = left_part + smoothed_nonhom_portion2 + right_part;
            if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
                //file.GetCountersForRead(left_right, w);
                get_local_counts(kmer_hood, left_right, w, k);
            } else {
                                file.GetCountersForRead(left_right, w);
                        }
            //file.GetCountersForRead(left_part + smoothed_nonhom_portion2 + right_part, w);
        }
        //std::sort(w.begin(), w.end());
        float average2 = std::accumulate(w.begin(), w.end(), 0.0)/w.size();
        return average1 > average2;
    };
    pdqsort(smoothed_nonhom_portions.begin(), smoothed_nonhom_portions.end(), is_greater_than);
    std::ofstream* hetwrite_output_file;
    //we finished the search and presumably we have found two heterozygous paths
    //actually, let's relax the !queue_broken requirement, jk restoring this requirement
    //and when strict==1 we only smooth when there are exactly two paths
    //and when strict==0 we smooth when there are at least two paths
    //if ((!queue_broken) && (smoothed_nonhom_portions.size() == 2))
    bool end_condition;
    bool was_smoothed = false;
    std::size_t found;
    if (strict==1)
    {
        end_condition = ((!queue_broken) && (smoothed_nonhom_portions.size() == 2));
    }
    else
    {
        end_condition = ((!queue_broken) && (smoothed_nonhom_portions.size() >= 2));
    }
    if (end_condition)
    {
        found = smoothed_nonhom_portions[0].find("-");
        if (found!=std::string::npos)
        {
            smoothed_read += left_part.substr(0, left_part.size() - smoothed_nonhom_portions[0].size());
        }
        else
        {
            smoothed_read += left_part + smoothed_nonhom_portions[0];
        }
        //We have checked that the condition for smoothing was met
        //Now we check whether the path chosen actually differs from the original
        //If so, then there was a true smooth and we add this to hetedits_output_file
        if (smoothed_nonhom_portions[0] != original_nonhom_portion)
        {
            was_smoothed = true;
        }
    }
    //there is not exactly two paths, or the search wasn't finished.
    //We are currently not smoothing.
    else
    {
        found = original_nonhom_portion.find("-");
        if (found!=std::string::npos)
        {
            smoothed_read += left_part.substr(0, left_part.size() - original_nonhom_portion.size());
        }
        else
        {
            smoothed_read += left_part + original_nonhom_portion;
        }
        
    }
    if ((!queue_broken) && (smoothed_nonhom_portions.size() == 2))
    {
        hetwrite_output_file = &hetpaths_queuecomplete1_numpaths2_output_file;
    }
    if ((!queue_broken) && (smoothed_nonhom_portions.size() < 2))
    {
        hetwrite_output_file = &hetpaths_queuecomplete1_numpaths0to1_output_file;
    }
    if ((!queue_broken) && (smoothed_nonhom_portions.size() > 2))
    {
        hetwrite_output_file = &hetpaths_queuecomplete1_numpaths3plus_output_file;
    }
    if ((queue_broken) && (smoothed_nonhom_portions.size() == 2))
    {
        hetwrite_output_file = &hetpaths_queuecomplete0_numpaths2_output_file;
    }
    if ((queue_broken) && (smoothed_nonhom_portions.size() < 2))
    {
        hetwrite_output_file = &hetpaths_queuecomplete0_numpaths0to1_output_file;
    }
    if ((queue_broken) && (smoothed_nonhom_portions.size() > 2))
    {
        hetwrite_output_file = &hetpaths_queuecomplete0_numpaths3plus_output_file;
    }
    //outputfileMutex.lock();
    std::string original_nonhom_block;
    found = original_nonhom_portion.find("-");
    if (found!=std::string::npos)
    {
        original_nonhom_block = before_first_nonhom_kmer + after_last_nonhom_kmer.substr(original_nonhom_portion.size());
    }
    else
    {
        original_nonhom_block = before_first_nonhom_kmer + original_nonhom_portion + after_last_nonhom_kmer;
    }
    *hetwrite_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_original" << '\n';
    *hetwrite_output_file << before_first_nonhom_kmer << " " << original_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
    if (was_smoothed)
    {
        hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_original" << '\n';
        hetedits_output_file << before_first_nonhom_kmer << " " << original_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
    }
    std::vector<uint32_t> w;
    if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
        //file.GetCountersForRead(original_nonhom_block, w);
        get_local_counts(kmer_hood, original_nonhom_block, w, k);
    } else {
        file.GetCountersForRead(original_nonhom_block, w);
        }
    //file.GetCountersForRead(original_nonhom_block, w);
    for (int j=0; j < w.size(); j++)
    {
        *hetwrite_output_file << w.at(j) << " ";
        if (was_smoothed)
        {
            hetedits_output_file << w.at(j) << " ";
        }
    }
    *hetwrite_output_file << '\n';
    if (was_smoothed)
    {
        hetedits_output_file << '\n';
    }
    for (int l = 0; l < smoothed_nonhom_portions.size(); l++)
    {
        std::string smoothed_nonhom_portion = smoothed_nonhom_portions[l];
        std::string smoothed_nonhom_block;
        found = smoothed_nonhom_portion.find("-");
        if (found!=std::string::npos)
        {
            smoothed_nonhom_block = before_first_nonhom_kmer + after_last_nonhom_kmer.substr(smoothed_nonhom_portion.size());
        }
        else
        {
            smoothed_nonhom_block = before_first_nonhom_kmer + smoothed_nonhom_portion + after_last_nonhom_kmer;
        }
        *hetwrite_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_path" << l << '\n';
        *hetwrite_output_file << before_first_nonhom_kmer << " " << smoothed_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
        if ((l==0) && (was_smoothed))
        {
            hetedits_output_file << ">read" << read_number << "_firstnonhomkmer" << first_nonhom_idx << "_lastnonhomkmer" << last_nonhom_idx << "_smoothed" << '\n';
            hetedits_output_file << before_first_nonhom_kmer << " " << smoothed_nonhom_portion << " " << after_last_nonhom_kmer << '\n';
        }
        w.clear();
        if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
            //file.GetCountersForRead(smoothed_nonhom_block, w);
            get_local_counts(kmer_hood, smoothed_nonhom_block, w, k);
        } else {
                    file.GetCountersForRead(smoothed_nonhom_block, w);
            }
        //file.GetCountersForRead(smoothed_nonhom_block, w);
        for (int j=0; j < w.size(); j++)
        {
            *hetwrite_output_file << w.at(j) << " ";
            if ((l==0) && (was_smoothed))
            {
                hetedits_output_file << w.at(j) << " ";
            }
        }
        *hetwrite_output_file << '\n';
        if ((l==0) && (was_smoothed))
        {
            hetedits_output_file << '\n';
        }
    }
    //outputfileMutex.unlock();
    //#endif
}

std::string remove_err (int& k, counter_kmer &kmer_hood, int shared_reads_l, CKMCFile& file, std::vector<uint32_t>& v, std::string& read, int& read_number, std::ofstream& errpaths_queuecomplete0_numpaths0_output_file, std::ofstream& errpaths_queuecomplete0_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& errpaths_queuecomplete1_numpaths0_output_file, std::ofstream& errpaths_queuecomplete1_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& erredits_output_file, int& error_threshold, int& het_threshold, int& unique_threshold, int& max_nodes_to_search, double& distance_multiplier)
{
    
    //initialize variables
    std::string edited_read;
    //int k = 21;
#ifdef DEBUG_OUT        
        auto start_time = chrono::high_resolution_clock::now();
#endif
    //iterate over counts to edit errors
    int previous_type = -1;
    int first_nonerror_idx;
    int last_nonerror_idx;
    int first_error_idx;
    int last_error_idx;
    std::string before_first_error_kmer;
    std::string after_last_error_kmer;
    
    for (int i = 0; i < v.size(); i++)
    {


#ifdef DEBUG_OUT        
        auto start_time_loop = chrono::high_resolution_clock::now();
#endif

        int current_type = get_type(v[i], error_threshold, het_threshold, unique_threshold);
        //std::cout << "Current type: " << current_type << "\n";
        //std::cout.flush();
#ifdef DEBUG_OUT
        auto end_time_loop = chrono::high_resolution_clock::now();
        auto total_time = chrono::duration_cast<chrono::nanoseconds>(end_time_loop - start_time_loop);
        section10 = section10 + total_time;
        start_time_loop = chrono::high_resolution_clock::now();
#endif        
        //if kmer is error
        if (current_type == 0)
        {    
            
            //if this is the first kmer of the read
            if (previous_type == -1)
            {
                first_error_idx = i;
                last_nonerror_idx = i-1;
            }
            //if previous kmer was an error, we are continuing the error block
            
            if (previous_type == 0)
            {
                ;
            }
            //if previous kmer was not an error, we are leaving the nonerror block
            if (previous_type > 0)
            {
                //get kmer that is right before the first error kmer of error block
                first_error_idx = i;
                last_nonerror_idx = i-1;
                before_first_error_kmer = read.substr(i-1, k);
                
                std::string nonerror_portion = read.substr(first_nonerror_idx, last_nonerror_idx - first_nonerror_idx + 1);
                edited_read += nonerror_portion;
            }
            previous_type = current_type;
        }
        //if kmer is nonerror
        
        if (current_type > 0)
        {
            
            //if this is the first kmer of the read
            if (previous_type == -1)
            {
                first_nonerror_idx = i;
                last_error_idx = i-1;
            }
            
            //if previous kmer was error, and we are at the beginning of the read
            if (previous_type == 0 && before_first_error_kmer.empty())
            {
                
                //The very beginning of the read is an error portion
                after_last_error_kmer = read.substr(i, k);
                int min_distance_of_path = 0;
                int max_distance_of_path = i;
                bool queue_broken = false;
                std::vector<std::string> edited_error_portions = get_paths(k, kmer_hood, shared_reads_l, file, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, queue_broken);
                last_error_idx = i-1;
                first_nonerror_idx = i;
                std::string original_error_portion = read.substr(0, i);
                
                write_error_paths(k, kmer_hood, shared_reads_l, file, queue_broken, edited_error_portions, errpaths_queuecomplete0_numpaths0_output_file, errpaths_queuecomplete0_numpaths1to2_output_file, errpaths_queuecomplete0_numpaths3plus_output_file, errpaths_queuecomplete1_numpaths0_output_file, errpaths_queuecomplete1_numpaths1to2_output_file, errpaths_queuecomplete1_numpaths3plus_output_file, erredits_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer);
            }
            
            //if previous kmer was error, we have left the error block
            
            if (previous_type == 0 && !before_first_error_kmer.empty())
            {
                
                int number_of_error_kmers = i - first_error_idx;
                //If the position of after_last_error_kmer overlaps before_first_error_kmer
                //we keep progressing as if nothing has happened, waiting to find another non_error kmer
                //if (number_of_error_kmers < k)
                //{
                //    current_type = previous_type;
                //    continue;
                //}
                //get kmer that is right after the last error kmer of block
                after_last_error_kmer = read.substr(i, k);
                //int min_distance_of_path = k;
                int min_distance_of_path = std::min(number_of_error_kmers, k);
                int max_distance_of_path = ceil(distance_multiplier * number_of_error_kmers);
                bool queue_broken;
                
                std::vector<std::string> edited_error_portions = get_paths(k, kmer_hood, shared_reads_l, file, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, queue_broken);
                
                last_error_idx = i-1;
                first_nonerror_idx = i;
                std::string original_error_portion;
                //If the after_last_error_kmer overlaps before_first_error_kmer
                if (last_error_idx - first_error_idx + 2 - k < 0)
                {
                    
                    for (int number_overlaps=0; number_overlaps < first_error_idx - last_error_idx + k - 2; number_overlaps++)
                    {
                        original_error_portion += "-";
                    }
                }
                else
                {
                    original_error_portion = read.substr(first_error_idx+k-1, last_error_idx - first_error_idx + 2 - k);
                }
                
                write_error_paths(k, kmer_hood,shared_reads_l, file, queue_broken, edited_error_portions, errpaths_queuecomplete0_numpaths0_output_file, errpaths_queuecomplete0_numpaths1to2_output_file, errpaths_queuecomplete0_numpaths3plus_output_file, errpaths_queuecomplete1_numpaths0_output_file, errpaths_queuecomplete1_numpaths1to2_output_file, errpaths_queuecomplete1_numpaths3plus_output_file, erredits_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer);
            }
            //if previous kmer is nonerror, we are continuing a non error block
            if (previous_type > 0)
            {
                ;
                
            }
            previous_type = current_type;
        }
    }
#ifdef DEBUG_OUT
        auto end_time_loop = chrono::high_resolution_clock::now();
        auto total_time = chrono::duration_cast<chrono::nanoseconds>(end_time_loop - start_time);
        section11 = section11 + total_time;
        auto start_time_loop = chrono::high_resolution_clock::now();
#endif
    //We have reached the end of the read, let's make sure we have added the last bit of the read
    if (previous_type == 0)
    {
        //std::cout << previous_type << "\n";
                //std::cout.flush();
        //We have "left" the error portion of the read
        after_last_error_kmer = "";
        int min_distance_of_path = 0;
        int max_distance_of_path = v.size()-first_error_idx;
        bool queue_broken = false;

        std::vector<std::string> edited_error_portions = get_paths(k, kmer_hood, shared_reads_l, file, error_threshold, het_threshold, unique_threshold, before_first_error_kmer, after_last_error_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, queue_broken);
        
        last_error_idx = v.size()-1;
        first_nonerror_idx = v.size();
        
        std::string original_error_portion = read.substr(first_error_idx+k-1);
        write_error_paths(k, kmer_hood, shared_reads_l, file, queue_broken, edited_error_portions, errpaths_queuecomplete0_numpaths0_output_file, errpaths_queuecomplete0_numpaths1to2_output_file, errpaths_queuecomplete0_numpaths3plus_output_file, errpaths_queuecomplete1_numpaths0_output_file, errpaths_queuecomplete1_numpaths1to2_output_file, errpaths_queuecomplete1_numpaths3plus_output_file, erredits_output_file, edited_read, read_number, first_error_idx, last_error_idx, before_first_error_kmer, original_error_portion, after_last_error_kmer);

        //std::cout << previous_type << "\n";
                //std::cout.flush();
    }
    
    if (previous_type > 0)
    {
        //std::cout << previous_type << "\n";
        //std::cout.flush();
        //We have "left" the nonerror portion of the read
        first_error_idx = v.size();
        last_nonerror_idx = v.size()-1;
        std::string nonerror_portion = read.substr(first_nonerror_idx, last_nonerror_idx - first_nonerror_idx + k);
        edited_read += nonerror_portion;
        
    }
#ifdef DEBUG_OUT
        end_time_loop = chrono::high_resolution_clock::now();
        total_time = chrono::duration_cast<chrono::nanoseconds>(end_time_loop - start_time_loop);
        section12 = section12 + total_time;
        start_time_loop = chrono::high_resolution_clock::now();
#endif
    edited_read += '\n';
    return edited_read;
}

std::string smooth_het (int& k, counter_kmer &kmer_hood, int shared_reads_l, CKMCFile& file, std::vector<uint32_t>& v, std::string& read, int& read_number, std::ofstream& hetpaths_queuecomplete0_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths3plus_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& hetedits_output_file, int& error_threshold, int& het_threshold, int& unique_threshold, int& anchor_threshold, int& max_nodes_to_search, double& distance_multiplier, int& strict)
{
#ifdef DEBUG_OUT    
        auto start_time = chrono::high_resolution_clock::now();
#endif
    //initialize variables
    std::string smoothed_read;
    //int k = 21;

    //iterate over counts to smoothe het
    int previous_type = -1;
    int first_hom_idx;
    int last_hom_idx;
    int first_nonhom_idx;
    int last_nonhom_idx;
    std::string before_first_nonhom_kmer;
    std::string after_last_nonhom_kmer;
    for (int i = 0; i < v.size(); i++)
    {
                auto start_time_loop = chrono::high_resolution_clock::now();
        std::string previous_kmer;
        int previous_count;
        if (i == 0)
        {
            previous_kmer = "";
            previous_count = 0;
        }
        else
        {
            previous_kmer = read.substr(i-1, k);
            std::vector<uint32_t> previous_count_vector;
            //file.GetCountersForRead(previous_kmer, previous_count_vector);
            if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
		get_local_counts(kmer_hood, previous_kmer, previous_count_vector, k);
            } else {
                file.GetCountersForRead(previous_kmer, previous_count_vector);
            }
            previous_count = previous_count_vector[0];
        }
#ifdef DEBUG_OUT        
        auto end_time_loop = chrono::high_resolution_clock::now();
        auto total_time = chrono::duration_cast<chrono::nanoseconds>(end_time_loop - start_time_loop);
                
        section15 = section15 + total_time;
        start_time_loop = chrono::high_resolution_clock::now();
#endif
        std::string current_kmer = read.substr(i, k);
        std::vector<uint32_t> current_count_vector;
        if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
            //file.GetCountersForRead(current_kmer, current_count_vector);
            get_local_counts(kmer_hood, current_kmer, current_count_vector , k);
        } else {
            file.GetCountersForRead(current_kmer, current_count_vector);
        }
        //file.GetCountersForRead(current_kmer, current_count_vector);
        int current_count = current_count_vector[0];
	//TODO        
	//std::cout << current_count << "\t" ; 
	//file.GetCountersForRead(current_kmer, current_count_vector);
	//if (current_count != current_count_vector[0]) {
	//	std::cout << current_kmer << "\t" << current_kmer.size() << "\t" << current_count << "\t" << current_count_vector[0] << "\n" ;
	//	std::cout << "\n";
	//}

	std::string anchor_found;
        int current_type = get_type_het(k, kmer_hood,  shared_reads_l, file, previous_type, previous_kmer, current_kmer, previous_count, current_count, error_threshold, het_threshold, unique_threshold, anchor_threshold, anchor_found);
#ifdef DEBUG_OUT        
        end_time_loop = chrono::high_resolution_clock::now();
        total_time = chrono::duration_cast<chrono::nanoseconds>(end_time_loop - start_time_loop);
        section16 = section16 + total_time;
        start_time_loop = chrono::high_resolution_clock::now();
#endif
        //if kmer is nonhom
        if ((current_type == 0) || (current_type == 1) || (current_type == 3))
        {
            //if this is the first kmer of the read
            if (previous_type == -1)
            {
                first_nonhom_idx = i;
                last_hom_idx = i-1;
            }
            //if previous kmer was nonhom, we are continuing the nonhom block
            if ((previous_type == 0) || (previous_type == 1) || (previous_type == 3))
            {
                ;
            }
            //if previous kmer was hom, we are leaving the hom block
            if (previous_type == 2)
            {
                //get kmer that is right before the first nonhom kmer of nonhom block
                first_nonhom_idx = i;
                last_hom_idx = i-1;
                before_first_nonhom_kmer = read.substr(i-1, k);
                std::string hom_portion = read.substr(first_hom_idx, last_hom_idx - first_hom_idx + 1);
                smoothed_read += hom_portion;
                
            }
            previous_type = current_type;
        }
#ifdef DEBUG_OUT        
        end_time_loop = chrono::high_resolution_clock::now();
        total_time = chrono::duration_cast<chrono::nanoseconds>(end_time_loop - start_time_loop);
        section17 = section17 + total_time;
        start_time_loop = chrono::high_resolution_clock::now();
#endif
        //if kmer is hom
        if (current_type == 2)
        {
            
                    start_time_loop = chrono::high_resolution_clock::now();
            //if this is the first kmer of the read
            if (previous_type == -1)
            {
                first_hom_idx = i;
                last_nonhom_idx = i-1;
            }
            //if previous kmer was nonhom, and we are at the beginning of the read
            if (((previous_type == 0) || (previous_type == 1) || (previous_type == 3)) && before_first_nonhom_kmer.empty())
            {
                //The very beginning of the read is an nonhom portion
                after_last_nonhom_kmer = read.substr(i, k);
                int min_distance_of_path = 0;
                int max_distance_of_path = i;
                bool queue_broken = false;
                std::vector<std::string> smoothed_nonhom_portions = get_paths(k, kmer_hood, shared_reads_l, file, error_threshold, het_threshold, unique_threshold, before_first_nonhom_kmer, after_last_nonhom_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, queue_broken);
                last_nonhom_idx = i-1;
                first_hom_idx = i;
                std::string original_nonhom_portion = read.substr(0, i);
                
                write_nonhom_paths(k, kmer_hood, shared_reads_l, file, queue_broken, smoothed_nonhom_portions, hetpaths_queuecomplete0_numpaths0to1_output_file, hetpaths_queuecomplete0_numpaths2_output_file, hetpaths_queuecomplete0_numpaths3plus_output_file, hetpaths_queuecomplete1_numpaths0to1_output_file, hetpaths_queuecomplete1_numpaths2_output_file, hetpaths_queuecomplete1_numpaths3plus_output_file, hetedits_output_file, smoothed_read, read_number, first_nonhom_idx, last_nonhom_idx, before_first_nonhom_kmer, original_nonhom_portion, after_last_nonhom_kmer, strict);
            
            }
#ifdef DEBUG_OUT            
            end_time_loop = chrono::high_resolution_clock::now();
            total_time = chrono::duration_cast<chrono::nanoseconds>(end_time_loop - start_time_loop);
            section18 = section18 + total_time;
            start_time_loop = chrono::high_resolution_clock::now();
#endif    
            //if previous kmer was nonhom, we have left the nonhom block
            if (((previous_type == 0) || (previous_type == 1) || (previous_type == 3)) && !before_first_nonhom_kmer.empty())
            {
                if (anchor_found == "left")
                {
                    after_last_nonhom_kmer = read.substr(i, k);
                    std::string hom_portion = read.substr(first_hom_idx+1, i-first_hom_idx-1);
                    
                    last_nonhom_idx = i-1;
                    first_hom_idx = i;
                    smoothed_read += hom_portion;
                }
                else
                {
                int number_of_nonhom_kmers = i - first_nonhom_idx;
                //If the position of after_last_nonhom_kmer overlaps before_first_nonhom_kmer
                //we keep progressing as if nothing has happened, waiting to find another hom kmer
                //if ((anchor_found == "right") && (number_of_nonhom_kmers < k))
                //{
                //    
                //    current_type = previous_type;
                //    continue;
                //}
                //get kmer that is right after the last nonhom kmer of block
                after_last_nonhom_kmer = read.substr(i, k);
                //int min_distance_of_path = k;
                int min_distance_of_path = std::min(number_of_nonhom_kmers, k);
                int max_distance_of_path = ceil(distance_multiplier * number_of_nonhom_kmers);
                bool queue_broken;
                std::vector<std::string> smoothed_nonhom_portions = get_paths(k, kmer_hood, shared_reads_l, file, error_threshold, het_threshold, unique_threshold, before_first_nonhom_kmer, after_last_nonhom_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, queue_broken);
                last_nonhom_idx = i-1;
                first_hom_idx = i;
                std::string original_nonhom_portion;
                //If the after_last_nonhom_kmer overlaps before_first_nonhom_kmer
                if (last_nonhom_idx - first_nonhom_idx + 2 - k < 0)
                {
                    for (int number_overlaps=0; number_overlaps < first_nonhom_idx - last_nonhom_idx + k - 2; number_overlaps++)
                    {
                        original_nonhom_portion += "-";
                    }
                }
                else
                {
                    original_nonhom_portion = read.substr(first_nonhom_idx+k-1, last_nonhom_idx - first_nonhom_idx + 2 - k);
                }
                write_nonhom_paths(k, kmer_hood, shared_reads_l, file, queue_broken, smoothed_nonhom_portions, hetpaths_queuecomplete0_numpaths0to1_output_file, hetpaths_queuecomplete0_numpaths2_output_file, hetpaths_queuecomplete0_numpaths3plus_output_file, hetpaths_queuecomplete1_numpaths0to1_output_file, hetpaths_queuecomplete1_numpaths2_output_file, hetpaths_queuecomplete1_numpaths3plus_output_file, hetedits_output_file, smoothed_read, read_number, first_nonhom_idx, last_nonhom_idx, before_first_nonhom_kmer, original_nonhom_portion, after_last_nonhom_kmer, strict);
                }
            }
#ifdef DEBUG_OUT            
            end_time_loop = chrono::high_resolution_clock::now();
            total_time = chrono::duration_cast<chrono::nanoseconds>(end_time_loop - start_time_loop);
            section19 = section19 + total_time;
            start_time_loop = chrono::high_resolution_clock::now();
#endif    
            //if previous kmer is hom, we are continuing a hom block
            if (previous_type == 2)
            {
                if (anchor_found == "left")
                {
                    std::string hom_portion = read.substr(first_hom_idx, i-first_hom_idx);
                    smoothed_read += hom_portion;
                    
                    first_hom_idx=i;
                    last_nonhom_idx=i-1;
                }
            }
            previous_type = current_type;
#ifdef DEBUG_OUT            
            end_time_loop = chrono::high_resolution_clock::now();
            total_time = chrono::duration_cast<chrono::nanoseconds>(end_time_loop - start_time_loop);
            section20 = section20 + total_time;
	    start_time_loop = chrono::high_resolution_clock::now();
#endif    
        }
    }
#ifdef DEBUG_OUT
//    auto end_time = chrono::high_resolution_clock::now();
//    auto total_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
//    section21 = section21 + total_time;
    start_time = chrono::high_resolution_clock::now();
#endif
    //We have reached the end of the read, let's make sure we have added the last bit of the read
    if ((previous_type == 0) || (previous_type == 1) || (previous_type == 3))
    {
        //We have "left" the nonhom portion of the read
        //If we have a homozygous on left with which to anchor
        if (first_nonhom_idx > 0)
        {
            after_last_nonhom_kmer = "";
            int min_distance_of_path = 0;
            int max_distance_of_path = v.size()-first_nonhom_idx;
            bool queue_broken = false;
            std::vector<std::string> smoothed_nonhom_portions = get_paths(k, kmer_hood, shared_reads_l, file, error_threshold, het_threshold, unique_threshold, before_first_nonhom_kmer, after_last_nonhom_kmer, min_distance_of_path, max_distance_of_path, max_nodes_to_search, queue_broken);
            last_nonhom_idx = v.size()-1;
            first_hom_idx = v.size();
            std::string original_nonhom_portion = read.substr(first_nonhom_idx+k-1);
            write_nonhom_paths(k, kmer_hood, shared_reads_l, file, queue_broken, smoothed_nonhom_portions, hetpaths_queuecomplete0_numpaths0to1_output_file, hetpaths_queuecomplete0_numpaths2_output_file, hetpaths_queuecomplete0_numpaths3plus_output_file, hetpaths_queuecomplete1_numpaths0to1_output_file, hetpaths_queuecomplete1_numpaths2_output_file, hetpaths_queuecomplete1_numpaths3plus_output_file, hetedits_output_file, smoothed_read, read_number, first_nonhom_idx, last_nonhom_idx, before_first_nonhom_kmer, original_nonhom_portion, after_last_nonhom_kmer, strict);
        }
        else
        {
            smoothed_read += read;
        }
    }

    if (previous_type == 2)
    {
        //We have "left" the hom portion of the read
        first_nonhom_idx = v.size();
        last_hom_idx = v.size()-1;
        std::string hom_portion = read.substr(first_hom_idx, last_hom_idx - first_hom_idx + k);
        smoothed_read += hom_portion;
        
    }

     smoothed_read += '\n';
#ifdef DEBUG_OUT
    auto end_time = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
    section21 = section21 + total_time;
#endif
    //coutMutex.unlock();
              

    return smoothed_read;
}

std::string getFileExt (const std::string &s)
{
    size_t i = s.find_last_of('.');
    if (i != std::string::npos)
    {
        return(s.substr(i+1, s.length() - i));
    }
    return("");
}

void processRead(std::mutex& coutMutex, int& main_k, std::mutex& inputfileMutex,
		int& cov_t, int& min_t, int& m, std::ifstream& input_file,
		int& line_num, int& num_lines_per_read, CKMCFile& file, int& error_threshold, 
		std::ofstream& err_output_file, std::ofstream& errpaths_queuecomplete0_numpaths0_output_file, 
		std::ofstream& errpaths_queuecomplete0_numpaths1to2_output_file, std::ofstream& errpaths_queuecomplete0_numpaths3plus_output_file, 
		std::ofstream& errpaths_queuecomplete1_numpaths0_output_file, std::ofstream& errpaths_queuecomplete1_numpaths1to2_output_file, 
		std::ofstream& errpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& erredits_output_file, 
		int& het_threshold, int& unique_threshold, int& max_nodes_to_search, double& distance_multiplier, int& polish, double& allowed_rep_fraction, 
		std::ofstream& het_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths0to1_output_file, 
		std::ofstream& hetpaths_queuecomplete0_numpaths2_output_file, std::ofstream& hetpaths_queuecomplete0_numpaths3plus_output_file, 
		std::ofstream& hetpaths_queuecomplete1_numpaths0to1_output_file, std::ofstream& hetpaths_queuecomplete1_numpaths2_output_file, 
		std::ofstream& hetpaths_queuecomplete1_numpaths3plus_output_file, std::ofstream& hetedits_output_file, 
		int& anchor_threshold, int& strict, int& verbose)
{

    int k = main_k;
    std::string line;
    int thread_line_num;
    std::string header;
    inputfileMutex.lock();
    //for (uint32_t rdnum=first_seq; rdnum < (last_seq); rdnum++) 
    while (getline(input_file, line))
    {
#ifdef DEBUG_OUT
        auto start_time = chrono::high_resolution_clock::now();
#endif
        line_num++;
	thread_line_num = line_num;
	if (thread_line_num % num_lines_per_read == 0){
		header = ">" + line.substr(1,50);
	}
        //header = headers[rdnum];
        if (thread_line_num % num_lines_per_read == 1){    
            inputfileMutex.unlock();
            std::string read = line;
            //seqs[rdnum] = "";
	    int read_number = (thread_line_num+num_lines_per_read-1)/num_lines_per_read - 1;
            
            //std::time_t thread_result = std::time(nullptr);
            //if (verbose == 1)
            //{
            //    coutMutex.lock();
            //    std::cout << "Analyzing read/contig/scaffold number  " << read_number << " at " << std::asctime(std::localtime(&thread_result)) << '\n';
            //    coutMutex.unlock();
            //}
            //std::cout << "Read number: " << read_number << "\n";
            //std::cout.flush();
            
	    //Find the shared reads TODO
            //std::vector <long long> kmer_arr[1]; 
            robin_hood::unordered_set<uint32_t> shared_reads;
            //robin_hood::unordered_set<uint32_t> shared_reads_all;
	    //Now we want to go through the all_v_all for that read. 
	    //Will need to loop through all of the all v alls, since we didn't want to merge them earlier (space saving) 
	    std::vector <uint32_t> these_reads;
	    for (uint32_t indx=0; indx < all_v_all.size(); indx++) { // Go through each all v all index 
	    	
		if (all_v_all[indx].find(read_number) != all_v_all[indx].end()) { // If this read number is in this particular all_v_all index
			these_reads = all_v_all[indx][read_number];
			all_v_all[indx][read_number].clear();
			all_v_all[indx][read_number].shrink_to_fit();
			break;
	    	}
	    }
	    //these_reads.reserve(10000);
	    //Go through the minimizers for this read
	    //std::vector <long long> these_read_mins = read_mins_all[read_number];
	    //"these_reads" are a list of all of the reads that have overlapped with this particular read 
	    //std::cout << "About to find other reads\n";
	    //std::cout.flush();
	    //for (long long min: read_mins_all[read_number]) {
	    //	these_reads.insert(these_reads.end(), min_db_all[min].begin(), min_db_all[min].end());
	    //}

	    //std::cout << these_reads.size() << "\n";
            //std::cout.flush();

	    //all_v_all[read_number].clear();
	    //all_v_all[read_number].shrink_to_fit();
	    //std::vector <uint32_t> these_reads = all_v_all[read_number];
            int shared_reads_l = 0;
	    //std::cout << min_t << "\n";
	    if (these_reads.size() != 0) {
           	 
	    	pdqsort(these_reads.begin(), these_reads.end()) ;

	    	int overlap_cntr = 0;

            	uint32_t last_read = these_reads.front();
            
	    	for (int y=0; y< these_reads.size(); y++) {
                	if (these_reads[y] == last_read) {
                    		++overlap_cntr;
                	} else {
                    		if (overlap_cntr > min_t ) { // && last_read != read_number) { 
                        		shared_reads.insert(last_read);
                    	}
                    	last_read = these_reads[y];
                    	overlap_cntr = 1;
                	}
            	}
            	if (overlap_cntr > min_t ) { // && last_read != read_number) {
			
                	shared_reads.insert(last_read);
            	} 

            	if (shared_reads.size() > shared_read_t){
			if (shared_reads.size() < 200){
                		shared_reads_l = 1;
			}
			
            	}
		
	    }  
	    //coutMutex.lock();
	    //std::cout << shared_reads.size() << "\n";
	    //std::cout.flush();
	    //coutMutex.unlock();
	    these_reads.clear();
	    //uint32_t t = 0;
	    //all_v_all[read_number] = {t};
	    //all_v_all[read_number].shrink_to_fit();

	    //these_reads.shrink_to_fit();	    
	    std::vector<uint32_t> v;

#ifdef DEBUG_OUT            
            auto end_time = chrono::high_resolution_clock::now();
            auto total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
            section7 = section7 + total_time;
            start_time = chrono::high_resolution_clock::now();
#endif    
	    //std::cout << "Found the overlaps\n";
            //std::cout.flush();

            counter_kmer kmer_hood;
            //std::vector <unsigned long long> cnt_vec;
	    //std::cout << read_number << "\t" << shared_reads.size() << "\n";
	    //std::cout.flush();
	    if (shared_reads_l == 1) {
                //std::vector<unsigned long long> kmer_test_kmer;
                //kmerize_all(this_read, k, kmer_test_kmer);
		//std::cout << shared_reads.size() << "\n";
		//std::cout.flush();
		make_local_kmer_db(read.length(), kmer_hood, seqs, shared_reads, k);
		//make_local_kmer_db(read.length(), kmer_hood, seq_kmers, shared_reads);
		v.reserve(read.length());
		int error_rate = get_local_counts_with_error_rate(kmer_hood, seqs[read_number], v, k); 
		
		//Add a case for re-making the kmer db and re-setting k to something lower if there is too much error and het. 
		//int error_rate = error_rate_finder(v);
		//std::cout << error_rate << "\n";
		//std::cout.flush();
		/*if (error_rate > 100000) {
			std::cout << header << "\n";
			k = 17;
			make_local_kmer_db(read.length(), kmer_hood, seqs, shared_reads, k);
			v.clear();
			get_local_counts(kmer_hood, seqs[read_number], v, k);
		}*/
		//get_local_counts_precomputed_kmers(seq_kmers[read_number], kmer_hood, v); 
		
		//seq_kmers[read_number].clear();
		//seq_kmers[read_number].shrink_to_fit();
		//get_local_counts_precomputed_kmers(kmer_hood, kmer_arr, cnt_vec, tmp_vec, read, v, k);
            } else {
                file.GetCountersForRead(read, v);
            }
#ifdef DEBUG_OUT            
            end_time = chrono::high_resolution_clock::now();
            total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
            section8 = section8 + total_time;
            start_time = chrono::high_resolution_clock::now();
#endif    
            int num_kmers = v.size();
            
            //Calculate the lambda value and re-evaluate the thresholds 
            int error_threshold_local = error_threshold ;
            int het_threshold_local = het_threshold;
            int anchor_threshold_local = anchor_threshold;
            int unique_threshold_local = unique_threshold;
	    
	    //std::cout << read_number << "\t" << shared_reads.size() << "\n" ;
	    //std::cout.flush();
	    if (shared_reads.size() > shared_read_t) {
                //if (read.length() > 2000) {
                    if (shared_reads.size() < 90) {
			float local_lambda = (shared_reads.size()/5.96) + 0.936;
                    	error_threshold_local = round(ceil(0.25 * local_lambda)) ; //round(ceil(0.25 * l))
                    	het_threshold_local = round(ceil(1.5 * local_lambda));
                    	anchor_threshold_local = round(ceil(2.5 * local_lambda));
                    	unique_threshold_local = round(ceil(3.5 * local_lambda));
		    }
                //}
            } 

#ifdef DEBUG_OUT            
            end_time = chrono::high_resolution_clock::now();
            total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
            section9 = section9 + total_time;
            start_time = chrono::high_resolution_clock::now();
#endif    
            //std::cout << "About to edit errors";
            //std::cout.flush();
	    //remove errors from the read to get edited read
            std::string edited_read = remove_err(k, kmer_hood, shared_reads_l, file, v, read, read_number, errpaths_queuecomplete0_numpaths0_output_file, errpaths_queuecomplete0_numpaths1to2_output_file, errpaths_queuecomplete0_numpaths3plus_output_file, errpaths_queuecomplete1_numpaths0_output_file, errpaths_queuecomplete1_numpaths1to2_output_file, errpaths_queuecomplete1_numpaths3plus_output_file, erredits_output_file, error_threshold_local, het_threshold_local, unique_threshold_local, max_nodes_to_search, distance_multiplier);
//#ifdef DEBUG_OUT
//            end_time = chrono::high_resolution_clock::now();
//            total_time = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
//            section6 = section6 + total_time;
//            start_time = chrono::high_resolution_clock::now();
//#endif    
            //write header and edited read to err_output_file
            //outputfileMutex.lock(#ifdef DEBUG_OUT
#ifdef DEBUG_OUT
	    start_time = chrono::high_resolution_clock::now();
#endif            
	    err_output_file << header << '\n';
            err_output_file << edited_read;
            //outputfileMutex.unlock();
            edited_read.pop_back();

#ifdef DEBUG_OUT            
            end_time = chrono::high_resolution_clock::now();
            total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
            section13 = section13 + total_time;
            start_time = chrono::high_resolution_clock::now();
#endif      
      	    if (polish == 1)
            {
                continue;
            }

            //get counters of kmers in edited read
            v.clear();
            if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
                v.reserve(read.length());
                get_local_counts(kmer_hood, edited_read, v, k);
                
            } else {
                file.GetCountersForRead(edited_read, v);
            }

	    num_kmers = v.size();

#ifdef DEBUG_OUT
            end_time = chrono::high_resolution_clock::now();
            total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
            section14 = section14 + total_time;
            start_time = chrono::high_resolution_clock::now();
#endif    
	    //std::cout << "About to edit het\n";
            //std::cout.flush();
            std::string smoothed_read = smooth_het(k, kmer_hood, shared_reads_l, file, v, edited_read, read_number, hetpaths_queuecomplete0_numpaths0to1_output_file, hetpaths_queuecomplete0_numpaths2_output_file, hetpaths_queuecomplete0_numpaths3plus_output_file, hetpaths_queuecomplete1_numpaths0to1_output_file, hetpaths_queuecomplete1_numpaths2_output_file, hetpaths_queuecomplete1_numpaths3plus_output_file, hetedits_output_file, error_threshold_local, het_threshold_local, unique_threshold_local, anchor_threshold_local, max_nodes_to_search, distance_multiplier, strict);

#ifdef DEBUG_OUT
            start_time = chrono::high_resolution_clock::now();
#endif

	    //v.clear();
            //if (shared_reads_l == 1){ //(shared_reads.size() > shared_read_t) {
            //    v.reserve(smoothed_read.length());
            //    get_local_counts(kmer_hood, smoothed_read, v, k);

            //} else {
            //    file.GetCountersForRead(smoothed_read, v);
            //}
	    //coutMutex.lock();
	    //std::cout << header << "\t";
	    //for (auto thisthing: v) {
	//	std::cout << thisthing << "," ;
	    //}
	    //std::cout << "\n";
	    //coutMutex.unlock();
	    //std::cout << "About to close the files\n";
            //std::cout.flush();
            //write header and smoothed read to het_output_file
            //outputfileMutex.lock();
            het_output_file << header << '\n';
            het_output_file << smoothed_read;
            //outputfileMutex.unlock();
            inputfileMutex.lock();
#ifdef DEBUG_OUT
            end_time = chrono::high_resolution_clock::now();
            total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
            section22 = section22 + total_time;
            start_time = chrono::high_resolution_clock::now();
#endif

            
        }
    }
    inputfileMutex.unlock();
}

void read_fasta_slow(std::ifstream& fasta_file_que){
    int seq_cntr; seq_cntr=0;
    int name_cntr; name_cntr=0;
    unsigned long long kmer_f;
    unsigned long long kmer_r;
    //std::string line;
    
    std::ifstream &file(fasta_file_que);
        if (file.is_open()) {
        std::string line;
                while (std::getline(file, line)) {
                        //cout << line << endl;
                        if (name_cntr > seq_cntr)
                        {
                		//seq_kmers.push_back({});
                                seqs[seq_cntr] = line; //seq_cntr] = line;
                
                		//string r_line = line;
                		//reverse_c(r_line);
                		//seqs_r.push_back(r_line);
                		
				//std::vector<unsigned long long> kmer_test;
                		//kmerize_all(line, k, kmer_test);
				
				/*for (size_t x=0; x < (line.length() - k+1) ; x++) {
                    			//std::string tmp = line.substr(x, k);
                    			kmer_f = kmerize(line, x, k);
                   			kmer_r = kmerize(r_line, (line.length()-k-x) , k);

                    			if (kmer_f < kmer_r) {
                        			seq_kmers[seq_cntr].push_back(kmer_f);
                    			} else {
                        			seq_kmers[seq_cntr].push_back(kmer_r);
                    			}
                		}*/

				//seq_kmers[seq_cntr] = kmer_test;
                		
				//pdqsort(seq_kmers[seq_cntr].begin(), seq_kmers[seq_cntr].end());
                		seq_cntr++;

                        }
                        else
                        {      
				//headers.push_back(line.substr(0, 50));
                                name_cntr++;
            		}
                }
        }
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
//find_similar_reads, std::ref(file), std::ref(min_dbs[i]), std::ref(read_mins[i]), std::ref(first_seqs[i]), std::ref(last_seqs[i]), std::ref(m), std::ref(cov_t)
//std::ref(line_num), std::ref(num_lines_per_read), std::ref(inputfileMutex1), input_file, std::ref(seq_kmers_threads[i]) , std::ref(read_mins[i]), std::ref(all_v_all_threads[i]) , std::ref(file), std::ref(min_dbs[i]), std::ref(m), std::ref(cov_t)
void find_similar_reads(int& k, robin_hood::unordered_flat_map<uint32_t, std::string>& seqs_threads, int& line_num, 
		int& num_lines_per_read, std::mutex& inputfileMutex, std::ifstream& fasta_file_que, 
		all_v_all_type& all_v_all_thread ,CKMCFile& file, 
		counter_t& min_db, int& m, int& cov_t)
{
#ifdef DEBUG_OUT
	auto start_time = chrono::high_resolution_clock::now();
#endif

    
    unsigned long long this_min;
    //std::string forward;
    //std::string reverse;
    //unsigned long long kmer_for;
    //unsigned long long kmer_rev;
    vector<uint32_t> tmp_v;
    //vector<uint32_t> another_tmp;
    unsigned long long prev_min;
    //std::string this_read;
    std::string line;
    int thread_line_num;
    //kmer_test;
    //std::cout << "In find similar reads\n";
    //std::cout.flush();
    inputfileMutex.lock();
    
    while (getline(fasta_file_que, line))
    {
	
	line_num++;
	//std::cout << line_num << "\n";
        //std::cout.flush();
        thread_line_num = line_num;
        if (thread_line_num % num_lines_per_read != 0){
                
		//header = ">" + line.substr(1);
        	inputfileMutex.unlock();
            	//this_read = line;
            	//seqs[rdnum] = "";
            	uint32_t read_number = (thread_line_num+num_lines_per_read-1)/num_lines_per_read - 1;
        	std::vector<uint32_t> cnts_q;
        	file.GetCountersForRead(line, cnts_q);
		
		//std::cout << "Got read counts\n";
    		//std::cout.flush();	
        	
		std::vector<unsigned long long> kmer_test;
        	kmerize_all(line, m, kmer_test);
		//std::cout << "kmerized\n";
                //std::cout.flush();
		//std::cout << this_min << "\n";
		
		for (size_t xm=0; xm <= (line.length() - k+ 1) ; xm+=15) {
             			
            			if (cnts_q[xm] > 6 && cnts_q[xm] < cov_t) {
                			
					this_min = 9223372036854775806;
					
					//What if we just took the smallest value in the kmer test in this range? 
                			for (int jm=0; jm < k-m ; jm++) {
                    				if (kmer_test[xm+jm] < this_min) {
                        				this_min = kmer_test[xm+jm];
                    				}
                			}

                			if (this_min != 9223372036854775806 && this_min != prev_min) {
                    				bool in_db = false;
                    				//Maybe we should save the minimizzers that are found for each read 

                    				for (auto rds1: min_db[this_min]) {
                        				if (rds1 == read_number) {
                            					in_db = true;
								//break();
                        				} 
							//all_v_all_thread[x].push_back(rds1);
							//all_v_all_thread[rds1].push_back(x);
                    				}

                    				if (!in_db) {
                        				min_db[this_min].push_back(read_number);
                        				//read_mins[read_number].push_back(this_min);
                    				}

                   				prev_min = this_min;
                			}
            			}
        	}
		
		//std::vector<unsigned long long> kmer_test_kmer;
		//kmerize_all(this_read, k, kmer_test_kmer);
		//seq_kmers_threads[read_number] = kmer_test_kmer;
		seqs_threads[read_number] = line;
      		//}
    		inputfileMutex.lock();
	}
    }
        
#ifdef DEBUG_OUT
	auto end_time = chrono::high_resolution_clock::now();
        auto total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
        section2 = section2 + total_time;
	//std::cout << "Section 2.5: " << total_time.count() << "\n";
	//std::cout.flush();
#endif
	
	//}
	
	inputfileMutex.unlock();
	
}

void make_all_v_all (uint32_t& first_seq, uint32_t& last_seq, int& k, int& m, all_v_all_type& all_v_all_threads , const counter_t& min_db_all)
{
#ifdef DEBUG_OUT
        auto start_time = chrono::high_resolution_clock::now();
#endif
	//uint32_t rdnum;
	//std::vector<unsigned long long> rdmins;
	
	//std::cout << read_mins.size() << "\n";
	std::cout << "Making all v all\n";
	//std::cout << "Start: " << first_seq << " Last: " << last_seq << "\n";
	std::cout.flush();
	//for ( robin_hood::pair<uint32_t, std::vector<unsigned long long>> element : read_mins ) {
	for (uint32_t i = first_seq; i < (last_seq); i++) {

		std::vector<unsigned long long> rdmins;
		//std::cout << i << "\n";
		//std::cout.flush();
		kmerize_all(seqs[i], m, rdmins);
		//std::cout << rdmins.size() << ",";
		//std::cout.flush();
		//rdnum = element.first;
		//rdmins = element.second;
		//for (auto this_min: rdmins) {
		for (int y = 0; y < rdmins.size(); y+=5) {
			auto this_min = rdmins[y];	
			//if (min_db_all.find() == m.) {//(rdmins[y] in min_db_all) {
			//auto tmpvec = min_db_all.find(rdmins[y]);
			//all_v_all_threads[i].insert(all_v_all_threads[i].end(), tmpvec.begin(), tmpvec.end());	
			if (min_db_all.find(this_min) != min_db_all.end()){
				all_v_all_threads[i].insert(all_v_all_threads[i].end(), min_db_all.at(this_min).begin(), min_db_all.at(this_min).end());
			}
				/*for (auto ele: min_db_all[this_min]) {
				all_v_all_threads[ele].push_back(this_min);
				all_v_all_threads[this_min].push_back(ele);
			}*/
		}
		//read_mins[rdnum].clear();
		//read_mins[rdnum].shrink_to_fit();
	}
	
	//std::cout << all_v_all_threads"\n";
	//std::cout.flush();
#ifdef DEBUG_OUT
        auto end_time = chrono::high_resolution_clock::now();
        auto total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
        section4 = section4 + total_time;
        std::cout << "Section 2.5: " << total_time.count() << "\n";
        std::cout.flush();
#endif

}

int main(int argc, char* argv[])
{
    //parse arguments
    std::time_t result = std::time(nullptr);
        auto start_time = chrono::high_resolution_clock::now();

    int c;
    std::ifstream input_file;
    int num_lines_per_read;
    std::string kmcdb;
    std::string outdir;
    int main_k = 0;
    double l = 0;
    int error_threshold = 0;
    int het_threshold = 0;
    int unique_threshold = 0;
    int anchor_threshold = 0;
    double allowed_err_fraction = 1;
    double allowed_rep_fraction = 1;
    int max_nodes_to_search = 1000;
    double distance_multiplier = 1.2;
    int strict = 1;
    int polish = 0;
    int num_threads = 1;
    int verbose = 0;
    int min_t = -1;
    int cov_t = 35;
    int m = 21;
    //std::cout << "Have not made it very far at all\n" ;    
    while ((c = getopt(argc, argv, "hi:j:o:k:M:C:z:l:g:m:a:u:t:pn:d:e:r:s:v")) != -1)
    {
        switch (c)
        {
            case 'h':
                fprintf(stderr, "Usage: %s -i input.fa/fq -j kmcdb -o outdir -k kmersize -l lambda [-t num_threads (default 1)] [-p (run in polish mode, i.e. run only error correction and not het smoothing)] [-n max_nodes_to_search (default 1000)] [-d distance_multiplier (default 1.2)] [-e allowed_err_fraction (default 1.0)] [-r allowed_rep_fraction (default 1.0)] [-s strict (0 or 1, default 1)] [-v (run in verbose mode, i.e. print time stamps for analyses)]\n", argv[0]);
                exit(EXIT_FAILURE);
            case 'i':
                //is the input fasta or fastq?
                //note: the output will be fasta format, since quality values
                //will not match once the read is edited and smoothed
                if ((getFileExt(optarg) == "fasta") || (getFileExt(optarg) == "fa"))
                {
                    num_lines_per_read = 2;
                }
                else if ((getFileExt(optarg) == "fastq") || (getFileExt(optarg) == "fq"))
                {
                    num_lines_per_read = 4;
                }
                else
                {
                    fprintf(stderr, "Input filename must end in .fa .fasta .fq or .fastq.\n");
                    exit(EXIT_FAILURE);
                }
                input_file.open(optarg);
                if (!input_file.is_open())
                {
                    fprintf(stderr, "Please ensure %s exists.\n", optarg);
                    exit(EXIT_FAILURE);
                }
                break;
            case 'j':
                kmcdb = optarg;
                break;
            case 'o':
                outdir = optarg;
		break;
            case 'k':
                main_k = atoi(optarg);
                break;
            case 'M':
                //similarity threshold (min_t)
                                min_t = atoi(optarg);
                                break;
            case 'C':
                                //similarity threshold (min_t)
                                cov_t = atoi(optarg);
                                break;
            case 'z':
                                //similarity threshold (min_t)
                                m = atoi(optarg);
                                break;
            case 'l':
                //only set values based on lambda if not already set by hidden parameters g, m, a, or u
                l = std::stod(optarg);
                if (error_threshold == 0)
                {
                    error_threshold = round(ceil(0.25 * l));
                }
                if (het_threshold == 0)
                {
                    het_threshold = round(ceil(1.5 * l));
                }
                if (unique_threshold == 0)
                {
                    unique_threshold = round(ceil(3.5 * l));
                }
                if (anchor_threshold == 0)
                {
                    anchor_threshold = round(ceil(2.5 * l));
                }
                break;
            case 'g':
                error_threshold = atoi(optarg);
                break;
            case 'm':
                het_threshold = atoi(optarg);
                break;
            case 'a':
                anchor_threshold = atoi(optarg);
                break;
            case 'u':
                unique_threshold = atoi(optarg);
                break;
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'p':
                polish = 1;
                break;
            case 'n':
                max_nodes_to_search = atoi(optarg);
                break;
            case 'd':
                distance_multiplier = std::stod(optarg);
                break;
            case 'e':
                allowed_err_fraction = std::stod(optarg);
                break;
            case 'r':
                allowed_rep_fraction = std::stod(optarg);
                break;
            case 's':
                strict = atoi(optarg);
                break;
            case 'v':
                verbose = 1;
                break;
            case '?':
                fprintf(stderr, "Option -%c is invalid or requires an argument.\n", optopt);
            default:
                fprintf(stderr, "Usage: %s -i input.fa/fq -j kmcdb -o outdir -k kmersize -l lambda [-t num_threads (default 1)] [-p (run in polish mode, i.e. run only error correction and not het smoothing)] [-n max_nodes_to_search (default 1000)] [-d distance_multiplier (default 1.2)] [-e allowed_err_fraction (default 1.0)] [-r allowed_rep_fraction (default 1.0)] [-s strict (0 or 1, default 1)] [-v (run in verbose mode, i.e. print time stamps for analyses)]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }
    
    //check that required arguments are given
    if (!input_file.is_open())
    {
        fprintf(stderr, "Please provide input file with -i argument.\n");
        exit(EXIT_FAILURE);
    }
    if (kmcdb.empty())
    {
        fprintf(stderr, "Please provide kmcdb with -j argument.\n");
        exit(EXIT_FAILURE);
    }
    if (outdir.empty())
    {
        fprintf(stderr, "Please provide output directory with -o argument.\n");
        exit(EXIT_FAILURE);
    }
    if (main_k==0)
    {
        fprintf(stderr, "Please provide kmer size with -k argument.\n");
        exit(EXIT_FAILURE);
    }
    if ((error_threshold == 0) || (het_threshold == 0) || (anchor_threshold == 0) || (unique_threshold == 0))
    {
        fprintf(stderr, "Please provide average kmer coverage with -l argument. (Or specify the thresholds manually).\n");
        exit(EXIT_FAILURE);
    }
    //std::cout << "We parsed the args\n";
    //std::cout.flush();
    //load KMC database
    always_call_this();
        CKMCFile file;
    
    file.OpenForRA(kmcdb);
    result = std::time(nullptr);
    //std::cout << "Opened the kmcdb\n";
    //std::cout.flush();
    //Now we need to make the minimizer DB. This can be a global variable?  TODO
    //First we will have to go sequence by sequence, find the kmers in each sequence,
    //use the getcounts function to get the counts, then send this off to the minimizer 
    //db maker. 
    //exit(1);
    //read_fasta_slow(input_file);
    //int seq_cntr=seqs.size();
    //std::cout << seq_cntr << "\n";
    /*
     * Make all of these maps and vectors to feed to the threads 
     */
    //exit(1); 
    //std::vector<read_min_map> read_mins(num_threads);
    std::vector<counter_t> min_dbs(num_threads);
    //std::vector<uint32_t> first_seqs;
    //std::vector<uint32_t> last_seqs;
    //std::vector<all_v_all_type> all_v_all_threads(num_threads);
    all_v_all.resize(num_threads);
    //std::vector<robin_hood::unordered_map<uint32_t, std::vector<unsigned long long>>> seq_kmers_threads(num_threads);
    std::vector<robin_hood::unordered_flat_map<uint32_t, std::string>> seqs_threads(num_threads);
    //std::cout << "Set up seq_threads\n";
    //std::cout.flush();
    
    //    int reads_per_thread = floor(seq_cntr/num_threads) + 1;
//    for (int i = 0; i < num_threads; i++) {
	/*read_min_map read_min;
	read_mins.push_back(read_min);

	counter_t min_db;
	min_dbs.push_back(min_db);

	robin_hood::unordered_flat_map<uint32_t, std::vector<unsigned long long>> skt;
        seq_kmers_threads.push_back(skt);
	
	all_v_all_type all_v_all_t;
	all_v_all_threads.push_back(all_v_all_t);
	*/
	
    //exit(1);
#ifdef DEBUG_OUT
    auto end_time = chrono::high_resolution_clock::now();
    auto total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
    std::cout << "Section 1:  " << total_time.count() << " NANOSECONDS" << '\n';
    std::cout.flush();
    start_time = chrono::high_resolution_clock::now();
#endif

    //Read in the sequences at the same time as finding the minimizers 
    //Need mutex to corridinate
    int line_num=-1; 
    std::mutex inputfileMutex1;
    std::vector<std::thread> threads1;
    //std::cout << "About to start finding similar reads\n";
    //std::cout.flush();
    for (int i = 0; i < num_threads; i++)
    {
            threads1.push_back(std::thread(find_similar_reads, std::ref(main_k), std::ref(seqs_threads[i]), std::ref(line_num), 
				    std::ref(num_lines_per_read), std::ref(inputfileMutex1), std::ref(input_file), 
				    std::ref(all_v_all[i]) , std::ref(file), std::ref(min_dbs[i]), std::ref(m), std::ref(cov_t)));
    
    }
    //std::cout << "First multithreading started\n";
    //std::cout.flush();
    for (auto &th : threads1)
    {
        th.join();
    }
    //std::cout << "First multithreading stopped\n";
    //std::cout.flush();
    input_file.clear();
    input_file.seekg(0, std::ios::beg);
    //int seq_cntr=seqs.size();
    //#ifdef DEBUG_OUT 
//    end_time = chrono::high_resolution_clock::now();
//    total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
//    std::cout << "Section 2:  " << total_time.count() << " NANOSECONDS" << '\n';
//    std::cout.flush();
//    auto start_time_final = chrono::high_resolution_clock::now();
//#endif

    /* 
     * Now we have to go through and merge all of these maps. A bit annoying, but hopefully will save time overall? 
     */
    //exit(1);
    //seqs.clear();
    //seqs.shrink_to_fit();
    //read_min_map read_mins_all;
    counter_t min_db_all;
    //all_v_all_type all_v_all;
    std::string read_tmp_vec;
    std::vector<uint32_t> min_tmp_vec;
    uint32_t read_key;
    unsigned long long min_key;
    std::vector<uint32_t> first_seqs(num_threads);
    std::vector<uint32_t> last_seqs(num_threads);
    //int reads_per_thread = floor(seq_cntr/num_threads) + 1;
    for (int i = 0; i < num_threads; i++) {
        
	for (robin_hood::pair<unsigned long long, std::vector<uint32_t>> element : min_dbs[i]) {
                min_key = element.first;
                min_tmp_vec = element.second;
		//min_db_all[min_key].insert(min_db_all[min_key].begin(), element.second.begin(), element.second.end());
		for (auto j: min_tmp_vec){
                        min_db_all[min_key].push_back(j);
			
                }
        }
	robin_hood::unordered_flat_map<unsigned long long, std::vector<uint32_t>> empty; 
	min_dbs[i] = empty;
	
	//uint32_t first_seq = (i*reads_per_thread);
	//first_seqs[i] = first_seq;
	//last_seqs[i] = (first_seq + reads_per_thread);
	/*for (robin_hood::pair<uint32_t, std::vector<unsigned long long>> element : seq_kmers_threads[i]) {
                read_key = element.first;
                read_tmp_vec = element.second;
                for (auto j: read_tmp_vec){
                        seq_kmers[read_key].push_back(j);

                }
        	
	}
	robin_hood::unordered_map<uint32_t, std::vector<unsigned long long>> empty2;
	seq_kmers_threads[i] = empty2;
    	*/

	for (robin_hood::pair<uint32_t, std::string> element : seqs_threads[i]) {
                //read_key = element.first;
                //read_tmp_vec = element.second;
                seqs[element.first] = element.second;
		/*for (auto j: read_tmp_vec){
                        seq_kmers[read_key].push_back(j);

                }*/

        }
        robin_hood::unordered_flat_map<uint32_t, std::string> empty2;
        seqs_threads[i] = empty2;

    }
    int reads_per_thread = floor(seqs.size()/num_threads) + 1;
    for (int i = 0; i < num_threads; i++) {

    	uint32_t first_seq = (i*reads_per_thread);
        first_seqs[i] = first_seq;
        last_seqs[i] = (first_seq + reads_per_thread);
    }
    //seq_kmers_threads.clear();
    min_dbs.clear();
    //seqs_threads.clear();
    //exit(1);
#ifdef DEBUG_OUT
    end_time = chrono::high_resolution_clock::now();
    total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
    std::cout << "Section 3:  " << total_time.count() << " NANOSECONDS" << '\n';
    std::cout.flush();
    start_time = chrono::high_resolution_clock::now();
#endif

    //Send this out to all of the threads (have them figure out the all_v_all
    std::vector<std::thread> threads2;
    for (int i = 0; i < num_threads; i++)
    {
	    
            threads2.push_back(std::thread(make_all_v_all, std::ref(first_seqs[i]), std::ref(last_seqs[i]), std::ref(main_k), std::ref(m) , std::ref(all_v_all[i]) , std::ref(min_db_all)));
    }

    //Combine the threads 
    for (auto &th : threads2)
    {
        th.join();
    }
    //exit(1);

    std::cout << "Finished making the all v all\n";
    std::cout.flush();
    uint32_t all_key ;
    std::vector<uint32_t> all_tmp_vec ;
    //for (int i = 0; i < num_threads; i++) {
        //Deal with read mins 
        //Deal with the min db
	/*for (robin_hood::pair<uint32_t, std::vector<uint32_t>> element : all_v_all_threads[i]) {
        //for (robin_hood::pair<unsigned long long, std::vector<uint32_t>> element : all_v_all_threads[i]) {
                all_key = element.first;
		//std::cout << all_key << "\n";
		//std::cout.flush();
		all_tmp_vec = element.second;
                for (auto j: all_tmp_vec){
                        all_v_all[all_key].push_back(j);

                }
        }*/
	//all_v_all.push_back(all_v_all_threads);
	//robin_hood::unordered_flat_map<uint32_t, std::vector<uint32_t>> empty;
	//all_v_all_threads[i] = empty;
    //}
    //exit(1);
    //std::cout << "About to clean up\n";
    //std::cout.flush();
    //clean up the min_db_all
    min_db_all.clear();
    counter_t min_db_empty;
    min_db_all = min_db_empty;
    //min_db_all.shrink_to_fit();
    //clean up the min_dbs
    min_dbs.clear();
    min_dbs.shrink_to_fit();
    //clean up the read_mins
    //read_mins.clear();
    //read_mins.shrink_to_fit();
    
    //min_dbs.clear();
    //std::cout << "Clean up done\n";
    //std::cout.flush();

    //delete min_db_all;
    //Now go through and build up the all_v_all system. This can maybe be done with multi-threads? 
    //This could look like 
    line_num = -1;    
    //make locks for the threads
    std::mutex inputfileMutex;
    //std::mutex outputfileMutex;
    std::mutex coutMutex;
    
//#ifdef DEBUG_OUT 
//    end_time = chrono::high_resolution_clock::now();
//    total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
//    std::cout << "Section 4:  " << total_time.count() << " NANOSECONDS" << '\n';
//    std::cout.flush();
//    start_time_final = chrono::high_resolution_clock::now();
//#endif

    //std::vector<std::vector <long long>> tmp_vec_vec(num_threads, std::vector<long long>(70000));
    int thread_reads ;

    std::vector<std::ofstream> err_output_file(num_threads);
    std::vector<std::ofstream> errpaths_queuecomplete0_numpaths0_output_file(num_threads);
    std::vector<std::ofstream> errpaths_queuecomplete0_numpaths1to2_output_file(num_threads);
    std::vector<std::ofstream> errpaths_queuecomplete0_numpaths3plus_output_file(num_threads);
    std::vector<std::ofstream> errpaths_queuecomplete1_numpaths0_output_file(num_threads);
    std::vector<std::ofstream> errpaths_queuecomplete1_numpaths1to2_output_file(num_threads);
    std::vector<std::ofstream> errpaths_queuecomplete1_numpaths3plus_output_file(num_threads);
    std::vector<std::ofstream> erredits_output_file(num_threads);
    std::vector<std::ofstream> het_output_file(num_threads);
    std::vector<std::ofstream> hetpaths_queuecomplete0_numpaths0to1_output_file(num_threads);
    std::vector<std::ofstream> hetpaths_queuecomplete0_numpaths2_output_file(num_threads);
    std::vector<std::ofstream> hetpaths_queuecomplete0_numpaths3plus_output_file(num_threads);
    std::vector<std::ofstream> hetpaths_queuecomplete1_numpaths0to1_output_file(num_threads);
    std::vector<std::ofstream> hetpaths_queuecomplete1_numpaths2_output_file(num_threads);
    std::vector<std::ofstream> hetpaths_queuecomplete1_numpaths3plus_output_file(num_threads);
    std::vector<std::ofstream> hetedits_output_file(num_threads);

    
    for (int i = 0; i < num_threads; i++) {
	//    std::ofstream a; a.open(outdir+"/errremoved.fasta");
	err_output_file[i].open(outdir+"/errremoved."+std::to_string(i)+".fasta");
        errpaths_queuecomplete0_numpaths0_output_file[i].open(outdir+"/errpaths_queuecomplete0_numpaths0."+std::to_string(i)+".fasta");
        errpaths_queuecomplete0_numpaths1to2_output_file[i].open(outdir+"/errpaths_queuecomplete0_numpaths1to2."+std::to_string(i)+".fasta");
        errpaths_queuecomplete0_numpaths3plus_output_file[i].open(outdir+"/errpaths_queuecomplete0_numpaths3plus."+std::to_string(i)+".fasta");
        errpaths_queuecomplete1_numpaths0_output_file[i].open(outdir+"/errpaths_queuecomplete1_numpaths0."+std::to_string(i)+".fasta");
        errpaths_queuecomplete1_numpaths1to2_output_file[i].open(outdir+"/errpaths_queuecomplete1_numpaths1to2."+std::to_string(i)+".fasta");
        errpaths_queuecomplete1_numpaths3plus_output_file[i].open(outdir+"/errpaths_queuecomplete1_numpaths3plus."+std::to_string(i)+".fasta");
        erredits_output_file[i].open(outdir+"/erredits."+std::to_string(i)+".fasta");
        het_output_file[i].open(outdir+"/hetremoved."+std::to_string(i)+".fasta");
        hetpaths_queuecomplete0_numpaths0to1_output_file[i].open(outdir+"/hetpaths_queuecomplete0_numpaths0to1."+std::to_string(i)+".fasta");
        hetpaths_queuecomplete0_numpaths2_output_file[i].open(outdir+"/hetpaths_queuecomplete0_numpaths2."+std::to_string(i)+".fasta");
        hetpaths_queuecomplete0_numpaths3plus_output_file[i].open(outdir+"/hetpaths_queuecomplete0_numpaths3plus."+std::to_string(i)+".fasta");
        hetpaths_queuecomplete1_numpaths0to1_output_file[i].open(outdir+"/hetpaths_queuecomplete1_numpaths0to1."+std::to_string(i)+".fasta");
        hetpaths_queuecomplete1_numpaths2_output_file[i].open(outdir+"/hetpaths_queuecomplete1_numpaths2."+std::to_string(i)+".fasta");
        hetpaths_queuecomplete1_numpaths3plus_output_file[i].open(outdir+"/hetpaths_queuecomplete1_numpaths3plus."+std::to_string(i)+".fasta");
        hetedits_output_file[i].open(outdir+"/hetedits."+std::to_string(i)+".fasta");

    
    }

#ifdef DEBUG_OUT 
    end_time = chrono::high_resolution_clock::now();
    total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
    std::cout << "Section 5:  " << total_time.count() << " NANOSECONDS" << '\n';
    std::cout.flush();
    start_time = chrono::high_resolution_clock::now();
#endif

    //std::vector<std::string>& seqs, std::vector<std::string>& headers, std::vector<uint32_t>& our_reads
    std::vector<std::thread> threads;
    for (int i = 0; i < num_threads; i++)
    {
	    //std::cout <<  i << "\n";
	    //std::cout.flush();
	    threads.push_back(std::thread(processRead, std::ref(coutMutex) , std::ref(main_k), std::ref(inputfileMutex),
				    std::ref(cov_t), std::ref(min_t), std::ref(m), std::ref(input_file),
				    std::ref(line_num), std::ref(num_lines_per_read), std::ref(file), std::ref(error_threshold), 
				    std::ref(err_output_file[i]), std::ref(errpaths_queuecomplete0_numpaths0_output_file[i]), std::ref(errpaths_queuecomplete0_numpaths1to2_output_file[i]), std::ref(errpaths_queuecomplete0_numpaths3plus_output_file[i]), 
				    std::ref(errpaths_queuecomplete1_numpaths0_output_file[i]), std::ref(errpaths_queuecomplete1_numpaths1to2_output_file[i]), std::ref(errpaths_queuecomplete1_numpaths3plus_output_file[i]), std::ref(erredits_output_file[i]), 
				    std::ref(het_threshold), std::ref(unique_threshold), std::ref(max_nodes_to_search), std::ref(distance_multiplier), std::ref(polish), std::ref(allowed_rep_fraction),
				    std::ref(het_output_file[i]), std::ref(hetpaths_queuecomplete0_numpaths0to1_output_file[i]), std::ref(hetpaths_queuecomplete0_numpaths2_output_file[i]), std::ref(hetpaths_queuecomplete0_numpaths3plus_output_file[i]), 
				    std::ref(hetpaths_queuecomplete1_numpaths0to1_output_file[i]), std::ref(hetpaths_queuecomplete1_numpaths2_output_file[i]), std::ref(hetpaths_queuecomplete1_numpaths3plus_output_file[i]), std::ref(hetedits_output_file[i]), 
				    std::ref(anchor_threshold), std::ref(strict), std::ref(verbose)));
    	
    }
//#ifdef DEBUG_OUT     
//    auto end_time_final = chrono::high_resolution_clock::now();
//        total_time = chrono::duration_cast<chrono::milliseconds>(end_time_final - start_time_final);

        //std::cout << "Processing the reads took " << total_time.count() << " MILLISECONDS" << '\n';
        //std::cout.flush();
    
//#endif
    for (auto &th : threads)
    {
        th.join();
    }
    //delete min_dbs; 
#ifdef DEBUG_OUT
    end_time = chrono::high_resolution_clock::now();
    total_time = chrono::duration_cast<chrono::nanoseconds>(end_time - start_time);
    std::cout << "Section 6:  " << total_time.count() << " NANOSECONDS" << '\n';
    std::cout.flush();
    start_time = chrono::high_resolution_clock::now();
#endif
    


    //close files
    input_file.close();
    for (int i = 0; i < num_threads; i++) {
    err_output_file[i].close();
    errpaths_queuecomplete0_numpaths0_output_file[i].close();
    errpaths_queuecomplete0_numpaths1to2_output_file[i].close();
    errpaths_queuecomplete0_numpaths3plus_output_file[i].close();
    errpaths_queuecomplete1_numpaths0_output_file[i].close();
    errpaths_queuecomplete1_numpaths1to2_output_file[i].close();
    errpaths_queuecomplete1_numpaths3plus_output_file[i].close();
    erredits_output_file[i].close();
    het_output_file[i].close();
    hetpaths_queuecomplete0_numpaths0to1_output_file[i].close();
    hetpaths_queuecomplete0_numpaths2_output_file[i].close();
    hetpaths_queuecomplete0_numpaths3plus_output_file[i].close();
    hetpaths_queuecomplete1_numpaths0to1_output_file[i].close();
    hetpaths_queuecomplete1_numpaths2_output_file[i].close();
    hetpaths_queuecomplete1_numpaths3plus_output_file[i].close();
    hetedits_output_file[i].close();
    }
#ifdef DEBUG_OUT 

    std::cout << "Section 2:  " << section2.count() << "\n";
    std::cout << "Section 4:  " << section4.count() << "\n";
    std::cout << "Section 7:  " << section7.count() << "\n";
        std::cout << "Section 8:  " << section8.count() << "\n";
        std::cout << "Section 9:  " << section9.count() << "\n";
        std::cout << "Section 10:  " << section10.count() << "\n";
        std::cout << "Section 11:  " << section11.count() << "\n";
        std::cout << "Section 12:  " << section12.count() << "\n";
        //std::cout << "Section 13:  " << section9.count() << "\n";
        //std::cout << "Section 10: " << section10.count() << "\n";
        std::cout << "Section 13:  " << section13.count() << "\n";
        std::cout << "Section 14: " << section14.count() << "\n";
        std::cout << "Section 15: " << section15.count() << "\n";
        std::cout << "Section 16: " << section16.count() << "\n";
        std::cout << "Section 17: " << section17.count() << "\n";
        std::cout << "Section 18: " << section18.count() << "\n";
        std::cout << "Section 19: " << section19.count() << "\n";
        std::cout << "Section 20: " << section20.count() << "\n";
        std::cout << "Section 21: " << section21.count() << "\n";
	std::cout << "Section 22: " << section22.count() << "\n";
#endif
    return 0;
}
