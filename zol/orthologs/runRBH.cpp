/*
AUTHOR - Rauf Salamzade
DATE - 09/14/22
PROGRAM NAME - runRBH.cpp
DESCRIPTION - Find orthologs using RBH after normalization with reflexive query bitscore and phylogenetic normalization
              using expected differences based on universal single-copy genes from Hug et al. 2016. In-paralogs
              identified if they are not part of an RBH pairing but exhibit a higher bitscore than what was observed for
              the lowest RBH pairing (discounting for phylogenetic bitscore). Result is a pairing listing orthologs and
              in-paralogs with weights assigned as normalized bitscores in the ABC format needed to run Markov Clustering.

Usage:
intraClusterRBH <BLAST/DIAMOND output format 6 Result file>
""
*/

#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <math.h>
#include <tuple>

using namespace std;

string delim = "\t";
string pipe_delim = "|";

/*
Sort vector of tuples in descending order. Solution taken/adapted from:
https://www.geeksforgeeks.org/sorting-vector-tuple-c-descending-order/
*/
bool sortdesc(const tuple<double, string>& a, const tuple<double, string>& b) {
    return (get<0>(a) > get<0>(b));
}

/*
Driver function to sort the vector elements by second element of pairs. Solution taken from:
https://www.geeksforgeeks.org/sorting-vector-of-pairs-in-c-set-1-sort-by-first-and-second/
*/
bool sortbysec(const pair<int,int> &a, const pair<int,int> &b) {
    return (a.second < b.second);
}

/*
Taken from Arafat Hasan's response on StackOverflow:
https://stackoverflow.com/questions/14265581/parse-split-a-string-in-c-using-string-delimiter-standard-c
*/
// for string delimiter
vector<string> split (string s, string delimiter) {
    size_t pos_start = 0, pos_end, delim_len = delimiter.length();
    string token;
    vector<string> res;

    while ((pos_end = s.find (delimiter, pos_start)) != string::npos) {
        token = s.substr (pos_start, pos_end - pos_start);
        pos_start = pos_end + delim_len;
        res.push_back (token);
    }

    res.push_back (s.substr (pos_start));
    return res;
}

string doubleToString(double val) {
	string result;
	ostringstream convertd;
	convertd << val;
	result = convertd.str();
	return result;
}

int main (int argc, char* argv[]) {
    if ( argv[1]==NULL || argv[2]==NULL || (argv[1][0]=='-' && argv[1][1]=='h') || (argv[1][0]=='-' && argv[1][1]=='-' && argv[1][2]=='h') ) {
	    cout << "Usage:" << endl;
	    cout << "runRBH <BLAST/DIAMOND Alignment Results - format 6 with query coverage as an additional final column.> <identity cutoff> <coverage cutoff>" << endl;
	    return 0;
    }
    else {
        int split_counter;
        string line;
        vector<string> v;
        string sA, sB;

        // read in arguments for coverage and identity cutoffs
        double coverage_cutoff, identity_cutoff;
        sscanf(argv[2],"%lf",&identity_cutoff);
        sscanf(argv[3],"%lf",&coverage_cutoff);

        /*
        Parse out all hits to perform gene length normalization.
        */
        ifstream input_file;
        input_file.open (argv[1]);
        string q, h;
        double bitscore, identity, coverage;
        map<string, double> reflexive_bitscores;
        if ( input_file.is_open() ) {
            while ( input_file.good()  ) {
                getline (input_file,line);
                if (!line.empty()) {
                    split_counter = 0;
                    v = split (line, delim);
                    for (auto i : v) {
                        if (split_counter == 0) {
                            q = i;
                        }
                        else if (split_counter == 1) {
                            h = i;
                        }
                        else if (split_counter == 11) {
                            bitscore = stod(i);
                        }
                        split_counter++;
                    }
                    if (q.compare(h) == 0) {
                        reflexive_bitscores[q] = bitscore;
                    }
                }
            }
        } else {
  	        cout << "ERROR: Unable to open file " + (string)argv[1] << endl;
        }

        ifstream input_file_reopen;
        input_file_reopen.open (argv[1]);

        map<pair<string, string>, double> paired_normalized_bitscores;
        map<pair<string, string>, double> query_max_bitscore;
        map<pair<string, string>, set<string>> query_max_hits;
        set<pair<string, string>>::iterator align_iter;
        pair<string, string> align, query_hit_sample_pair, qshs_pair;
        string query, hit, query_sample, hit_sample;
        set<string> all_samples;
        set<pair<string, string>> all_queries;
        double normalized_bitscore, reflexive_bitscore_factor;
        if ( input_file_reopen.is_open() ) {
            while ( input_file_reopen.good()  ) {
                getline (input_file_reopen,line);
                if (!line.empty()) {
                    split_counter = 0;
                    v = split (line, delim);
                    for (auto i : v) {
                        if (split_counter == 0) {
                            query = i;
                        }
                        else if (split_counter == 1) {
                            hit = i;
                        }
                        else if (split_counter == 11) {
                            bitscore = stod(i);
                        }
                        else if (split_counter == 2) {
                            identity = stod(i);
                        }
                        else if (split_counter == 12) {
                            coverage = stod(i);
                        }
                        split_counter++;
                    }
                    if (query.compare(hit) != 0 && identity >= identity_cutoff && coverage >= coverage_cutoff) {

                        query_sample = query.substr(0, query.find('|'));
                        hit_sample = hit.substr(0, hit.find('|'));
                        reflexive_bitscore_factor = reflexive_bitscores[query];
                        qshs_pair = make_pair(query_sample, hit_sample);
                        normalized_bitscore = (bitscore)/(reflexive_bitscore_factor);

                        /*
                        cout << hit_sample << endl;
                        cout << query_sample << endl;
                        cout << reflexive_bitscore_factor << endl;
                        cout << phylo_normalization_factor << endl;
                        cout << qh_bitscore_prime << endl;
                        cout << normalized_bitscore << endl;
                        cout << "--------------------" << endl;
                        */

                        query_hit_sample_pair = make_pair(query, hit_sample);
                        if (query_sample.compare(hit_sample) != 0) {
                            if (query_max_bitscore.find(query_hit_sample_pair) == query_max_bitscore.end()) {
                                query_max_bitscore[query_hit_sample_pair] = normalized_bitscore;
                                query_max_hits[query_hit_sample_pair].insert(hit);
                            } else if (query_max_bitscore.at(query_hit_sample_pair) == normalized_bitscore) {
                                query_max_hits[query_hit_sample_pair].insert(hit);
                            } else if (query_max_bitscore.at(query_hit_sample_pair) < normalized_bitscore) {
                                query_max_bitscore[query_hit_sample_pair] = normalized_bitscore;
                                query_max_hits[query_hit_sample_pair].clear();
                                query_max_hits[query_hit_sample_pair].insert(hit);
                            }
                        }
                        paired_normalized_bitscores[make_pair(query, hit)] = normalized_bitscore;
                        all_queries.insert(make_pair(query, query_sample));
                        all_samples.insert(hit_sample);
                    }
                }
            }
        }
        else {
  	        cout << "ERROR: Unable to open file " + (string)argv[1] << endl;
        }

        // Look for high-quality ortholog prospects based on direct RBH
        string qs, qid, hs, hid;
        set<pair<string, string>>::iterator qid_itr, hid_itr;
        set<string>::iterator hs_itr;
        set<string> best_hits_for_qhs, best_hits_for_hqs;
        pair<string, string> pair_qhs, pair_hqs, qh_pair, hq_pair, pair_qshs;
        double bidir_bs, threshold, qh_max_bitscore;
        set<pair<string, string>> accounted;
        map<string, double> minimum_query_rbh_bitscores;
        bool rbh_found;
        for (qid_itr = all_queries.begin(); qid_itr != all_queries.end(); qid_itr++) {
            qid = (*qid_itr).first;
            qs = (*qid_itr).second;
            for (hs_itr = all_samples.begin(); hs_itr != all_samples.end(); hs_itr++) {
                hs = *hs_itr;
                rbh_found = false;
                if (qs.compare(hs) != 0) {
                    pair_qhs = make_pair(qid, hs);
                    pair_qshs = make_pair(qs, hs);
                    qh_max_bitscore = query_max_bitscore[pair_qhs];
                    best_hits_for_qhs = query_max_hits[pair_qhs];
                    for (auto hid : best_hits_for_qhs) {
                        qh_pair = make_pair(qid, hid);
                        hq_pair = make_pair(hid, qid);
                        pair_hqs = make_pair(hid, qs);
                        best_hits_for_hqs = query_max_hits[pair_hqs];
                        if (best_hits_for_hqs.find(qid) != best_hits_for_hqs.end()) {
                            if (accounted.find(qh_pair) == accounted.end() && accounted.find(hq_pair) == accounted.end()) {
                                bidir_bs = 50.0*(paired_normalized_bitscores[qh_pair] + paired_normalized_bitscores[hq_pair]);
                                cout <<  qid + '\t' + qs + '\t' + hid + '\t' + hs + '\t' + doubleToString(bidir_bs) << endl;
                                if (bidir_bs >= 0.0) {
                                    rbh_found = true;
                                    accounted.insert(qh_pair);
                                }
                            }
                        }
                    }
                    if (rbh_found) {
                        if (minimum_query_rbh_bitscores.find(qid) == minimum_query_rbh_bitscores.end()) {
                            minimum_query_rbh_bitscores[qid] = qh_max_bitscore;
                        } else if (minimum_query_rbh_bitscores.at(qid) > qh_max_bitscore) {
                            minimum_query_rbh_bitscores[qid] = qh_max_bitscore;
                        }
                    }
                }
            }
        }

        // Look for non-RBH pairs to consider as part of orthogroups
        for (qid_itr = all_queries.begin(); qid_itr != all_queries.end(); qid_itr++) {
            qid = (*qid_itr).first;
            qs = (*qid_itr).second;
            if (minimum_query_rbh_bitscores.find(qid) != minimum_query_rbh_bitscores.end()) {
                threshold = minimum_query_rbh_bitscores[qid];
            } else {
                threshold = 10000.0;
            }
            //cout << "-------------------" << endl;
            //cout << qid << endl;
            //cout << threshold << endl;
            for (hid_itr = all_queries.begin(); hid_itr != all_queries.end(); hid_itr++) {
                hid = (*hid_itr).first;
                hs = (*hid_itr).second;
                if (hs.compare(qs) == 0) {
                    qh_pair = make_pair(qid, hid);
                    hq_pair = make_pair(hid, qid);
                    pair_qshs = make_pair(qs, hs);
                    if (accounted.find(qh_pair) == accounted.end() && accounted.find(hq_pair) == accounted.end() && paired_normalized_bitscores.find(qh_pair) != paired_normalized_bitscores.end()) {
                        if (paired_normalized_bitscores.at(qh_pair) >= threshold && paired_normalized_bitscores.find(hq_pair) != paired_normalized_bitscores.end()) {
                            bidir_bs = 50.0*(paired_normalized_bitscores.at(qh_pair) + paired_normalized_bitscores.at(hq_pair));
                            cout << qid + '\t' + qs + '\t' + hid + '\t' + hs + '\t' + doubleToString(bidir_bs) << endl;
                            accounted.insert(qh_pair);
                        }
                    }
                }
            }
        }
    }
    return 0;
}