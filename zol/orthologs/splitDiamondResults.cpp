/*
AUTHOR - Rauf Salamzade
DATE - 12/24/22
PROGRAM NAME - splitDiamondResults.cpp
DESCRIPTION - Split Diamond alignment results to different files based on queries.

Large

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
using namespace std;

string delim = "\t";

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

int main (int argc, char* argv[]) {
    if ( argv[1]==NULL || (argv[1][0]=='-' && argv[1][1]=='h') || (argv[1][0]=='-' && argv[1][1]=='-' && argv[1][2]=='h') ) {
	    cout << "Usage:" << endl;
	    cout << "./splitDiamondResults <BLAST/DIAMOND output> <Sample Listing>" << endl;
	    return 0;
    }
    else {
        /*
        Read in sample listing and create map of sample names to out-files
        */
        map<string, ofstream*> sample_to_outfile;
        string line, sample, outfile_path, query, subject, query_sample, subject_sample;
        int split_counter;
        vector<string> v;
        ifstream input_file;
        input_file.open (argv[2]);
        set<string> all_samples;
        if (input_file.is_open()) {
            while (input_file.good()) {
                getline (input_file,line);
                if (!line.empty()) {
                    split_counter = 0;
                    v = split (line, delim);
                    for (auto i : v) {
                        if (split_counter == 0) {
                            sample = i;
                        }
                        else if (split_counter == 1) {
                            outfile_path = i;
                        }
                        split_counter++;
                    }
                    sample_to_outfile[sample] = new ofstream(outfile_path);
                    all_samples.insert(sample);
                }
            }
        } else {
  	        cout << "ERROR: Unable to open file " + (string)argv[2] << endl;
        }
        input_file.close();

        /*
        Parse DIAMOND results and start writing to individual files.
        */
        input_file.open (argv[1]);
        string line_with_newline;
        if ( input_file.is_open() ) {
            while ( input_file.good()  ) {
                getline (input_file,line);
                if (!line.empty()) {
                    split_counter = 0;
                    v = split(line, delim);
                    for (auto i: v) {
                         if (split_counter == 0) {
                             query = i;
			     query_sample = i.substr(0, i.find('|'));

                         } else if (split_counter == 1) {
		             subject = i;
                             subject_sample = i.substr(0, i.find('|'));
                         }
                         split_counter++;
                    }
                    line_with_newline = line + '\n';
		    if (query.compare(subject) == 0) {
			for (auto& pair : sample_to_outfile) {
		             *sample_to_outfile[pair.first] << line_with_newline;
			}
		    } else if (query_sample.compare(subject_sample) == 0) {
	                *sample_to_outfile[query_sample] << line_with_newline;
	            } else {
                        *sample_to_outfile[query_sample] << line_with_newline;
		        *sample_to_outfile[subject_sample] << line_with_newline;
		    }
		}
            }
        } else {
  	        cout << "ERROR: Unable to open file " + (string)argv[1] << endl;
        }

        // Close all of the files
        for(auto& pair : sample_to_outfile) {
            delete pair.second;
            pair.second = 0;
        }
        return 0;
    }
}