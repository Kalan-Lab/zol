/*
AUTHOR - Rauf Salamzade
DATE - 09/14/22 (rewritten 11/20/26)
PROGRAM NAME - splitDiamondResults.cpp
DESCRIPTION - Split a reflexive (all-vs-all) DIAMOND alignment into per-sample files for runRBH.

Input alignment columns (DIAMOND --outfmt 6): qseqid sseqid pident bitscore qcovhsp
Protein identifiers are expected as "<sample>|<protein>".

For a focal sample A, runRBH needs: hits within A, hits between A and every other sample (both directions) and,
for each query, its reflexive (self-hit) bitscore used for normalization. So:
  - pass 1 records the reflexive bitscore of every protein (one entry per protein),
  - pass 2 writes each non-self hit with the query's reflexive bitscore appended as a 6th column:
      * hits between two different samples go to both samples' files,
      * hits within a sample go only to that sample's file,
      * self hits are not written (their bitscore is carried in the 6th column).
Output lines are buffered in memory and flushed to the per-sample files once the buffered total exceeds a cap, so
the program never keeps more than a handful of files open and memory stays bounded regardless of sample count.

Usage:
splitDiamondResults <DIAMOND alignment> <Sample Listing: sample<TAB>output file> [buffer cap in MB, default 256]
*/

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
using namespace std;

// Extract the first n tab-delimited fields of line into fields; returns false if the line has fewer.
static bool get_fields(const string& line, size_t n, vector<string>& fields) {
    fields.clear();
    size_t start = 0;
    while (fields.size() < n) {
        size_t end = line.find('\t', start);
        if (end == string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, end - start));
        start = end + 1;
    }
    return fields.size() >= n;
}

static string sample_of(const string& protein) {
    return protein.substr(0, protein.find('|'));
}

int main(int argc, char* argv[]) {
    if (argc < 3 || (argv[1][0] == '-' && argv[1][1] == 'h') || (argv[1][0] == '-' && argv[1][1] == '-' && argv[1][2] == 'h')) {
        cout << "Usage:" << endl;
        cout << "./splitDiamondResults <DIAMOND alignment: qseqid sseqid pident bitscore qcovhsp> <Sample Listing> [buffer cap in MB, default 256]" << endl;
        return 0;
    }
    double buffer_cap_mb = 256.0;
    if (argc > 3) {
        buffer_cap_mb = atof(argv[3]);
        if (buffer_cap_mb <= 0.0) buffer_cap_mb = 256.0;
    }
    const size_t buffer_cap = (size_t)(buffer_cap_mb * 1024.0 * 1024.0);

    // Sample listing: sample -> output file
    unordered_map<string, string> sample_to_outfile_path;
    vector<string> samples;
    {
        ifstream listing(argv[2]);
        if (!listing.is_open()) {
            cerr << "ERROR: Unable to open file " << argv[2] << endl;
            return 1;
        }
        string line;
        vector<string> fields;
        while (getline(listing, line)) {
            if (line.empty()) continue;
            if (!get_fields(line, 2, fields)) continue;
            if (sample_to_outfile_path.find(fields[0]) == sample_to_outfile_path.end()) samples.push_back(fields[0]);
            sample_to_outfile_path[fields[0]] = fields[1];
        }
    }

    // Pass 1: reflexive (self-hit) bitscores, kept as the original text so runRBH parses exactly the same value.
    unordered_map<string, string> reflexive_bitscores;
    {
        ifstream alignment(argv[1]);
        if (!alignment.is_open()) {
            cerr << "ERROR: Unable to open file " << argv[1] << endl;
            return 1;
        }
        string line;
        vector<string> fields;
        while (getline(alignment, line)) {
            if (line.empty() || !get_fields(line, 4, fields)) continue;
            if (fields[0] == fields[1]) reflexive_bitscores[fields[0]] = fields[3];
        }
    }

    // Pass 2: route non-self hits (with the query's reflexive bitscore appended) to per-sample buffers.
    unordered_map<string, string> buffers;
    unordered_set<string> created;
    size_t buffered = 0;

    auto flush_sample = [&](const string& sample, string& buffer) -> bool {
        ios_base::openmode mode = ios_base::out | (created.count(sample) ? ios_base::app : ios_base::trunc);
        ofstream out(sample_to_outfile_path[sample], mode);
        if (!out.is_open()) {
            cerr << "ERROR: Unable to write to file " << sample_to_outfile_path[sample] << endl;
            return false;
        }
        out << buffer;
        created.insert(sample);
        string().swap(buffer);  // release the memory (clear() would keep the capacity)
        return true;
    };
    auto flush_all = [&]() -> bool {
        for (auto& kv : buffers) {
            if (!kv.second.empty() && !flush_sample(kv.first, kv.second)) return false;
        }
        buffered = 0;
        return true;
    };

    {
        ifstream alignment(argv[1]);
        string line, out_line;
        vector<string> fields;
        while (getline(alignment, line)) {
            if (line.empty() || !get_fields(line, 2, fields)) continue;
            const string& query = fields[0];
            const string& hit = fields[1];
            if (query == hit) continue;

            auto refl = reflexive_bitscores.find(query);
            out_line = line;
            out_line += '\t';
            out_line += (refl == reflexive_bitscores.end() ? string("0") : refl->second);
            out_line += '\n';

            string query_sample = sample_of(query);
            string hit_sample = sample_of(hit);
            if (sample_to_outfile_path.count(query_sample)) {
                buffers[query_sample] += out_line;
                buffered += out_line.size();
            }
            if (hit_sample != query_sample && sample_to_outfile_path.count(hit_sample)) {
                buffers[hit_sample] += out_line;
                buffered += out_line.size();
            }
            if (buffered > buffer_cap && !flush_all()) return 1;
        }
    }
    if (!flush_all()) return 1;

    // Every listed sample gets a (possibly empty) file for runRBH.
    for (const auto& sample : samples) {
        if (!created.count(sample)) {
            string empty;
            if (!flush_sample(sample, empty)) return 1;
        }
    }
    return 0;
}
