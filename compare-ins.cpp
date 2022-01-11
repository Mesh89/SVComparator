#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <set>
#include <sstream>
#include <stdio.h>
#include <unistd.h>
#include "lib/IntervalTree.h"
#include "lib/ssw_cpp.h"
#include "lib/kseq.h"
KSEQ_INIT(int, read)

#include "common.h"

std::unordered_map<std::string, std::string> chrs;
int maxdist;
StripedSmithWaterman::Aligner aligner(2,2,4,1);
StripedSmithWaterman::Filter filter(false, false, 0, 32767);
StripedSmithWaterman::Alignment alignment;

bool check_seq = true;

int dist(bp_t& bp1, bp_t& bp2) {
    if (bp1.chr != bp2.chr) return INT32_MAX;
    return abs(bp1.pos - bp2.pos);
}
int dist(sv_t& bsv, sv_t& csv) {
	int d1 = std::min(dist(bsv.bp1, csv.bp1), dist(bsv.bp2, csv.bp1));
    int d2 = std::min(dist(bsv.bp1, csv.bp2), dist(bsv.bp2, csv.bp2));
	return std::min(d1, d2);
}

std::string rev_comp(std::string& s) {
    std::string rc_s(s.rbegin(), s.rend());
    for (int i = 0; i < rc_s.length(); i++) {
        rc_s[i] = toupper(rc_s[i]);
        if (rc_s[i] == 'A') rc_s[i] = 'T';
        else if (rc_s[i] == 'C') rc_s[i] = 'G';
        else if (rc_s[i] == 'G') rc_s[i] = 'C';
        else if (rc_s[i] == 'T') rc_s[i] = 'A';
    }
    return rc_s;
}

bool check_ins_seq(sv_t& bsv, sv_t& csv) {
    if (!check_seq) return true;
    if (bsv.size() > 50000 || csv.size() > 50000) return true;

    std::string cseq;
    if (csv.type == "DUP") {
        cseq = chrs[csv.bp1.chr].substr(csv.bp1.pos, csv.bp2.pos-csv.bp1.pos);
        std::string extended_cseq;
        while (extended_cseq.length() <= csv.seq.length()) extended_cseq += cseq;
        cseq = extended_cseq;
    } else if (csv.type == "TRA") {
        // the stable bp in a TRA is the one near the insertion point, unstable the other one
        bp_t unstable_bp = dist(bsv.bp1, csv.bp1) > dist(bsv.bp1, csv.bp2) ? csv.bp1 : csv.bp2;
        int a = unstable_bp.pos - maxdist - bsv.seq.length();
        a = std::max(a, 0);
        int b = unstable_bp.pos + maxdist + bsv.seq.length();
        cseq = chrs[unstable_bp.chr].substr(a, b-a);
    } else if (csv.type == "INS") {
        if (abs(bsv.size()-csv.size()) > maxdist) return false;
        cseq = csv.seq;
    }

    aligner.Align(bsv.seq.data(), cseq.data(), cseq.length(), filter, &alignment);
    if (alignment.sw_score >= std::min(bsv.seq.length(), cseq.length())) return true;

//    cseq = rev_comp(cseq);
//    aligner.Align(bsv.seq.data(), cseq.data(), cseq.length(), filter, &alignment);
//    return alignment.sw_score >= std::min(bsv.seq.length(), cseq.length());
    return false;
}

int main(int argc, char* argv[]) {
    
	if (argc < 6) {
        std::cout << "Given a VCF or SV file with benchmark insertions and one with the called ones, reports "
        		"for each benchmark insertion if it is being called." << std::endl;
        std::cout << "Usage: exec benchmark_file called_file ucsc_repeats ref_genome maxdist [--ignore-seq]" << std::endl;
		return 0;
	}

	std::vector<sv_t> benchmark_svs = read_sv_list(argv[1]);
	std::vector<sv_t> called_svs = read_sv_list(argv[2]);

	auto is_ins_func = [](const sv_t& sv) {return sv.type != "INS" && sv.type != "DUP";};
	benchmark_svs.erase(std::remove_if(benchmark_svs.begin(), benchmark_svs.end(), is_ins_func), benchmark_svs.end());
	called_svs.erase(std::remove_if(called_svs.begin(), called_svs.end(), is_ins_func), called_svs.end());

    std::ifstream rep_f(argv[3]);
	maxdist = std::stoi(argv[5]);

	bool report = false, print_fp = false;
	for (int i = 6; i < argc; i++) {
		std::string argv_str = std::string(argv[i]);
		if (argv_str == "--ignore-seq") {
			check_seq = false;
		}
		if (argv_str == "--report") {
			report = true;
		}
		if (argv_str == "--print-fp") {
			print_fp = true;
		}
		if (argv_str == "--force-ids") {
			for (int j = 0; j < benchmark_svs.size(); j++) benchmark_svs[j].id = "INS_" + std::to_string(j);
			for (int j = 0; j < called_svs.size(); j++) called_svs[j].id = "INS_" + std::to_string(j);
		}
	}

    FILE* fastaf = fopen(argv[4], "r");
    kseq_t *seq = kseq_init(fileno(fastaf));
    while (kseq_read(seq) >= 0) {
        chrs[std::string(seq->name.s)] = seq->seq.s;
    }
    kseq_destroy(seq);

    std::unordered_map<std::string, std::vector<sv_t> > called_svs_by_chr;
    for (sv_t& sv : called_svs) {
    	if (check_seq && !sv.has_seq() && sv.type == "INS") continue;
		called_svs_by_chr[sv.bp1.chr].push_back(sv);
	}

    std::unordered_map<std::string, std::vector<repeat_t>> reps;
    std::unordered_map<std::string, std::vector<Interval<repeat_t>>> reps_iv;
    std::unordered_map<std::string, IntervalTree<repeat_t>*> reps_i;
    std::string line;
	while (getline(rep_f, line)) {
		if (line[0] == '#') continue;
		repeat_t r(line);
		reps[r.chr].push_back(r);
        reps_iv[r.chr].push_back(Interval<repeat_t>(r.start, r.end, r));
	}
    
    for (auto it = reps.begin(); it != reps.end(); it++) {
        reps_i[it->first] = new IntervalTree<repeat_t>(reps_iv[it->first]);
    }

    /* == Compare benchmark svs with called svs == */
    std::set<std::string> b_tps, c_tps;

	for (int i = 0; i < benchmark_svs.size(); i++) {
		sv_t bsv = benchmark_svs[i];

		std::vector<repeat_t> reps_containing_bsv;
        if (reps_i[bsv.bp1.chr] != NULL) {
            std::vector<Interval<repeat_t>> intervals_temp = reps_i[bsv.bp1.chr]->findOverlapping(bsv.bp1.pos, bsv.bp2.pos);
            for (auto& iv : intervals_temp) {
                repeat_t rep = iv.value;
                if (rep.contains(bsv)) {
                    reps_containing_bsv.push_back(rep);
                }
            }
        }

        bool matched = false;
        for (int j = 0; j < called_svs_by_chr[bsv.bp1.chr].size(); j++) {
            sv_t csv = called_svs_by_chr[bsv.bp1.chr][j];
            if (dist(bsv, csv) <= maxdist && check_ins_seq(bsv, csv)) {
                if (!report) std::cout << bsv.id << " " << csv.id << std::endl;
                b_tps.insert(bsv.id);
                c_tps.insert(csv.id);
                matched = true;
            }

            for (repeat_t& rep : reps_containing_bsv) {
                if (rep.intersects(csv) && check_ins_seq(bsv, csv)) {
                	if (!report) std::cout << bsv.id << " " << csv.id << " REP" << std::endl;
                    b_tps.insert(bsv.id);
                    c_tps.insert(csv.id);
                    matched = true;
                }
            }
        }

        if (!matched) {
        	if (!report) std::cout << bsv.id << " NONE" << std::endl;
        }
	}

	if (report) {
		std::cout.precision(2);
		std::cout << "RECALL: " << b_tps.size() << "/" << benchmark_svs.size() << " = " << double(b_tps.size())/benchmark_svs.size() << std::endl;
		std::cout << "PRECISION: " << c_tps.size() << "/" << called_svs.size() << " = " << double(c_tps.size())/called_svs.size() << std::endl;
	}

	if (print_fp)
	for (sv_t& sv : called_svs) {
		if (!c_tps.count(sv.id)) {
			std::cerr << sv.id << std::endl;
		}
	}
}
