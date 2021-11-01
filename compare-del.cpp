#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <set>
#include <stdio.h>
#include <unistd.h>

#include "htslib/vcf.h"
#include "lib/IntervalTree.h"
#include "lib/ssw_cpp.h"
#include "lib/kseq.h"

KSEQ_INIT(int, read)

#include "common.h"

int dist(sv_t& sv1, sv_t& sv2) {
    if (sv1.bp1.chr != sv2.bp1.chr) return INT32_MAX;
    return abs(sv1.bp1.pos-sv2.bp1.pos) + abs(sv1.bp2.pos-sv2.bp2.pos);
}

int get_overlap(sv_t& sv1, sv_t& sv2) {
    return std::max(0, std::min(sv1.bp2.pos, sv2.bp2.pos)-std::max(sv1.bp1.pos, sv2.bp1.pos));
}

int main(int argc, char* argv[]) {

    if (argc < 6) {
        std::cout << "Given a VCF or SV file with benchmark deletions and one with the called ones, reports for "
                     "each benchmark deletions if it is being called." << std::endl;
        std::cout << "Usage: exec benchmark_file called_file ucsc_repeats ref_genome maxdist [--overlap min_overlap_frac "
        		" --report --force-ids]" << std::endl;
        return 0;
    }

    std::vector<sv_t> benchmark_svs = read_sv_list(argv[1]);
	std::vector<sv_t> called_svs = read_sv_list(argv[2]);

	auto is_del_func = [](const sv_t& sv) {return sv.type != "DEL";};
	benchmark_svs.erase(std::remove_if(benchmark_svs.begin(), benchmark_svs.end(), is_del_func));
	called_svs.erase(std::remove_if(called_svs.begin(), called_svs.end(), is_del_func));

    std::ifstream rep_f(argv[3]);
    bool report = false, print_fp = false;
    int maxdist = std::stoi(argv[5]);

    double overlap_frac = 0.0;
	for (int i = 6; i < argc; i++) {
		std::string argv_str = std::string(argv[i]);
		if (argv_str == "--overlap") {
			overlap_frac = std::stod(argv[i+1]);
		}
		if (argv_str == "--report") {
			report = true;
		}
		if (argv_str == "--print-fp") {
			print_fp = true;
		}
		if (argv_str == "--force-ids") {
			for (int j = 0; j < benchmark_svs.size(); j++) benchmark_svs[j].id = "DEL_" + std::to_string(j);
			for (int j = 0; j < called_svs.size(); j++) called_svs[j].id = "DEL_" + std::to_string(j);
		}
	}

	std::unordered_map<std::string, std::vector<sv_t> > called_svs_by_chr;
	for (sv_t& sv : called_svs) {
		called_svs_by_chr[sv.bp1.chr].push_back(sv);
	}

    std::unordered_map<std::string, std::string> chrs;
    
    FILE* fastaf = fopen(argv[4], "r");
    kseq_t* seq = kseq_init(fileno(fastaf));
    while (kseq_read(seq) >= 0) {
        chrs[std::string(seq->name.s)] = seq->seq.s;
    }
    kseq_destroy(seq);

    std::string line;
    std::unordered_map<std::string, std::vector<repeat_t>> reps;
    std::unordered_map<std::string, std::vector<Interval<repeat_t>>> reps_iv;
    std::unordered_map<std::string, IntervalTree<repeat_t>*> reps_i;
    while (getline(rep_f, line)) {
        if (line[0] == '#') continue;
        repeat_t r(line);
        reps[r.chr].push_back(r);
        reps_iv[r.chr].push_back(Interval<repeat_t>(r.start, r.end, r));
    }
    
    for (auto it = reps.begin(); it != reps.end(); it++) {
        reps_i[it->first] = new IntervalTree<repeat_t>(reps_iv[it->first]);
    }
    
    StripedSmithWaterman::Aligner aligner(2,2,4,1);
    StripedSmithWaterman::Filter filter(false, false, 0, 32767);
    StripedSmithWaterman::Alignment alignment;

    /* == Compare benchmark svs with called svs == */
    std::set<std::string> b_tps, c_tps;

    for (sv_t& bsv : benchmark_svs) {
        std::vector<repeat_t> reps_containing_bsv;
        if (reps_i[bsv.bp1.chr] != NULL) {
            std::vector<Interval<repeat_t>> intervals_temp = reps_i[bsv.bp1.chr]->findOverlapping(bsv.bp1.pos, bsv.bp2.pos);
            for (auto &iv : intervals_temp) {
                repeat_t rep = iv.value;
                if (rep.contains(bsv)) {
                    reps_containing_bsv.push_back(rep);
                }
            }
        }

        bool matched = false;
        for (sv_t& csv : called_svs_by_chr[bsv.bp1.chr]) {
        	bool good_overlap = get_overlap(bsv,csv) >= std::min(bsv.size(), csv.size()) * overlap_frac;
            if (dist(bsv, csv) <= maxdist && good_overlap) {
                std::cout << bsv.id << " " << csv.id << std::endl;
                b_tps.insert(bsv.id);
                c_tps.insert(csv.id);
                matched = true;
            }

            for (repeat_t& rep : reps_containing_bsv) {
                if (rep.intersects(csv) && good_overlap) {
                    std::string bseq = chrs[bsv.bp1.chr].substr(bsv.bp1.pos, bsv.bp2.pos-bsv.bp1.pos);
                    std::string cseq = chrs[csv.bp1.chr].substr(csv.bp1.pos, csv.bp2.pos-csv.bp1.pos);
                    if (bseq.length() > cseq.length()) {
                        bseq = bseq + bseq;
                    } else {
                        cseq = cseq + cseq;
                    }

                    if (cseq.length() < 50000 && bseq.length() < 50000) {
                        aligner.Align(bseq.data(), cseq.data(), cseq.length(), filter, &alignment);
                        if (alignment.sw_score >= std::min(bseq.length(), cseq.length())) {
                            std::cout << bsv.id << " " << csv.id << " REP" << std::endl;
                            b_tps.insert(bsv.id);
							c_tps.insert(csv.id);
                            matched = true;
                        }
                    }
                }
            }
        }

        if (!matched) {
            std::cout << bsv.id << " NONE" << std::endl;
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
