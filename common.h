#ifndef SVCOMPARE_COMMON_H
#define SVCOMPARE_COMMON_H

#include "htslib/vcf.h"

std::string get_sv_type(bcf_hdr_t* hdr, bcf1_t* sv) {
    char* data = NULL;
    int len = 0;
    if (bcf_get_info_string(hdr, sv, "SVTYPE", &data, &len) < 0) {
        throw std::runtime_error("Failed to determine SVTYPE for sv " + std::string(sv->d.id));
    }
    std::string svtype = data;
    delete[] data;
    return svtype;
}

int get_sv_end(bcf_hdr_t* hdr, bcf1_t* sv) {
    int* data = NULL;
    int size = 0;
    bcf_get_info_int32(hdr, sv, "END", &data, &size);
    if (size > 0) {
        int end = data[0];
        delete[] data;
        return end-1; // return 0-based
    }

    bcf_get_info_int32(hdr, sv, "SVLEN", &data, &size);
    if (size > 0) {
        int svlen = data[0];
        delete[] data;
        return sv->pos + abs(svlen);
    }

    std::cerr << "Warning: SV (ID=" << sv->d.id << ") has no END or SVLEN annotation. Assuming END == POS." << std::endl;
    return sv->pos;
}

std::string get_ins_seq(bcf_hdr_t* hdr, bcf1_t* sv) {
	// priority to the ALT allele, if it is not symbolic and longer than just the padding base
	char c = toupper(sv->d.allele[1][0]);
	if ((c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'N') && strlen(sv->d.allele[1]) > 1) {
		return sv->d.allele[1];
	}

	// otherwise, look for SVINSSEQ (compliant with Manta)
	char* data = NULL;
	int size = 0;
	bcf_get_info_string(hdr, sv, "SVINSSEQ", (void**) &data, &size);
	if (data) return data;

	return "";
}


struct bp_t {
    std::string chr;
    int pos;
    char strand;
};

struct sv_t {
    std::string id;
    bp_t bp1, bp2;
    std::string type, seq;

    sv_t(std::string& line) {
        std::stringstream ss(line);
        ss >> id >> bp1.chr >> bp1.pos >> bp1.strand >> bp2.chr >> bp2.pos >> bp2.strand >> type >> seq;
    }
    sv_t(bcf_hdr_t* hdr, bcf1_t* line) {
    	bcf_unpack(line, BCF_UN_STR);
    	id = line->d.id;
    	std::string chr = bcf_hdr_id2name(hdr, line->rid);
    	bp1.chr = bp2.chr = chr;
    	bp1.pos = line->pos;
    	bp2.pos = get_sv_end(hdr, line);
    	type = get_sv_type(hdr, line);
    	seq = get_ins_seq(hdr, line);
    }

    int size() {
        if (type == "INS") return seq.size();
        else if (type == "TRA") return 0;
        else return bp2.pos - bp1.pos;
    }

    bool has_seq() {
        return !seq.empty() && seq != "NA";
    }
};

bool ends_with(const char* str, const char* suffix) {
	int l_string = strlen(str);
	int l_suffix = strlen(suffix);
	return strcmp(str+(l_string-l_suffix), suffix) == 0;
}
std::vector<sv_t> read_sv_list(const char* filename) {
	std::vector<sv_t> svs;
	if (ends_with(filename, ".vcf.gz") || ends_with(filename, ".vcf") || ends_with(filename, ".bcf")) { // vcf format
		htsFile* file = bcf_open(filename, "r");
		bcf_hdr_t* hdr = bcf_hdr_read(file);
		bcf1_t* line = bcf_init();
		while (bcf_read(file, hdr, line) == 0) {
			svs.push_back(sv_t(hdr, line));
		}
		bcf_destroy(line);
		bcf_hdr_destroy(hdr);
		bcf_close(file);
	} else {
		std::ifstream file(filename);
		std::string line;
		while (getline(file, line)) {
			sv_t sv(line);
			svs.push_back(sv);
		}
	}
	return svs;
}

struct repeat_t {
    std::string chr;
    int start, end;

    repeat_t() {}
    repeat_t(std::string& line) {
        std::stringstream ss(line);
        ss >> chr >> start >> end;
    }

    bool contains(sv_t& sv) {
        return sv.bp1.chr == chr && start <= sv.bp1.pos && end >= sv.bp2.pos;
    }
    bool intersects(sv_t& sv) {
        return sv.bp1.chr == chr && sv.bp1.pos <= end && sv.bp2.pos >= start;
    }
};
bool operator == (const repeat_t& r1, const repeat_t& r2) {
    return r1.chr == r2.chr && r1.start == r2.start && r1.end == r2.end;
}

#endif //SVCOMPARE_COMMON_H
