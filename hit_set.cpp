/*
 * hit_set.cpp
 *
 *  Created on: Jul 31, 2009
 *      Author: Ben Langmead
 */

#include <iostream>
#include "alphabet.h"
#include "hit_set.h"
#include "sstring.h"

using namespace std;

/**
 * Report up to 'khits' hits from this HitSet.
 */
void HitSet::reportUpTo(ostream& os, int khits) {
	khits = min(khits, (int)size());
	BTDnaString seqrc;
	BTString qualr;
	for(int i = 0; i < khits; i++) {
		const HitSetEnt& h = ents[i];
		if(!h.fw && seqrc.empty()) {
			// Lazily initialize seqrc and qualr
			seqrc = seq;
			seqrc.reverseComp(color);
			assert_eq(seqrc.length(), seq.length());
			qualr = qual;
			qualr.reverse();
			assert_eq(qualr.length(), qual.length());
		}
		os << name << '\t'
		   << (h.fw ? '+' : '-') << '\t'
		   << h.h.first << '\t'
		   << h.h.second << '\t'
		   << (h.fw ? seq : seqrc) << '\t'
		   << (h.fw ? qual : qualr) << '\t'
		   << h.oms << '\t';
		for(size_t i = 0; i < h.edits.size(); i++) {
			const Edit& e = h.edits[i];
			os << e.pos;
			if(e.type == EDIT_TYPE_SNP) os << "S";
			os << ":" << (char)e.chr << ">" << (e.qchr != 0 ? (char)e.qchr : (char)seq[e.pos]);
			if(i < h.edits.size()-1 || !h.cedits.empty()) os << ",";
		}
		for(size_t i = 0; i < h.cedits.size(); i++) {
			const Edit& e = h.cedits[i];
			os << e.pos;
			if(e.type == EDIT_TYPE_SNP) os << "S";
			os << ":" << (char)e.chr << ">" << (e.qchr != 0 ? (char)e.qchr : (char)seq[e.pos]);
			if(i < h.cedits.size()-1) os << ",";
		}
		os << endl;
	}
}

ostream& operator << (ostream& os, const HitSetEnt& hs) {
	os << "\t" << hs.h.first << ":" << hs.h.second;
	return os;
}

ostream& operator << (ostream& os, const HitSet& hs) {
	os << hs.name << ":" << hs.seq << ":" << hs.qual << endl;
	for(size_t i = 0; i < hs.ents.size(); i++) {
		os << hs.ents[i];
	}
	return os;
}
