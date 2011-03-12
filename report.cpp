/*
 *  report.cpp
 *  bowtie-beta1
 *
 *  Created by Benjamin Langmead on 5/26/10.
 *  Copyright 2010 Johns Hopkins University. All rights reserved.
 *
 */

#include <utility>
#include <iostream>
#include "sstring.h"
#include "edit.h"
#include "ebwt.h"
#include "search_globals.h"

using namespace std;

/**
 * Report a hit.  Returns true iff caller can call off the search.
 */
void EbwtSearchParams::reportHitImpl(
	const BTDnaString& seq, // read sequence
	BTString* quals, // read quality values
	BTString* name,  // read name
	const BitPairReference* ref, // reference (= NULL if not necessary)
	const ReferenceMap* rmap, // map to another reference coordinate system
	ColorspaceDecoder * dec, // colorspace decoder
	RandomSource& rand,
	bool fw,
	bool ebwtFw,         // whether index is forward (true) or mirror (false)
	const EList<Edit>& edits, // edits
	U32Pair h,          // ref coords
	U32Pair a,          // arrow pair
	uint32_t tlen,      // length of text
	uint32_t qlen,      // length of query
	int stratum,        // alignment stratum
	uint16_t cost,      // cost of alignment
	uint32_t oms,       // approx. # other valid alignments
	uint32_t patid,
	uint32_t seed,
	uint8_t mate,
	Hit& hit)
{
	// If ebwtFw is true, then 'seq' and 'quals' are reversed
	// If fw is false, then 'seq' and 'quals' are reverse complemented
	assert(!gColor || ref != NULL);
	assert(quals != NULL);
	assert(name != NULL);
	assert_gt(qlen, 0);
	hit.stratum = stratum;
	hit.cost = cost;
	hit.seq = seq;
	hit.quals = *quals;
	hit.readGaps = 0;
	hit.refGaps = 0;
	// Count gaps
	for(size_t i = 0; i < edits.size(); i++) {
		if     (edits[i].isInsert()) hit.readGaps++;
		else if(edits[i].isDelete()) hit.refGaps++;
	}
	if(!ebwtFw) {
		// Re-reverse the pattern and the quals back to how they
		// appeared in the read file
		hit.seq.reverse();
		hit.quals.reverse();
	}
#ifndef NDEBUG
	// Sanity check incoming edits
	for(size_t i = 0; i < edits.size(); i++) {
		const Edit& e = edits[i];
		assert(e.isGap() || e.isMismatch());
		assert_neq(e.chr, e.qchr);
		assert_in(e.chr, "ACGT-");
		assert_in(e.qchr, "ACGTN-");
	}
#endif
	hit.cseq.clear();
	hit.cquals.clear();
	hit.cedits.clear();
	hit.ccedits.clear();
	if(gColor) {
		// Keep pre-decoding versions of the sequence and qualities
		hit.cseq = hit.seq;
		hit.cquals = hit.quals;
		hit.cstratum = hit.stratum;
		hit.ccost = hit.cost;
		hit.seq.clear();
		hit.quals.clear();
		assert(ref != NULL);
		// Get reference chars involved in the alignment
		size_t readi = 0;
		size_t readf = hit.cseq.length();
		size_t refi = 0;
		assert_gt((int)readf + hit.readGaps - hit.refGaps + 1, 0);
		size_t reff = readf + hit.readGaps - hit.refGaps + 1;
		{
			uint32_t rfbuf[(1024+16)/4];
			ASSERT_ONLY(char rfbuf2[1024]);
#ifndef NDEBUG
			for(size_t i = 0; i < qlen + hit.readGaps + 1; i++) {
				rfbuf2[i] = ref->getBase(h.first, h.second + i);
			}
#endif
			int offset = ref->getStretch(rfbuf, h.first, h.second, qlen + hit.readGaps + 1);
			char *rfc = (char*)rfbuf + offset;
			for(size_t i = 0; i < qlen + hit.readGaps + 1; i++) {
				// If gaps are allowed, it's possible for the reference
				// portion of the alignment to include an N.  I don't
				// fully understand why we have to check hit.refGaps
				// though.
				assert(hit.readGaps > 0 || hit.refGaps > 0 || (int)rfc[i] < 4);
				assert_eq(rfc[i], rfbuf2[i]);
				rfc[i] = (1 << rfc[i]);
			}
			BTString rf(rfc, qlen + hit.readGaps + 1);
			edits_.clear();
			for(size_t i = 0; i < edits.size(); i++) {
				if(edits[i].isGap()) {
#ifndef NDEBUG
					int fromLeft = (int)edits[i].pos;
					int fromRight = (int)(hit.cseq.length() - edits[i].pos - 1);
					assert_geq(fromLeft, gGapBarrier);
					if(edits[i].isInsert()) {
						assert_geq(fromRight, gGapBarrier-1);
					} else {
						assert_geq(fromRight, gGapBarrier);
					}
#endif
					edits_.push_back(edits[i]);
				}
			}
			if(!fw) Edit::invertPoss(edits_, hit.cseq.length());
			rand.init(seed);
			dec->decode(
				hit.cseq,     // ASCII colors, '0', '1', '2', '3', '.'
				hit.cquals,   // ASCII quals, Phred+33 encoded
				readi,        // offset of first character within 'read' to consider
				readf,        // offset of last char (exclusive) in 'read' to consider
				rf,           // reference sequence, as masks
				refi,         // offset of first character within 'ref' to consider
				reff,         // offset of last char (exclusive) in 'ref' to consider
				gSnpPhred,    // penalty incurred by a SNP
				MAX_SCORE,    // max cost
				hit.readGaps, // # gaps on read side of color-to-color alignment
				hit.refGaps,  // # gaps on ref side of color-to-color alignment
				edits_,       // all the colorspace-to-colorspace edits
				gInsOpen,     // penalty for opening a new gap in the read
				gInsExtend,   // penalty for extending a gap in the read
				gDelOpen,     // penalty for opening a new gap in the reference
				gDelExtend,   // penalty for extending a gap in the reference
				gGapBarrier,  // # of alignment chars on either side that can't have gaps
				gColorExEnds, // true -> trim either end of nucleotide string
				!gNoMaqRound, // true -> use Maq-like rounding
				hit.seq,      // decoded nucleotides appended here
				hit.quals,    // decoded qualities appended here
				hit.edits,    // destination for decoded nucleotide edits
				hit.aedits,   // destination for decoded ambiguous nucleotides
				hit.cedits,   // destination for decoded color miscalls
				hit.ccedits,  // destination for decoded color edits
				rand);
			assert_eq(hit.seq.length(), hit.cseq.length() + (gColorExEnds ? -1 : 1));
			assert_eq(hit.seq.length(), hit.quals.length());
			//assert_geq((int)hit.ccedits.size(), hit.cstratum);
			// rf is used as scratch space in the decode, so we let
			// it fall out of scope here
		}
#ifndef NDEBUG
		// Make sure decoded alignment has same type and number of
		// gaps as colorspace-to-colorspace alignment
		{
			int barrier = gGapBarrier;
			if(gColor && gColorExEnds) {
				barrier--;
			}
			int ins = 0, del = 0;
			for(size_t i = 0; i < hit.edits.size(); i++) {
				int fromLeft = (int)hit.edits[i].pos;
				int fromRight = (int)(hit.seq.length() - hit.edits[i].pos - 1);
				if(hit.edits[i].isGap()) {
					assert_geq(fromLeft, barrier-1);
					if(hit.edits[i].isInsert()) {
						assert_geq(fromRight, barrier-2);
					} else {
						assert_eq((int)hit.edits[i].qchr, hit.seq.toChar(hit.edits[i].pos));
						assert_geq(fromRight, barrier-1);
					}
				} else if(hit.edits[i].isMismatch()) {
					assert_eq((int)hit.edits[i].qchr, hit.seq.toChar(hit.edits[i].pos));
				}
				if(hit.edits[i].isInsert()) ins++;
				else if(hit.edits[i].isDelete()) del++;
			}
			assert_eq(ins, hit.readGaps);
			assert_eq(del, hit.refGaps);
		}
#endif
		if(!fw) {
			// Invert hit.edits and hit.cedits so that low values
			// for 'pos' refer to the 5' end
			Edit::invertPoss(hit.edits, hit.seq.length());
			Edit::invertPoss(hit.cedits, hit.cseq.length());
			Edit::invertPoss(hit.ccedits, hit.cseq.length());
		}
		if(gColorExEnds) {
			// Extreme bases have been removed; that makes the
			// nucleotide alignment one character shorter than the
			// color alignment
			qlen--;
			// It also shifts the alignment's offset up by 1
			h.second++;
		} else {
			// Extreme bases are included; that makes the
			// nucleotide alignment one character longer than the
			// color alignment
			qlen++;
		}
	} else {
		hit.edits = edits;
	}
	// Check the hit against the original text, if it's available
	if(_texts.size() > 0) {
		assert_lt(h.first, _texts.size());
		BTDnaString editRef;
		edits_ = hit.edits;
		if(!fw) {
			Edit::invertPoss(edits_, hit.seq.length());
		}
		Edit::toRef(hit.seq, edits_, editRef);
		assert_leq(h.second + qlen + hit.readGaps - hit.refGaps, length(_texts[h.first]));
		BTDnaString trueRef;
		for(size_t i = 0; i < qlen + hit.readGaps - hit.refGaps; i++) {
			int refc = (int)_texts[h.first][h.second + i];
			assert_neq(4, refc);
			trueRef.append(refc);
		}
		bool error = false;
		if(editRef != trueRef) {
			cerr << endl;
			cerr << "Edited nucleotide read didn't equal nucleotide reference" << endl;
			Edit::printQAlign(cerr, "          ", hit.seq, edits_);
			cerr << endl;
			cerr << " Obs ref: " << editRef << endl;
			cerr << "True ref: " << trueRef << endl;
			cerr << "  FW: " << fw << endl;
			cerr << "  Ebwt FW: " << ebwtFw << endl;
			cerr << "  Edits: ";
			for(size_t i = 0; i < hit.edits.size(); i++) {
				cerr << hit.edits[i];
				if(i < hit.edits.size()-1) cerr << ",";
			}
			if(hit.edits.empty()) cerr << "-";
			cerr << endl;
			error = true;
		}
		if(gColor) {
			// Construct the colorspace reference portion
			trueRef.clear();
			uint32_t off = h.second + (gColorExEnds ? -1 : 0);
			for(size_t i = 0; i < hit.cseq.length() + hit.readGaps - hit.refGaps; i++) {
				int refc1 = (int)_texts[h.first][off + i];
				int refc2 = (int)_texts[h.first][off + i + 1];
				assert_neq(4, refc1);
				assert_neq(4, refc2);
				int c = dinuc2color[refc1][refc2];
				trueRef.append(c);
			}
			// Transform the color alignment using ccedits
			edits_ = hit.ccedits;
			if(!fw) {
				Edit::invertPoss(edits_, hit.cseq.length());
			}
			editRef.clear();
			Edit::toRef(hit.cseq, edits_, editRef);
			assert_leq(h.second + qlen + hit.readGaps - hit.refGaps, length(_texts[h.first]));
			if(editRef != trueRef) {
				cerr << endl;
				cerr << "Edited color read didn't equal color reference" << endl;
				Edit::printQAlign(cerr, "          ", hit.cseq, edits_);
				cerr << endl;
				cerr << " Obs ref: " << editRef << endl;
				cerr << "True ref: " << trueRef << endl;
				cerr << "  FW: " << fw << endl;
				cerr << "  Ebwt FW: " << ebwtFw << endl;
				cerr << "  Edits: ";
				for(size_t i = 0; i < hit.ccedits.size(); i++) {
					cerr << hit.ccedits[i];
					if(i < hit.ccedits.size()-1) cerr << ",";
				}
				if(hit.ccedits.empty()) cerr << "-";
				cerr << endl;
				error = true;
			}
		}
		assert(!error);
	}
	hit.h = h;
	if(rmap != NULL) rmap->map(hit.h);
	hit.patid = ((patid == 0xffffffff) ? _patid : patid);
	hit.name = *name;
	hit.fw = fw;
	hit.oms = oms;
	hit.mate = mate;
	hit.color = gColor;
	hit.seed = seed;
	hit.pmate = NULL;
}
