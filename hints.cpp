/*
 * Copyright 2016, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "hints.h"

/** Return true iff read name indicates hints are present */
bool has_hint(const Read& r) {
	return r.name[0] == '!' &&
	       r.name[1] == 'h' &&
	       r.name[2] == '!';
}

/**
 * Parse hints out of the read name and into the given list of SeedHits.
 */
void parse_hints(const Read& r,
                 EList<SeedHit>& hints,
                 const EMap<std::string, TRefId>& refidMap)
{
	const char *origbuf = r.name.buf();
	const char *buf = origbuf;
	const size_t namelen = r.name.length();
	char refname[1024];
	assert_eq(buf[0], '!');
	assert_eq(buf[1], 'h');
	assert_eq(buf[2], '!');
	buf += 2;
	while(*buf == '!' && buf - origbuf < namelen) {
		SeedHit hint;
		assert_eq(*buf, '!');
		// parse reference name
		char *refbuf = refname;
		buf++;
		while(*buf != '!') {
			*refbuf++ = *buf++;
		}
		*refbuf = '\0';
		assert_eq(*buf, '!');
		buf++;
		TRefId refid = -1;
		bool found = refidMap.find(refname, refid);
		if(!found) {
			std::cerr << "Hint parsing error: Bad reference name: "
			          << refname << endl;
			throw 1;
		}
		// parse reference offset
		TRefOff refoff = 0;
		while(*buf != '!') {
			assert(isdigit((int)*buf));
			refoff *= 10;
			refoff += (*buf++) - '0';
		}
		assert_eq(*buf, '!');
		buf++;
		// parse orientation
		bool fw = (*buf++) == '+';
		assert_eq(*buf, '!');
		buf++;
		// parse length of seed hit
		size_t len = 0;
		while(*buf != '!') {
			assert(isdigit((int)*buf));
			len *= 10;
			len += (*buf++) - '0';
		}
		hint.refival.init(refid, refoff, fw, len);
		buf++;
		// parse read 5' off of seed hit
		size_t fivepoff = 0;
		while(isdigit(*buf) && buf - origbuf < namelen) {
			fivepoff *= 10;
			fivepoff += (*buf++) - '0';
		}
		hint.rd5primeOff = fivepoff;
		hints.push_back(hint);
	}
}

/**
 * Parse hints out of the read name and into the given list of SeedHits.
 */
void parse_interval_hints(const Read& r,
                          EList<IntervalHit>& hints,
                          const EMap<std::string, TRefId>& refidMap)
{
	const char *origbuf = r.name.buf();
	const char *buf = origbuf;
	const size_t namelen = r.name.length();
	char refname[1024];
	assert_eq(buf[0], '!');
	assert_eq(buf[1], 'h');
	assert_eq(buf[2], '!');
	buf += 2;
	while(*buf == '!' && buf - origbuf < namelen) {
		IntervalHit hint;
		assert_eq(*buf, '!');
		// parse reference name
		char *refbuf = refname;
		buf++;
		while(*buf != '!') {
			*refbuf++ = *buf++;
		}
		*refbuf = '\0';
		assert_eq(*buf, '!');
		buf++;
		TRefId refid = -1;
		bool found = refidMap.find(refname, refid);
		if(!found) {
			std::cerr << "Hint parsing error: Bad reference name: "
			          << refname << endl;
			throw 1;
		}
		// parse left reference offset
		TRefOff refoff_l = 0;
		while(*buf != '!') {
			assert(isdigit((int)*buf));
			refoff_l *= 10;
			refoff_l += (*buf++) - '0';
		}
		assert_eq(*buf, '!');
		buf++;
		// parse right reference offset
		TRefOff refoff_r = 0;
		while(*buf != '!') {
			assert(isdigit((int)*buf));
			refoff_r *= 10;
			refoff_r += (*buf++) - '0';
		}
		assert_eq(*buf, '!');
		assert_gt(refoff_r, refoff_l);
		buf++;
		// parse length of seed hit
		size_t len = 0;
		while(*buf != '!') {
			assert(isdigit((int)*buf));
			len *= 10;
			len += (*buf++) - '0';
		}
		// note: orientation is always true, but that is ignored
		hint.refival.init(refid, refoff_l, true, refoff_r - refoff_l + 1);
		buf++;
		// parse read 5' off of seed hit
		size_t fivepoff = 0;
		while(isdigit(*buf) && buf - origbuf < namelen) {
			fivepoff *= 10;
			fivepoff += (*buf++) - '0';
		}
		hint.rd5primeOff = fivepoff;
		hint.hitlen = len;
		hints.push_back(hint);
	}
}
