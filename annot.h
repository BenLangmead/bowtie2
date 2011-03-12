/*
 * annot.h
 *
 *  Created on: Aug 3, 2009
 *      Author: Ben Langmead
 */

#ifndef ANNOT_H_
#define ANNOT_H_

#include <stdint.h>
#include <map>
#include <iostream>
#include <fstream>

/**
 * Encapsulates a sorted list of reference positions that are annotated
 * somehow (e.g. as a SNP).
 */
class AnnotationMap {
public:
	typedef std::pair<uint32_t, uint32_t> U32Pair;
	typedef std::pair<char, char> CharPair;
	typedef std::map<U32Pair, CharPair> AnnotMap;
	typedef std::map<U32Pair, CharPair>::const_iterator Iter;

	AnnotationMap(const char *fname) {
		fname_ = fname;
		parse();
	}

	/**
	 * Give a reference coordinate in the index, translate it into a
	 * new reference coordinate via the reference map supplied by the
	 * user.
	 */
	Iter lower_bound(const U32Pair& h) const {
		return map_.lower_bound(h);
	}

	Iter begin() const {
		return map_.begin();
	}

	Iter end() const {
		return map_.end();
	}

protected:

	/**
	 * Parse an annotation-map file.
	 */
	void parse();

	/// filename of file containing the annotation map
	const char *fname_;
	/// maps reference positions to character annotations
	AnnotMap map_;
};

#endif /* ANNOT_H_ */
