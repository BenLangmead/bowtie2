/*
 * refmap.h
 *
 *  Created on: Aug 3, 2009
 *      Author: Ben Langmead
 */

#ifndef REFMAP_H_
#define REFMAP_H_

#include <stdint.h>
#include <cassert>
#include <iostream>
#include <fstream>
#include "ds.h"

class ReferenceMap {
	typedef std::pair<uint32_t, uint32_t> U32Pair;

public:
	ReferenceMap(const char *fname, bool parseNames) {
		fname_ = fname;
		parseNames_ = parseNames;
		parse();
	}

	/**
	 * Give a reference coordinate in the index, translate it into a
	 * new reference coordinate via the reference map supplied by the
	 * user.
	 */
	void map(U32Pair& h) const;

	/**
	 * Return true iff we have a name for reference with id 'i'.
	 */
	bool hasName(size_t i) const {
		if(!parseNames_) return false;
		return !names_[i].empty();
	}

	/**
	 * Get the name for reference with id 'i'.
	 */
	const std::string& getName(size_t i) const {
		assert(parseNames_);
		assert(hasName(i));
		return names_[i];
	}

protected:

	/**
	 * Parse a reference-map file.
	 */
	void parse();

	const char *fname_;
	EList<U32Pair> map_;
	bool parseNames_;
	EList<std::string> names_;
};

#endif /* REFMAP_H_ */
