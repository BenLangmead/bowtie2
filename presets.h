/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
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

/**
 * presets.h
 *
 * Maps simple command-line options to more complicated combinations of
 * options for ease-of-use.
 */

#ifndef PRESETS_H_
#define PRESETS_H_

#include <string>
#include <utility>
#include "ds.h"

class Presets {
public:
	
	Presets() { }
	
	virtual ~Presets() { }
	
	virtual void apply(
		const std::string& preset,
		std::string& policy,
		EList<std::pair<int, std::string> >& opts) = 0;
	
	virtual const char * name() = 0;
};

/**
 * Initial collection of presets: 8/14/2011 prior to first Bowtie 2 release.
 */
class PresetsV0 : public Presets {
public:
	
	PresetsV0() : Presets() { }
	
	virtual ~PresetsV0() { }
	
	virtual void apply(
		const std::string& preset,
		std::string& policy,
		EList<std::pair<int, std::string> >& opts);

	virtual const char * name() { return "V0"; }
};

#endif /*ndef PRESETS_H_*/
