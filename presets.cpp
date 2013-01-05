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

#include <iostream>
#include "presets.h"
#include "opts.h"

using namespace std;

void PresetsV0::apply(
	const std::string& preset,
	std::string& policy,
	EList<std::pair<int, std::string> >& opts)
{
	// Presets:                 Same as:
	//  For --end-to-end:
	//   --very-fast            -M 5 -R 1 -N 0 -L 22 -i S,1,2.50
	//   --fast                 -M 10 -R 2 -N 0 -L 22 -i S,1,2.50
	//   --sensitive            -M 15 -R 2 -N 0 -L 22 -i S,1,1.15
	//   --very-sensitive       -M 25 -R 3 -N 0 -L 19 -i S,1,0.50
	if(preset == "very-fast") {
		policy += ";SEED=0,22";
		policy += ";DPS=5";
		policy += ";ROUNDS=1";
		policy += ";IVAL=S,0,2.50";
	} else if(preset == "fast") {
		policy += ";SEED=0,22";
		policy += ";DPS=10";
		policy += ";ROUNDS=2";
		policy += ";IVAL=S,0,2.50";
	} else if(preset == "sensitive") {
		policy += ";SEED=0,22";
		policy += ";DPS=15";
		policy += ";ROUNDS=2";
		policy += ";IVAL=S,1,1.15";
	} else if(preset == "very-sensitive") {
		policy += ";SEED=0,20";
		policy += ";DPS=20";
		policy += ";ROUNDS=3";
		policy += ";IVAL=S,1,0.50";
	}
	//  For --local:
	//   --very-fast-local      -M 1 -N 0 -L 25 -i S,1,2.00
	//   --fast-local           -M 2 -N 0 -L 22 -i S,1,1.75
	//   --sensitive-local      -M 2 -N 0 -L 20 -i S,1,0.75 (default)
	//   --very-sensitive-local -M 3 -N 0 -L 20 -i S,1,0.50
	else if(preset == "very-fast-local") {
		policy += ";SEED=0,25";
		policy += ";DPS=5";
		policy += ";ROUNDS=1";
		policy += ";IVAL=S,1,2.00";
	} else if(preset == "fast-local") {
		policy += ";SEED=0,22";
		policy += ";DPS=10";
		policy += ";ROUNDS=2";
		policy += ";IVAL=S,1,1.75";
	} else if(preset == "sensitive-local") {
		policy += ";SEED=0,20";
		policy += ";DPS=15";
		policy += ";ROUNDS=2";
		policy += ";IVAL=S,1,0.75";
	} else if(preset == "very-sensitive-local") {
		policy += ";SEED=0,20";
		policy += ";DPS=20";
		policy += ";ROUNDS=3";
		policy += ";IVAL=S,1,0.50";
	}
	else {
		cerr << "Unknown preset: " << preset.c_str() << endl;
	}
}
