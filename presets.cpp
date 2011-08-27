/**
 * presets.cpp
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
	// --very-fast:      -M 1 --multiseed=0,22,S,1,2.50,G,1.00,1.00,C,2,0 --score-min=L,-0.6,-0.9
	// --fast:           -M 5 --multiseed=0,22,S,1,2.50,G,1.00,1.00,C,2,0 --score-min=L,-0.6,-0.9
	// --sensitive:      -M 5 --multiseed=0,22,S,1,1.25,G,1.00,1.00,C,2,0 --score-min=L,-0.6,-0.9
	// --very-sensitive: -M 5 --multiseed=0,20,S,1,0.50,G,1.00,1.00,C,2,0 --score-min=L,-0.6,-0.9
	if(preset == "very-fast") {
		policy += ";SEED=0,22";
		policy += ";IVAL=S,1,2.50";
		policy += ";POSF=G,1.00,1.00";
		policy += ";ROWM=C,2,0";
		policy += ";MHITS=1";
	} else if(preset == "fast") {
		policy += ";SEED=0,22";
		policy += ";IVAL=S,1,2.50";
		policy += ";POSF=G,1.00,1.00";
		policy += ";ROWM=C,2,0";
		policy += ";MHITS=5";
	} else if(preset == "sensitive") {
		policy += ";SEED=0,22";
		policy += ";IVAL=S,1,1.25";
		policy += ";POSF=G,1.00,1.00";
		policy += ";ROWM=C,2,0";
		policy += ";MHITS=5";
	} else if(preset == "very-sensitive") {
		policy += ";SEED=0,20";
		policy += ";IVAL=S,1,0.50";
		policy += ";POSF=G,1.00,1.00";
		policy += ";ROWM=C,2,0";
		policy += ";MHITS=5";
	} else if(preset == "very-fast-local") {
		policy += ";SEED=0,25";
		policy += ";IVAL=S,1,2.00";
		policy += ";POSF=G,0.75,1.75";
		policy += ";ROWM=C,3,0";
	} else if(preset == "fast-local") {
		policy += ";SEED=0,22";
		policy += ";IVAL=S,1,1.75";
		policy += ";POSF=G,2.00,2.50";
		policy += ";ROWM=C,4,0";
	} else if(preset == "sensitive-local") {
		policy += ";SEED=0,20";
		policy += ";IVAL=S,0,1.25";
		policy += ";POSF=G,2.00,2.50";
		policy += ";ROWM=C,4,0";
	} else if(preset == "very-sensitive-local") {
		policy += ";SEED=0,20";
		policy += ";IVAL=S,1,0.50";
		policy += ";POSF=G,3.00,2.75";
		policy += ";ROWM=C,7,0";
	} else {
		cerr << "Unknown preset: " << preset << endl;
	}
}
