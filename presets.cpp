/**
 * presets.cpp
 */

#include <iostream>
#include "presets.h"

using namespace std;

void PresetsV0::apply(
	const std::string& preset,
	std::string& policy,
	std::string& args)
{
	if(preset == "very-fast") {
		policy += ";SEED=0,31";
		policy += ";IVAL=S,1,3.25";
		policy += ";POSF=G,0.50,0.35";
		policy += ";ROWM=C,1,0";
	} else if(preset == "fast") {
		policy += ";SEED=0,25";
		policy += ";IVAL=S,0,2.50";
		policy += ";POSF=G,0.50,0.50";
		policy += ";ROWM=C,1,0";
	} else if(preset == "sensitive") {
		policy += ";SEED=0,22";
		policy += ";IVAL=S,1,1.75";
		policy += ";POSF=G,2.00,2.50";
		policy += ";ROWM=C,3,0";
	} else if(preset == "very-sensitive") {
		policy += ";SEED=0,20";
		policy += ";IVAL=S,0,1.00";
		policy += ";POSF=G,2.00,2.75";
		policy += ";ROWM=C,5,0";
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
