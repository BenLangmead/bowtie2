/**
 * presets.h
 *
 * Maps simple command-line options to more complicated combinations of
 * options for ease-of-use.
 */

#ifndef PRESETS_H_
#define PRESETS_H_

#include <string>

class Presets {
public:
	
	Presets() { }
	
	virtual ~Presets() { }
	
	virtual void apply(
		const std::string& preset,
		std::string& policy,
		std::string& args) = 0;
	
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
		std::string& args);

	virtual const char * name() { return "V0"; }
};

#endif /*ndef PRESETS_H_*/
