//
//  read_sink.h
//

#ifndef READ_SINK_H_
#define READ_SINK_H_

#include <string>
#include "read.h"

/**
 * Writes reads and read pairs to one of a set of files (or file pairs)
 * depending on whether it aligned, failed to align, or was found to be
 * repetitive.
 */
class ReadSink {

public:

	explicit ReadSink(
		const std::string& dumpAl,   // filename to dump aligned reads to
		const std::string& dumpUnal, // filename to dump unaligned reads to
		const std::string& dumpMax,  // filename to dump repetitive reads to
		bool onePairFile) :          // true -> write both mates to same file
		dumpAlBase_(dumpAl),
		dumpUnalBase_(dumpUnal),
		dumpMaxBase_(dumpMax)
	{
		initDumps();
	}
	
	/**
	 * Destructor: close all files and destroy stream objects.
	 */
	~ReadSink() { destroyDumps(); }

	/**
	 * Dump an aligned read to all of the appropriate output streams.
	 * Be careful to synchronize correctly - there may be multiple
	 * simultaneous writers.
	 */
	void dumpAlign(
		const Read* m1,
		const Read* m2,
		TReadId rdid);

	/**
	 * Dump an unaligned read to all of the appropriate output streams.
	 * Be careful to synchronize correctly - there may be multiple
	 * simultaneous writers.
	 */
	void dumpUnal(
		const Read* m1,
		const Read* m2,
		TReadId rdid);

	/**
	 * Dump a maxed-out read to all of the appropriate output streams.
	 * Be careful to synchronize correctly - there may be multiple
	 * simultaneous writers.
	 */
	void dumpMaxed(
		const Read* m1,
		const Read* m2,
		TReadId rdid);

protected:

	/**
	 * Initialize all the locks for dumping.
	 */
	void initDumps();

	/**
	 * Close all file handles and destroy dynamically allocated stream
	 * objects.
	 */
	void destroyDumps();

	/**
	 * Open an ofstream with given name; output error message and quit
	 * if it fails.
	 */
	std::ofstream* openOf(
		const std::string& name,
		int mateType,
		const std::string& suffix);

	// Output filenames for dumping
	std::string dumpAlBase_;
	std::string dumpUnalBase_;
	std::string dumpMaxBase_;

	// Output streams for dumping
	std::ofstream *dumpAl_;   // for single-ended reads
	std::ofstream *dumpAl_1_; // for first mates
	std::ofstream *dumpAl_2_; // for second mates
	std::ofstream *dumpUnal_;   // for single-ended reads
	std::ofstream *dumpUnal_1_; // for first mates
	std::ofstream *dumpUnal_2_; // for second mates
	std::ofstream *dumpMax_;     // for single-ended reads
	std::ofstream *dumpMax_1_;   // for first mates
	std::ofstream *dumpMax_2_;   // for second mates

	// Locks for dumping
	MUTEX_T dumpAlignLock_;
	MUTEX_T dumpAlignLockPE_; // _1 and _2
	MUTEX_T dumpUnalLock_;
	MUTEX_T dumpUnalLockPE_; // _1 and _2
	MUTEX_T dumpMaxLock_;
	MUTEX_T dumpMaxLockPE_;   // _1 and _2

	// false -> no dumping
	bool dumpAlignFlag_;
	bool dumpUnalignFlag_;
	bool dumpMaxedFlag_;
	
	bool onePairFile_;

};

#endif /*ndef READ_SINK_H_*/
