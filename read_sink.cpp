//
//  read_sink.cpp
//

#include "read_sink.h"

/**
 * Dump an aligned read to all of the appropriate output streams.
 * Be careful to synchronize correctly - there may be multiple
 * simultaneous writers.
 */
void ReadSink::dumpAlign(
	const Read* m1,
	const Read* m2,
	TReadId rdid)
{
	if(!dumpAlignFlag_) return;
	bool paired = (m2 != NULL && !m2->empty());
	if(!paired || onePairFile_) {
		// Dump unpaired read to an aligned-read file of the same format
		if(!dumpAlBase_.empty()) {
			ThreadSafe ts(&dumpAlignLock_);
			if(dumpAl_ == NULL) {
				dumpAl_ = openOf(dumpAlBase_, 0, "");
				assert(dumpAl_ != NULL);
			}
			dumpAl_->write(m1->readOrigBuf.buf(), m1->readOrigBuf.length());
		}
	} else {
		// Dump paired-end read to an aligned-read file (or pair of
		// files) of the same format
		if(!dumpAlBase_.empty()) {
			ThreadSafe ts(&dumpAlignLockPE_);
			if(dumpAl_1_ == NULL) {
				dumpAl_1_ = openOf(dumpAlBase_, 1, "");
				dumpAl_2_ = openOf(dumpAlBase_, 2, "");
				assert(dumpAl_1_ != NULL);
				assert(dumpAl_2_ != NULL);
			}
			dumpAl_1_->write(m1->readOrigBuf.buf(), m1->readOrigBuf.length());
			dumpAl_2_->write(m2->readOrigBuf.buf(), m2->readOrigBuf.length());
		}
	}
}

/**
 * Dump an unaligned read to all of the appropriate output streams.
 * Be careful to synchronize correctly - there may be multiple
 * simultaneous writers.
 */
void ReadSink::dumpUnal(
	const Read* m1,
	const Read* m2,
	TReadId rdid)
{
	if(!dumpUnalignFlag_) return;
	bool paired = (m2 != NULL && !m2->empty());
	if(!paired || onePairFile_) {
		// Dump unpaired read to an unaligned-read file of the same format
		if(!dumpUnalBase_.empty()) {
			ThreadSafe ts(&dumpUnalLock_);
			if(dumpUnal_ == NULL) {
				dumpUnal_ = openOf(dumpUnalBase_, 0, "");
				assert(dumpUnal_ != NULL);
			}
			dumpUnal_->write(m1->readOrigBuf.buf(), m1->readOrigBuf.length());
		}
	} else {
		// Dump paired-end read to an unaligned-read file (or pair
		// of files) of the same format
		if(!dumpUnalBase_.empty()) {
			ThreadSafe ts(&dumpUnalLockPE_);
			if(dumpUnal_1_ == NULL) {
				dumpUnal_1_ = openOf(dumpUnalBase_, 1, "");
				dumpUnal_2_ = openOf(dumpUnalBase_, 2, "");
				assert(dumpUnal_1_ != NULL);
				assert(dumpUnal_2_ != NULL);
			}
			dumpUnal_1_->write(m1->readOrigBuf.buf(), m1->readOrigBuf.length());
			dumpUnal_2_->write(m2->readOrigBuf.buf(), m2->readOrigBuf.length());
		}
	}
}

/**
 * Dump a maxed-out read to all of the appropriate output streams.
 * Be careful to synchronize correctly - there may be multiple
 * simultaneous writers.
 */
void ReadSink::dumpMaxed(
	const Read* m1,
	const Read* m2,
	TReadId rdid)
{
	if(!dumpMaxedFlag_) {
		if(dumpUnalignFlag_) dumpUnal(m1, m2, rdid);
		return;
	}
	bool paired = (m2 != NULL && !m2->empty());
	if(paired || onePairFile_) {
		// Dump unpaired read to an maxed-out-read file of the same format
		if(!dumpMaxBase_.empty()) {
			ThreadSafe ts(&dumpMaxLock_);
			if(dumpMax_ == NULL) {
				dumpMax_ = openOf(dumpMaxBase_, 0, "");
				assert(dumpMax_ != NULL);
			}
			dumpMax_->write(m1->readOrigBuf.buf(), m1->readOrigBuf.length());
		}
	} else {
		// Dump paired-end read to a maxed-out-read file (or pair
		// of files) of the same format
		if(!dumpMaxBase_.empty()) {
			ThreadSafe ts(&dumpMaxLockPE_);
			if(dumpMax_1_ == NULL) {
				dumpMax_1_ = openOf(dumpMaxBase_, 1, "");
				dumpMax_2_ = openOf(dumpMaxBase_, 2, "");
				assert(dumpMax_1_ != NULL);
				assert(dumpMax_2_ != NULL);
			}
			dumpMax_1_->write(m1->readOrigBuf.buf(), m1->readOrigBuf.length());
			dumpMax_2_->write(m2->readOrigBuf.buf(), m2->readOrigBuf.length());
		}
	}
}

/**
 * Initialize all the locks for dumping.
 */
void ReadSink::initDumps() {
	dumpAl_       = dumpAl_1_     = dumpAl_2_     = NULL;
	dumpUnal_     = dumpUnal_1_   = dumpUnal_2_   = NULL;
	dumpMax_      = dumpMax_1_    = dumpMax_2_    = NULL;
	dumpAlignFlag_   = !dumpAlBase_.empty();
	dumpUnalignFlag_ = !dumpUnalBase_.empty();
	dumpMaxedFlag_   = !dumpMaxBase_.empty();
	MUTEX_INIT(dumpAlignLock_);
	MUTEX_INIT(dumpAlignLockPE_);
	MUTEX_INIT(dumpUnalLock_);
	MUTEX_INIT(dumpUnalLockPE_);
	MUTEX_INIT(dumpMaxLock_);
	MUTEX_INIT(dumpMaxLockPE_);
}

/**
 * Close all file handles and destroy dynamically allocated stream
 * objects.
 */
void ReadSink::destroyDumps() {
	if(dumpAl_       != NULL) { dumpAl_->close();       delete dumpAl_; }
	if(dumpAl_1_     != NULL) { dumpAl_1_->close();     delete dumpAl_1_; }
	if(dumpAl_2_     != NULL) { dumpAl_2_->close();     delete dumpAl_2_; }

	if(dumpUnal_     != NULL) { dumpUnal_->close();     delete dumpUnal_; }
	if(dumpUnal_1_   != NULL) { dumpUnal_1_->close();   delete dumpUnal_1_; }
	if(dumpUnal_2_   != NULL) { dumpUnal_2_->close();   delete dumpUnal_2_; }

	if(dumpMax_      != NULL) { dumpMax_->close();      delete dumpMax_; }
	if(dumpMax_1_    != NULL) { dumpMax_1_->close();    delete dumpMax_1_; }
	if(dumpMax_2_    != NULL) { dumpMax_2_->close();    delete dumpMax_2_; }
}

/**
 * Open an ofstream with given name; output error message and quit
 * if it fails.
 */
std::ofstream* ReadSink::openOf(
	const std::string& name,
	int mateType,
	const std::string& suffix)
{
	std::string s = name;
	size_t dotoff = name.find_last_of(".");
	if(mateType == 1) {
		if(dotoff == string::npos) {
			s += "_1"; s += suffix;
		} else {
			s = name.substr(0, dotoff) + "_1" + s.substr(dotoff);
		}
	} else if(mateType == 2) {
		if(dotoff == string::npos) {
			s += "_2"; s += suffix;
		} else {
			s = name.substr(0, dotoff) + "_2" + s.substr(dotoff);
		}
	} else if(mateType != 0) {
		cerr << "Bad mate type " << mateType << endl; throw 1;
	}
	std::ofstream* tmp = new ofstream(s.c_str(), ios::out);
	if(tmp->fail()) {
		if(mateType == 0) {
			cerr << "Could not open single-ended aligned/unaligned-read file for writing: " << name << endl;
		} else {
			cerr << "Could not open paired-end aligned/unaligned-read file for writing: " << name << endl;
		}
		throw 1;
	}
	return tmp;
}
