 /*==========================================================================
                SeqAn - The Library for Sequence Analysis
                          http://www.seqan.de 
 ============================================================================
  Copyright (C) 2007

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public
  License as published by the Free Software Foundation; either
  version 3 of the License, or (at your option) any later version.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  Lesser General Public License for more details.

 ============================================================================
  $Id: basic_counted_ptr.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

// THIS FILE IS CURRENTLY NOT USED IN SEQAN

#ifndef SEQAN_HEADER_BASIC_COUNTED_PTR_H
#define SEQAN_HEADER_BASIC_COUNTED_PTR_H

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// counted pointer

	template < typename Type >
	struct CountedPtr
	{
		typedef CountedPtr		_Self;
		typedef CountedPtr*	    _SelfPtr;
		typedef CountedPtr&	    _SelfRef;

		typedef Type&			reference;
		typedef const Type&		const_reference;
		typedef Type*			pointer;

        explicit CountedPtr(pointer p = 0): // allocate a new counter
            itsCounter(0)
        {
            if (p) itsCounter = new counter(p);
        }

        CountedPtr(const _Self& r) throw() {
            acquire(r.itsCounter);
        }

        ~CountedPtr() {
            release();
        }

        CountedPtr& operator=(const _Self& r)
        {
            if (this != &r) {
                release();
                acquire(r.itsCounter);
            }
            return *this;
        }

        reference operator*() const throw() {
            return *itsCounter->ptr;
        }

        pointer operator->() const throw() {
            return itsCounter->ptr;
        }

        pointer get() const throw() {
            return itsCounter ? itsCounter->ptr : 0;
        }

        bool unique() const throw() {
            return (itsCounter ? itsCounter->count == 1 : true);
        }

		inline operator pointer () const {
            return get();
		}

    private:

        struct counter {
            pointer     ptr;
            unsigned    count;
            counter(pointer p = 0, unsigned c = 1):
                ptr(p),
                count(c) { }
        }* itsCounter;

        void acquire(counter* c) throw()
        { // increment the count
            itsCounter = c;
            if (c) ++c->count;
        }

        void release()
        { // decrement the count, delete if it is 0
            if (itsCounter) {
                if (--itsCounter->count == 0) {
                    delete itsCounter->ptr;
                    delete itsCounter;
                }
                itsCounter = 0;
            }
        }
    };

}

#endif
