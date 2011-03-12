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
  $Id: basic_volatile_ptr.h,v 1.1 2008/08/25 16:20:01 langmead Exp $
 ==========================================================================*/

#ifndef SEQAN_HEADER_BASIC_VOLATILE_PTR_H
#define SEQAN_HEADER_BASIC_VOLATILE_PTR_H

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

	//////////////////////////////////////////////////////////////////////////////
	// volatile pointer
	// allows you to handle volatile data (used by ext. string during swapping)
    //
	// imagine volatile pointers as nodes in an undirected graph
    // when you assign one to another then they are connected
    // all pointers in a connection component points to the same value
    // by calling nukeCopies you can destroy the component and set all pointers to NULL

	template < typename Type >
	struct VolatilePtr
	{
		typedef VolatilePtr		_Self;
		typedef VolatilePtr*	_SelfPtr;
		typedef VolatilePtr&	_SelfRef;

		typedef Type&			reference;
		typedef const Type&		const_reference;
		typedef Type*			pointer;

		pointer			ptr;
		_SelfPtr		next;			// prev == NULL means this is the master node
		_SelfPtr		prev;			// prev == NULL means this is the master node

        VolatilePtr() {	    // volatile pinters behave like normal pointers
            prev = this;    // and are not initialized (ptr) per default
            next = this;
        };

        VolatilePtr(const pointer _p) {
			ptr = _p;
			prev = this;
			next = this;
        }

        VolatilePtr(const _Self& _vp) {
			ptr = _vp.ptr;
			prev = this;
			next = this;
        }

        VolatilePtr(_SelfRef _vp) {
			ptr = _vp.ptr;
			prev = this;
			next = this;
			if (ptr) hangOn(_vp);
		}

		~VolatilePtr() {
			hangOff();
		}
        
        template <typename size_type>
		inline reference operator[] (size_type offset) {
			return ptr[offset];
		}

        template <typename size_type>
		inline const_reference operator[] (size_type offset) const {
			return ptr[offset];
		}

		inline _Self& operator=(_Self const &_Right) {
			hangOff();
			ptr = _Right.ptr;
            if (ptr) hangOn(const_cast<_Self&>(_Right));
			return *this;
		}

		inline _Self& operator=(pointer const _Right) {
			hangOff();
			ptr = _Right;
			return *this;
		}

        inline bool isLonely() {
            return next == this;
        }

		inline void nukeCopies() {
			_SelfPtr p = next;
			while (p != this) {
				_SelfPtr tmp = p->next;
				p->ptr = NULL;
				p->prev = p;
				p->next = p;
				p = tmp;
			}
            prev = this;
			next = this;
		}

		inline bool operator== (const _Self &I) const {
			return ptr == I.ptr;
		}

		inline bool operator!= (const _Self &I) const {
			return ptr != I.ptr;
		}

		inline operator pointer () const {
			return ptr;
		}

	private:

		inline void hangOn(_SelfRef _prev) {
			// hang on between _prev and _prev.next
			prev = &_prev;
			next = _prev.next;
			_prev.next = this;
			next->prev = this;
		}

		inline void hangOff() {
			next->prev = prev;
			prev->next = next;
            next = this;
            prev = this;
		}
	};

    template <typename TValue>
    inline void nukeCopies(TValue* &) {}

    template <typename TValue>
    inline void nukeCopies(VolatilePtr<TValue> &ptr) { ptr.nukeCopies(); }

}

#endif
