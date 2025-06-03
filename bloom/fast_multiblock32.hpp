/* Copyright 2025 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See https://www.boost.org/libs/bloom for library home page.
 */

#ifndef BOOST_BLOOM_FAST_MULTIBLOCK32_HPP
#define BOOST_BLOOM_FAST_MULTIBLOCK32_HPP

#include "detail/avx2.hpp"
#include "detail/neon.hpp"
#include "detail/sse2.hpp"

#if defined(BOOST_BLOOM_AVX2)
#include "detail/fast_multiblock32_avx2.hpp"
#elif defined(BOOST_BLOOM_SSE2) /* important that this comes after AVX2 */
#include "detail/fast_multiblock32_sse2.hpp"
#elif defined(BOOST_BLOOM_LITTLE_ENDIAN_NEON)
#include "detail/fast_multiblock32_neon.hpp"
#else /* fallback */
#include "multiblock.hpp"
#include <boost/cstdint.hpp>
#include <cstddef>

namespace boost{
namespace bloom{

template<std::size_t K>
using fast_multiblock32=multiblock<boost::uint32_t,K>;

} /* namespace bloom */
} /* namespace boost */
#endif

#endif
