/* Copyright 2022-2023 Christian Mazakas.
 * Copyright 2024 Braden Ganetsky.
 * Copyright 2025 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * See https://www.boost.org/libs/bloom for library home page.
 */

#ifndef BOOST_BLOOM_DETAIL_TYPE_TRAITS_HPP
#define BOOST_BLOOM_DETAIL_TYPE_TRAITS_HPP

#include <boost/type_traits/make_void.hpp>
#include <type_traits>
#include <utility>

namespace boost{
namespace bloom{
namespace detail{
namespace is_nothrow_swappable_helper_detail{

using std::swap;

template<typename T,typename=void>
struct is_nothrow_swappable_helper
{
  constexpr static bool value=false;
};

template <typename T>
struct is_nothrow_swappable_helper<
  T,
  boost::void_t<decltype(swap(std::declval<T&>(),std::declval<T&>()))>
>
{
  constexpr static bool value=
    noexcept(swap(std::declval<T&>(),std::declval<T&>()));
};

} /* namespace is_nothrow_swappable_helper_detail */

template <class T>
struct is_nothrow_swappable:std::integral_constant<
  bool,
  is_nothrow_swappable_helper_detail::is_nothrow_swappable_helper<T>::value
>{};

#define BOOST_BLOOM_STATIC_ASSERT_IS_NOTHROW_SWAPPABLE(T) \
static_assert(                                            \
  boost::bloom::detail::is_nothrow_swappable< T >::value, \
  #T " must be nothrow swappable")

template<typename T>
struct is_cv_unqualified_object:std::integral_constant<
  bool,
  !std::is_const<T>::value&&
  !std::is_volatile<T>::value&&
  !std::is_function<T>::value&&
  !std::is_reference<T>::value&&
  !std::is_void<T>::value
>{};

#define BOOST_BLOOM_STATIC_ASSERT_IS_CV_UNQUALIFIED_OBJECT(T) \
static_assert(                                                \
  boost::bloom::detail::is_cv_unqualified_object< T >::value, \
  #T " must be a cv-unqualified object type")

template<typename T>
struct remove_cvref
{
  using type=
    typename std::remove_cv<typename std::remove_reference<T>::type>::type;
};

template<typename T>
using remove_cvref_t=typename remove_cvref<T>::type;

template<typename T,typename=void>
struct is_transparent:std::false_type{};

template<typename T>
struct is_transparent<T,void_t<typename T::is_transparent>>:std::true_type{};

template<typename T,class Q=void>
using enable_if_transparent_t=
  typename std::enable_if<is_transparent<T>::value,Q>::type;

} /* namespace detail */
} /* namespace bloom */
} /* namespace boost */

#endif
