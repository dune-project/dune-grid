// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GRID_IO_FILE_VTK_B64ENC_HH
#define DUNE_GRID_IO_FILE_VTK_B64ENC_HH

#include <assert.h>

namespace Dune {

  /** @file
      @author Christian Engwer
      @brief Simple base64 encode

      We implement the base64 encoding (c.f. RFC 4648 https://tools.ietf.org/html/rfc4648).

      @{
   */

  /** @brief endoing table */
  const char base64table[] =
  {
    'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M',
    'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z',
    'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm',
    'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z',
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '/'
  };

  /** @brief struct representing the three byte text as well as the four 6 bit chunks */
  struct b64chunk
  {
    using size_type = unsigned char;
    size_type size;
    char txt[3];

    void reset()
    {
      size = 0;
      txt[0] = txt[1] = txt[2] = 0;
    }

    int read(const char* t, size_type s)
    {
      size = s>=3 ? 3 : s;
      txt[0] = s>0 ? t[0] : 0;
      txt[1] = s>1 ? t[1] : 0;
      txt[2] = s>2 ? t[2] : 0;
      return size;
    }

    void put(const char c)
    {
      assert (size < 3);
      txt[size++] = c;
    }

    void write(char* t)
    {
      const unsigned A = (txt[0] & 0b1111'1100) >> 2;
      const unsigned B = (txt[0] & 0b0000'0011) << 4 | (txt[1] & 0b1111'0000) >> 4;
      const unsigned C = (txt[1] & 0b0000'1111) << 2 | (txt[2] & 0b1100'0000) >> 6;
      const unsigned D = txt[2] & 0b0011'1111;
      t[0] = size>0 ? base64table[A] : '=';
      t[1] = size>0 ? base64table[B] : '=';
      t[2] = size>1 ? base64table[C] : '=';
      t[3] = size>2 ? base64table[D] : '=';
      size = 0;
    }
  };

  /** @} */

} // namespace Dune

#endif // DUNE_GRID_IO_FILE_VTK_B64ENC_HH
