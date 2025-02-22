/*
 * exception.h - exception class definitions
 *
 * Copyright (C) 2004 Stefan Jahn <stefan@lkcc.org>
 *
 * This is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * This software is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this package; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __EXCEPTION_H__
#define __EXCEPTION_H__

namespace qucs {

enum exception_type {
  EXCEPTION_UNKNOWN = -1,
  EXCEPTION_PIVOT,
  EXCEPTION_NA_FAILED,
  EXCEPTION_NO_CONVERGENCE,
  EXCEPTION_ZERO_DIVISION,
  EXCEPTION_WRONG_VOLTAGE,
  EXCEPTION_SINGULAR,
  EXCEPTION_MATH,
  EXCEPTION_UNKNOWN_ETR_MODE,
};

class exception final {
public:
  exception();
  explicit exception(int);
  exception(const exception &) = delete;
  ~exception();
  [[nodiscard]] int getCode() const { return this->code; }
  void setCode(const int code) { this->code = code; }
  [[nodiscard]] int getData() const { return this->data; }
  void setData(const int data) { this->data = data; }
  [[nodiscard]] char *getText() const { return this->text; }
  void setText(const char *, ...);
  [[nodiscard]] exception *getNext() const { return this->next; }
  void setNext(exception *next) { this->next = next; }

private:
  int code;
  int data;
  char *text;
  exception *next;
};

} /* namespace qucs */

#endif /* __EXCEPTION_H__ */
