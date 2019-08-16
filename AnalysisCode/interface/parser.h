#ifndef PARSER_H
#define PARSER_H

#include <iterator>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "UserCode/AnalysisCode/interface/types.h"

class CSVRow
{
 public:
  std::string const& operator[](std::size_t index) const
    {
      return m_data[index];
    }
  std::size_t size() const
    {
      return m_data.size();
    }
  void read_next_row(std::istream&);
  void bad_row();

 private:
  std::vector<std::string> m_data;
};

std::istream& operator>>(std::istream&, CSVRow&);

class CSVIterator
{   
 public:
  typedef std::input_iterator_tag     iterator_category;
  typedef CSVRow                      value_type;
  typedef std::size_t                 difference_type;
  typedef CSVRow*                     pointer;
  typedef CSVRow&                     reference;

 CSVIterator(std::istream& str)  :m_str(str.good()?&str:nullptr) { ++(*this); }
 CSVIterator()                   :m_str(nullptr) {}

  // Pre Increment
  CSVIterator& operator++();
  // Post increment
  CSVIterator operator++(int);
  CSVRow const& operator*()   const;
  CSVRow const* operator->()  const;
  bool operator==(CSVIterator const&);
  bool operator!=(CSVIterator const&);

 private:
  std::istream*       m_str;
  CSVRow              m_row;
};

#endif //PARSER_H


/////////////////////////////////
////How to use this header
/*
#include "interface/parser.h"

int main()
{
    std::ifstream       file("plop.csv");

    for(CSVIterator loop(file); loop != CSVIterator(); ++loop)
    {
      std::cout << "2nd Element(" << (*loop)[1] << ", " << (*loop).size() << ")\n";
    }
}
*/
///////////////////////////////////
