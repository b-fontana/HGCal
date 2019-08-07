#include "../interface/parser.h"

void CSVRow::read_next_row(std::istream& str)
{
  std::string line;
  std::getline(str, line);

  std::stringstream lineStream(line);
  std::string cell;

  m_data.clear();
  while(std::getline(lineStream, cell, ','))
    {
      m_data.push_back(cell);
    }
  // This checks for a trailing comma with no data after it.
  if (!lineStream && cell.empty())
    {
      // If there was a trailing comma then add an empty element.
      m_data.push_back("");
    }
}

void CSVRow::bad_row()
{
  std::cout << "Row contents: " << std::endl;
  for(uint_ i=0; i<this->size(); ++i)
    std::cout << (*this)[i] << "  ";
  throw std::length_error("The parsed row has the wrong number of elements.");
}

std::istream& operator>>(std::istream& str, CSVRow& data)
{
  data.read_next_row(str);
  return str;
}

CSVIterator& CSVIterator::operator++()
{
  if (m_str) 
    { 
      if (!((*m_str) >> m_row))
	{
	  m_str = NULL;
	}
    }
  return *this;
}

CSVIterator CSVIterator::operator++(int) 
{
  CSVIterator tmp(*this);
  ++(*this);
  return tmp;
}

CSVRow const& CSVIterator::operator*() const       
{
  return m_row;
}

CSVRow const* CSVIterator::operator->()  const {
  return &m_row;
}

bool CSVIterator::operator==(CSVIterator const& rhs) 
{
  return ((this == &rhs) || ((this->m_str == NULL) && (rhs.m_str == NULL)));
}

bool CSVIterator::operator!=(CSVIterator const& rhs) 
{
  return !((*this) == rhs);
}
