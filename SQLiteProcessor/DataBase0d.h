#ifndef __DataBasse0d_H__
#define __DataBasse0d_H__

#include <vector>
#include <map>
#include <string>
#include "sqlite3.h"
#include <iostream>

using std::vector;
using std::map;
using std::string;

class DataBase0d
{
public:
	DataBase0d(void);
	~DataBase0d(void);
public:
	void Open(const string &dbFileName);
	void LoadData(unsigned int firstTag, unsigned int lastTag, vector<string> &fields); 
	void ShowTable();
	double GetData(unsigned int tag, unsigned int fieldNumber);
private:
	typedef map< unsigned int, vector<double> > tableIntDouble;
	tableIntDouble table;
	vector<string> fieldNames;
	sqlite3* dataBase;
	sqlite3_stmt *statement;
};

//Functor
class outputTabDelim
{
public:
	outputTabDelim(std::ostream& outStream,size_t sizeOfvec)
		:out(outStream),n(1),size(sizeOfvec){}
	template <typename T>
	void operator()(T buf)
	{ 
		out << buf;
		if (n < size) out << "\t";
		else out << std::endl;
		n++;
	}
private:
	std::ostream& out;
	const size_t size;
	size_t n;
};
#endif