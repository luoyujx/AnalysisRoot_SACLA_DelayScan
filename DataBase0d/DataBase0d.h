#ifndef __DataBasse0d_H__
#define __DataBasse0d_H__

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include "sqlite3.h"
#include "mysql.h"

using std::vector;
using std::map;
using std::string;

class DataBase0d
{
public:
	DataBase0d(void);
	~DataBase0d(void);
public:
	void Connect(const char* DBHOST, const char* DBUSER, const char* DBPASS, const char* DBNAME);
	void LoadDataM(int firstTag, int lastTag, vector<string> &fields, const char* TABLE);
	void ShowTable();
	double GetData(unsigned int tag, unsigned int fieldNumber);
	double GetDataSize()					{ return table.size(); }
	std::pair<int, double> DataBase0d::GetStatusAndData(unsigned int tag, unsigned int fieldNumber);
	int GetLatestTag(const char* TABLE);
	void DataBase0d::CloseMySQL();
private:
	typedef map< unsigned int, vector<double> > tableIntDouble;
	tableIntDouble table;
	vector<string> fieldNames;
	MYSQL* dataBaseM;
	MYSQL_RES *res;
	MYSQL_ROW row;
};
//Functor
class outputTabDelim
{
public:
	outputTabDelim(std::ostream& outStream, size_t sizeOfvec)
		:out(outStream), n(1), size(sizeOfvec){}
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