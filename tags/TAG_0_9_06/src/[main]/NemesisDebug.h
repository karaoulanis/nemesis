
#ifndef _NEMESIS_DEBUG_
#define _NEMESIS_DEBUG_

#include <fstream>
#include <iostream>
#include <string>
#include <Matrix.h>
#include <Vector.h>
#include <string>
using namespace std;

class LogFile
{
private:
	ofstream log;
public:
	LogFile(const char* name)						{log.open(name,ios_base::out);}
	~LogFile()										{log.close();}
	LogFile& operator<<(char c)						{log<<c; return *this;}
	LogFile& operator<<(unsigned char c)			{log<<c; return *this;}
	LogFile& operator<<(signed char c)				{log<<c; return *this;}
	LogFile& operator<<(const char *s)				{log<<s; return *this;}
	LogFile& operator<<(const unsigned char *s)		{log<<s; return *this;}
	LogFile& operator<<(const signed char *s)		{log<<s; return *this;}
	LogFile& operator<<(const void *p)				{log<<p; return *this;}
	LogFile& operator<<(int n)						{log<<n; return *this;}
	LogFile& operator<<(unsigned int n)				{log<<n; return *this;}
	LogFile& operator<<(long n)						{log<<n; return *this;}
	LogFile& operator<<(unsigned long n)			{log<<n; return *this;}
	LogFile& operator<<(short n)					{log<<n; return *this;}
	LogFile& operator<<(unsigned short n)			{log<<n; return *this;}
	LogFile& operator<<(double d)					{log<<d; return *this;}
	LogFile& operator<<(float d)					{log<<d; return *this;}
	LogFile& operator<<(bool b)						{log<<b; return *this;}
	LogFile& operator<<(string s)					{log<<s; return *this;}
	LogFile& operator<<(LogFile& (*f)(LogFile&))	{f(*this);return *this;}
	LogFile& operator>>(LogFile& (*f)(LogFile&))	{f(*this);return *this;}
	LogFile& flush()								{log.flush();return *this;}
	void width(int n)								{log.width(n);}
	void fill(int n)								{log.fill(n);}
	void write(const Vector& v)						{log<<v<<endl;}
	void write(const Matrix& m)						{log<<m<<endl;}
};

namespace Counters
{
	int c1=0;
	int c2=0;
};
#endif
