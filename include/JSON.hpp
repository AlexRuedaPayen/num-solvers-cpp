#pragma once
#include<fstream>
#include<string>

#ifndef JSON_
#define JSON_

class JSON{
private:
std::ofstream file;
bool comma_add;
bool comma_start;
public:
JSON(std::string);
void start(std::string);
void add(std::string, int);
void add(std::string, double);
void add(std::string, std::vector<double> &);
void add(std::string, std::vector<int> &);
void end();
void close();
};
#endif