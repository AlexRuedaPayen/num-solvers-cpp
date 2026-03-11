#include "JSON.hpp"
#include<vector>

JSON::JSON (std::string filename) {
	file.open(filename,std::ofstream::out);
	file<<"{";
	comma_add=false;
	comma_start=false;
}
/*
JSON::JSON (const JSON & j) {
	file=j.file;
	comma_add=j.comma_add;
	comma_start=j.comma_start;
	}
*/
void JSON::close() {
	if (file.is_open()) {
		file<<"\n}";file.close();
	}
}

void JSON::start(std::string key) {
	if (file.is_open()) {
		if (comma_start) {file<<",";}
		file<<"\n \""<<key<<"\" : {";
		comma_add=false;
		comma_start=false;
	}
}

void JSON::end() {
	if (file.is_open()) {
		file<<"\n }";
		comma_add =false;
		comma_start=true;
	}
}

void JSON::add(std::string key, int val) {
	if (file.is_open()) {
		if(comma_add) {file<<",";}
		file<<"\n \""<<key<<"\" : "<<val;
		comma_add=true;
	}
}
void JSON::add(std::string key, double val) {
	if (file.is_open()) {
		if(comma_add) {file<<",";}
		file<<"\n \""<<key<<"\" : "<<val;
		comma_add=true;
	}
}

void JSON::add(std::string key, std::vector<double> & val) {
	if (file.is_open()) {
		if(comma_add) {file<<",";}
		file<<"\n \""<<key<<"\" : [";
		for (int i=0;i<val.size()-1;i++) {
			file<<val[i]<<",";
		}
		file<<val[val.size()-1]<<"]";
		comma_add=true;
	}
}

void JSON::add(std::string key, std::vector<int> & val) {
	if (file.is_open()) {
		if(comma_add) {file<<",";}
		file<<"\n \""<<key<<"\" : [";
		for (int i=0;i<val.size()-1;i++) {
			file<<val[i]<<",";
		}
		file<<val[val.size()-1]<<"]";
		comma_add=true;
	}
}