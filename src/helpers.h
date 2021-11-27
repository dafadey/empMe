#pragma once
#include <vector>
#include <sstream>
#include <string>
#include <istream>
#include <ostream>

template <typename T>
std::istream& operator>>(std::istream& ss, std::vector<T>& v) {
	char c{};
	std::stringstream vss;
	while(ss.good())
	{
		ss >> c;
		if(c == ',' || c == ' ' || !ss.good())
		{
			T val;
			if(vss >> val)
				v.push_back(val);
			vss = std::stringstream();
		} else
			vss << c;
	}
	
	if(v.size() == 0) // try single value
	{
		ss.clear();
		ss.seekg(0);
		T val;
		if(ss >> val)
			v.push_back(val);
	}
	return ss;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T>& v) {
	os << "{";
	for(int i = 0; i < static_cast<int>(v.size()) - 1; i++)
		os << v[i] << ", ";
	if(v.size()>0)
		os << v[v.size()-1];
	os << "}";
	return os;
}
