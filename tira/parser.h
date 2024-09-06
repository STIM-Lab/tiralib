#ifndef TIRA_PARSER_H
#define TIRA_PARSER_H

#include <fstream>
#include <vector>
#include <sstream>
#include <unordered_map>

namespace tira{

class parser{
	typedef std::vector< std::string > Arguments;		// an array of Arguments in string format
	typedef std::vector< Arguments > Entries;			// an Entry is a set of arguments for a specific option

private:
	// the arguments will be stored in an unordered map using the argument name as the key
	std::unordered_map<std::string, Entries> _options;

public:

	static std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
		std::stringstream ss(s);
		std::string item;
		while (std::getline(ss, item, delim)) {
			if(item.size() != 0)
				elems.push_back(item);
		}
		return elems;
	}

	static std::vector<std::string> split(const std::string &s, char delim) {
		std::vector<std::string> elems;
		split(s, delim, elems);
		return elems;
	}

	// constructor parses an input file
	parser(std::string filename) {
		std::ifstream infile(filename);										// open the option file for parsing
		if(!infile.is_open()) {												// if the file doesn't exist, return w/ error
			std::cerr << "Error opening file " << filename << std::endl;
			return;
		}

		std::string line;									// stores a single line from the file
		while(std::getline(infile, line)) {					// for each line in the file
			std::vector<std::string> tokens;
			tokens = split(line, ' ');						// split the line into a series of tokens
			if (tokens.size() == 0) continue;				// if the line has no tokens, ignore it
			if (tokens[0][0] == '#') continue;				// if the first character of the first token is #, ignore the line

			// the first token is the option name, so create a list of arguments from the following tokens
			Arguments args = std::vector<std::string>(tokens.begin() + 1, tokens.end());
			_options[tokens[0]].push_back(args);			// add the tokens to the specified option
		}
	}

	// Get an argument from the first
	template <typename T>
	T get(std::string name, size_t arg) {
		Entries e = _options[name];
		std::stringstream ss;
		ss << e[0][arg];
		T data;
		ss >> data;
		return data;
	}

	template <typename T>
	T get(std::string name, size_t entry, size_t arg) {
		Entries e = _options[name];

		std::stringstream ss;
		ss << e[entry][arg];
		T data;
		ss >> data;
		return data;
	}

	template <typename T>
	std::vector< std::vector<T> > get(std::string name) {
		Entries e = _options[name];

		std::vector< std::vector<T> > data;
		size_t E = e.size();						// number of entries
		data.resize(E);

		for(size_t ei = 0; ei < E; ei++) {			// for each entry
			size_t A = e[ei].size();				// number of arguments
			data[ei].resize(A);						// allocate space for all of the arguments
			for(size_t ai = 0; ai < A; ai++) {
				data[ei][ai] = (T) atof(e[ei][ai].c_str());
			}
		}

		return data;
	}

	// returns the number of entries for an option name
	size_t count(std::string name) {
		return _options[name].size();
	}

	std::string str() {
		std::stringstream ss;
		for(auto opt = _options.begin(); opt != _options.end(); ++opt) {

			if(opt->second.size() == 1) {
				ss<<opt->first<<":  ";
				for(size_t i = 0; i < opt->second[0].size(); ++i) {
					ss<<opt->second[0][i]<<" ";
				}
				ss<<std::endl;
			}
			else {
				ss<<"[ "<<opt->first<<" ]"<<std::endl;
				for(size_t ei = 0; ei < opt->second.size(); ++ei) {
					for(size_t ai = 0; ai < opt->second[ei].size(); ++ai) {
						ss<<opt->second[ei][ai]<<" ";
					}
					ss<<std::endl;
				}
			}
		}
		return ss.str();
	}

};

}	//end namespace tira

#endif