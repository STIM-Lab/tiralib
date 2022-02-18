#ifndef TIRA_FILENAME_H
#define TIRA_FILENAME_H

#include <stdio.h>  /* defines FILENAME_MAX */
#include <sstream>
#include <vector>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include <tira/parser.h>

// set OS dependent defines, including:
//	1) path division character ('/' or '\') used for output only
//	2) command for getting the current working directory
//	3) header files for getting the current working directory
//	4) include BOOST filesystem class if compiling on GCC
#define STIM_DIV '/'
#ifdef _WIN32
	#include <windows.h>
	#include <direct.h>
	#define GetCurrentDir _getcwd
	#define STIM_FILENAME_DIV '\\'
#else
	#ifdef BOOST_PRECOMPILED
		#include <boost/filesystem.hpp>
	#endif
	#include <unistd.h>
	#define GetCurrentDir getcwd
	#define STIM_FILENAME_DIV '/'
 #endif

namespace tira{

class filepath{
protected:
	std::string _drive;					//drive on which the file is located (used for Windows)
	std::vector<std::string> _path;		//path for the specified file (list of hierarchical directories)

	/// replaces win32 dividers with the Linux standard (used internally by default)
	std::string unix_div(std::string s) {
		std::replace( s.begin(), s.end(), '\\', '/');
		return s;
	}

	/// gets the directory hierarchy for the current working directory
	static std::string cwd(){
		char cCurrentPath[FILENAME_MAX];

		 if (!GetCurrentDir(cCurrentPath, sizeof(cCurrentPath))){
			 std::cout<<"ERROR getting current working directory."<<std::endl;
			 exit(1);
		 }

		 std::stringstream ss;
		 ss<<cCurrentPath;
		 return ss.str();
	}

	/// convert a relative path to an absolute path
	void get_absolute(std::string &drive, std::vector<std::string> &absolute, std::vector<std::string> relative){

		std::string current = cwd();						//get the current directory as a string

		std::string current_drive;
		std::vector<std::string> current_dir;
		parse_path(current_drive, current_dir, current);			//get the current drive and directories
		drive = current_drive;										//all relative paths have to be relative to the current drive
		if (current_dir.size() > 0) {

			// step through each directory in the relative path, adjusting the current directory
			//		index depending on whether the relative path is specified with '.' or '..'
			int current_i = (int)current_dir.size() - 1;
			int relative_i;
			for (relative_i = 0; relative_i < relative.size(); relative_i++) {
				if (relative[relative_i] == "..") current_i--;
				else if (relative[relative_i] != ".") break;
			}
			if (current_i < 0) {
				std::cerr << "ERROR tira::filepath - relative path is incompatible with working directory" << std::endl;
				exit(1);
			}

			absolute.clear();
			for (size_t i = 0; i <= current_i; i++) {
				absolute.push_back(current_dir[i]);
			}
			for (size_t i = relative_i; i < relative.size(); i++) {
				absolute.push_back(relative[i]);
			}
		}
		else {
			if (relative[0] == ".")
				relative = std::vector<std::string>(relative.begin() + 1, relative.end());
			absolute = relative;
		}
	}

	/// Parses a directory string into a drive (NULL if not Windows) and list of hierarchical directories
	void parse_path(std::string &drive, std::vector<std::string> &absolute, std::string dir){
		drive = "";											//initialize the drive to NULL (it will stay that way for Windows)
		std::vector<std::string> path;
		bool relative = true;								//the default path is relative

		if(dir.length() == 0) return;						//return if the file locator is empty
		std::string unix_dir = unix_div(dir);				//replace any Windows division characters with Unix

		if(unix_dir.length() > 1 && unix_dir[1] == ':'){	//if a drive identifier is given
			if(unix_dir[0] > 64 && unix_dir[0] < 91)		//if the drive letter is upper case
				drive = unix_dir[0] + 32;					//convert it to lower case
			else if(unix_dir[0] > 96 && unix_dir[0] < 123)	//if the drive character is lower case
				drive = unix_dir[0];					//store it as-is
			else{											//otherwise throw an error
				std::cerr<<"ERROR tira::filename - drive letter is invalid: "<<unix_dir[0]<<std::endl;
				exit(1);
			}
			unix_dir = unix_dir.substr(2, unix_dir.length()-2);	//extract the directory structure
		}

		if(drive.size() != 0){
			relative = false;
		}
		if (unix_dir.size() > 0) {										//if there is a directory specified, remove surrounding slashes
			if (unix_dir[0] == '/') {						//if there is a leading slash
				relative = false;								//the path is not relative
				unix_dir = unix_dir.substr(1, unix_dir.length() - 1);	//remove the slash
			}
		}
		if(unix_dir.size() > 0){
			if(unix_dir[unix_dir.size()-1] == '/')
				unix_dir = unix_dir.substr(0, unix_dir.length() - 1);
		}

		path = parser::split(unix_dir, '/');					//split up the directory structure
		
		if(relative)
			get_absolute(drive, absolute, path);
		else
			absolute = path;
	}

	

public:
	filepath(){
		_drive = "";
	}

	filepath(const filepath& p){
		*this = p;
	}

	filepath(const std::string s){
		parse_path(_drive, _path, s);
	}

	filepath& operator=(const std::string s){
		parse_path(_drive, _path, s);		//parse the string to get the drive and relative path
		return *this;
	}

	std::string str(){
		std::stringstream ss;
		if(_drive != "")							//if a drive letter is specified
			ss<<_drive<<":";						//add it to the string
		for(size_t i = 0; i < _path.size(); i++)
			ss<<STIM_FILENAME_DIV<<_path[i];
		ss<<STIM_FILENAME_DIV;
		return ss.str();
	}	
};				//end filepath

class filename : public filepath{
protected:
	std::string _prefix;				//filename prefix
	std::string _extension;				//filename extension

	void set_filename(std::string fname){
		size_t ext_i = fname.find_last_of('.');								//calculate the index of the '.'
		if(ext_i != std::string::npos){											//if there is an extension
			_prefix = fname.substr(0, ext_i);				//store the prefix
			_extension = fname.substr(ext_i + 1, fname.size() - ext_i - 1);	//store the extension
		}
		else												//otherwise just store the prefix
			_prefix = fname;
	}

public:
	filename(){}

	filename(filepath p, std::string fname) : filepath(p){
		set_filename(fname);
	}

	filename& operator=(const std::string s){
		std::string unix_s = unix_div(s);					//convert dividers to unix
		size_t name_i = unix_s.find_last_of('/');		//find the index of the last divider

		if(name_i == std::string::npos){					//if there is no divider, this is just a filename
			unix_s = "./" + unix_s;							//append a ./ to the beginning so that the working directory is used
			name_i = 1;
		}

		name_i++;


		std::string filename = unix_s.substr(name_i, unix_s.size() - name_i);	//extract the filename
		std::string filepath = unix_s.substr(0, name_i-1);						//extract the filepath

		filepath::operator=(filepath);						//parse and store the file path

		set_filename(filename);
		return *this;
	}

	filename(std::string s){
		operator=(s);
	}

	bool is_relative(){
		return false;
	}


	std::string str(){
		std::stringstream ss;
		ss<<filepath::str()<<_prefix;
		if(_extension.size() != 0)
			ss<<"."<<_extension;
		return ss.str();
	}

	//return a string for the filename without an extension
	std::string str_noext(){
		std::stringstream ss;
		ss<<filepath::str()<<_prefix;
		return ss.str();
	}

	/// Create a matching file locator with a prefix s
	filename with_prefix(std::string s){
		filename result = *this;
		result._prefix = s;
		return result;
	}

	std::string prefix(){
		return _prefix;
	}

	std::string get_prefix(){
		return _prefix;
	}

	/// Create a matching file locator with the extension changed to s
	filename extension(std::string s){
		filename result = *this;
		result._extension = s;
		return result;
	}

	std::string extension(){
		return _extension;
	}

	filename fname(std::string s){
		filename result = *this;
		size_t ext_i = s.find_last_of('.');								//calculate the index of the '.'
		if(ext_i != std::string::npos){											//if there is an extension
			result._prefix = s.substr(0, ext_i);				//store the prefix
			result._extension = s.substr(ext_i + 1, s.size() - ext_i - 1);	//store the extension
		}
		else												//otherwise just store the prefix
			result._prefix = s;
		return result;
	}

	std::string fname(){
		std::string result = prefix();
		result += ".";
		result += extension();
		return result;
	}

	/// create a matching file locator with the path changed to s
	filename path(std::string s){
		filename result = *this;
		result.parse_path(result._drive, result._path, s);
		return result;
	}

	std::string path(){
		return filepath::str();
	}

	/// Casting operator, casts to a string
	operator std::string(){
		return str();
	}

	/// This function replaces a wildcard in the prefix with the specified string
	filename insert(std::string str){

		filename result = *this;				//initialize the result filename to the current filename
		size_t loc = result._prefix.find('*');		//find a wild card in the string
		if(loc == std::string::npos)						//if a wildcard isn't found
			result._prefix += str;							//append the value to the prefix
		else
			result._prefix.replace(loc, 1, str);			//replace the wildcard with the string
		return result;								//return the result
	}

	/// This function replaces a wildcard in the prefix with the specified integer (with a padding of n)
	filename insert(size_t i, size_t n = 2){

		std::stringstream ss;
		ss << std::setw(n) << std::setfill('0') << i;
		return insert(ss.str());
	}

	///This method returns true if any characters in the filename contain '*' or '?'
	bool wildcards() {
		if (_prefix.find('*') != std::string::npos) return true;
		if (_prefix.find('?') != std::string::npos) return true;
		if (_extension.find('*') != std::string::npos) return true;
		if (_extension.find('?') != std::string::npos) return true;
		return false;
	}

	
	/// Returns a list of files using the current filename as a template.
	/// For example:
	///			C:\this\is\a\path\file*.txt
	///		can be used as a template to find a series of files file001.txt, file002.txt, file003.txt, etc.
	std::vector<filename> get_list(){
		//this is OS dependent, so we're starting with Windows
		//the Unix version requires Boost

#ifdef _WIN32

		HANDLE            hFind = INVALID_HANDLE_VALUE;							//allocate data structures for looping through files
		WIN32_FIND_DATAA   FindFileData;
		std::vector<filename> file_list;									//initialize a list to hold all matching filenames

		std::string path_string = str();
		hFind = FindFirstFileA(path_string.c_str(), &FindFileData);		//find the first file that matches the specified file path

		if (hFind == INVALID_HANDLE_VALUE) { 									//if there are no matching files
			//printf ("Invalid file handle. Error is %u.\n", GetLastError());		//print an error
			return file_list;
		}
		else {
			std::string file_name = FindFileData.cFileName;						//save the file name
			std::string file_path = path();										//the file is in the specified directory, so save it
			filename current_file(file_path + file_name);					//create a filename structure representing this file
			if(!(FindFileData.cFileName[0] == '.' && FindFileData.cFileName[1] == '\0'))
				file_list.push_back(current_file);									//push the new filename to the file list


			// List all the other files in the directory.
			while (FindNextFileA(hFind, &FindFileData) != 0){ 					//iterate until there are no more matching files
				if(!(FindFileData.cFileName[0] == '.' && FindFileData.cFileName[1] == '.' && FindFileData.cFileName[2] == '\0')){	//account for the possibility of the '..' filename
					file_name = FindFileData.cFileName;								//save the next file
					current_file = fname(file_name);								//append the directory
					file_list.push_back(current_file);								//push it to the list
				}
			}
			FindClose(hFind);													//close the file data structure
		}
		return file_list;														//return the list of files

#elif BOOST_PRECOMPILED

		boost::filesystem::path p(path());	//create a path from the current filename
		std::vector<filename> file_list;
		if(boost::filesystem::exists(p)){
			if(boost::filesystem::is_directory(p)){
				typedef std::vector<boost::filesystem::path> vec;             // store paths,
				vec v;                                // so we can sort them later
				std::copy(boost::filesystem::directory_iterator(p), boost::filesystem::directory_iterator(), back_inserter(v));
				std::sort(v.begin(), v.end());             // sort, since directory iteration
				  // is not ordered on some file systems
				//compare file names to the current template (look for wild cards)
				for (vec::const_iterator it(v.begin()), it_end(v.end()); it != it_end; ++it)
				{
					//if the filename is a wild card *or* it matches the read file name
					if( _prefix == "*" || _prefix == (*it).filename().stem().string()){
						//if the extension is a wild card *or* it matches the read file extension
						if( _extension == "*" || "." + _extension == (*it).filename().extension().string()){
							file_list.push_back((*it).string());	//include it in the final list
						}

					}



				}

			}

		}
		return file_list;
#else
		std::cout<<"ERROR: UNIX systems don't support file lists without the Boost precompiled libraries."<<std::endl;
		exit(1);
#endif
		
	}
};				//end filename
}				//end namespace tira
#endif
