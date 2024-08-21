#pragma once
#include <cstdio>
#include <stdexcept>
#include <string>

template<class T>
inline void save_vector(std::FILE*file, const std::vector<T>&v){
  int s = v.size();
  if(std::fwrite(&s, sizeof(s), 1, file) != 1)
    throw std::runtime_error("std::fwrite failed");
  if(std::fwrite(&v[0], sizeof(T), v.size(), file) != v.size())
    throw std::runtime_error("std::fwrite failed");
}
template<class T>
inline void save_cpd_vector(const std::string& file_name, const std::vector<T>&v){
    FILE*file = fopen(file_name.c_str(), "wb");
    int s = v.size();
    if(std::fwrite(&s, sizeof(s), 1, file) != 1)
        throw std::runtime_error("std::fwrite failed");
    if(std::fwrite(&v[0], sizeof(T), v.size(), file) != v.size())
        throw std::runtime_error("std::fwrite failed");
    fclose(file);
}
template<class T>
inline std::vector<T>load_vector(std::FILE*file){
  int s;
  if(std::fread(&s, sizeof(s), 1, file) != 1)
    throw std::runtime_error("std::fread failed");
  std::vector<T>v(s);

  if((int)std::fread(&v[0], sizeof(T), s, file) != s)
    throw std::runtime_error("std::fread failed");

  return v; // NVRO
}
template<class T>
inline std::vector<T>load_cpd_vector(const std::string& file_name){
    FILE*file = fopen(file_name.c_str(), "r");
    int s;
    if(std::fread(&s, sizeof(s), 1, file) != 1)
        throw std::runtime_error("std::fread failed");
    std::vector<T>v(s);

    if((int)std::fread(&v[0], sizeof(T), s, file) != s)
        throw std::runtime_error("std::fread failed");
    fclose(file);
    return v; // NVRO
}
inline void save_string(std::FILE* f, const std::string& s) {
  size_t len = s.size();
  if (std::fwrite(&len, sizeof(size_t), 1, f) != 1)
    throw std::runtime_error("std::fwrite failed");
  if (std::fwrite(&s[0], sizeof(char), s.size(), f) != len)
    throw std::runtime_error("std::fwrite failed");
}


template<class T>
inline std::vector<T>load_vector(const std::string& file_name){
    FILE*file = fopen(file_name.c_str(), "r");
    int s;
    if(std::fread(&s, sizeof(s), 1, file) != 1)
        throw std::runtime_error("std::fread failed");
    std::vector<T>v(s);

    if((int)std::fread(&v[0], sizeof(T), s, file) != s)
        throw std::runtime_error("std::fread failed");
    fclose(file);
    return v; // NVRO
}

inline std::string load_string(std::FILE* f) {
  size_t len;
  if (std::fread(&len, sizeof(size_t),1, f) != 1)
    throw std::runtime_error("std::fread failed");

  char* temp = new char[len+1];
  if (std::fread(temp, sizeof(char), len, f) != len)
    throw std::runtime_error("std::fread failed");

  temp[len] = '\0';
  std::string res = std::string(temp);
  return res;
}
