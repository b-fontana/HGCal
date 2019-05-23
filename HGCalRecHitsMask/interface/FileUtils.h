#ifndef FileUtils_h
#define FileUtils_h

namespace FileUtils {
  class outData {
  public:
    std::ofstream sfile;
    std::unordered_map<int, std::string> names;
    std::string ext = ".txt";

    outData() {}
    outData(std::unordered_map<int, std::string> names_) {
      this->names=names_;
    }
    outData(std::unordered_map<int, std::string> names_, std::string ext_) {
      this->names=names_;
      this->ext = ext_;
    }

    void setName(std::string n_, int pos_) {
      this->names.insert(std::make_pair(pos_, n_));
    }    
    void setNames(std::unordered_map<int, std::string> names_) {
      this->names=names_;
    }
    void setExt(std::string e_) {this->ext=e_;}

    std::string getName(int pos_) {
      return this->names[pos_];
    }
    std::string getExt() {
      return this->ext;
    }
  };

  //templates
  template <typename Data>
    void create(Data& d, int pos) {
    d.sfile.open(d.names[pos]+d.ext, std::ios_base::out);
    d.sfile.close();
    d.sfile.clear();
  }

  template <typename Data>
    void write(Data& d, std::string t) {
    d.sfile << t << std::endl;
  }
  
  template <typename Data>
    void reopen(Data& d, int pos, std::ios_base::openmode m) {
    d.sfile.close();
    d.sfile.clear();
    d.sfile.open(d.names[pos]+d.ext, m);
  }

  template <typename Data>
    void close(Data& d) {
    d.sfile.close();
    d.sfile.clear();
  }
}

#endif
