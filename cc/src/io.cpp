#include <iostream>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <string>
#include <math.h>
#include <vector>
#include <stdexcept>

// read point particles from text file
std::vector<float> readPoints(std::string fname)
{
    int npar;
    std::vector<float> r;
    float xv, yv, zv;
    float mnx, mny, mnz;
    
    // read position data from a file (indices; no header)
    std::ifstream ifile;
    if (not std::filesystem::exists(fname)){
        throw std::invalid_argument("Requested file does not exist!");
    }
    ifile.open(fname);
  
    // skip 3 comment lines at beginning of file
    constexpr auto max_size = std::numeric_limits<std::streamsize>::max();
    ifile.ignore(max_size, '\n');
    ifile.ignore(max_size, '\n');
    ifile.ignore(max_size, '\n');
  
    // mean position
    mnx = 0.;
    mny = 0.;
    mnz = 0.;
  
    while (!ifile.eof()){
        ifile >> xv;
        ifile >> yv;
        ifile >> zv;
        r.push_back(xv);
        r.push_back(yv);
        r.push_back(zv);
    
        mnx += xv;
        mny += yv;
        mnz += zv;

    }
  
    ifile.close();
    npar = r.size()/3;
    mnx = mnx/npar;
    mny = mny/npar;
    mnz = mnz/npar;
  
    // center particle com at origin
    for (int i = 0; i < npar; i++){   
        r[3*i] = r[3*i]-mnx;
        r[3*i+1] = r[3*i+1]-mny;
        r[3*i+2] = r[3*i+2]-mnz;
    }
    return r;
}

// write vector to text file element by element
void writeVector(std::vector<float> vec, std::string outFile)
{
    std::ofstream file_vec(outFile);
    //std::ostream_iterator<float> fiter(file_vec, "\n");
    //std::copy(vec.begin(), vec.end(), fiter);

    for (auto const &elem : vec) {
        file_vec << std::setw(11) << std::setprecision(7) << std::fixed << elem << "\n";
    }
}
