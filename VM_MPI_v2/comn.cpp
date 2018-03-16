#include "comn.h"

using namespace std;

void mkdir(const char *folder) {
#ifdef _MSC_VER
  if (_access(folder, 0) != 0)
#else
  if (access(folder, 0) != 0)
#endif
  {
    char command[100];
    snprintf(command, 100, "mkdir %s", folder);
    if (system(command))
      cout << "create folder: " << folder << " successfully" << endl;
  } else
    cout << "folder: " << folder << " already exists" << endl;
}

vector<string> split(const string &str, const string &delim) {
  string::size_type pos;
  vector<string> res;
  int str_size = str.size();
  int dlm_size = delim.size();
  for (int i = 0; i < str.size(); i++) {
    pos = str.find(delim, i);
    if (pos == string::npos) {
      res.push_back(str.substr(i, str_size));
      break;
    } else {
      res.push_back(str.substr(i, pos - i));
      i = pos + dlm_size - 1;
    }
  }
  return res;
}

double cal_packing_fraction_2(int n, double Lx, double Ly, double sigma) {
  double a = sigma * 0.5;
  return PI * a * a * n / (Lx * Ly);
}

int cal_particle_number_2(double phi, double Lx, double Ly, double sigma) {
  double phi_max = PI / (2 * std::sqrt(3.0));
  if (phi > phi_max) {
    std::cout << "Input packing fraction phi = " << phi
      << " is larger than phi_max = " << phi_max << std::endl;
    exit(1);
  } else {
    double a = sigma * 0.5;
    return int(std::round(phi * Lx * Ly / (PI * a * a)));
  }
}