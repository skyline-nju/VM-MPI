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
