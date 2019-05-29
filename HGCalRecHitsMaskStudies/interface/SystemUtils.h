#include <dirent.h>
#include <errno.h>

namespace SystemUtils {
 
  int createDir(std::string path) {
    DIR* dir = opendir(path);
    if (dir)
      {
	closedir(dir);
	return 1;
      }
    else if (ENOENT == errno)
      {
	return 0;
      }
    else
      {
	std::cout << "Some error ocurred." << std::endl,
	  return -1;
      }
  }
}
