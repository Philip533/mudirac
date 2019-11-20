#include <exception>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <chrono>
#include "../lib/config.hpp"
#include "../lib/atom.hpp"
#include "../lib/output.hpp"
#include "../lib/utils.hpp"
#include "../lib/elements.hpp"
#include "../lib/constants.hpp"

#include "../vendor/aixlog/aixlog.hpp"

#define PROJECT_VERSION "0.1.0"

using namespace std;

struct TransLineSpec
{
    int n1, n2;
    int l1, l2;
    bool s1, s2;
};
