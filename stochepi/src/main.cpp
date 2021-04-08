#include "main.hpp"

#include <iostream>
#include <sstream>
#include <getopt.h>
#include <map>

#include "macros.hpp" // RIGHT_HERE
#include "D614G.hpp" // for the SARS-CoV-2 mutation model

enum class Mode {
  HELP=0,
  TEST,
  FILTER,
  SIMULATE
};

const std::map<std::string, Mode> modeMap = {
  std::make_pair("help",     Mode::HELP),
  std::make_pair("test",     Mode::TEST),
  std::make_pair("filter",   Mode::FILTER),
  std::make_pair("simulate", Mode::SIMULATE),
};


/** a parameter value can be given with the command line option -F.
 * the argument has the form "parname=parval", where parval is a number.
 */
std::pair<std::string, double> parse_cl_param(std::string s) {
  std::string delim = "=";
  size_t delim_pos = s.find(delim);
  if ( delim_pos == std::string::npos ) {
    throw std::runtime_error("invalid parameter option given" + RIGHT_HERE);
  } // else...
  std::string parname = s.substr(0, delim_pos);
  // remove the name
  s.erase(0, delim_pos + delim.length());
  // get the value
  std::stringstream ss(s);
  double parval;
  ss >> parval;
  return std::make_pair(parname, parval);
}

const std::string helpmessage =
"This is estavoir version 0.2                          \n"
"Copyright Christiaan H. van Dorp (2018)               \n"
"                                                      \n"
"Command line options:                                 \n"
"    -s --seed                                         \n"
"        a seed the the RNG                            \n"
"    -t --threads                                      \n"
"        the number of CPU threads used                \n"
"    -J --particles                                    \n"
"        the number of particles for SMC or simlations \n"
"    -M --iter                                         \n"
"        the number of iterations for IPF              \n"
"    -D --duplicates                                   \n"
"        the number of duplicate MC computations,      \n"
"        e.g. for comuting the log-likelihood          \n"
"    -G --traject                                      \n"
"        the number of filtered trajectories sampled   \n"
"    -h --help                                         \n"
"        print this help message                       \n"
"    -m --mode {help, test, filter, simulate}          \n"
"        choose what you want estavoir to do           \n"
"    -F --fix-par PARNAME=VALUE                        \n"
"        fix a parameter with name PARNAME to VALUE    \n"
"        e.g for profile likelihood computation        \n"
"    -d --data-file                                    \n"
"        name of the file with time series data        \n"
"    -p --param-file                                   \n"
"        name of the file with parameters              \n"
"    -x --addl-data-file                               \n"
"        name of file with additional data             \n"
"    -i --id                                           \n"
"        string used to identify output file names     \n";


int main(int argc, char* argv[]) {
  std::cout << "# starting estavoir..." << std::endl;

  unsigned long seed = 144169;
  int threads = 24; // use 24 threads by default
  Mode mode = Mode::FILTER;
  int particles = 2000; // the number of particles (J)
  int iter = 200; // the number of iterations for iterated filtering (M)
  int duplicates = 1; // number of duplicate PFs at the end of IPF (with no parameter RW)
  int traject = 25; // number of sampled trajecties
  std::map<std::string, double> fixedParamMap; // used for profile likelihood
  std::string dataFileName; // import data for filtering
  std::string paramFileName; // import parameter values
  std::string id = ""; // use default file names
  // variables for additional data
  std::string addlDataFileName; // additional data
  bool addlDataProvided = false; // unless filename is given

  // get command line options
  int option_char;

  while ( true ) {
    int option_index = 0;
    static struct option long_options[] = {
      {"seed",         required_argument, 0,  's' },
      {"threads",      required_argument, 0,  't' },
      {"particles",    required_argument, 0,  'J' },
      {"iter",         required_argument, 0,  'M' },
      {"duplicates",   required_argument, 0,  'D' },
      {"traject",      required_argument, 0,  'G' },
      {"help",         no_argument,       0,  'h' },
      {"mode",         required_argument, 0,  'm' },
      {"fix-par",      required_argument, 0,  'F' },
      {"data-file",    required_argument, 0,  'd' },
      {"param-file",   required_argument, 0,  'p' },
      {"id",           required_argument, 0,  'i' },
      {"addl-data-file",required_argument, 0,  'x' },
      {0,              0,                 0,   0  }
    };

    option_char = getopt_long(argc, argv, "hs:t:J:M:D:G:m:F:d:p:i:x:", long_options, &option_index);
    if ( option_char == -1) {
      break;
    } // else ...
    switch ( option_char ) {
      case 0: {
        std::cerr << "WARNING: option " << long_options[option_index].name
                  << " not implemented" << std::endl;
        mode = Mode::HELP;
        break;
      }
      case 's': {
        std::stringstream(optarg) >> seed;
        break;
      }
      case 't': {
        std::stringstream(optarg) >> threads;
        break;
      }
      case 'h': {
        mode = Mode::HELP;
        break;
      }
      case 'J': {
        std::stringstream(optarg) >> particles;
        break;
      }
      case 'M': {
        std::stringstream(optarg) >> iter;
        break;
      }
      case 'D': {
        std::stringstream(optarg) >> duplicates;
        break;
      }
      case 'G': {
        std::stringstream(optarg) >> traject;
        break;
      }
      case 'm': {
        try {
          mode = modeMap.at(optarg);
        } catch ( const std::out_of_range & e ) {
          std::cerr << "WARNING: invalid mode given (" << optarg << ")" << std::endl;
          mode = Mode::HELP;
        }
        break;
      }
      case 'F': {
        auto [parname, parval] = parse_cl_param(optarg);
        fixedParamMap[parname] = parval;
        break;
      }
      case 'd': {
        dataFileName = optarg;
        break;
      }
      case 'p': {
        paramFileName = optarg;
        break;
      }
      case 'x': {
        addlDataFileName = optarg;
        addlDataProvided = true;
        break;
      }
      case 'i': {
        id = optarg;
        break;
      }
      case '?': {
        std::cerr << "WARNING: invalid command-line argument" << std::endl;
        mode = Mode::HELP;
        break;
      }
      default: {
        throw std::logic_error("getopt error" + RIGHT_HERE);
      }
    } // switch
  } // while true

  switch ( mode ) {
    case Mode::HELP: {
      std::cout << helpmessage << std::endl;
      return 0;
      break; // redundant
    }
    case Mode::TEST: { // JUST FOR TESTING...
      break;
    }
    case Mode::SIMULATE: {
      simulate_sars2mut_model(seed, threads, particles, paramFileName);
      // todo: options for choosing a model
      return 0;
      break;
    }
    case Mode::FILTER: {
      filter_sars2mut_model(seed, threads, particles, traject, iter,
          duplicates, fixedParamMap, dataFileName, paramFileName, id,
          addlDataProvided, addlDataFileName);
      // todo: options for choosing a model
      return 0;
      break;
    }
    default: {
      throw std::logic_error("invalid mode" + RIGHT_HERE);
      break; // redundant
    }
  }

  // put tests here...

  std::cout << "\n# goodbye!" << std::endl;
  return 0;
}
