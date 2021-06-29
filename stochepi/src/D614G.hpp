#ifndef D614G_HPP_
#define D614G_HPP_

/** @file D164G.hpp
 * @short Fit model to COVID-19 D164G mutant data
 */

#include <map>
#include <string>

void simulate_sars2mut_model(unsigned long seed, int threads, int J,
    const std::string & paramFileName, double tmax);


void filter_sars2mut_model(unsigned long seed, int threads, int J, int G, int M,
    int D, const std::map<std::string, double> & fixedParamMap,
    const std::string & dataFileName, const std::string & paramFileName,
    const std::string & id); // no splines for migration

void filter_sars2mut_model(unsigned long seed, int threads, int J, int G, int M,
    int D, const std::map<std::string, double> & fixedParamMap,
    const std::string & dataFileName, const std::string & paramFileName,
    const std::string & id, bool timeDepMigr, const std::string & extFoiFileName,
    bool endTimeProvided, double endTime);



#endif
