/* mlev-hiv-model - A program for running multi-level HIV-1 simulations
 * Copyright (C) 2017 Christiaan H. van Dorp <chvandorp@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "eta.hpp"

#include <ctime>
#include <algorithm> // std::max
#include <string>
#include <sstream>


/* methods for TimeConstants */

std::string TimeConstants::mkTimeString(double t) const {
	std::stringstream ss;
	int days = floor(t / secperday);
    t -= days * secperday;
    int hours = floor(t / secperhour);
    t -= hours * secperhour;
    int minutes = floor(t / secperminute);
    t -= minutes * secperminute;
    int seconds = floor(t);
    ss << (days > 0 ? std::to_string(days) + " " : "")
       << hours << ":"
       << (minutes < 10 ? "0" : "") << minutes << ":"
       << (seconds < 10 ? "0" : "") << seconds;
    return ss.str();
}

/* methods for EtaEstimator */

EtaEstimator::EtaEstimator(int N, double alpha) :
        ct(0.0), wct(0.0), etl(0.0), wetl(0.0), n(0), N(N), alpha(alpha) {
    tic = time(nullptr); // the time the Eta was created
    toc = tic;
}

void EtaEstimator::update() {
	time_t now = time(nullptr);
	double dt = difftime(now, toc); // time between now and last update
    toc = now; // the current time
    ct = difftime(toc, tic); // time between start and now
    wct = alpha * wct + dt;
    ++n; // number of steps taken
    etl = std::max((ct/n) * (N-n), 0.0); // expected time left
    wetl = std::max((wct/geometric_sum(alpha, n)) * (N-n), 0.0); // expected time left
}

void EtaEstimator::print(std::ostream & os) const {
	// os << mkTimeString(etl);
    os << mkTimeString(wetl);
}

std::ostream & operator<<(std::ostream & os, const EtaEstimator & eta) {
    eta.print(os);
    return os;
}


/* methods for Stopwatch */

StopWatch::StopWatch() {
	tic = time(nullptr);
}

void StopWatch::print(std::ostream & os) const {
	time_t toc = time(nullptr);
	double dt = difftime(toc, tic);
	os << mkTimeString(dt);
}

std::ostream & operator<<(std::ostream & os, const StopWatch & sw) {
    sw.print(os);
    return os;
}
