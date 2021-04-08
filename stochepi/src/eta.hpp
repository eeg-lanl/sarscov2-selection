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

#ifndef ETA_HPP_
#define ETA_HPP_

#include <cmath>
#include <iostream>

/* first n terms of the geometric series:
 * 1 + alpha + alpha^2 + ... + alpha^{n-1}
 */
inline double geometric_sum(double alpha, int n) {
	if ( alpha == 1.0 ) {
		return n;
	}
	else {
		double an = pow(alpha, n);
		return (an-1) / (alpha-1);
	}
}

class TimeConstants {
	/** a super class that defines time constants and can convert between
	 * time formats. To be inherited by other time classes.
	 */
public:
	virtual ~TimeConstants() { /* empty */ }
	std::string mkTimeString(double t) const; // time in seconds
protected:
    static const int secperday = 86400;
    static const int secperhour = 3600;
    static const int secperminute = 60;
};

class EtaEstimator : public TimeConstants {
public:
    EtaEstimator(int N, double alpha=1.0);
    // constuction starts the clock. Pass the number of steps
    void update();
    void print(std::ostream & ) const;
protected:
    double ct, wct, etl, wetl; // TODO: wct: more weight on recent time
    // cumulative time, weighted cumulative time, estimated time left, historic weight
    int n, N; // steps taken, total amount of steps
    double alpha; // discounting history
    time_t tic;
    time_t toc;
};

std::ostream & operator<<(std::ostream & , const EtaEstimator & );

class StopWatch : public TimeConstants {
public:
	StopWatch();
	void print(std::ostream & ) const;
protected:
	time_t tic; // moment the stopwatch was created
};

std::ostream & operator<<(std::ostream & , const StopWatch & );


#endif
