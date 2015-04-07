/*
 * Progressbar.h
 *
 *  Created on: Feb 24, 2015
 *      Author: morpheus
 *
 *  Special Thanks to EdT (enrico.diteodoro@gmail.com)
 *  who allowed me to basically copy the following version
 */

#ifndef INCLUDE_PROGRESSBAR_H_
#define INCLUDE_PROGRESSBAR_H_


#include <iostream>
#include <string>
#include <iomanip>
#include <omp.h>
#include <sstream>


template <class T> std::string to_string (const T& t) {

	/// Convert another type T to string.

	std::stringstream ss;
	ss << t;
	return ss.str();
}
template std::string to_string (const int&);
template std::string to_string (const long&);
template std::string to_string (const float&);
template std::string to_string (const double&);



class ProgressBar			/// A class that prints out a progress
{							///  bar in the form: |\#\#\#\ \ \ \ \ \ \ \ \ \ |
public:
	ProgressBar(int nlength=20, std::string ss="#",
				std::string timeflag="time");		///< Constructor
	virtual ~ProgressBar() {};              		///< Destructor.
	enum POS {BEG=0,END};                   		///< So that we can record where we are.

	void init(int size);                    		///< Prints empty bar, defines increment.
	void update(int num);                   		///< Prints correct number of hashes
	void rewind();                          		///< Prints backspaces over bar.
	void remove();                          		///< Overwrites bar with blanks
	void fillSpace(std::string someString); 		///< Overwrites bar with a string.
	std::string getTimeLeft ();

private:
	POS loc;                                		///< Are we at the start or end?
	float stepSize;                         		///< What is the interval between hashes?
	int length;                             		///< What's the maximum number of hashes?
	int numVisible;                         		///< How many hashes are there visible?
    int stepMade;                                   ///< How many steps we've made so far
	std::string s;									///< String that fill the bar.
	std::string tflag;								///< Do you want the ETE output?

	void printBackSpace (int n, std::ostream &str=std::cout) {for(int i=0;i<n;i++) str<<'\b';};
	void printSpace     (int n, std::ostream &str=std::cout) {for(int i=0;i<n;i++) str<<' ';};
  	void printString    (int n, std::ostream &str=std::cout) {for(int i=0;i<n;i++) str<< s;};

};


inline ProgressBar::ProgressBar(int nlength, std::string ss, std::string timeflag) {

	/// This constructor enables the user to define how many string should
	/// appear. The number visible is set to 0 and the location to be at the beginning.
	///
	/// \param nlength 	The new number of char to appear in the bar.
	/// \param ss 		The character to appear in the bar.


	length=nlength;
	loc=BEG;
	numVisible = 0;
	s = ss.at(0);
	tflag = timeflag;
    stepMade=0;
    stepSize=0;
}


inline void ProgressBar::init(int size) {

	/// This initialises the bar to deal with a loop of a certain size.
	/// This size will imply a certain step size, dependent on the number
	/// of hashes that will be written.  A blank bar is written out as
	/// well, and we remain at the end.
	///
	/// \param size 	The maximum number of iterations to be covered by the
	/// 				progress bar.

	stepSize = float(size) / float(length);
	std::cout << "|";
	printSpace(length);
	std::cout << "|" << std::flush;
	loc = END;
}


inline void ProgressBar::update(int num) {

	/// This makes sure the correct number of hashes are drawn.
	///
	/// Based on the number provided, as well as the stepsize, we compare
	/// the number of hashes we expect to see with the number that are
	/// there, and if they differ, the correct number are drawn. Again,
	/// we remain at the end.
	///
	/// \param num 	The loop counter to be translated into the
	/// 				progress bar.

	std::string timestring;
	if (tflag=="notime") timestring = "no time info";
	else timestring = getTimeLeft();

	std::cout << std::fixed << std::setprecision(2);
	std::cout << std::setw(10) << float(num/(stepSize*length)*100)
			  << std::setw(2) <<" %"
			  << std::setw(20) << timestring;
	std::cout << std::flush;
	printBackSpace(32);

	int numNeeded = 0;
	for(int i=0;i<length;i++)
		if(num>(i*stepSize)) numNeeded++;

	if(numNeeded != numVisible){
		numVisible = numNeeded;
		if(loc==END) printBackSpace(length+2);
		std::cout << "|";
		printString(numNeeded);
		printSpace(length-numNeeded);
		std::cout << "|" << std::flush;
		loc=END;
	}
}


inline void ProgressBar::rewind() {

	/// If we are at the end, we print out enough backspaces to wipe out
	/// the entire bar.  If we are not, the erasing does not need to be
	/// done.

	if(loc==END) printBackSpace(length+2);
	loc=BEG;
	std::cout << std::flush;
}


inline void ProgressBar::remove() {

	/// We first rewind() to the beginning, overwrite the bar with blank spaces,
	/// and then rewind(). We end up at the beginning.

	rewind();
	printSpace(length+34);
	printBackSpace(32);
	loc=END;
	rewind();
	std::cout << std::flush;
}


inline void ProgressBar::fillSpace(std::string someString) {

	/// First remove() the bar and then write out the requested string.
	/// \param someString The string to be written over the bar area.

	remove();
	std::cout << someString;
	loc=END;
}


inline std::string ProgressBar::getTimeLeft () {

	/// Define the time needed to end the loop. Return the time as
	/// format hh mm ss

	std::string tstr;

	static clock_t start;
	clock_t stop;
	double timestep=0, lefttime;

#ifdef _OPENMP
	if (stepMade==0) {
#pragma omp parallel
{
		start = clock()/omp_get_num_threads();
}
		stepMade++;
		return tstr;
	}
#pragma omp parallel
{
	stop  = clock()/omp_get_num_threads();
	timestep = (double(stop-start)/CLOCKS_PER_SEC)/double(stepMade);
	lefttime = timestep*(stepSize*length-stepMade);
}
#else
    if (stepMade==0) {
		start = clock();
		stepMade++;
		return tstr;
	}
	stop  = clock();
	timestep = (double(stop-start)/CLOCKS_PER_SEC)/double(stepMade);
	lefttime = timestep*(stepSize*length-stepMade);
#endif
	stepMade++;

	if (lefttime/60.>1) {
		if (lefttime/3600.>1) {
			int hours = int(lefttime/3600);
			int min = (int(lefttime)%3600)/60;
			int sec = (int(lefttime)%3600)%60;
			tstr = to_string(hours)+"h";
			tstr = min<10 ? tstr+"0"+to_string(min)+"m": tstr+to_string(min)+"m";
			tstr = sec<10 ? tstr+"0"+to_string(sec)+"s": tstr+to_string(sec)+"s";
		}
		else {
			int min = int(lefttime/60);
			int sec = int(lefttime)%60;
			tstr = to_string(min)+"m";
			tstr = sec<10 ? tstr+"0"+to_string(sec)+"s":tstr+to_string(sec)+"s";
		}
	}
	else tstr = to_string(int(lefttime))+"s";

	return tstr;
}

#endif /* INCLUDE_PROGRESSBAR_H_ */
