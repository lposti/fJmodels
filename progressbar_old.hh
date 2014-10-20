// -----------------------------------------------------------------------
// progressbar.hh: Definitions & functions for the ProgressBar class.
// -----------------------------------------------------------------------


#ifndef PROGRESSBAR_HH_
#define PROGRESSBAR_HH_

#include <iostream>
#include <string>
#include <iomanip>

class ProgressBar			/// A class that prints out a progress 
{							///  bar in the form: |\#\#\#\ \ \ \ \ \ \ \ \ \ |		
public:
	ProgressBar(int nlength=20, std::string s="#"); ///< Constructor
	virtual ~ProgressBar() {};              		///< Destructor.
	enum POS {BEG=0,END};                   		///< So that we can record where we are.

	void init(int size);                    		///< Prints empty bar, defines increment.
	void update(int num);                   		///< Prints correct number of hashes
	void rewind();                          		///< Prints backspaces over bar.
	void remove();                          		///< Overwrites bar with blanks
	void fillSpace(std::string someString); 		///< Overwrites bar with a string.

private:
	POS loc;                                		///< Are we at the start or end?
	float stepSize;                         		///< What is the interval between hashes?
	int length;                             		///< What's the maximum number of hashes?
	int numVisible;                         		///< How many hashes are there visible?
	std::string s;									///< String that fill the bar.

	void printBackSpace (int n, std::ostream &str=std::cout) {for(int i=0;i<n;i++) str<<'\b';};	
	void printSpace     (int n, std::ostream &str=std::cout) {for(int i=0;i<n;i++) str<<' ';};
  	void printString    (int n, std::ostream &str=std::cout) {for(int i=0;i<n;i++) str<< s;};

};


ProgressBar::ProgressBar(int nlength, std::string ss) {
	
	/// This constructor enables the user to define how many string should
	/// appear. The number visible is set to 0 and the location to be at the beginning.  
	/// 
	/// \param nlength 	The new number of char to appear in the bar.
	/// \param ss 		The character to appear in the bar.
	
	
	length=nlength; 
	loc=BEG; 
	numVisible = 0;
	s = ss.at(0);
}


void ProgressBar::init(int size) { 
	
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


void ProgressBar::update(int num) {
	
	/// This makes sure the correct number of hashes are drawn.
	/// 
	/// Based on the number provided, as well as the stepsize, we compare
	/// the number of hashes we expect to see with the number that are
	/// there, and if they differ, the correct number are drawn. Again,
	/// we remain at the end.  
	/// 
	/// \param num 	The loop counter to be translated into the
	/// 				progress bar.
	
	std::cout << std::fixed << std::setprecision(2);	
	std::cout.setf(std::ios::right);
	std::cout << std::setw(10) << float(num/(stepSize*length)*100) << std::setw(2) <<" %";  
	std::cout.unsetf(std::ios::right);
	std::cout << std::flush;
	printBackSpace(12);	
		
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


void ProgressBar::rewind() {
	
	/// If we are at the end, we print out enough backspaces to wipe out
	/// the entire bar.  If we are not, the erasing does not need to be
	/// done.
	
	if(loc==END) printBackSpace(length+2); 
	loc=BEG;
	std::cout << std::flush;
}


void ProgressBar::remove() {
	
	/// We first rewind() to the beginning, overwrite the bar with blank spaces, 
	/// and then rewind(). We end up at the beginning.
	
	rewind();
	printSpace(length+14);
	printBackSpace(12);
	loc=END; 
	rewind(); 
	std::cout << std::flush;
}


void ProgressBar::fillSpace(std::string someString) {
	
	/// First remove() the bar and then write out the requested string.
	/// \param someString The string to be written over the bar area.
	
	remove();
	std::cout << someString;
	loc=END;
}

#endif
