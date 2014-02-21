#ifndef TIMER_H
#define TIMER_H

// Include
#include <string>
#include <memory>

using namespace std;

namespace xolotlPerf {

//Realizations of this interface are responsible for the collection of performance timing statistics.
class Timer {

private:

	/** The default constructor is declared private since all timers
	 *  must be initialized with a name.
	 */
	Timer();

protected:

	/**
	 * The name of this Timer.
	 */
	std::string name;

//	/**
//	 * The value of this Timer.
//	 */
//	double value;

public:

	/**
	 * Timer constructor that takes the argument name
	 * to distinguish specific timer.
	 *
	 * @param aname The Timer's name
	 */
	Timer(std::string aname);

	/**
	 * The copy constructor.
	 * @param other The Timer to copy
	 */
	Timer(const Timer &other);

	/** The destructor
	 */
	virtual ~Timer();

	/**
	 * This operation returns a Timer that is created using the copy
	 * constructor. If this Timer is actually a subclass of Timer, the
	 * clone will be of the same type and therefore carry all of the members
	 * and virtual functions of the subclass in addition to those of the
	 * Timer.
	 * @return A copy of this Timer.
	 */
//	virtual std::shared_ptr<Timer> clone();

    // This operations starts the Timer.
    virtual void start() = 0;

    // This operation stops the Timer.
    virtual void stop() = 0;

	/**
	 * This operation returns the name of the Timer.
	 *
	 * @return The name of this Timer
	 */
	virtual const std::string getName() const;

    // This operation returns the value of the Timer.
    virtual double getValue() = 0;
//    virtual double getValue() const;

    // This operation returns the units of the Timer.
    virtual long getUnits() const = 0;

};
//end class Timer

}//end namespace xolotlPerf

#endif
