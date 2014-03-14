#ifndef GPTLTIMER_H
#define GPTLTIMER_H

#include <string>
#include "ITimer.h"
#include "Identifiable.h"


namespace xolotlPerf {

/**
 * The GPTLTimer class is instantiated by the StandardHandlerRegistry class
 * and realizes the ITimer interface to access the timing statistics found
 * General Purpose Timing Library (GPTL).
 */
class GPTLTimer: public ITimer, public xolotlCore::Identifiable {

private:

	/**
	 * The default constructor is declared private since all Timers
	 *  must be initialized with a name.
	 */
    GPTLTimer()
      : xolotlCore::Identifiable("unused")
    { }

public:

	/**
	 * GPTLTimer constructor that takes the argument
	 * timerName
	 *
	 * @param name The GPTLTimer's name
	 */
	GPTLTimer(std::string name)
      : xolotlCore::Identifiable(name)
    { }

	/**
	 * The destructor
	 */
	virtual ~GPTLTimer() { }

    /**
     * This operations starts the ITimer.
     */
	virtual void start();

    /**
     * This operation stops the ITimer.
     */
	virtual void stop();

    /**
     * This operation returns the value of the GPTLTimer.
     */
	virtual double getValue() const;

	/**
	 * This operation returns the units of the GPTLTimer.
	 *
	 * NOTE:  wall -- wallclock time (seconds)
	 * 		  usr -- user CPU time (seconds)
	 * 		  sys -- system CPU time (seconds)
	 *
	 */
	virtual std::string getUnits() const;
};
//end class GPTLTimer

}//end namespace xolotlPerf

#endif
