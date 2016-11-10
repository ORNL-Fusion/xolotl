// Includes
#include "VideoPlot.h"

using namespace xolotlViz;

VideoPlot::VideoPlot(const std::string& name) :
		Plot(name), frameRate(0.0) {
}

VideoPlot::~VideoPlot() {
}

void VideoPlot::render(const std::string&) {
}

void VideoPlot::setFrameRate(double fRate) {
	frameRate = fRate;
	return;
}

double VideoPlot::getFrameRate() const {
	return frameRate;
}
