// Includes
#include <iostream>

#include <xolotl/viz/standard/plot/SeriesPlot.h>

// General VTKM includes
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CellSetStructured.h>
#include <vtkm/cont/CoordinateSystem.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderRectilinear.h>
#include <vtkm/cont/DeviceAdapter.h>
#include <vtkm/cont/DynamicCellSet.h>
#include <vtkm/cont/ErrorBadValue.h>
#include <vtkm/cont/Field.h>

// VTKM Dataset creator includes
#include <vtkm/cont/DeviceAdapterAlgorithm.h>
#include <vtkm/cont/testing/MakeTestDataSet.h>

// VTKM Rendering includes
#include <vtkm/rendering/Actor.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/Color.h>
#include <vtkm/rendering/MapperWireframer.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/rendering/TextAnnotation.h>
#include <vtkm/rendering/TextAnnotationScreen.h>
#include <vtkm/rendering/View1D.h>

namespace xolotl
{
namespace viz
{
namespace standard
{
namespace plot
{
#define W_WIDTH 1024
#define W_HEIGHT 1024
#define BG_COLOR vtkm::rendering::Color::white
#define FG_COLOR vtkm::rendering::Color::black

SeriesPlot::SeriesPlot(const std::string& name) :
	Plot(name),
	plotDataProviders(std::make_shared<
		std::vector<std::shared_ptr<dataprovider::IDataProvider>>>())
{
}

SeriesPlot::~SeriesPlot()
{
}

void
SeriesPlot::render(const std::string& fileName)
{
	// Check if the label provider is set
	if (!plotLabelProvider) {
		std::cout << "The LabelProvider is not set!!" << std::endl;
		return;
	}

	// Check if the data provider is set
	if (plotDataProviders->empty()) {
		std::cout << "No DataProvider!!" << std::endl;
		return;
	}

	// Define a series of line colors
	vtkm::rendering::Color lineColor[18] = {vtkm::rendering::Color::magenta,
		vtkm::rendering::Color::blue, vtkm::rendering::Color::red,
		vtkm::rendering::Color::cyan, vtkm::rendering::Color::green,
		vtkm::rendering::Color::yellow, vtkm::rendering::Color(0.5, 0.0, 1.0),
		vtkm::rendering::Color(1.0, 0.0, 0.5),
		vtkm::rendering::Color(0.5, 0.0, 0.5),
		vtkm::rendering::Color(0.0, 0.5, 0.5),
		vtkm::rendering::Color(1.0, 0.5, 0.0),
		vtkm::rendering::Color(0.5, 0.0, 0.0),
		vtkm::rendering::Color(0.0, 0.5, 0.0),
		vtkm::rendering::Color(0.0, 0.0, 0.5),
		vtkm::rendering::Color(0.5, 0.5, 0.5),
		vtkm::rendering::Color(0.2, 0.8, 0.2),
		vtkm::rendering::Color(0.8, 0.2, 0.2),
		vtkm::rendering::Color(0.2, 0.2, 0.1)};

	// Create the view
	vtkm::rendering::View1D* view = nullptr;
	vtkm::rendering::CanvasRayTracer canvas(W_WIDTH, W_HEIGHT);
	vtkm::rendering::MapperWireframer mapper;

	// Create the scene
	vtkm::rendering::Scene scene;

	// Create bounds objects to track bounds
	vtkm::Bounds bounds;
	vtkm::Bounds fieldBounds;

	// Create initial data set
	vtkm::cont::DataSetFieldAdd dsf;
	vtkm::cont::DataSetBuilderRectilinear dsb;
	auto xVector = plotDataProviders->at(0)->getAxis1Vector();
	vtkm::cont::DataSet dataSet = dsb.Create(xVector);

	// Loop on all the data providers to plot the different series
	for (int i = 0; i < getDataProviderNumber(); i++) {
		// Get the value that will be plotted on X and Y
		auto yVector = plotDataProviders->at(i)->getAxis2Vector();

		// Add the 1D value to plot
		dsf.AddPointField(
			dataSet, plotDataProviders->at(i)->getDataName(), yVector);

		// Accumulate the bounds of our data to focus camera
		fieldBounds.X = dataSet.GetCoordinateSystem().GetBounds().X;
		dataSet.GetField(i).GetRange(&fieldBounds.Y);
		bounds.Include(fieldBounds);

		// Add Plot to our scene for later rendering
		scene.AddActor(vtkm::rendering::Actor(dataSet.GetCellSet(),
			dataSet.GetCoordinateSystem(),
			dataSet.GetField(plotDataProviders->at(i)->getDataName()),
			lineColor[i % 18]));
	}

	// Set camera position
	vtkm::rendering::Camera camera =
		vtkm::rendering::Camera(vtkm::rendering::Camera::MODE_2D);
	camera.ResetToBounds(bounds, 0, 0, 0);
	camera.SetClippingRange(1.f, 100.f);
	camera.SetViewport(-0.43f, +0.85f, -0.7f, +0.7f);

	// Create vtkm rendering infrastructure
	view = new vtkm::rendering::View1D(
		scene, mapper, canvas, camera, BG_COLOR, FG_COLOR);

	// Set the log scale
	if (enableLogScale)
		view->SetLogY(true);

	// Print the title
	std::unique_ptr<vtkm::rendering::TextAnnotationScreen> titleAnnotation(
		new vtkm::rendering::TextAnnotationScreen(plotLabelProvider->titleLabel,
			FG_COLOR, .1, vtkm::Vec<vtkm::Float32, 2>(-.05, .8), 0));
	view->AddAnnotation(std::move(titleAnnotation));

	// Print x axis label
	std::unique_ptr<vtkm::rendering::TextAnnotationScreen> axis1Annotation(
		new vtkm::rendering::TextAnnotationScreen(plotLabelProvider->axis1Label,
			FG_COLOR, .065, vtkm::Vec<vtkm::Float32, 2>(-.1, -.87), 0));
	view->AddAnnotation(std::move(axis1Annotation));

	// Print y axis label
	std::unique_ptr<vtkm::rendering::TextAnnotationScreen> axis2Annotation(
		new vtkm::rendering::TextAnnotationScreen(plotLabelProvider->axis2Label,
			FG_COLOR, .065, vtkm::Vec<vtkm::Float32, 2>(-.82, -.15), 90));
	view->AddAnnotation(std::move(axis2Annotation));

	// Add the time information
	std::unique_ptr<vtkm::rendering::TextAnnotationScreen> timeAnnotation(
		new vtkm::rendering::TextAnnotationScreen(plotLabelProvider->timeLabel,
			FG_COLOR, .055, vtkm::Vec<vtkm::Float32, 2>(-.85, -.85), 0));
	view->AddAnnotation(std::move(timeAnnotation));

	std::unique_ptr<vtkm::rendering::TextAnnotationScreen> timeStepAnnotation(
		new vtkm::rendering::TextAnnotationScreen(
			plotLabelProvider->timeStepLabel, FG_COLOR, .055,
			vtkm::Vec<vtkm::Float32, 2>(-.85, -.90), 0));
	view->AddAnnotation(std::move(timeStepAnnotation));

	// Set the view
	view->Initialize();

	// Paint
	view->Paint();

	// Save the final buffer as an image
	view->SaveAs(fileName);

	return;
}

void
SeriesPlot::addDataProvider(
	std::shared_ptr<dataprovider::IDataProvider> dataProvider)
{
	plotDataProviders->push_back(dataProvider);
	return;
}

std::shared_ptr<dataprovider::IDataProvider>
SeriesPlot::getDataProvider(int i) const
{
	return plotDataProviders->at(i);
}

int
SeriesPlot::getDataProviderNumber() const
{
	return plotDataProviders->size();
}

} // namespace plot
} // namespace standard
} // namespace viz
} // namespace xolotl
