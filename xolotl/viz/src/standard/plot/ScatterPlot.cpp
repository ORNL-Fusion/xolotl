// Includes
#include <iostream>

#include <xolotl/viz/standard/plot/ScatterPlot.h>

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

ScatterPlot::ScatterPlot(const std::string& name) : Plot(name)
{
}

ScatterPlot::~ScatterPlot()
{
}

void
ScatterPlot::render(const std::string& fileName)
{
	// Check if the label provider is set
	if (!plotLabelProvider) {
		std::cout << "The LabelProvider is not set!!" << std::endl;
		return;
	}

	// Check if the data provider is set
	if (!plotDataProvider) {
		std::cout << "The DataProvider is not set!!" << std::endl;
		return;
	}

	// Get the value that will be plotted on X and Y
	std::vector<double> xVector = plotDataProvider->getAxis1Vector();
	std::vector<double> yVector = plotDataProvider->getAxis2Vector();

	// Create the vtk-m data set
	vtkm::cont::DataSetBuilderUniform dsb;
	vtkm::cont::DataSet dataSet = dsb.Create(xVector.size());

	// Add the 1D value to plot
	dataSet.AddPointField(plotDataProvider->getDataName(), yVector);

	// Create the view
	vtkm::rendering::View1D* view = nullptr;
	vtkm::rendering::CanvasRayTracer canvas(W_WIDTH, W_HEIGHT);
	vtkm::rendering::MapperWireframer mapper;

	// Create the scene
	vtkm::rendering::Scene scene;
	scene.AddActor(vtkm::rendering::Actor(dataSet.GetCellSet(),
		dataSet.GetCoordinateSystem(),
		dataSet.GetField(plotDataProvider->getDataName()),
		vtkm::rendering::Color::magenta));

	// Set camera position
	vtkm::Bounds bounds;
	vtkm::rendering::Camera camera =
		vtkm::rendering::Camera(vtkm::rendering::Camera::MODE_2D);
	bounds.X = dataSet.GetCoordinateSystem().GetBounds().X;
	dataSet.GetField(plotDataProvider->getDataName()).GetRange(&bounds.Y);
	camera.ResetToBounds(bounds, 0, .02, 0);
	camera.SetClippingRange(1.f, 100.f);
	camera.SetViewport(-0.7f, +0.7f, -0.7f, +0.7f);

	// Create vtkm rendering infrastructure
	view = new vtkm::rendering::View1D(
		scene, mapper, canvas, camera, BG_COLOR, FG_COLOR);

	// Set the log scale
	if (enableLogScale)
		view->SetLogY(true);

	// Print the title
	std::unique_ptr<vtkm::rendering::TextAnnotationScreen> titleAnnotation(
		new vtkm::rendering::TextAnnotationScreen(plotLabelProvider->titleLabel,
			FG_COLOR, .1, vtkm::Vec<vtkm::Float32, 2>(-.27, .87), 0));
	view->AddAnnotation(std::move(titleAnnotation));

	// Print the axis labels
	std::unique_ptr<vtkm::rendering::TextAnnotationScreen> axis1Annotation(
		new vtkm::rendering::TextAnnotationScreen(plotLabelProvider->axis1Label,
			FG_COLOR, .065, vtkm::Vec<vtkm::Float32, 2>(-.1, -.87), 0));
	view->AddAnnotation(std::move(axis1Annotation));

	std::unique_ptr<vtkm::rendering::TextAnnotationScreen> axis2Annotation(
		new vtkm::rendering::TextAnnotationScreen(plotLabelProvider->axis2Label,
			FG_COLOR, .065, vtkm::Vec<vtkm::Float32, 2>(-.92, -.15), 90));
	view->AddAnnotation(std::move(axis2Annotation));

	// Add the time information
	std::unique_ptr<vtkm::rendering::TextAnnotationScreen> timeAnnotation(
		new vtkm::rendering::TextAnnotationScreen(plotLabelProvider->timeLabel,
			FG_COLOR, .055, vtkm::Vec<vtkm::Float32, 2>(.55, -.94), 0));
	view->AddAnnotation(std::move(timeAnnotation));

	std::unique_ptr<vtkm::rendering::TextAnnotationScreen> timeStepAnnotation(
		new vtkm::rendering::TextAnnotationScreen(
			plotLabelProvider->timeStepLabel, FG_COLOR, .055,
			vtkm::Vec<vtkm::Float32, 2>(.55, -.99), 0));
	view->AddAnnotation(std::move(timeStepAnnotation));

	// Set the view
	view->Initialize();

	// Paint
	view->Paint();

	// Save the final buffer as an image
	view->SaveAs(fileName);
	return;
}

} // namespace plot
} // namespace standard
} // namespace viz
} // namespace xolotl
