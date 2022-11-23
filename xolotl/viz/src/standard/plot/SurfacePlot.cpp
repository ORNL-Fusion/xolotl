// Includes
#include <string>

#include <xolotl/util/Log.h>
#include <xolotl/viz/standard/plot/SurfacePlot.h>

// General VTKM includes
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CellSetStructured.h>
#include <vtkm/cont/CoordinateSystem.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderRectilinear.h>
#include <vtkm/cont/DynamicCellSet.h>
#include <vtkm/cont/ErrorBadValue.h>
#include <vtkm/cont/Field.h>

// VTKM Dataset creator includes
#include <vtkm/cont/DeviceAdapterAlgorithm.h>

// VTKM Rendering includes
#include <vtkm/cont/ColorTable.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/rendering/Actor.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/Color.h>
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/rendering/TextAnnotation.h>
#include <vtkm/rendering/TextAnnotationScreen.h>
#include <vtkm/rendering/View2D.h>

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

void
SurfacePlot::render(const std::string& fileName)
{
	// Check if the label provider is set
	if (!plotLabelProvider) {
		XOLOTL_LOG_WARN << "The LabelProvider is not set!!";
		return;
	}

	// Check if the data provider is set
	if (!plotDataProvider) {
		XOLOTL_LOG_WARN << "The DataProvider is not set!!";
		return;
	}

	// Get the value that will be plotted on X, Y, and Z
	auto xVector = plotDataProvider->getAxis1Vector();
	auto yVector = plotDataProvider->getAxis2Vector();
	auto zVector = plotDataProvider->getAxis3Vector();

	if (enableLogScale) {
		for (int i = 0; i < zVector.size(); i++) {
			if (zVector[i] > 0) {
				zVector[i] = std::log10(zVector[i]);
			}
			else {
				zVector[i] = 0;
			}
		}
	}

	// Create the vtk-m data set
	vtkm::cont::DataSetBuilderRectilinear dsb;
	vtkm::cont::DataSet dataSet = dsb.Create(xVector, yVector);
	dataSet.AddCellField(plotDataProvider->getDataName(), zVector);

	// Create the view
	vtkm::rendering::Color bg(0.5, 0.5, 0.5, 1.f);
	vtkm::rendering::CanvasRayTracer canvas(W_WIDTH, W_HEIGHT);
	vtkm::rendering::MapperRayTracer mapper;
	mapper.SetShadingOn(false);

	// Create the scene, default to log scale color map
	vtkm::rendering::Scene scene;
	vtkm::rendering::Actor actor(dataSet.GetCellSet(),
		dataSet.GetCoordinateSystem(),
		dataSet.GetField(plotDataProvider->getDataName()),
		vtkm::cont::ColorTable::Preset::Jet);
	scene.AddActor(actor);

	// Set camera position
	vtkm::rendering::Camera camera =
		vtkm::rendering::Camera(vtkm::rendering::Camera::MODE_2D);
	camera.ResetToBounds(dataSet.GetCoordinateSystem().GetBounds());
	camera.SetClippingRange(1.f, 100.f);
	camera.SetViewport(-0.75f, +0.75f, -0.75f, +0.75f);
	camera.SetXScale(1.01f);

	// Create vtkm rendering infrastructure
	vtkm::rendering::View2D view(scene, mapper, canvas, camera, bg);

	// Print the title
	std::string titleLabel;
	float labelLeftPos;
	if (enableLogScale) {
		titleLabel = "Log of " + plotLabelProvider->titleLabel;
		labelLeftPos = -.45;
	}
	else {
		titleLabel = plotLabelProvider->titleLabel;
		labelLeftPos = -.25;
	}
	view.AddAnnotation(std::make_unique<vtkm::rendering::TextAnnotationScreen>(
		titleLabel, vtkm::rendering::Color::white, .07,
		vtkm::Vec<vtkm::Float32, 2>(labelLeftPos, .93), 0));

	// Print the axis labels
	view.AddAnnotation(std::make_unique<vtkm::rendering::TextAnnotationScreen>(
		plotLabelProvider->axis1Label, vtkm::rendering::Color::white, .065,
		vtkm::Vec<vtkm::Float32, 2>(-.12, -.9), 0));

	view.AddAnnotation(std::make_unique<vtkm::rendering::TextAnnotationScreen>(
		plotLabelProvider->axis2Label, vtkm::rendering::Color::white, .065,
		vtkm::Vec<vtkm::Float32, 2>(-.85, -.15), 90));

	view.AddAnnotation(std::make_unique<vtkm::rendering::TextAnnotationScreen>(
		plotLabelProvider->axis3Label, vtkm::rendering::Color::white, .065,
		vtkm::Vec<vtkm::Float32, 2>(-.15, .73), 0));

	// Add the time information
	view.AddAnnotation(std::make_unique<vtkm::rendering::TextAnnotationScreen>(
		plotLabelProvider->timeLabel, vtkm::rendering::Color::white, .055,
		vtkm::Vec<vtkm::Float32, 2>(.6, -.91), 0));

	view.AddAnnotation(std::make_unique<vtkm::rendering::TextAnnotationScreen>(
		plotLabelProvider->timeStepLabel, vtkm::rendering::Color::white, .055,
		vtkm::Vec<vtkm::Float32, 2>(.6, -.96), 0));

	// Paint
	view.Paint();

	// Save the final buffer as an image
	view.SaveAs(fileName);

	return;
}

} // namespace plot
} // namespace standard
} // namespace viz
} // namespace xolotl
