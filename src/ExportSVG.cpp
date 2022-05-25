#include "ExportSVG.h"
#include "simple_svg_1.0.0.hpp"

void exportSVG(std::string filename, const std::vector<MyPolyline>& polys, double width, double height, std::string bgFilename)
{
	svg::Document doc(filename, svg::Layout(svg::Dimensions(width, height), svg::Layout::TopLeft));
	if (!bgFilename.empty())
	doc << "\n<g>";
	{
		svg::Image bgImg(bgFilename, width, height, -0.5, 0, 0.6);
		doc << bgImg;
	}
    doc << "\n</g>";
    doc << "\n<g>";
	for (int i = 0; i < polys.size(); ++i)
	{
		if (!polys[i].empty())
		{
			svg::Polyline poly(svg::Fill(), svg::Stroke(1, svg::Color(0,0,255)));
			for (auto& p : polys[i])
				poly << svg::Point(p.x() + 0.5, p.y() + 0.5);
			doc << poly;
		}
	}
    doc << "</g>";
	doc.save();
}
