// Visualization test for NURBS2D isoline intersection using sciplot.
// Run with: xmake run test_visual

#include <cmath>
#include <iostream>
#include <vector>

#include "nurbs2d_iso_intersect.h"
#include "sciplot/sciplot.hpp"
#include "tinynurbs/tinynurbs.h"

using namespace nurbs_solver;
using namespace sciplot;

// ============================================================================
// Helper Functions: Test Curve Creation
// ============================================================================

tinynurbs::Curve<double> create_line_curve() {
  tinynurbs::Curve<double> curve;
  curve.degree = 1;
  curve.knots = {0.0, 0.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(1.0, 1.0, 0.0)};
  return curve;
}

tinynurbs::Curve<double> create_quadratic_curve() {
  tinynurbs::Curve<double> curve;
  curve.degree = 2;
  curve.knots = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(0.5, 1.0, 0.0),
                          glm::dvec3(1.0, 0.0, 0.0)};
  return curve;
}

tinynurbs::Curve<double> create_cubic_s_curve() {
  tinynurbs::Curve<double> curve;
  curve.degree = 3;
  curve.knots = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(0.0, 1.0, 0.0),
                          glm::dvec3(1.0, 0.0, 0.0), glm::dvec3(1.0, 1.0, 0.0)};
  return curve;
}

tinynurbs::RationalCurve<double> create_quarter_circle() {
  tinynurbs::RationalCurve<double> curve;
  curve.degree = 2;
  curve.knots = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(1.0, 0.0, 0.0), glm::dvec3(1.0, 1.0, 0.0),
                          glm::dvec3(0.0, 1.0, 0.0)};
  curve.weights = {1.0, std::sqrt(2.0) / 2.0, 1.0};
  return curve;
}

/// Creates a horizontal line at y = 0.5 (for coincident region test)
tinynurbs::Curve<double> create_horizontal_line() {
  tinynurbs::Curve<double> curve;
  curve.degree = 1;
  curve.knots = {0.0, 0.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(0.0, 0.5, 0.0), glm::dvec3(1.0, 0.5, 0.0)};
  return curve;
}

/// Creates a curve with partial coincident segment
/// The curve has a horizontal segment at y = 0.5 from x = 0.3 to x = 0.7
tinynurbs::Curve<double> create_partial_horizontal_curve() {
  tinynurbs::Curve<double> curve;
  curve.degree = 2;
  curve.knots = {0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0};
  curve.control_points = {
      glm::dvec3(0.0, 0.0, 0.0),  // Start at origin
      glm::dvec3(0.3, 0.5, 0.0),  // Rise to y=0.5
      glm::dvec3(0.7, 0.5, 0.0),  // Horizontal at y=0.5
      glm::dvec3(1.0, 1.0, 0.0)   // Rise to top
  };
  return curve;
}

// ============================================================================
// Helper Functions: Curve Sampling for Plotting
// ============================================================================

template <typename CurveType>
void sample_curve(const CurveType& curve, int num_samples,
                  std::vector<double>& x, std::vector<double>& y) {
  x.clear();
  y.clear();
  x.reserve(num_samples);
  y.reserve(num_samples);

  ParameterInterval domain = GetCurveDomain(curve);
  double step = domain.Length() / (num_samples - 1);

  for (int i = 0; i < num_samples; ++i) {
    double t = domain.start + i * step;
    Vec2 pt = EvaluateCurvePoint(curve, t);
    x.push_back(pt.x);
    y.push_back(pt.y);
  }
}

template <typename CurveType>
void get_intersection_points(const CurveType& curve,
                         const std::vector<double>& params,
                         std::vector<double>& x, std::vector<double>& y) {
  x.clear();
  y.clear();
  x.reserve(params.size());
  y.reserve(params.size());

  for (double t : params) {
    Vec2 pt = EvaluateCurvePoint(curve, t);
    x.push_back(pt.x);
    y.push_back(pt.y);
  }
}

// ============================================================================
// Plotting Functions
// ============================================================================

void plot_all_curves_multiplot() {
  auto line_curve = create_line_curve();
  std::vector<double> line_x, line_y;
  sample_curve(line_curve, 100, line_x, line_y);

  std::vector<double> line_x_lines = {0.3, 0.7};
  std::vector<double> line_y_lines = {0.5};
  std::vector<double> line_params;
  IntersectCurveWithIsolines(line_curve, nullptr, &line_x_lines, &line_y_lines,
                     line_params);
  std::vector<double> line_intersect_x, line_intersect_y;
  get_intersection_points(line_curve, line_params, line_intersect_x, line_intersect_y);

  Plot2D plot1;
  plot1.xlabel("x");
  plot1.ylabel("y");
  plot1.xrange(-0.1, 1.1);
  plot1.yrange(-0.1, 1.1);
  plot1.drawCurve(line_x, line_y).label("Line").lineWidth(2);
  for (double xv : line_x_lines) {
    std::vector<double> vx = {xv, xv};
    std::vector<double> vy = {-0.1, 1.1};
    plot1.drawCurve(vx, vy).label("").lineColor("gray").lineWidth(1);
  }
  for (double yv : line_y_lines) {
    std::vector<double> hx = {-0.1, 1.1};
    std::vector<double> hy = {yv, yv};
    plot1.drawCurve(hx, hy).label("").lineColor("gray").lineWidth(1);
  }
  plot1.drawPoints(line_intersect_x, line_intersect_y)
      .label("Intersections")
      .pointSize(2)
      .pointType(7);

  auto quad_curve = create_quadratic_curve();
  std::vector<double> quad_x, quad_y;
  sample_curve(quad_curve, 100, quad_x, quad_y);

  std::vector<double> quad_x_lines = {0.5};
  std::vector<double> quad_y_lines = {0.25, 0.4};
  std::vector<double> quad_params;
  IntersectCurveWithIsolines(quad_curve, nullptr, &quad_x_lines, &quad_y_lines,
                     quad_params);
  std::vector<double> quad_intersect_x, quad_intersect_y;
  get_intersection_points(quad_curve, quad_params, quad_intersect_x, quad_intersect_y);

  Plot2D plot2;
  plot2.xlabel("x");
  plot2.ylabel("y");
  plot2.xrange(-0.1, 1.1);
  plot2.yrange(-0.1, 0.6);
  plot2.drawCurve(quad_x, quad_y).label("Quadratic").lineWidth(2);
  for (double xv : quad_x_lines) {
    std::vector<double> vx = {xv, xv};
    std::vector<double> vy = {-0.1, 0.6};
    plot2.drawCurve(vx, vy).label("").lineColor("gray").lineWidth(1);
  }
  for (double yv : quad_y_lines) {
    std::vector<double> hx = {-0.1, 1.1};
    std::vector<double> hy = {yv, yv};
    plot2.drawCurve(hx, hy).label("").lineColor("gray").lineWidth(1);
  }
  plot2.drawPoints(quad_intersect_x, quad_intersect_y)
      .label("Intersections")
      .pointSize(2)
      .pointType(7);

  auto s_curve = create_cubic_s_curve();
  std::vector<double> s_x, s_y;
  sample_curve(s_curve, 100, s_x, s_y);

  std::vector<double> s_x_lines = {0.5};
  std::vector<double> s_y_lines = {0.3, 0.5, 0.7};
  std::vector<double> s_params;
  IntersectCurveWithIsolines(s_curve, nullptr, &s_x_lines, &s_y_lines, s_params);
  std::vector<double> s_intersect_x, s_intersect_y;
  get_intersection_points(s_curve, s_params, s_intersect_x, s_intersect_y);

  Plot2D plot3;
  plot3.xlabel("x");
  plot3.ylabel("y");
  plot3.xrange(-0.1, 1.1);
  plot3.yrange(-0.1, 1.1);
  plot3.drawCurve(s_x, s_y).label("S-curve").lineWidth(2);
  for (double xv : s_x_lines) {
    std::vector<double> vx = {xv, xv};
    std::vector<double> vy = {-0.1, 1.1};
    plot3.drawCurve(vx, vy).label("").lineColor("gray").lineWidth(1);
  }
  for (double yv : s_y_lines) {
    std::vector<double> hx = {-0.1, 1.1};
    std::vector<double> hy = {yv, yv};
    plot3.drawCurve(hx, hy).label("").lineColor("gray").lineWidth(1);
  }
  plot3.drawPoints(s_intersect_x, s_intersect_y)
      .label("Intersections")
      .pointSize(2)
      .pointType(7);

  auto circle_curve = create_quarter_circle();
  std::vector<double> circle_x, circle_y;
  sample_curve(circle_curve, 100, circle_x, circle_y);

  std::vector<double> circle_x_lines = {0.5, 0.8};
  std::vector<double> circle_y_lines = {0.5, 0.8};
  std::vector<double> circle_params;
  IntersectCurveWithIsolines(circle_curve, nullptr, &circle_x_lines, &circle_y_lines,
                     circle_params);
  std::vector<double> circle_intersect_x, circle_intersect_y;
  get_intersection_points(circle_curve, circle_params, circle_intersect_x,
                      circle_intersect_y);

  Plot2D plot4;
  plot4.xlabel("x");
  plot4.ylabel("y");
  plot4.xrange(-0.1, 1.2);
  plot4.yrange(-0.1, 1.2);
  plot4.drawCurve(circle_x, circle_y).label("Circle").lineWidth(2);
  for (double xv : circle_x_lines) {
    std::vector<double> vx = {xv, xv};
    std::vector<double> vy = {-0.1, 1.2};
    plot4.drawCurve(vx, vy).label("").lineColor("gray").lineWidth(1);
  }
  for (double yv : circle_y_lines) {
    std::vector<double> hx = {-0.1, 1.2};
    std::vector<double> hy = {yv, yv};
    plot4.drawCurve(hx, hy).label("").lineColor("gray").lineWidth(1);
  }
  plot4.drawPoints(circle_intersect_x, circle_intersect_y)
      .label("Intersections")
      .pointSize(2)
      .pointType(7);

  Figure fig = {{plot1, plot2}, {plot3, plot4}};
  Canvas canvas = {{fig}};
  canvas.title("NURBS Curves with Iso-line Intersections");
  canvas.size(1000, 1000);
  canvas.save("iso_intersections.svg");

  std::cout << "Saved: iso_intersections.svg" << std::endl;
}

// ============================================================================
// Coincident Region Visualization
// ============================================================================

void plot_coincident_test() {
  // Test 1: Fully coincident horizontal line
  auto h_line = create_horizontal_line();
  std::vector<double> h_x, h_y;
  sample_curve(h_line, 100, h_x, h_y);

  std::vector<double> h_y_lines = {0.5};
  std::vector<double> h_params;
  IntersectCurveWithIsolines(h_line, nullptr, nullptr, &h_y_lines, h_params);
  std::vector<double> h_intersect_x, h_intersect_y;
  get_intersection_points(h_line, h_params, h_intersect_x, h_intersect_y);

  Plot2D plot1;
  plot1.xlabel("x");
  plot1.ylabel("y");
  plot1.xrange(-0.1, 1.1);
  plot1.yrange(-0.1, 1.1);
  plot1.drawCurve(h_x, h_y)
      .label("Horizontal Line (y=0.5)")
      .lineWidth(3)
      .lineColor("blue");

  // Draw the coincident iso-line
  std::vector<double> iso_x = {-0.1, 1.1};
  std::vector<double> iso_y = {0.5, 0.5};
  plot1.drawCurve(iso_x, iso_y)
      .label("Iso-line y=0.5")
      .lineColor("red")
      .lineWidth(2)
      .lineType(2);

  // Mark the coincident endpoints
  plot1.drawPoints(h_intersect_x, h_intersect_y)
      .label("Coincident endpoints")
      .pointSize(3)
      .pointType(7);

  std::cout << "Horizontal line coincident test:" << std::endl;
  std::cout << "  Found " << h_params.size() << " intersection params: ";
  for (double t : h_params) {
    std::cout << t << " ";
  }
  std::cout << std::endl;

  // Test 2: Partial coincident curve
  auto partial_curve = create_partial_horizontal_curve();
  std::vector<double> p_x, p_y;
  sample_curve(partial_curve, 100, p_x, p_y);

  std::vector<double> p_y_lines = {0.5};
  std::vector<double> p_params;
  IntersectCurveWithIsolines(partial_curve, nullptr, nullptr, &p_y_lines, p_params);
  std::vector<double> p_intersect_x, p_intersect_y;
  get_intersection_points(partial_curve, p_params, p_intersect_x, p_intersect_y);

  Plot2D plot2;
  plot2.xlabel("x");
  plot2.ylabel("y");
  plot2.xrange(-0.1, 1.1);
  plot2.yrange(-0.1, 1.1);
  plot2.drawCurve(p_x, p_y)
      .label("Curve with horizontal segment")
      .lineWidth(2)
      .lineColor("blue");

  // Draw the iso-line
  plot2.drawCurve(iso_x, iso_y)
      .label("Iso-line y=0.5")
      .lineColor("red")
      .lineWidth(2)
      .lineType(2);

  // Mark intersections
  plot2.drawPoints(p_intersect_x, p_intersect_y)
      .label("Intersections")
      .pointSize(3)
      .pointType(7);

  std::cout << "Partial coincident curve test:" << std::endl;
  std::cout << "  Found " << p_params.size() << " intersection params: ";
  for (double t : p_params) {
    std::cout << t << " ";
  }
  std::cout << std::endl;

  // Test 3: Multiple iso-lines including non-coincident
  auto quad_curve = create_quadratic_curve();
  std::vector<double> q_x, q_y;
  sample_curve(quad_curve, 100, q_x, q_y);

  std::vector<double> q_y_lines = {0.25, 0.5};  // 0.5 is tangent point
  std::vector<double> q_params;
  IntersectCurveWithIsolines(quad_curve, nullptr, nullptr, &q_y_lines, q_params);
  std::vector<double> q_intersect_x, q_intersect_y;
  get_intersection_points(quad_curve, q_params, q_intersect_x, q_intersect_y);

  Plot2D plot3;
  plot3.xlabel("x");
  plot3.ylabel("y");
  plot3.xrange(-0.1, 1.1);
  plot3.yrange(-0.1, 0.6);
  plot3.drawCurve(q_x, q_y)
      .label("Quadratic curve")
      .lineWidth(2)
      .lineColor("blue");

  // Draw iso-lines
  for (double yv : q_y_lines) {
    std::vector<double> hx = {-0.1, 1.1};
    std::vector<double> hy = {yv, yv};
    plot3.drawCurve(hx, hy).label("").lineColor("red").lineWidth(1).lineType(2);
  }

  plot3.drawPoints(q_intersect_x, q_intersect_y)
      .label("Intersections (incl. tangent)")
      .pointSize(3)
      .pointType(7);

  std::cout << "Quadratic with tangent test:" << std::endl;
  std::cout << "  Found " << q_params.size() << " intersection params: ";
  for (double t : q_params) {
    std::cout << t << " ";
  }
  std::cout << std::endl;

  // Test 4: Vertical line coincident test
  tinynurbs::Curve<double> v_line;
  v_line.degree = 1;
  v_line.knots = {0.0, 0.0, 1.0, 1.0};
  v_line.control_points = {glm::dvec3(0.5, 0.0, 0.0),
                           glm::dvec3(0.5, 1.0, 0.0)};

  std::vector<double> v_x, v_y;
  sample_curve(v_line, 100, v_x, v_y);

  std::vector<double> v_x_lines = {0.5};
  std::vector<double> v_params;
  IntersectCurveWithIsolines(v_line, nullptr, &v_x_lines, nullptr, v_params);
  std::vector<double> v_intersect_x, v_intersect_y;
  get_intersection_points(v_line, v_params, v_intersect_x, v_intersect_y);

  Plot2D plot4;
  plot4.xlabel("x");
  plot4.ylabel("y");
  plot4.xrange(-0.1, 1.1);
  plot4.yrange(-0.1, 1.1);
  plot4.drawCurve(v_x, v_y)
      .label("Vertical Line (x=0.5)")
      .lineWidth(3)
      .lineColor("blue");

  // Draw the coincident iso-line
  std::vector<double> viso_x = {0.5, 0.5};
  std::vector<double> viso_y = {-0.1, 1.1};
  plot4.drawCurve(viso_x, viso_y)
      .label("Iso-line x=0.5")
      .lineColor("red")
      .lineWidth(2)
      .lineType(2);

  plot4.drawPoints(v_intersect_x, v_intersect_y)
      .label("Coincident endpoints")
      .pointSize(3)
      .pointType(7);

  std::cout << "Vertical line coincident test:" << std::endl;
  std::cout << "  Found " << v_params.size() << " intersection params: ";
  for (double t : v_params) {
    std::cout << t << " ";
  }
  std::cout << std::endl;

  // Create figure with all plots
  Figure fig = {{plot1, plot2}, {plot3, plot4}};
  Canvas canvas = {{fig}};
  canvas.title("Coincident Region Tests");
  canvas.size(1200, 1000);
  canvas.save("iso_coincident.svg");

  std::cout << "\nSaved: iso_coincident.svg" << std::endl;
}

int main(int argc, char** argv) {
  std::cout << "Generating iso-intersection visualization..." << std::endl;
  plot_all_curves_multiplot();

  std::cout << "\nGenerating coincident region visualization..." << std::endl;
  plot_coincident_test();

  std::cout << "\nDone!" << std::endl;
  return 0;
}
