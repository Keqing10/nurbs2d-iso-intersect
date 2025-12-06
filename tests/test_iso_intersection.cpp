// test_iso_intersection.cpp
//
// Unit tests for NURBS2D isoline intersection solver

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "nurbs2d_iso_intersect.h"

using namespace nurbs_solver;
using namespace Catch::Matchers;

// ============================================================================
// Helper Functions: Test Curve Creation
// ============================================================================

/// Creates a simple line from (0,0) to (1,1)
tinynurbs::Curve<double> create_line_curve() {
  tinynurbs::Curve<double> curve;
  curve.degree = 1;
  curve.knots = {0.0, 0.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(1.0, 1.0, 0.0)};
  return curve;
}

/// Creates a quadratic B-spline curve (parabola shape)
tinynurbs::Curve<double> create_quadratic_curve() {
  tinynurbs::Curve<double> curve;
  curve.degree = 2;
  curve.knots = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(0.0, 0.0, 0.0),
                          glm::dvec3(0.5, 1.0, 0.0),  // apex
                          glm::dvec3(1.0, 0.0, 0.0)};
  return curve;
}

/// Creates a cubic B-spline S-curve
tinynurbs::Curve<double> create_cubic_s_curve() {
  tinynurbs::Curve<double> curve;
  curve.degree = 3;
  curve.knots = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(0.0, 1.0, 0.0),
                          glm::dvec3(1.0, 0.0, 0.0), glm::dvec3(1.0, 1.0, 0.0)};
  return curve;
}

/// Creates a quarter circle arc using rational B-spline.
/// Goes from (1,0) to (0,1).
tinynurbs::RationalCurve<double> create_quarter_circle() {
  tinynurbs::RationalCurve<double> curve;
  curve.degree = 2;
  curve.knots = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(1.0, 0.0, 0.0), glm::dvec3(1.0, 1.0, 0.0),
                          glm::dvec3(0.0, 1.0, 0.0)};
  // Arc weight: middle point = cos(45 deg) = sqrt(2)/2
  curve.weights = {1.0, std::sqrt(2.0) / 2.0, 1.0};
  return curve;
}

// ============================================================================
// ParameterParameterInterval Class Tests
// ============================================================================

TEST_CASE("ParameterInterval basic operations", "[interval]") {
  ParameterInterval interval(0.0, 1.0);

  SECTION("length") { REQUIRE_THAT(interval.Length(), WithinAbs(1.0, 1e-10)); }

  SECTION("mid") { REQUIRE_THAT(interval.Midpoint(), WithinAbs(0.5, 1e-10)); }

  SECTION("contains") {
    REQUIRE(interval.Contains(0.0));
    REQUIRE(interval.Contains(0.5));
    REQUIRE(interval.Contains(1.0));
    REQUIRE_FALSE(interval.Contains(-0.1));
    REQUIRE_FALSE(interval.Contains(1.1));
  }

  SECTION("clamp") {
    REQUIRE_THAT(interval.Clamp(-0.5), WithinAbs(0.0, 1e-10));
    REQUIRE_THAT(interval.Clamp(0.5), WithinAbs(0.5, 1e-10));
    REQUIRE_THAT(interval.Clamp(1.5), WithinAbs(1.0, 1e-10));
  }
}

// ============================================================================
// Curve Evaluation Tests
// ============================================================================

TEST_CASE("2D curve point evaluation", "[curve]") {
  SECTION("Line curve") {
    auto curve = create_line_curve();

    Vec2 pt0 = EvaluateCurvePoint(curve, 0.0);
    REQUIRE_THAT(pt0.x, WithinAbs(0.0, 1e-10));
    REQUIRE_THAT(pt0.y, WithinAbs(0.0, 1e-10));

    Vec2 pt_mid = EvaluateCurvePoint(curve, 0.5);
    REQUIRE_THAT(pt_mid.x, WithinAbs(0.5, 1e-10));
    REQUIRE_THAT(pt_mid.y, WithinAbs(0.5, 1e-10));

    Vec2 pt1 = EvaluateCurvePoint(curve, 1.0);
    REQUIRE_THAT(pt1.x, WithinAbs(1.0, 1e-10));
    REQUIRE_THAT(pt1.y, WithinAbs(1.0, 1e-10));
  }

  SECTION("Quadratic curve") {
    auto curve = create_quadratic_curve();

    // Start point
    Vec2 pt0 = EvaluateCurvePoint(curve, 0.0);
    REQUIRE_THAT(pt0.x, WithinAbs(0.0, 1e-10));
    REQUIRE_THAT(pt0.y, WithinAbs(0.0, 1e-10));

    // Midpoint (should be parabola vertex)
    Vec2 pt_mid = EvaluateCurvePoint(curve, 0.5);
    REQUIRE_THAT(pt_mid.x, WithinAbs(0.5, 1e-10));
    // Y value at midpoint for quadratic B-spline
    REQUIRE(pt_mid.y > 0.0);

    // End point
    Vec2 pt1 = EvaluateCurvePoint(curve, 1.0);
    REQUIRE_THAT(pt1.x, WithinAbs(1.0, 1e-10));
    REQUIRE_THAT(pt1.y, WithinAbs(0.0, 1e-10));
  }
}

TEST_CASE("2D curve derivatives evaluation", "[curve]") {
  auto curve = create_line_curve();

  auto ders = EvaluateCurveDerivatives(curve, 0.5, 1);

  REQUIRE(ders.size() == 2);

  // Position
  REQUIRE_THAT(ders[0].x, WithinAbs(0.5, 1e-10));
  REQUIRE_THAT(ders[0].y, WithinAbs(0.5, 1e-10));

  // First derivative (constant for a line)
  REQUIRE_THAT(ders[1].x, WithinAbs(1.0, 1e-10));
  REQUIRE_THAT(ders[1].y, WithinAbs(1.0, 1e-10));
}

// ============================================================================
// Iso-line Intersection Tests
// ============================================================================

TEST_CASE("Line curve iso intersection", "[iso_intersection]") {
  auto curve = create_line_curve();

  SECTION("Vertical line x = 0.5") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(curve, 0.5, IsoLineType::kVertical, roots);

    REQUIRE(roots.size() == 1);
    REQUIRE_THAT(roots[0], WithinAbs(0.5, 1e-6));

    // Verify intersection point
    Vec2 pt = EvaluateCurvePoint(curve, roots[0]);
    REQUIRE_THAT(pt.x, WithinAbs(0.5, 1e-6));
  }

  SECTION("Horizontal line y = 0.3") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(curve, 0.3, IsoLineType::kHorizontal, roots);

    REQUIRE(roots.size() == 1);
    REQUIRE_THAT(roots[0], WithinAbs(0.3, 1e-6));

    // Verify intersection point
    Vec2 pt = EvaluateCurvePoint(curve, roots[0]);
    REQUIRE_THAT(pt.y, WithinAbs(0.3, 1e-6));
  }
}

TEST_CASE("Quadratic curve iso intersection", "[iso_intersection]") {
  auto curve = create_quadratic_curve();

  SECTION("Horizontal line y = 0.25") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(curve, 0.25, IsoLineType::kHorizontal, roots);

    // Parabola should have two intersections with y = 0.25
    REQUIRE(roots.size() == 2);

    // Verify intersection points
    for (double t : roots) {
      Vec2 pt = EvaluateCurvePoint(curve, t);
      REQUIRE_THAT(pt.y, WithinAbs(0.25, 1e-5));
    }
  }

  SECTION("Vertical line x = 0.5") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(curve, 0.5, IsoLineType::kVertical, roots);

    // Should have one intersection (curve symmetric about x = 0.5)
    REQUIRE(roots.size() == 1);

    Vec2 pt = EvaluateCurvePoint(curve, roots[0]);
    REQUIRE_THAT(pt.x, WithinAbs(0.5, 1e-6));
  }
}

TEST_CASE("Cubic S-curve iso intersection", "[iso_intersection]") {
  auto curve = create_cubic_s_curve();

  SECTION("Horizontal line y = 0.5") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(curve, 0.5, IsoLineType::kHorizontal, roots);

    // S-curve may have multiple intersections with y = 0.5
    REQUIRE(roots.size() >= 1);

    // Verify all intersection points
    for (double t : roots) {
      Vec2 pt = EvaluateCurvePoint(curve, t);
      REQUIRE_THAT(pt.y, WithinAbs(0.5, 1e-5));
    }
  }
}

TEST_CASE("Rational curve (quarter circle) iso intersection",
          "[iso_intersection]") {
  auto curve = create_quarter_circle();

  SECTION("Vertical line x = 0.5") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(curve, 0.5, IsoLineType::kVertical, roots);

    REQUIRE(roots.size() == 1);

    Vec2 pt = EvaluateCurvePoint(curve, roots[0]);
    REQUIRE_THAT(pt.x, WithinAbs(0.5, 1e-5));

    // Point on circle should satisfy x^2 + y^2 = 1
    double radius = std::sqrt(pt.x * pt.x + pt.y * pt.y);
    REQUIRE_THAT(radius, WithinAbs(1.0, 1e-4));
  }

  SECTION("Horizontal line y = 0.5") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(curve, 0.5, IsoLineType::kHorizontal, roots);

    REQUIRE(roots.size() == 1);

    Vec2 pt = EvaluateCurvePoint(curve, roots[0]);
    REQUIRE_THAT(pt.y, WithinAbs(0.5, 1e-5));

    // Point on circle should satisfy x^2 + y^2 = 1
    double radius = std::sqrt(pt.x * pt.x + pt.y * pt.y);
    REQUIRE_THAT(radius, WithinAbs(1.0, 1e-4));
  }
}

// ============================================================================
// Multiple Iso-lines Intersection Tests
// ============================================================================

TEST_CASE("Multiple iso lines intersection", "[iso_intersection]") {
  auto curve = create_line_curve();

  std::vector<double> x_lines = {0.2, 0.4, 0.6, 0.8};
  std::vector<double> y_lines = {0.3, 0.5, 0.7};
  std::vector<double> out_params;

  IntersectCurveWithIsolines(curve, nullptr, &x_lines, &y_lines, out_params);

  // Line intersects 4 vertical + 3 horizontal lines = 7 intersections
  // Since y = x, some may coincide but we get unique params
  REQUIRE(out_params.size() == 7);

  // Verify parameters are sorted
  for (size_t i = 1; i < out_params.size(); ++i) {
    REQUIRE(out_params[i] > out_params[i - 1]);
  }
}

TEST_CASE("Iso intersection with interval restriction", "[iso_intersection]") {
  auto curve = create_line_curve();

  ParameterInterval restricted(0.3, 0.7);
  std::vector<double> x_lines = {0.1, 0.5, 0.9};
  std::vector<double> out_params;

  IntersectCurveWithIsolines(curve, &restricted, &x_lines, nullptr, out_params);

  // Only x = 0.5 is within interval [0.3, 0.7]
  REQUIRE(out_params.size() == 1);
  REQUIRE_THAT(out_params[0], WithinAbs(0.5, 1e-6));
}

// ============================================================================
// Edge Cases Tests
// ============================================================================

TEST_CASE("Edge cases", "[iso_intersection]") {
  auto curve = create_line_curve();

  SECTION("Empty candidates") {
    std::vector<double> empty_lines;
    std::vector<double> out_params;

    IntersectCurveWithIsolines(curve, nullptr, &empty_lines, nullptr, out_params);

    REQUIRE(out_params.empty());
  }

  SECTION("No intersection (line outside range)") {
    std::vector<double> roots;
    // Line goes from (0,0) to (1,1), x = 2 has no intersection
    IntersectCurveWithIsoline(curve, 2.0, IsoLineType::kVertical, roots);

    REQUIRE(roots.empty());
  }

  SECTION("Intersection at endpoint") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(curve, 0.0, IsoLineType::kVertical, roots);

    REQUIRE(roots.size() == 1);
    REQUIRE_THAT(roots[0], WithinAbs(0.0, 1e-6));
  }
}

// ============================================================================
// Utility Functions Tests
// ============================================================================

TEST_CASE("RemoveDuplicateParameters", "[utility]") {
  std::vector<double> params = {0.5, 0.1, 0.5, 0.3, 0.100001, 0.9};
  RemoveDuplicateParameters(params, 1e-4);

  // Should remove duplicate 0.5 and near-duplicate 0.1/0.100001
  REQUIRE(params.size() == 4);  // 0.1, 0.3, 0.5, 0.9

  // Verify sorted
  for (size_t i = 1; i < params.size(); ++i) {
    REQUIRE(params[i] > params[i - 1]);
  }
}

TEST_CASE("GetCurveDomain", "[utility]") {
  auto curve = create_quadratic_curve();

  ParameterInterval domain = GetCurveDomain(curve);

  REQUIRE_THAT(domain.start, WithinAbs(0.0, 1e-10));
  REQUIRE_THAT(domain.end, WithinAbs(1.0, 1e-10));
}

TEST_CASE("EvaluateCurveAtParameters", "[utility]") {
  auto curve = create_line_curve();

  std::vector<double> params = {0.0, 0.5, 1.0};
  auto points = EvaluateCurveAtParameters(curve, params);

  REQUIRE(points.size() == 3);

  REQUIRE_THAT(points[0].x, WithinAbs(0.0, 1e-10));
  REQUIRE_THAT(points[1].x, WithinAbs(0.5, 1e-10));
  REQUIRE_THAT(points[2].x, WithinAbs(1.0, 1e-10));
}

// ============================================================================
// Coincident Region Tests
// ============================================================================

/// Creates a horizontal line segment at y = 0.5 from x=0 to x=1
tinynurbs::Curve<double> create_horizontal_line() {
  tinynurbs::Curve<double> curve;
  curve.degree = 1;
  curve.knots = {0.0, 0.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(0.0, 0.5, 0.0), glm::dvec3(1.0, 0.5, 0.0)};
  return curve;
}

/// Creates a vertical line segment at x = 0.5 from y=0 to y=1
tinynurbs::Curve<double> create_vertical_line() {
  tinynurbs::Curve<double> curve;
  curve.degree = 1;
  curve.knots = {0.0, 0.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(0.5, 0.0, 0.0), glm::dvec3(0.5, 1.0, 0.0)};
  return curve;
}

TEST_CASE("Coincident region - horizontal line on y iso-line", "[coincident]") {
  auto curve = create_horizontal_line();  // y = 0.5 everywhere

  std::vector<double> params;

  SECTION("Curve coincident with y = 0.5 iso-line") {
    std::vector<double> y_lines = {0.5};
    IntersectCurveWithIsolines(curve, nullptr, nullptr, &y_lines, params);

    // Should return both endpoints of the coincident segment
    REQUIRE(params.size() == 2);
    REQUIRE_THAT(params[0], WithinAbs(0.0, 1e-6));
    REQUIRE_THAT(params[1], WithinAbs(1.0, 1e-6));
  }

  SECTION("Curve not coincident with y = 0.3 iso-line") {
    std::vector<double> y_lines = {0.3};
    IntersectCurveWithIsolines(curve, nullptr, nullptr, &y_lines, params);

    // No intersection
    REQUIRE(params.empty());
  }
}

TEST_CASE("Coincident region - vertical line on x iso-line", "[coincident]") {
  auto curve = create_vertical_line();  // x = 0.5 everywhere

  std::vector<double> params;

  SECTION("Curve coincident with x = 0.5 iso-line") {
    std::vector<double> x_lines = {0.5};
    IntersectCurveWithIsolines(curve, nullptr, &x_lines, nullptr, params);

    // Should return both endpoints of the coincident segment
    REQUIRE(params.size() == 2);
    REQUIRE_THAT(params[0], WithinAbs(0.0, 1e-6));
    REQUIRE_THAT(params[1], WithinAbs(1.0, 1e-6));
  }

  SECTION("Curve not coincident with x = 0.7 iso-line") {
    std::vector<double> x_lines = {0.7};
    IntersectCurveWithIsolines(curve, nullptr, &x_lines, nullptr, params);

    // No intersection
    REQUIRE(params.empty());
  }
}
