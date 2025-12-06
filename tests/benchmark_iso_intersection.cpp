// Copyright 2024 The NurbsDemo Authors
// SPDX-License-Identifier: MIT
//
// benchmark_iso_intersection.cpp
//
// Performance benchmarks for NURBS2D isoline intersection solver.
// Uses Catch2 BENCHMARK feature.
// Run with: xmake run benchmark

#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_test_macros.hpp>
#include <cmath>
#include <random>

#include "nurbs2d_iso_intersect.h"

using namespace nurbs_solver;

// ============================================================================
// Test Curve Creation Functions
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

tinynurbs::Curve<double> create_cubic_curve() {
  tinynurbs::Curve<double> curve;
  curve.degree = 3;
  curve.knots = {0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(0.0, 1.0, 0.0),
                          glm::dvec3(1.0, 0.0, 0.0), glm::dvec3(1.0, 1.0, 0.0)};
  return curve;
}

tinynurbs::RationalCurve<double> create_rational_curve() {
  tinynurbs::RationalCurve<double> curve;
  curve.degree = 2;
  curve.knots = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(1.0, 0.0, 0.0), glm::dvec3(1.0, 1.0, 0.0),
                          glm::dvec3(0.0, 1.0, 0.0)};
  curve.weights = {1.0, std::sqrt(2.0) / 2.0, 1.0};
  return curve;
}

// 创建复杂的高阶曲线
tinynurbs::Curve<double> create_complex_curve() {
  tinynurbs::Curve<double> curve;
  curve.degree = 5;
  // 6 阶需要 n+degree+1 个节点，n=10 个控制点
  curve.knots = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.2, 0.4,
                 0.6, 0.8, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  curve.control_points = {glm::dvec3(0.0, 0.0, 0.0), glm::dvec3(0.1, 0.5, 0.0),
                          glm::dvec3(0.2, 0.2, 0.0), glm::dvec3(0.4, 0.8, 0.0),
                          glm::dvec3(0.5, 0.3, 0.0), glm::dvec3(0.6, 0.7, 0.0),
                          glm::dvec3(0.7, 0.1, 0.0), glm::dvec3(0.8, 0.6, 0.0),
                          glm::dvec3(0.9, 0.4, 0.0), glm::dvec3(1.0, 0.5, 0.0)};
  return curve;
}

// ============================================================================
// Benchmark: Single Iso-line Intersection
// ============================================================================

TEST_CASE("Benchmark: Single iso-line intersection", "[benchmark][single]") {
  auto line_curve = create_line_curve();
  auto quad_curve = create_quadratic_curve();
  auto cubic_curve = create_cubic_curve();
  auto rational_curve = create_rational_curve();
  auto complex_curve = create_complex_curve();

  BENCHMARK("Line curve - single x intersection") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(line_curve, 0.5, IsoLineType::kVertical, roots);
    return roots.size();
  };

  BENCHMARK("Quadratic curve - single x intersection") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(quad_curve, 0.5, IsoLineType::kVertical, roots);
    return roots.size();
  };

  BENCHMARK("Quadratic curve - single y intersection (2 roots)") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(quad_curve, 0.25, IsoLineType::kHorizontal,
                              roots);
    return roots.size();
  };

  BENCHMARK("Cubic curve - single y intersection") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(cubic_curve, 0.5, IsoLineType::kHorizontal,
                              roots);
    return roots.size();
  };

  BENCHMARK("Rational curve (circle) - single x intersection") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(rational_curve, 0.5, IsoLineType::kVertical,
                              roots);
    return roots.size();
  };

  BENCHMARK("Complex degree-5 curve - single y intersection") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(complex_curve, 0.5, IsoLineType::kHorizontal,
                              roots);
    return roots.size();
  };
}

// ============================================================================
// Benchmark: Multiple Iso-lines Intersection
// ============================================================================

TEST_CASE("Benchmark: Multiple iso-lines intersection",
          "[benchmark][multiple]") {
  auto quad_curve = create_quadratic_curve();
  auto cubic_curve = create_cubic_curve();
  auto complex_curve = create_complex_curve();

  // Few iso-lines
  std::vector<double> few_x_lines = {0.25, 0.5, 0.75};
  std::vector<double> few_y_lines = {0.2, 0.4};

  // Medium number of iso-lines
  std::vector<double> medium_x_lines, medium_y_lines;
  for (int i = 1; i < 10; ++i) {
    medium_x_lines.push_back(i * 0.1);
    medium_y_lines.push_back(i * 0.1);
  }

  // Large number of iso-lines
  std::vector<double> many_x_lines, many_y_lines;
  for (int i = 1; i < 50; ++i) {
    many_x_lines.push_back(i * 0.02);
    many_y_lines.push_back(i * 0.02);
  }

  BENCHMARK("Quadratic - few iso-lines (3x + 2y)") {
    std::vector<double> params;
    IntersectCurveWithIsolines(quad_curve, nullptr, &few_x_lines, &few_y_lines,
                               params);
    return params.size();
  };

  BENCHMARK("Quadratic - medium iso-lines (9x + 9y)") {
    std::vector<double> params;
    IntersectCurveWithIsolines(quad_curve, nullptr, &medium_x_lines,
                               &medium_y_lines, params);
    return params.size();
  };

  BENCHMARK("Quadratic - many iso-lines (49x + 49y)") {
    std::vector<double> params;
    IntersectCurveWithIsolines(quad_curve, nullptr, &many_x_lines,
                               &many_y_lines, params);
    return params.size();
  };

  BENCHMARK("Cubic - medium iso-lines (9x + 9y)") {
    std::vector<double> params;
    IntersectCurveWithIsolines(cubic_curve, nullptr, &medium_x_lines,
                               &medium_y_lines, params);
    return params.size();
  };

  BENCHMARK("Complex degree-5 - medium iso-lines (9x + 9y)") {
    std::vector<double> params;
    IntersectCurveWithIsolines(complex_curve, nullptr, &medium_x_lines,
                               &medium_y_lines, params);
    return params.size();
  };
}

// ============================================================================
// Benchmark: Root Finding Methods Comparison
// ============================================================================

TEST_CASE("Benchmark: Root finding methods", "[benchmark][methods]") {
  auto quad_curve = create_quadratic_curve();

  ParameterInterval domain = GetCurveDomain(quad_curve);

  BENCHMARK("Newton method - quadratic curve") {
    IsoIntersectionSolver<tinynurbs::Curve<double>> solver(
        quad_curve, domain, IsoLineType::kHorizontal, 0.25, kDefaultTolerance);
    auto root = solver.FindRootNewton(0.5);
    return root.value_or(-1.0);
  };

  BENCHMARK("Bisection method - quadratic curve") {
    IsoIntersectionSolver<tinynurbs::Curve<double>> solver(
        quad_curve, domain, IsoLineType::kHorizontal, 0.25, kDefaultTolerance);
    auto root = solver.FindRootBisection(0.0, 0.5);
    return root.value_or(-1.0);
  };

  BENCHMARK("Hybrid method - quadratic curve") {
    IsoIntersectionSolver<tinynurbs::Curve<double>> solver(
        quad_curve, domain, IsoLineType::kHorizontal, 0.25, kDefaultTolerance);
    auto root = solver.FindRootHybrid();
    return root.value_or(-1.0);
  };

  BENCHMARK("Find all roots - quadratic curve") {
    IsoIntersectionSolver<tinynurbs::Curve<double>> solver(
        quad_curve, domain, IsoLineType::kHorizontal, 0.25, kDefaultTolerance);
    std::vector<double> roots;
    solver.FindAllRoots(roots);
    return roots.size();
  };
}

// ============================================================================
// Benchmark: Curve Evaluation (Reference)
// ============================================================================

TEST_CASE("Benchmark: Curve evaluation (reference)",
          "[benchmark][evaluation]") {
  auto quad_curve = create_quadratic_curve();
  auto cubic_curve = create_cubic_curve();
  auto rational_curve = create_rational_curve();
  auto complex_curve = create_complex_curve();

  BENCHMARK("EvaluateCurvePoint - quadratic") {
    return EvaluateCurvePoint(quad_curve, 0.5);
  };

  BENCHMARK("EvaluateCurvePoint - cubic") {
    return EvaluateCurvePoint(cubic_curve, 0.5);
  };

  BENCHMARK("EvaluateCurvePoint - rational") {
    return EvaluateCurvePoint(rational_curve, 0.5);
  };

  BENCHMARK("EvaluateCurvePoint - complex degree-5") {
    return EvaluateCurvePoint(complex_curve, 0.5);
  };

  BENCHMARK("EvaluateCurveDerivatives - quadratic") {
    return EvaluateCurveDerivatives(quad_curve, 0.5, 1);
  };

  BENCHMARK("EvaluateCurveDerivatives - cubic") {
    return EvaluateCurveDerivatives(cubic_curve, 0.5, 1);
  };

  BENCHMARK("EvaluateCurveDerivatives - rational") {
    return EvaluateCurveDerivatives(rational_curve, 0.5, 1);
  };
}

// ============================================================================
// Benchmark: Tolerance Impact
// ============================================================================

TEST_CASE("Benchmark: Tolerance impact", "[benchmark][tolerance]") {
  auto quad_curve = create_quadratic_curve();

  BENCHMARK("Tolerance 1e-4") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(quad_curve, 0.25, IsoLineType::kHorizontal, roots,
                              1e-4);
    return roots.size();
  };

  BENCHMARK("Tolerance 1e-6 (default)") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(quad_curve, 0.25, IsoLineType::kHorizontal, roots,
                              1e-6);
    return roots.size();
  };

  BENCHMARK("Tolerance 1e-8") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(quad_curve, 0.25, IsoLineType::kHorizontal, roots,
                              1e-8);
    return roots.size();
  };

  BENCHMARK("Tolerance 1e-10") {
    std::vector<double> roots;
    IntersectCurveWithIsoline(quad_curve, 0.25, IsoLineType::kHorizontal, roots,
                              1e-10);
    return roots.size();
  };
}

// ============================================================================
// Benchmark: Random Batch Test
// ============================================================================

TEST_CASE("Benchmark: Random batch test", "[benchmark][random]") {
  auto quad_curve = create_quadratic_curve();
  auto cubic_curve = create_cubic_curve();

  std::mt19937 rng(42);
  std::uniform_real_distribution<double> dist(0.1, 0.9);

  std::vector<double> random_x_lines(20), random_y_lines(20);
  for (int i = 0; i < 20; ++i) {
    random_x_lines[i] = dist(rng);
    random_y_lines[i] = dist(rng);
  }

  BENCHMARK("Quadratic - 20 random x + 20 random y lines") {
    std::vector<double> params;
    IntersectCurveWithIsolines(quad_curve, nullptr, &random_x_lines,
                               &random_y_lines, params);
    return params.size();
  };

  BENCHMARK("Cubic - 20 random x + 20 random y lines") {
    std::vector<double> params;
    IntersectCurveWithIsolines(cubic_curve, nullptr, &random_x_lines,
                               &random_y_lines, params);
    return params.size();
  };
}
