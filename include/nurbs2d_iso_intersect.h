#ifndef NURBS2D_ISO_INTERSECT_H_
#define NURBS2D_ISO_INTERSECT_H_

#include <tinynurbs/tinynurbs.h>

#include <algorithm>
#include <cmath>
#include <optional>
#include <vector>

namespace nurbs_solver {

// ============================================================================
// Type Definitions
// ============================================================================

/// 2D vector type (using glm)
using Vec2 = glm::dvec2;

/// 2D non-rational curve type
template <typename T>
using Curve2D = tinynurbs::Curve<T>;

/// 2D rational curve type
template <typename T>
using RationalCurve2D = tinynurbs::RationalCurve<T>;

// ============================================================================
// Constants
// ============================================================================

/// Maximum recursion depth for interval subdivision
inline constexpr int kMaxRecursionDepth = 64;

/// Machine epsilon tolerance for floating point comparisons
inline constexpr double kMachineEpsilon = 1e-12;

/// Default tolerance for UV parameter space
inline constexpr double kDefaultTolerance = 1e-6;

/// Maximum iterations for Newton-Raphson method
inline constexpr int kMaxNewtonIterations = 50;

/// Maximum iterations for bisection method
inline constexpr int kMaxBisectionIterations = 100;

// ============================================================================
// ParameterInterval
// ============================================================================

/// Represents a closed interval [start, end] in parameter space.
struct ParameterInterval {
  double start = 0.0;
  double end = 1.0;

  constexpr ParameterInterval() = default;
  constexpr ParameterInterval(double start_val, double end_val)
      : start(start_val), end(end_val) {}

  /// Returns the length of the interval.
  [[nodiscard]] constexpr double Length() const { return end - start; }

  /// Returns the midpoint of the interval.
  [[nodiscard]] constexpr double Midpoint() const {
    return (start + end) * 0.5;
  }

  /// Returns true if param is within the interval (inclusive).
  [[nodiscard]] constexpr bool Contains(double param) const {
    return param >= start && param <= end;
  }

  /// Returns true if param is within the interval with tolerance.
  [[nodiscard]] constexpr bool Contains(double param, double tolerance) const {
    return param >= (start - tolerance) && param <= (end + tolerance);
  }

  /// Clamps value to the interval bounds.
  [[nodiscard]] constexpr double Clamp(double param) const {
    return std::clamp(param, start, end);
  }
};

// ============================================================================
// Curve Evaluation Functions
// ============================================================================

/// Result of evaluating a curve point and its first derivative.
struct CurveEvaluation {
  Vec2 position;  ///< Curve point C(t)
  Vec2 tangent;   ///< First derivative C'(t)
};

/// Evaluates a 2D point on a non-rational curve (ignoring z-component).
template <typename T>
[[nodiscard]] inline Vec2 EvaluateCurvePoint(const Curve2D<T>& curve, T param) {
  const auto point3d = tinynurbs::curvePoint(curve, param);
  return Vec2(point3d.x, point3d.y);
}

/// Evaluates a 2D point on a rational curve (ignoring z-component).
template <typename T>
[[nodiscard]] inline Vec2 EvaluateCurvePoint(const RationalCurve2D<T>& curve,
                                             T param) {
  const auto point3d = tinynurbs::curvePoint(curve, param);
  return Vec2(point3d.x, point3d.y);
}

/// Evaluates curve point and first derivative (optimized, no heap allocation).
template <typename T>
[[nodiscard]] inline CurveEvaluation EvaluateCurveWithTangent(
    const Curve2D<T>& curve, T param) {
  const auto derivatives = tinynurbs::curveDerivatives(curve, 1, param);
  return {Vec2(derivatives[0].x, derivatives[0].y),
          Vec2(derivatives[1].x, derivatives[1].y)};
}

/// Evaluates rational curve point and first derivative.
template <typename T>
[[nodiscard]] inline CurveEvaluation EvaluateCurveWithTangent(
    const RationalCurve2D<T>& curve, T param) {
  const auto derivatives = tinynurbs::curveDerivatives(curve, 1, param);
  return {Vec2(derivatives[0].x, derivatives[0].y),
          Vec2(derivatives[1].x, derivatives[1].y)};
}

/// Evaluates curve derivatives up to specified order.
template <typename T>
std::vector<Vec2> EvaluateCurveDerivatives(const Curve2D<T>& curve, T t,
                                           int num_derivatives = 1) {
  const auto ders3d = tinynurbs::curveDerivatives(curve, num_derivatives, t);
  std::vector<Vec2> result;
  result.reserve(ders3d.size());
  for (const auto& d : ders3d) {
    result.emplace_back(d.x, d.y);
  }
  return result;
}

/// Evaluates rational curve derivatives up to specified order.
template <typename T>
std::vector<Vec2> EvaluateCurveDerivatives(const RationalCurve2D<T>& curve, T t,
                                           int num_derivatives = 1) {
  const auto ders3d = tinynurbs::curveDerivatives(curve, num_derivatives, t);
  std::vector<Vec2> result;
  result.reserve(ders3d.size());
  for (const auto& d : ders3d) {
    result.emplace_back(d.x, d.y);
  }
  return result;
}

/// Returns the valid parameter domain for a curve.
template <typename T>
ParameterInterval GetCurveDomain(const Curve2D<T>& curve) {
  if (curve.knots.empty()) {
    return ParameterInterval(0.0, 1.0);
  }
  return ParameterInterval(curve.knots[curve.degree],
                           curve.knots[curve.knots.size() - curve.degree - 1]);
}

/// Returns the valid parameter domain for a rational curve.
template <typename T>
ParameterInterval GetCurveDomain(const RationalCurve2D<T>& curve) {
  if (curve.knots.empty()) {
    return ParameterInterval(0.0, 1.0);
  }
  return ParameterInterval(curve.knots[curve.degree],
                           curve.knots[curve.knots.size() - curve.degree - 1]);
}

// ============================================================================
// IsoIntersectionSolver
// ============================================================================

/// Enum specifying the type of iso-line.
enum class IsoLineType {
  kVertical,   ///< x = constant (vertical line)
  kHorizontal  ///< y = constant (horizontal line)
};

/// Represents the equation f(t) = Component(C(t)) - iso_value = 0.
///
/// This class encapsulates the root-finding problem for curve-isoline
/// intersections and provides multiple solving strategies.
template <typename CurveType>
class IsoIntersectionSolver {
 public:
  /// Constructs an iso-intersection solver.
  ///
  /// @param curve The curve to intersect
  /// @param domain The parameter domain to search
  /// @param iso_type Type of iso-line (vertical or horizontal)
  /// @param iso_value The constant value of the iso-line
  /// @param tolerance Convergence tolerance for root finding
  IsoIntersectionSolver(const CurveType& curve, const ParameterInterval& domain,
                        IsoLineType iso_type, double iso_value,
                        double tolerance = kDefaultTolerance)
      : curve_(curve),
        domain_(domain),
        iso_type_(iso_type),
        iso_value_(iso_value),
        tolerance_(tolerance) {}

  /// Result of function evaluation: {f(t), f'(t)}.
  struct FunctionEvaluation {
    double value;       ///< Function value f(t)
    double derivative;  ///< Derivative f'(t)
  };

  /// Evaluates f(t) and f'(t) simultaneously.
  [[nodiscard]] FunctionEvaluation Evaluate(double param) const {
    const auto [position, tangent] = EvaluateCurveWithTangent(curve_, param);
    if (iso_type_ == IsoLineType::kVertical) {
      return {position.x - iso_value_, tangent.x};
    }
    return {position.y - iso_value_, tangent.y};
  }

  /// Evaluates only f(t).
  [[nodiscard]] double EvaluateValue(double param) const {
    const Vec2 position = EvaluateCurvePoint(curve_, param);
    return (iso_type_ == IsoLineType::kVertical) ? (position.x - iso_value_)
                                                 : (position.y - iso_value_);
  }

  /// Finds a root using Newton-Raphson method with damping.
  ///
  /// @param initial_guess Starting parameter value
  /// @return The found root, or std::nullopt if not found
  [[nodiscard]] std::optional<double> FindRootNewton(
      double initial_guess) const {
    double current_param = domain_.Clamp(initial_guess);

    for (int iter = 0; iter < kMaxNewtonIterations; ++iter) {
      const auto [func_value, derivative] = Evaluate(current_param);

      if (std::abs(func_value) < tolerance_) {
        return current_param;
      }

      if (std::abs(derivative) < kMachineEpsilon) {
        break;  // Derivative too small, Newton may fail
      }

      const double step = func_value / derivative;
      double next_param = current_param - step;

      // Apply damping if Newton step leaves domain
      if (!domain_.Contains(next_param)) {
        double damping_factor = 0.5;
        while (!domain_.Contains(next_param) && damping_factor > 0.01) {
          next_param = current_param - damping_factor * step;
          damping_factor *= 0.5;
        }
        if (!domain_.Contains(next_param)) {
          break;
        }
      }

      if (std::abs(next_param - current_param) < kMachineEpsilon) {
        return std::abs(EvaluateValue(next_param)) < tolerance_
                   ? std::optional{next_param}
                   : std::nullopt;
      }

      current_param = next_param;
    }

    return std::nullopt;
  }

  /// Finds a root using bisection method.
  ///
  /// @param bracket_left Left endpoint of bracketing interval
  /// @param bracket_right Right endpoint of bracketing interval
  /// @return The found root, or std::nullopt if not found
  [[nodiscard]] std::optional<double> FindRootBisection(
      double bracket_left, double bracket_right) const {
    double value_left = EvaluateValue(bracket_left);
    double value_right = EvaluateValue(bracket_right);

    if (std::abs(value_left) < tolerance_) {
      return bracket_left;
    }
    if (std::abs(value_right) < tolerance_) {
      return bracket_right;
    }

    if (value_left * value_right > 0) {
      return std::nullopt;  // No sign change, no root guaranteed
    }

    for (int iter = 0; iter < kMaxBisectionIterations; ++iter) {
      const double midpoint = (bracket_left + bracket_right) * 0.5;
      const double value_mid = EvaluateValue(midpoint);

      if (std::abs(value_mid) < tolerance_ ||
          (bracket_right - bracket_left) < kMachineEpsilon) {
        return midpoint;
      }

      if (value_left * value_mid < 0) {
        bracket_right = midpoint;
        value_right = value_mid;
      } else {
        bracket_left = midpoint;
        value_left = value_mid;
      }
    }

    const double final_midpoint = (bracket_left + bracket_right) * 0.5;
    return std::abs(EvaluateValue(final_midpoint)) < tolerance_
               ? std::optional{final_midpoint}
               : std::nullopt;
  }

  /// Finds a root using hybrid Newton-bisection method (recommended).
  ///
  /// This method combines the speed of Newton's method with the
  /// robustness of bisection, automatically falling back when needed.
  ///
  /// @return The found root, or std::nullopt if not found
  [[nodiscard]] std::optional<double> FindRootHybrid() const {
    double bracket_left = domain_.start;
    double bracket_right = domain_.end;
    double value_left = EvaluateValue(bracket_left);
    double value_right = EvaluateValue(bracket_right);

    if (std::abs(value_left) < tolerance_) {
      return bracket_left;
    }
    if (std::abs(value_right) < tolerance_) {
      return bracket_right;
    }

    // No sign change - try Newton from midpoint
    if (value_left * value_right > 0) {
      auto newton_result = FindRootNewton(domain_.Midpoint());
      if (newton_result && domain_.Contains(*newton_result)) {
        return newton_result;
      }
      return std::nullopt;
    }

    double current_param = (bracket_left + bracket_right) * 0.5;

    for (int iter = 0; iter < kMaxNewtonIterations; ++iter) {
      const auto [func_value, derivative] = Evaluate(current_param);

      if (std::abs(func_value) < tolerance_) {
        return current_param;
      }

      double next_param;
      if (std::abs(derivative) > kMachineEpsilon) {
        next_param = current_param - func_value / derivative;
        if (next_param <= bracket_left || next_param >= bracket_right) {
          next_param =
              (bracket_left + bracket_right) * 0.5;  // Bisection fallback
        }
      } else {
        next_param = (bracket_left + bracket_right) * 0.5;
      }

      // Update bracket
      const double value_next = EvaluateValue(next_param);
      if (value_left * value_next < 0) {
        bracket_right = next_param;
        value_right = value_next;
      } else {
        bracket_left = next_param;
        value_left = value_next;
      }

      if (std::abs(bracket_right - bracket_left) < kMachineEpsilon) {
        const double final_param = (bracket_left + bracket_right) * 0.5;
        return std::abs(EvaluateValue(final_param)) < tolerance_
                   ? std::optional{final_param}
                   : std::nullopt;
      }

      current_param = next_param;
    }

    return FindRootBisection(bracket_left, bracket_right);
  }

  /// Finds all roots in the domain using recursive subdivision.
  ///
  /// Uses adaptive interval splitting with derivative analysis
  /// to reliably find all intersections including tangent points.
  /// Also handles coincident regions where curve lies on the iso-line.
  ///
  /// @param[out] roots Vector to store found roots
  void FindAllRoots(std::vector<double>& roots) const {
    roots.clear();

    const double param_start = domain_.start;
    const double param_end = domain_.end;
    const auto [value_start, deriv_start] = Evaluate(param_start);
    const auto [value_end, deriv_end] = Evaluate(param_end);

    // Check for coincident region (curve lies on iso-line)
    if (IsCoincidentRegion(param_start, param_end, value_start, deriv_start,
                           value_end, deriv_end)) {
      // Return both endpoints for coincident segment
      roots.push_back(param_start);
      if (std::abs(param_end - param_start) > tolerance_) {
        roots.push_back(param_end);
      }
      return;
    }

    // Check start endpoint
    if (std::abs(value_start) < tolerance_) {
      roots.push_back(param_start);
    }

    // Recursively subdivide interval
    SubdivideInterval(param_start, param_end, value_start, deriv_start,
                      value_end, deriv_end, kMaxRecursionDepth, roots);

    // Check end point
    if (std::abs(value_end) < tolerance_) {
      if (roots.empty() || std::abs(param_end - roots.back()) > tolerance_) {
        roots.push_back(param_end);
      }
    }
  }

 private:
  /// Check if the entire interval is a coincident region.
  ///
  /// A coincident region means the curve segment lies entirely on the iso-line.
  /// This occurs when function values and derivatives are all near zero.
  ///
  /// @param param_start Start parameter
  /// @param param_end End parameter
  /// @param value_start Function value at param_start
  /// @param deriv_start Derivative at param_start
  /// @param value_end Function value at param_end
  /// @param deriv_end Derivative at param_end
  /// @return true if the interval is coincident with the iso-line
  [[nodiscard]] bool IsCoincidentRegion(double param_start, double param_end,
                                        double value_start, double deriv_start,
                                        double value_end,
                                        double deriv_end) const {
    // Both endpoints must be on the iso-line
    if (std::abs(value_start) >= tolerance_ ||
        std::abs(value_end) >= tolerance_) {
      return false;
    }

    // Both derivatives must be near zero (curve tangent to iso-line)
    if (std::abs(deriv_start) >= tolerance_ ||
        std::abs(deriv_end) >= tolerance_) {
      return false;
    }

    // Check midpoint to confirm it's not just a tangent touch
    const double param_mid = (param_start + param_end) * 0.5;
    const auto [value_mid, deriv_mid] = Evaluate(param_mid);

    // Midpoint should also be on iso-line with near-zero derivative
    if (std::abs(value_mid) >= tolerance_ ||
        std::abs(deriv_mid) >= tolerance_) {
      return false;
    }

    // Check third points for extra confidence
    const double param_third1 = param_start + (param_end - param_start) / 3.0;
    const double param_third2 =
        param_start + (param_end - param_start) * 2.0 / 3.0;
    const auto [value_third1, deriv_third1] = Evaluate(param_third1);
    const auto [value_third2, deriv_third2] = Evaluate(param_third2);

    return std::abs(value_third1) < tolerance_ &&
           std::abs(value_third2) < tolerance_;
  }

  /// Recursively subdivide interval to find all roots.
  void SubdivideInterval(double param_start, double param_end,
                         double value_start, double deriv_start,
                         double value_end, double deriv_end, int depth,
                         std::vector<double>& roots) const {
    if (depth <= 0) return;

    const double interval_length = param_end - param_start;
    if (interval_length < kMachineEpsilon * 10) return;

    // Case 0: Check for local coincident region
    if (IsCoincidentRegion(param_start, param_end, value_start, deriv_start,
                           value_end, deriv_end)) {
      // Add both endpoints of coincident region
      if (roots.empty() || std::abs(param_start - roots.back()) > tolerance_) {
        roots.push_back(param_start);
      }
      if (roots.empty() || std::abs(param_end - roots.back()) > tolerance_) {
        roots.push_back(param_end);
      }
      return;
    }

    // Case 1: Sign change - definite root
    if (value_start * value_end < 0) {
      auto root =
          SolveInBracket(param_start, param_end, value_start, value_end);
      if (root) {
        if (roots.empty() || std::abs(*root - roots.back()) > tolerance_) {
          roots.push_back(*root);
        }
      }
      return;
    }

    // Case 2: Check if monotonic and moving away from zero
    if (CanSkipInterval(value_start, deriv_start, value_end, deriv_end,
                        interval_length)) {
      return;
    }

    // Case 3: Potential extremum - split and recurse
    if (HasPotentialExtremum(value_start, deriv_start, value_end, deriv_end)) {
      SplitAtExtremum(param_start, param_end, value_start, deriv_start,
                      value_end, deriv_end, depth, roots);
      return;
    }

    // Case 4: Default subdivision at midpoint
    const double param_mid = (param_start + param_end) * 0.5;
    const auto [value_mid, deriv_mid] = Evaluate(param_mid);

    SubdivideInterval(param_start, param_mid, value_start, deriv_start,
                      value_mid, deriv_mid, depth - 1, roots);
    SubdivideInterval(param_mid, param_end, value_mid, deriv_mid, value_end,
                      deriv_end, depth - 1, roots);
  }

  /// Check if interval can be skipped (no roots possible).
  [[nodiscard]] bool CanSkipInterval(double value_start, double deriv_start,
                                     double value_end, double deriv_end,
                                     double interval_length) const {
    // Both values same sign and far from zero
    if (value_start * value_end <= 0) return false;

    const double min_abs_value =
        std::min(std::abs(value_start), std::abs(value_end));
    const double max_abs_deriv =
        std::max(std::abs(deriv_start), std::abs(deriv_end));

    // Function can't reach zero based on derivative bounds
    return min_abs_value > max_abs_deriv * interval_length * 1.2;
  }

  /// Check if interval might contain an extremum.
  [[nodiscard]] bool HasPotentialExtremum(double value_start,
                                          double deriv_start, double value_end,
                                          double deriv_end) const {
    // Derivatives have opposite signs
    if (deriv_start * deriv_end < 0) return true;

    // One derivative near zero
    if (std::abs(deriv_start) < tolerance_ ||
        std::abs(deriv_end) < tolerance_) {
      return true;
    }

    // Function value trending toward zero
    const bool approaching_zero =
        (value_start * deriv_start < 0) && (value_end * deriv_end > 0);
    return approaching_zero;
  }

  /// Split interval at estimated extremum location.
  void SplitAtExtremum(double param_start, double param_end, double value_start,
                       double deriv_start, double value_end, double deriv_end,
                       int depth, std::vector<double>& roots) const {
    // Estimate extremum location using linear interpolation of derivative
    double param_extremum = (param_start + param_end) * 0.5;
    if (std::abs(deriv_start - deriv_end) > kMachineEpsilon) {
      // Parameter where derivative ≈ 0 by linear interpolation
      param_extremum = param_start - deriv_start * (param_end - param_start) /
                                         (deriv_end - deriv_start);
      param_extremum = std::clamp(param_extremum,
                                  param_start + (param_end - param_start) * 0.1,
                                  param_end - (param_end - param_start) * 0.1);
    }

    const auto [value_extremum, deriv_extremum] = Evaluate(param_extremum);

    // Check if extremum is a root (tangent intersection)
    if (std::abs(value_extremum) < tolerance_ &&
        std::abs(deriv_extremum) < tolerance_ * 10) {
      if (roots.empty() ||
          std::abs(param_extremum - roots.back()) > tolerance_) {
        roots.push_back(param_extremum);
      }
    }

    // Recurse on both sides
    SubdivideInterval(param_start, param_extremum, value_start, deriv_start,
                      value_extremum, deriv_extremum, depth - 1, roots);
    SubdivideInterval(param_extremum, param_end, value_extremum, deriv_extremum,
                      value_end, deriv_end, depth - 1, roots);
  }

  /// Solve for root in bracketed interval using hybrid method.
  [[nodiscard]] std::optional<double> SolveInBracket(double bracket_left,
                                                     double bracket_right,
                                                     double value_left,
                                                     double value_right) const {
    double current_param = (bracket_left + bracket_right) * 0.5;

    for (int iter = 0; iter < kMaxNewtonIterations; ++iter) {
      const auto [func_value, derivative] = Evaluate(current_param);

      if (std::abs(func_value) < tolerance_) {
        return current_param;
      }

      double next_param;
      if (std::abs(derivative) > kMachineEpsilon) {
        next_param = current_param - func_value / derivative;
        if (next_param <= bracket_left || next_param >= bracket_right) {
          next_param = (bracket_left + bracket_right) * 0.5;
        }
      } else {
        next_param = (bracket_left + bracket_right) * 0.5;
      }

      const double value_next = EvaluateValue(next_param);
      if (value_left * value_next < 0) {
        bracket_right = next_param;
        value_right = value_next;
      } else {
        bracket_left = next_param;
        value_left = value_next;
      }

      if (std::abs(bracket_right - bracket_left) < kMachineEpsilon) {
        const double final_param = (bracket_left + bracket_right) * 0.5;
        return std::abs(EvaluateValue(final_param)) < tolerance_
                   ? std::optional{final_param}
                   : std::nullopt;
      }

      current_param = next_param;
    }

    const double final_param = (bracket_left + bracket_right) * 0.5;
    return std::abs(EvaluateValue(final_param)) < tolerance_
               ? std::optional{final_param}
               : std::nullopt;
  }

  const CurveType& curve_;
  ParameterInterval domain_;
  IsoLineType iso_type_;
  double iso_value_;
  double tolerance_;
};

// ============================================================================
// Utility Functions
// ============================================================================

/// Removes duplicate parameter values within tolerance.
///
/// @param[in,out] params Parameter vector to deduplicate (will be sorted)
/// @param tolerance Tolerance for considering values as duplicates
inline void RemoveDuplicateParameters(std::vector<double>& params,
                                      double tolerance = kDefaultTolerance) {
  if (params.empty()) return;

  std::sort(params.begin(), params.end());

  auto is_duplicate = [tolerance](double a, double b) {
    return std::abs(a - b) <= tolerance;
  };

  params.erase(std::unique(params.begin(), params.end(), is_duplicate),
               params.end());
}

// ============================================================================
// Main API Functions
// ============================================================================

/// Computes intersections between a curve and multiple iso-lines.
///
/// Internal helper function that processes a single axis direction.
///
/// @tparam CurveType The curve type (Curve2D or RationalCurve2D)
/// @param curve The curve to intersect
/// @param search_domain The parameter domain to search
/// @param iso_values List of iso-line constant values
/// @param iso_type Type of iso-lines (vertical or horizontal)
/// @param[out] result_params Output parameter values at intersections
/// @param tolerance Root finding tolerance
template <typename CurveType>
void ComputeIsolineIntersections(const CurveType& curve,
                                 const ParameterInterval& search_domain,
                                 const std::vector<double>& iso_values,
                                 IsoLineType iso_type,
                                 std::vector<double>& result_params,
                                 double tolerance = kDefaultTolerance) {
  for (const double iso_value : iso_values) {
    IsoIntersectionSolver<CurveType> solver(curve, search_domain, iso_type,
                                            iso_value, tolerance);
    std::vector<double> roots;
    solver.FindAllRoots(roots);

    for (const double root : roots) {
      if (search_domain.Contains(root, tolerance)) {
        result_params.push_back(root);
      }
    }
  }
}

/// Computes all intersections between a 2D NURBS curve and iso-parametric grid.
///
/// Main API function for curve-isoline intersection computation.
/// Supports both vertical (x=const) and horizontal (y=const) iso-lines.
///
/// @tparam CurveType The curve type (Curve2D or RationalCurve2D)
/// @param curve The NURBS curve to intersect
/// @param search_domain Optional search domain (nullptr = use full curve
/// domain)
/// @param vertical_lines Vertical iso-lines x = const (can be nullptr)
/// @param horizontal_lines Horizontal iso-lines y = const (can be nullptr)
/// @param[out] result_params Output parameters (automatically sorted and
/// deduplicated)
/// @param tolerance Convergence tolerance for root finding
template <typename CurveType>
void IntersectCurveWithIsolines(const CurveType& curve,
                                const ParameterInterval* search_domain,
                                const std::vector<double>* vertical_lines,
                                const std::vector<double>* horizontal_lines,
                                std::vector<double>& result_params,
                                double tolerance = kDefaultTolerance) {
  result_params.clear();

  const ParameterInterval domain =
      search_domain ? *search_domain : GetCurveDomain(curve);

  if (vertical_lines && !vertical_lines->empty()) {
    ComputeIsolineIntersections(curve, domain, *vertical_lines,
                                IsoLineType::kVertical, result_params,
                                tolerance);
  }

  if (horizontal_lines && !horizontal_lines->empty()) {
    ComputeIsolineIntersections(curve, domain, *horizontal_lines,
                                IsoLineType::kHorizontal, result_params,
                                tolerance);
  }

  RemoveDuplicateParameters(result_params, tolerance);
}

/// Computes intersections between a curve and a single iso-line.
///
/// Convenience wrapper for single iso-line intersection queries.
///
/// @tparam CurveType The curve type (Curve2D or RationalCurve2D)
/// @param curve The NURBS curve to intersect
/// @param iso_value The constant coordinate value of the iso-line
/// @param iso_type Direction: kVertical (x=const) or kHorizontal (y=const)
/// @param[out] result_params Output parameter values at intersections
/// @param tolerance Convergence tolerance
template <typename CurveType>
void IntersectCurveWithIsoline(const CurveType& curve, double iso_value,
                               IsoLineType iso_type,
                               std::vector<double>& result_params,
                               double tolerance = kDefaultTolerance) {
  result_params.clear();
  const ParameterInterval domain = GetCurveDomain(curve);
  IsoIntersectionSolver<CurveType> solver(curve, domain, iso_type, iso_value,
                                          tolerance);
  solver.FindAllRoots(result_params);
}

/// Evaluates curve points at multiple parameter values.
///
/// @tparam CurveType The curve type
/// @param curve The curve to evaluate
/// @param params Parameter values to evaluate at
/// @return Vector of 2D points
template <typename CurveType>
std::vector<Vec2> EvaluateCurveAtParameters(const CurveType& curve,
                                            const std::vector<double>& params) {
  std::vector<Vec2> points;
  points.reserve(params.size());
  for (double t : params) {
    points.push_back(EvaluateCurvePoint(curve, t));
  }
  return points;
}

}  // namespace nurbs_solver

#endif  // NURBS2D_ISO_INTERSECT_H_
