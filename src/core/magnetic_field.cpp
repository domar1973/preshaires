#include "preshaires/magnetic_field.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace preshaires {
namespace {

void validate_field_vector(const MagneticFieldVector& field) {
  if (!is_finite(field.x) || !is_finite(field.y) || !is_finite(field.z)) {
    throw std::invalid_argument("Magnetic field components must be finite");
  }
}

void validate_node(const TabulatedMagneticFieldNode& node) {
  if (!is_finite(node.path_length) || node.path_length.meter < 0.0) {
    throw std::invalid_argument(
        "Tabulated field path lengths must be finite and nonnegative");
  }
  validate_field_vector(node.field);
}

}  // namespace

UniformMagneticField::UniformMagneticField(MagneticFieldVector field)
    : field_(field) {
  validate_field_vector(field_);
}

MagneticFieldVector UniformMagneticField::field_at(
    const Position&, Length path_length) const {
  if (!is_finite(path_length) || path_length.meter < 0.0) {
    throw std::invalid_argument("Path length must be finite and nonnegative");
  }
  return field_;
}

TabulatedMagneticField::TabulatedMagneticField(
    std::vector<TabulatedMagneticFieldNode> nodes)
    : nodes_(std::move(nodes)) {
  if (nodes_.size() < 2) {
    throw std::invalid_argument("Tabulated field requires at least two nodes");
  }
  validate_node(nodes_.front());
  for (std::size_t i = 1; i < nodes_.size(); ++i) {
    validate_node(nodes_[i]);
    if (nodes_[i].path_length.meter <= nodes_[i - 1].path_length.meter) {
      throw std::invalid_argument(
          "Tabulated field path lengths must be strictly increasing");
    }
  }
}

MagneticFieldVector TabulatedMagneticField::field_at(
    const Position&, Length path_length) const {
  if (!is_finite(path_length) || path_length.meter < 0.0) {
    throw std::invalid_argument("Path length must be finite and nonnegative");
  }
  if (path_length.meter < nodes_.front().path_length.meter ||
      path_length.meter > nodes_.back().path_length.meter) {
    throw std::out_of_range("Tabulated field does not extrapolate");
  }

  const auto exact = std::find_if(
      nodes_.begin(), nodes_.end(),
      [path_length](const TabulatedMagneticFieldNode& node) {
        return node.path_length.meter == path_length.meter;
      });
  if (exact != nodes_.end()) {
    return exact->field;
  }

  const auto upper = std::upper_bound(
      nodes_.begin(), nodes_.end(), path_length.meter,
      [](double value, const TabulatedMagneticFieldNode& node) {
        return value < node.path_length.meter;
      });
  const auto lower = upper - 1;

  const double span = upper->path_length.meter - lower->path_length.meter;
  const double weight = (path_length.meter - lower->path_length.meter) / span;
  const auto lerp = [weight](MagneticField a, MagneticField b) {
    return tesla(a.tesla + weight * (b.tesla - a.tesla));
  };

  return MagneticFieldVector{lerp(lower->field.x, upper->field.x),
                             lerp(lower->field.y, upper->field.y),
                             lerp(lower->field.z, upper->field.z)};
}

}  // namespace preshaires
