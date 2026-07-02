#pragma once

#include <vector>

#include "preshaires/geometry.hpp"

namespace preshaires {

class MagneticFieldModel {
 public:
  virtual ~MagneticFieldModel() = default;

  [[nodiscard]] virtual MagneticFieldVector field_at(
      const Position& position, Length path_length) const = 0;
};

class UniformMagneticField final : public MagneticFieldModel {
 public:
  explicit UniformMagneticField(MagneticFieldVector field);

  [[nodiscard]] MagneticFieldVector field_at(
      const Position& position, Length path_length) const override;

 private:
  MagneticFieldVector field_;
};

struct TabulatedMagneticFieldNode {
  Length path_length;
  MagneticFieldVector field;
};

class TabulatedMagneticField final : public MagneticFieldModel {
 public:
  // Nodes are defined on finite, nonnegative, strictly increasing path lengths.
  // Exact node queries return the stored vector, including both endpoints.
  // Interior queries are linearly interpolated; outside queries throw.
  explicit TabulatedMagneticField(
      std::vector<TabulatedMagneticFieldNode> nodes);

  [[nodiscard]] MagneticFieldVector field_at(
      const Position& position, Length path_length) const override;

  [[nodiscard]] const std::vector<TabulatedMagneticFieldNode>& nodes() const {
    return nodes_;
  }

 private:
  std::vector<TabulatedMagneticFieldNode> nodes_;
};

}  // namespace preshaires
