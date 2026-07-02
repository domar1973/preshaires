#pragma once

namespace preshaires {

class RandomEngine {
 public:
  virtual ~RandomEngine() = default;

  // Contract: return a finite value U with 0 < U < 1.
  [[nodiscard]] virtual double uniform_open01() = 0;
};

}  // namespace preshaires
