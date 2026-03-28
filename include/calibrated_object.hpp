#ifndef GM_PHD_INCLUDE_CALIBRATED_OBJECT_HPP_
#define GM_PHD_INCLUDE_CALIBRATED_OBJECT_HPP_

#include <any>
#include <concepts>
#include <expected>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
#include <utility>

namespace mot {
// Second element can not be reference or pointer
template <typename T>
struct is_reference_wrapper : std::false_type {};
template <typename U>
struct is_reference_wrapper<std::reference_wrapper<U>> : std::true_type {};
template <typename T>
concept ValueType =
    !std::is_reference_v<T> && !std::is_pointer_v<T> && !is_reference_wrapper<T>::value;

class CalibratedObject {
 public:
  CalibratedObject(void) = default;
  CalibratedObject(const CalibratedObject&) = default;
  CalibratedObject& operator=(const CalibratedObject&) = default;
  CalibratedObject(CalibratedObject&&) = default;
  CalibratedObject& operator=(CalibratedObject&&) = default;
  virtual ~CalibratedObject(void) = default;

  template <typename... Args>
  std::expected<void, std::string> SetCalibrations(Args&&... args) {
    static_assert(sizeof...(Args) % 2 == 0, "Arguments should be in pairs of name and value.");

    if (initialization_finished_) {
      return std::unexpected("Cannot set calibrations after initialization is finished.");
    }

    return SetCalibrationsImpl(std::make_index_sequence<sizeof...(Args) / 2>{},
                               std::forward<Args>(args)...);
  }

  template <typename T>
  std::expected<T, std::string> GetCalibration(const std::string& name) const {
    if (!calibrations_.contains(name)) {
      return std::unexpected("Invalid calibration name: " + name);
    }

    try {
      auto ref = std::any_cast<std::reference_wrapper<T>>(calibrations_.at(name));
      return ref.get();
    } catch (const std::bad_any_cast&) {
      return std::unexpected("Type mismatch for calibration: " + name);
    }
  }

 private:
  template <size_t... Is, typename... Args>
  std::expected<void, std::string> SetCalibrationsImpl(std::index_sequence<Is...>, Args&&... args) {
    auto tup = std::forward_as_tuple(std::forward<Args>(args)...);
    std::expected<void, std::string> result{};

    auto process_pair = [&](auto&& name, auto&& value) -> std::expected<void, std::string> {
      if (!result) {
        return result;
      }

      using NameT = std::decay_t<decltype(name)>;
      using ValueT = std::decay_t<decltype(value)>;

      static_assert(std::convertible_to<NameT, std::string_view>,
                    "Name argument must be convertible to std::string_view");
      static_assert(ValueType<ValueT>, "Value argument must satisfy ValueType concept");

      const auto key = std::string(std::string_view{name});
      if (!calibrations_.contains(key)) {
        return std::unexpected("Invalid calibration name: " + key);
      }

      try {
        auto ref = std::any_cast<std::reference_wrapper<ValueT>>(calibrations_[key]);
        ref.get() = static_cast<ValueT>(value);
      } catch (const std::bad_any_cast&) {
        return std::unexpected("Type mismatch for calibration: " + key);
      }

      return {};
    };

    (void)std::initializer_list<int>{
        (result = process_pair(std::get<2 * Is>(tup), std::get<2 * Is + 1>(tup)), 0)...};

    return result;
  }

 protected:
  void FinishInitialization(void) { initialization_finished_ = true; }
  bool CalibrationInitialized(void) const { return initialization_finished_; }
  void ResetCalibrations(void) { initialization_finished_ = false; }

  // Map to store calibration parameters, where the key is the parameter name and the value is a
  // pointer to the parameter value wrapped in std::any
  std::map<std::string, std::any> calibrations_{};
  bool initialization_finished_ = false;
};
}  // namespace mot

#endif  //  GM_PHD_INCLUDE_CALIBRATED_OBJECT_HPP_