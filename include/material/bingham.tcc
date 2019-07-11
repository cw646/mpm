//! Constructor with material properties
template <unsigned Tdim>
mpm::Bingham<Tdim>::Bingham(unsigned id, const Json& material_properties)
    : Material<Tdim>(id, material_properties) {
  try {
    density_ = material_properties["density"].template get<double>();
    youngs_modulus_ =
        material_properties["youngs_modulus"].template get<double>();
    poisson_ratio_ =
        material_properties["poisson_ratio"].template get<double>();
    tau0_ = material_properties["tau0"].template get<double>();
    mu_ = material_properties["mu"].template get<double>();
    critical_shear_rate_ =
        material_properties["critical_shear_rate"].template get<double>();
    // Calculate bulk modulus
    bulk_modulus_ = youngs_modulus_ / (3.0 * (1. - 2. * poisson_ratio_));

    properties_ = material_properties;
  } catch (std::exception& except) {
    console_->error("Material parameter not set: {}\n", except.what());
  }
}

//! Compute pressure
template <unsigned Tdim>
double mpm::Bingham<Tdim>::thermodynamic_pressure(
    double volumetric_strain) const {
  // Bulk modulus
  return (-bulk_modulus_ * volumetric_strain);
}

//! Compute stress
template <unsigned Tdim>
Eigen::Matrix<double, 6, 1> mpm::Bingham<Tdim>::compute_stress(
    const Vector6d& stress, const Vector6d& dstrain,
    const ParticleBase<Tdim>* ptr, mpm::dense_map* state_vars, Eigdoub& floc,
    Eigdoub& rest_t, Eigdoub& alpha) {

  const unsigned phase = 0;

  // Get strain rate
  auto strain_rate = ptr->strain_rate(phase);

  // Convert strain rate to rate of deformation tensor
  strain_rate.tail(3) *= 0.5;

  // Set threshold for minimum critical shear rate
  const double shear_rate_threshold = 1.0E-15;
  if (critical_shear_rate_ < shear_rate_threshold)
    critical_shear_rate_ = shear_rate_threshold;

  // Rate of shear = sqrt(2 * D_ij * D_ij)
  // Since D (D_ij) is in Voigt notation (D_i), and the definition above is in
  // matrix, the last 3 components have to be doubled D_ij * D_ij = D_0^2 +
  // D_1^2 + D_2^2 + 2*D_3^2 + 2*D_4^2 + 2*D_5^2 Yielding is defined: rate of
  // shear > critical_shear_rate_^2 Checking yielding from strain rate vs
  // critical yielding shear rate
  double shear_rate = std::sqrt(2. * (strain_rate.dot(strain_rate)));

  double apparent_viscosity = 0.;

  Eigen::Matrix<double, 6, 1> tau;

  double a_thix = 155. / 30.;

  double tollerance = 0.00001;

  double dt = 10E-4;

  if (shear_rate * shear_rate > critical_shear_rate_ * critical_shear_rate_) {

    if (floc(phase) > tollerance) {  // is lambda greater than 0

      auto floc_prev = floc(phase);  // assing previous floculation state
      double dldt = -alpha(phase) * floc(phase) *
                    shear_rate;             // calculate change in lambda
      auto floccheck = floc(phase) + dldt;  // check to see if lambda positive

      if (floccheck > tollerance) {      // if lambda > 0
        floc(phase) = floc_prev + dldt;  // update floculation state
        alpha(phase) = (log(floc_prev) - log(floc(phase))) /
                       shear_rate;  // calculate alpha using previous and
                                    // current floculation state
        double tau0_temp =
            (1 + floc(phase)) * tau0_;  // caluclate new apparent yield stress
        apparent_viscosity = 2. * ((tau0_temp / shear_rate) + mu_);

        // console_->error("Lambda {}, alpha {}, tau0_temp {}, shear_rate{}",
        //  floc(phase), alpha(phase), tau0_temp, shear_rate);

      } else {
        floc(phase) = 0.;
        apparent_viscosity = 2. * ((tau0_ / shear_rate) + mu_);
      }

    } else {

      apparent_viscosity = 2. * ((tau0_ / shear_rate) + mu_);
    }

    tau = apparent_viscosity * strain_rate;

    rest_t(phase) = 0.;

  } else {

    apparent_viscosity = 5000.;

    tau = (2 * apparent_viscosity) * strain_rate;

    rest_t(phase) += dt;

    double tau_t = tau0_ + (a_thix * rest_t(phase));
    floc(phase) = (tau_t / tau0_) - 1;

    //  console_->error("Lambda {}, rest_t {}, tau_t {}, a_thix {}",
    //  floc(phase),
    //      rest_t(phase), tau_t, a_thix);
  }

  const Eigen::Matrix<double, 6, 1> updated_stress =
      -ptr->pressure(phase) * this->dirac_delta() + tau;

  return updated_stress;
}

//! Dirac delta 2D
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Bingham<2>::dirac_delta() const {

  return (Eigen::Matrix<double, 6, 1>() << 1.f, 1.f, 0.f, 0.f, 0.f, 0.f)
      .finished();
}

//! Dirac delta 3D
template <>
inline Eigen::Matrix<double, 6, 1> mpm::Bingham<3>::dirac_delta() const {

  return (Eigen::Matrix<double, 6, 1>() << 1.f, 1.f, 1.f, 0.f, 0.f, 0.f)
      .finished();
}
